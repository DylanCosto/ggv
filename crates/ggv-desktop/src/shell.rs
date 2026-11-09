//! The window: one or more panes side by side (split view). Each pane has its own position, zoom,
//! history, outlined base and windows; the open files and display options are shared, so a file
//! opened or a colour chosen in the active pane (the last one clicked) applies to every pane.

use std::path::PathBuf;

use eframe::egui::{self, Stroke};

use crate::{
    app::{App, Args},
    gpu, session,
    session::Session,
    tracks::Style,
};

/// At most this many panes; beyond that each is too narrow to read.
const MAX_PANES: usize = 4;
/// The active pane's outline, drawn just inside its edge.
const PANE_OUTLINE: f32 = 2.0;
/// A pane's header (its locus, Close) starts this far in from its edge, clear of the outline.
const PANE_LABEL_INSET: f32 = 6.0;

type Files = (Vec<PathBuf>, Vec<PathBuf>, Vec<PathBuf>, Vec<PathBuf>, Option<PathBuf>);

/// The workspace's files, which every pane shows (signal tracks too: a session saved from a
/// pane lacking one lost it).
fn files_of(s: &Session) -> Files {
    (s.bams.clone(), s.vcfs.clone(), s.genes.clone(), s.signals.clone(), s.reference.clone())
}

pub struct Shell {
    panes: Vec<App>,
    active: usize,
    /// What the panes were last synchronised to.
    files: Files,
    style: (Style, bool),
    /// Per-file filters, variant colours and track arrangement, which every pane shares too.
    file_settings: crate::app::FileSettings,
    bookmarks: crate::bookmarks::Shared,
    overlays: crate::overlays::Shared,
    /// Bookmarks seen last frame, to tell the person adding one.
    seen_bookmarks: usize,
    /// What the person did this frame, for agents waiting on `events`.
    pending_events: Vec<serde_json::Value>,
    /// Set once every pane has loaded: what happens before (a session restoring its outlined
    /// base, say) is the window starting up, not the person.
    armed: bool,
    /// Per pane: the selection last reported as an event (reported once the drag ends).
    reported_selection: Vec<serde_json::Value>,
    /// Panes scroll and zoom together, keeping their offsets.
    linked: bool,
    /// The active pane's view last frame (pane, contig, centre, points per base).
    last_view: Option<(usize, String, f64, f64)>,
    /// Where `ggv mcp --attach` reaches this window (none for scripted runs).
    attach: Option<crate::attach::Server>,
    /// Per pane: its outlined base and selection last frame, and which of them changed last
    /// (what the person pointed at most recently).
    marks: Vec<(serde_json::Value, serde_json::Value, &'static str)>,
    /// A reopened workspace's linking, applied once its panes have loaded (linking earlier made
    /// every pane follow the first one to load).
    relink: Option<bool>,
    /// When the workspace was last autosaved, and what was saved.
    autosaved: (std::time::Instant, String),
}

/// How often the whole window is saved to the last session, so a crash loses little.
const AUTOSAVE: std::time::Duration = std::time::Duration::from_secs(20);

impl Shell {
    pub fn new(cc: &eframe::CreationContext<'_>, args: Args, split: Vec<String>, linked: bool) -> anyhow::Result<Self> {
        let rs = cc
            .wgpu_render_state
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("ggv-desktop needs the wgpu renderer"))?;
        gpu::init(rs);
        crate::app::set_gpu_info(&rs.adapter.get_info());
        cc.egui_ctx.set_theme(if args.style.dark { egui::Theme::Dark } else { egui::Theme::Light });
        let scripted = args.automation.is_some();
        let attach = if scripted { None } else { crate::attach::Server::start(&cc.egui_ctx) };
        let mut shell = Self::with_pane(&cc.egui_ctx, args, linked, attach)?;
        for locus in split.into_iter().take(MAX_PANES - 1) {
            shell.split(&cc.egui_ctx, Some(locus));
        }
        // The first pane (the command line's locus) starts active.
        shell.set_active(0);
        Ok(shell)
    }

    /// The window with its first pane (no GPU needed: `new` sets that up).
    fn with_pane(ctx: &egui::Context, mut args: Args, linked: bool, attach: Option<crate::attach::Server>) -> anyhow::Result<Self> {
        let bookmarks = crate::bookmarks::Shared::default();
        args.bookmarks = Some(bookmarks.clone());
        let overlays = crate::overlays::Shared::default();
        args.overlays = Some(overlays.clone());
        let first = App::new(ctx, args)?;
        let files = files_of(&first.current_session());
        let style = first.shared_style();
        let file_settings = first.file_settings();
        let seen_bookmarks = bookmarks.borrow().len();
        Ok(Self {
            panes: vec![first],
            active: 0,
            files,
            style,
            file_settings,
            bookmarks,
            overlays,
            seen_bookmarks,
            pending_events: Vec::new(),
            armed: false,
            reported_selection: Vec::new(),
            linked,
            last_view: None,
            attach,
            marks: Vec::new(),
            autosaved: (std::time::Instant::now(), String::new()),
            relink: None,
        })
    }

    /// Adds a pane showing the active pane's files, at `locus` (or the active pane's view).
    fn split(&mut self, ctx: &egui::Context, locus: Option<String>) {
        self.split_at(ctx, locus.map(|l| session::PaneView { locus: Some(l), highlight: None }));
    }

    /// Adds a pane at `view` (or the active pane's view).
    fn split_at(&mut self, ctx: &egui::Context, view: Option<session::PaneView>) {
        if self.panes.len() >= MAX_PANES {
            return;
        }
        let mut s = self.panes[self.active].current_session();
        if let Some(v) = view {
            s.locus = v.locus;
            s.highlight = v.highlight;
        }
        match App::new(ctx, Args { session: Some(s), bookmarks: Some(self.bookmarks.clone()), overlays: Some(self.overlays.clone()), ..Default::default() }) {
            Ok(mut pane) => {
                pane.set_shared_style(self.style);
                self.panes.push(pane);
                self.set_active(self.panes.len() - 1);
            }
            Err(e) => eprintln!("ggv: new pane: {e:#}"),
        }
    }

    fn close(&mut self, i: usize) {
        if self.panes.len() < 2 {
            return;
        }
        let mut gone = self.panes.remove(i);
        // Its GPU buffers are freed by the next paint of a pane that remains.
        let slots = gone.take_slots();
        self.panes[0].retire(slots);
        let active = if self.active > i || self.active >= self.panes.len() { self.active.saturating_sub(1) } else { self.active };
        self.set_active(active.min(self.panes.len() - 1));
        // The closed pane's view isn't a move of the new active pane, nor are its marks the
        // next pane's: comparing them made linked panes jump and phantom person events.
        self.last_view = None;
        if i < self.reported_selection.len() {
            self.reported_selection.remove(i);
        }
        if i < self.marks.len() {
            self.marks.remove(i);
        }
    }

    fn shell_action(&mut self, a: crate::app::Action, ctx: &egui::Context) {
        use crate::app::Action;
        match a {
            Action::Split => self.split(ctx, None),
            Action::ClosePane => self.close(self.active),
            Action::LinkPanes => self.linked = !self.linked,
            _ => {}
        }
    }

    /// `panes`: list (default), split (a new pane at `locus`, else the active pane's view),
    /// close (`pane`, else the active one), activate (`pane`), link (`linked`: true or false).
    fn attach_panes(&mut self, p: &serde_json::Value, pane: usize, ctx: &egui::Context) -> Result<serde_json::Value, String> {
        let action = p.get("action").and_then(serde_json::Value::as_str).unwrap_or("list");
        match action {
            "list" => {}
            "split" => {
                if self.panes.len() >= MAX_PANES {
                    return Err(format!("the window already has {MAX_PANES} panes, the most it shows"));
                }
                let locus = p.get("locus").and_then(serde_json::Value::as_str).map(str::trim).filter(|l| !l.is_empty()).map(str::to_owned);
                // A place the pane can't show is refused, as `goto` refuses it, before splitting.
                if let Some(l) = &locus {
                    self.panes[pane].check_locus(l)?;
                }
                self.set_active(pane);
                self.split(ctx, locus);
            }
            "close" => {
                if self.panes.len() < 2 {
                    return Err("the window has one pane; it can't be closed".into());
                }
                self.close(pane);
            }
            "activate" => self.set_active(pane),
            "link" => {
                self.linked = p.get("linked").and_then(serde_json::Value::as_bool).ok_or("give `linked`: true or false")?;
                self.last_view = None;
            }
            other => return Err(format!("unknown panes action `{other}`: list, split, close, activate or link")),
        }
        let panes: Vec<serde_json::Value> = self.panes.iter().enumerate().map(|(i, p)| serde_json::json!({"pane": i + 1, "region": p.locus_label()})).collect();
        Ok(serde_json::json!({"count": self.panes.len(), "active_pane": self.active + 1, "linked": self.linked, "panes": panes}))
    }

    /// A session was opened in the active pane: that pane becomes the first, and the others
    /// are rebuilt from the session (none when it saved one pane).
    fn rebuild_panes(&mut self, ctx: &egui::Context, s: &Session) {
        let mut first = self.panes.remove(self.active);
        // Its place is fitted again at the width it ends up with.
        first.refit(s);
        let mut old: Vec<App> = std::mem::take(&mut self.panes);
        self.panes.push(first);
        for mut p in old.drain(..) {
            let slots = p.take_slots();
            self.panes[0].retire(slots);
        }
        self.set_active(0);
        for v in s.panes.iter().skip(1) {
            self.set_active(0);
            self.split_at(ctx, Some(v.clone()));
        }
        self.linked = false;
        self.relink = Some(s.linked);
        let active = s.active_pane.unwrap_or(1).clamp(1, self.panes.len()) - 1;
        self.set_active(active);
        self.last_view = None;
        self.marks.clear();
        self.reported_selection.clear();
    }

    /// The whole window as a session: the active pane's files and options, every pane's view.
    fn workspace(&self) -> Session {
        self.record_workspace();
        self.panes[self.active].current_session().with_workspace()
    }

    /// Tells sessions saved from a pane (Save session, agents) about the other panes.
    fn record_workspace(&self) {
        session::set_workspace(self.panes.iter().map(App::pane_view).collect(), self.relink.unwrap_or(self.linked), self.active);
    }

    /// Saves the window to the last session every so often (when it changed), so a crash or a
    /// power cut keeps nearly everything.
    fn autosave(&mut self) {
        if self.autosaved.0.elapsed() < AUTOSAVE || self.panes.iter().any(App::is_scripted) || !self.armed {
            return;
        }
        self.autosaved.0 = std::time::Instant::now();
        let s = self.workspace();
        if s.is_empty() {
            return;
        }
        let text = serde_json::to_string(&s).unwrap_or_default();
        if text != self.autosaved.1
            && let Some(path) = session::last_session_path()
            && s.save(&path).is_ok()
        {
            self.autosaved.1 = text;
        }
    }

    fn set_active(&mut self, i: usize) {
        self.active = i;
        for (j, p) in self.panes.iter_mut().enumerate() {
            p.is_active = j == i;
        }
    }

    /// With linked panes: moves every other pane by as much as the active one moved (pan and
    /// zoom), or to the same place when the active pane changed contig.
    fn follow(&mut self) {
        let Some((contig, center, ppb)) = self.panes[self.active].view_state() else { return };
        if self.linked
            && let Some((pane, prev_contig, prev_center, prev_ppb)) = self.last_view.clone()
            && pane == self.active
            && (prev_contig != contig || prev_center != center || prev_ppb != ppb)
        {
            for (j, p) in self.panes.iter_mut().enumerate() {
                if j == self.active {
                    continue;
                }
                if prev_contig != contig {
                    p.follow(&contig, center, ppb);
                } else if let Some((c, oc, op)) = p.view_state() {
                    p.follow(&c, oc + (center - prev_center), op * ppb / prev_ppb);
                }
            }
        }
        self.last_view = Some((self.active, contig, center, ppb));
    }

    /// The window's state for attached agents: every pane, which is active, whether linked.
    /// With `person`, changes since last time were the person's doing (agents' commands read
    /// the state right after acting), and become events: outlining a base, selecting or
    /// clearing a region, adding a bookmark, and triage decisions.
    fn attach_state(&mut self, person: bool, settled: bool) -> serde_json::Value {
        let mut events = Vec::new();
        self.reported_selection.resize(self.panes.len(), serde_json::Value::Null);
        let reported = &mut self.reported_selection;
        self.marks.resize(self.panes.len(), (serde_json::Value::Null, serde_json::Value::Null, ""));
        let panes: Vec<serde_json::Value> = self
            .panes
            .iter()
            .zip(&mut self.marks)
            .enumerate()
            .map(|(i, (p, mark))| {
                let mut v = p.attach_state();
                v["pane"] = serde_json::json!(i + 1);
                v["active"] = serde_json::json!(i == self.active);
                // Which the person set last: the outlined base or the selection.
                if v["selection"] != mark.1 && !v["selection"].is_null() {
                    mark.2 = "selection";
                } else if v["highlight"] != mark.0 && !v["highlight"].is_null() {
                    mark.2 = "highlight";
                }
                if person {
                    let at = |kind: &str, what: &str, value: &serde_json::Value| {
                        serde_json::json!({"kind": kind, "pane": i + 1, what: value, "view": v["view"]["region"]})
                    };
                    // A selection is reported once, when the drag that makes it ends.
                    if settled && v["selection"] != reported[i] {
                        events.push(if v["selection"].is_null() { at("unselect", "region", &reported[i]) } else { at("select", "region", &v["selection"]) });
                        reported[i] = v["selection"].clone();
                    }
                    if v["highlight"] != mark.0 && !v["highlight"].is_null() {
                        events.push(at("outline", "locus", &v["highlight"]));
                    }
                } else {
                    reported[i] = v["selection"].clone();
                }
                mark.0 = v["highlight"].clone();
                mark.1 = v["selection"].clone();
                v["latest"] = if mark.2.is_empty() { serde_json::Value::Null } else { serde_json::json!(mark.2) };
                v
            })
            .collect();
        self.pending_events.append(&mut events);
        serde_json::json!({
            "pid": std::process::id(),
            "version": env!("CARGO_PKG_VERSION"),
            "build": {"commit": env!("GGV_COMMIT"), "target": env!("GGV_TARGET")},
            "active_pane": self.active + 1,
            "linked": self.linked,
            "api": {"version": crate::app::API_VERSION, "capabilities": crate::app::CAPABILITIES},
            "panes": panes,
        })
    }

    /// Does what attached agents asked for (on this frame, so the view is consistent).
    fn attach_commands(&mut self, ctx: &egui::Context) {
        let Some(server) = &self.attach else { return };
        // Only commands still wanted: one the agent gave up on is dropped, never done late.
        let commands: Vec<crate::attach::Command> = server.commands.try_iter().filter(crate::attach::Command::take).collect();
        for c in commands {
            let pane = match c.params.get("pane").and_then(serde_json::Value::as_u64) {
                Some(n) if n >= 1 && (n as usize) <= self.panes.len() => n as usize - 1,
                Some(n) => {
                    let _ = c.reply.send(Err(format!("no pane {n}: this window has {}", self.panes.len())));
                    continue;
                }
                None => self.active,
            };
            // An export from this pane mustn't overwrite what another pane has open (its queue
            // list, decisions) any more than its own.
            let others: Vec<PathBuf> = self.panes.iter().enumerate().filter(|&(j, _)| j != pane).flat_map(|(_, p)| p.own_files()).collect();
            self.panes[pane].other_panes_files = others;
            let text = c.params.get("target").and_then(serde_json::Value::as_str).unwrap_or("").trim().to_owned();
            let before = self.panes[pane].attach_state()["display"].clone();
            // What the command returns besides the window's state.
            let mut result = None;
            // The panes themselves (and which is active) change: published like every other
            // command below, so the agent's next `state` agrees with the reply.
            let done = if c.method == "panes" {
                self.attach_panes(&c.params, pane, ctx).map(|v| result = Some(v))
            } else {
                let p = &mut self.panes[pane];
                let mut with = |r: Result<serde_json::Value, String>| {
                    r.map(|v| {
                        result = Some(v);
                    })
                };
                let done = match c.method.as_str() {
                    "goto" if !text.is_empty() => p.attach_goto(&text),
                    "select" if !text.is_empty() => p.attach_select(&text),
                    "unselect" => {
                        p.attach_unselect();
                        Ok(())
                    }
                    "open" => {
                        let files: Vec<PathBuf> = c.params.get("files").and_then(serde_json::Value::as_array).map(|a| a.iter().filter_map(|f| f.as_str()).map(PathBuf::from).collect()).unwrap_or_default();
                        let replace = c.params.get("replace").and_then(serde_json::Value::as_bool) == Some(true);
                        if files.is_empty() { Err("no files given".to_owned()) } else { p.attach_open(files, replace, ctx) }
                    }
                    "display" => p.attach_display(&c.params),
                    "bookmarks" => with(p.attach_bookmarks(&c.params)),
                    "overlays" => with(p.attach_overlays(&c.params)),
                    "save_session" => with(p.attach_save_session(&c.params)),
                    "tracks" => with(p.attach_tracks(&c.params)),
                    "review" => with(p.attach_review_command(&c.params)),
                    "queue" => with(p.attach_queue(&c.params)),
                    "search" => with(p.attach_search(&c.params, ctx)),
                    "compare" => with(p.attach_compare(&c.params, ctx)),
                    "rows" => with(p.attach_rows(&c.params)),
                    "action" => match c.params.get("name").and_then(serde_json::Value::as_str).map(str::trim) {
                        None | Some("" | "list") => with(Ok(App::action_list())),
                        Some(id) => with(p.attach_action(id, ctx)),
                    },
                    "goto" | "select" => Err("give a locus, region or gene".to_owned()),
                    other => Err(format!("unknown command `{other}`")),
                };
                // Navigating a pane makes it the one the user and the agent look at.
                if done.is_ok() {
                    self.set_active(pane);
                }
                done
            };
            // Window-level actions an agent ran (split, close pane, link): now, not next frame.
            if let Some(a) = self.panes[self.active].shell_request.take() {
                self.shell_action(a, ctx);
            }
            // Display changes spread to the other panes as a toolbar change would.
            self.sync(ctx);
            // Publish before replying, so the agent's next read sees the change.
            let mut state = self.attach_state(false, true);
            self.seen_bookmarks = self.bookmarks.borrow().len();
            if c.method == "display" {
                state["before"] = before;
            }
            if let Some(r) = result {
                state["result"] = r;
            }
            if let Some(server) = &self.attach {
                server.publish(state.clone(), false, ctx.input(|i| i.viewport().focused.unwrap_or(true)));
            }
            let _ = c.reply.send(done.map(|()| state));
        }
    }

    /// Spreads file and display changes made in the active pane to the others.
    fn sync(&mut self, ctx: &egui::Context) {
        let src = self.panes[self.active].current_session();
        let files = files_of(&src);
        if files != self.files {
            self.files = files;
            for (j, p) in self.panes.iter_mut().enumerate() {
                if j != self.active {
                    p.sync_files(&src, ctx);
                }
            }
        }
        let style = self.panes[self.active].shared_style();
        if style != self.style {
            self.style = style;
            for (j, p) in self.panes.iter_mut().enumerate() {
                if j != self.active {
                    p.set_shared_style(style);
                }
            }
        }
        let settings = self.panes[self.active].file_settings();
        if settings != self.file_settings {
            for (j, p) in self.panes.iter_mut().enumerate() {
                if j != self.active {
                    p.set_file_settings(&settings);
                }
            }
            self.file_settings = settings;
        }
        // A VCF or list reviewed in two panes: each shows the decisions made in the other.
        for p in &mut self.panes {
            p.refresh_decisions();
        }
    }
}

impl eframe::App for Shell {
    fn on_exit(&mut self) {
        // Scripted runs (--screenshot, --bench) must not replace the user's last session; they
        // report where each pane ended up instead.
        if self.panes.iter().any(App::is_scripted) {
            if self.panes.len() > 1 {
                for (i, p) in self.panes.iter().enumerate() {
                    println!("pane {}: {}", i + 1, p.locus_label().unwrap_or_default());
                }
            }
            return;
        }
        let s = self.workspace();
        if !s.is_empty()
            && let Some(path) = session::last_session_path()
        {
            let _ = s.save(&path);
        }
    }

    /// Agents' commands run here: eframe calls this even while the window is hidden or
    /// minimised (when `ui` doesn't run), so a command never waits for someone to click the
    /// window.
    fn logic(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.attach_commands(ctx);
        let hidden = ctx.input(|i| i.viewport().minimized == Some(true) || i.viewport().occluded == Some(true));
        if hidden && self.attach.is_some() {
            let ppp = ctx.pixels_per_point();
            for p in &mut self.panes {
                p.load_hidden(ppp);
            }
            // Keep publishing, so `ready` and the view are current while nothing is drawn.
            let state = self.attach_state(false, true);
            if let Some(server) = &self.attach {
                server.publish(state, false, false);
            }
            // Loads finish on other threads, which wake the window as they do.
        }
    }

    fn ui(&mut self, ui: &mut egui::Ui, frame: &mut eframe::Frame) {
        let ctx = ui.ctx().clone();
        // The controls follow the tracks' theme, however it was set: --dark, a session (at
        // startup or loaded later), the toolbar, the palette.
        let dark = self.panes[self.active].shared_style().0.dark;
        if ctx.theme() != (if dark { egui::Theme::Dark } else { egui::Theme::Light }) {
            ctx.set_theme(if dark { egui::Theme::Dark } else { egui::Theme::Light });
        }
        let (dropped, hovering): (Vec<PathBuf>, bool) = ctx.input(|i| {
            (i.raw.dropped_files.iter().map(|f| f.path().to_path_buf()).collect(), !i.raw.hovered_files.is_empty())
        });
        if !dropped.is_empty() {
            self.panes[self.active].open_files(dropped, &ctx);
        }
        if hovering {
            egui::Area::new(egui::Id::unique("drop-hint")).fixed_pos(egui::pos2(20.0, 60.0)).show(&ctx, |ui| {
                egui::Frame::popup(ui.style()).show(ui, |ui| {
                    ui.heading("Drop to open");
                    ui.label("BAM/CRAM, VCF/BCF, GTF/GFF/BED, FASTA reference, or a session (.json)");
                });
            });
        }

        let n = self.panes.len();
        let split = false;
        egui::Panel::top("toolbar").show(ui, |ui| {
            // Split is in the pane's toolbar; the pane controls get their own line only when
            // there is more than one pane.
            if n > 1 {
                ui.horizontal(|ui| {
                    ui.checkbox(&mut self.linked, "Link panes").on_hover_text("Scroll and zoom every pane together, keeping their offsets");
                    ui.label(egui::RichText::new(format!("pane {} of {n} (click a pane to make it active)", self.active + 1)).weak());
                });
            }
            self.panes[self.active].toolbar_ui(ui, &ctx);
        });
        egui::Panel::bottom("status").show(ui, |ui| self.panes[self.active].status_ui(ui, frame));
        if split {
            self.split(&ctx, None);
        }
        self.sync(&ctx);
        for p in &mut self.panes {
            p.windows(&ctx);
        }

        let n = self.panes.len();
        let mut clicked = None;
        let mut close = None;
        let accent = ui.visuals().selection.bg_fill;
        let fill = ui.visuals().panel_fill;
        egui::CentralPanel::default().frame(egui::Frame::NONE.fill(fill)).show(ui, |ui| {
            if n == 1 {
                self.panes[0].body(ui);
                return;
            }
            ui.columns(n, |cols| {
                for (i, (col, pane)) in cols.iter_mut().zip(&mut self.panes).enumerate() {
                    let rect = col.max_rect();
                    col.horizontal(|ui| {
                        // Clear of the active pane's outline (drawn inside the pane's edge).
                        ui.add_space(PANE_LABEL_INSET);
                        let label = pane.locus_label().unwrap_or_else(|| "loading…".to_owned());
                        let text = egui::RichText::new(label).monospace();
                        ui.label(if i == self.active { text.strong() } else { text.weak() });
                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                            ui.add_space(PANE_LABEL_INSET);
                            if ui.small_button("Close").on_hover_text("Close this pane").clicked() {
                                close = Some(i);
                            }
                        });
                    });
                    pane.body(col);
                    if i == self.active {
                        col.painter().rect_stroke(rect.shrink(1.0), 0.0, Stroke::new(PANE_OUTLINE, accent), egui::StrokeKind::Inside);
                    }
                    let pressed = col.input(|inp| inp.pointer.any_pressed().then(|| inp.pointer.interact_pos()).flatten());
                    if pressed.is_some_and(|p| rect.contains(p)) {
                        clicked = Some(i);
                    }
                }
            });
        });
        // Popups (context menus, combo boxes) take presses outside the pane they belong to.
        if let Some(i) = clicked
            && !ctx.any_popup_open()
            && i != self.active
        {
            self.set_active(i);
        }
        if let Some(i) = close {
            self.close(i);
        }
        // A session opened in a pane: rebuild the panes to its.
        if let Some((i, s)) = self.panes.iter_mut().enumerate().find_map(|(i, p)| p.take_pending_panes().map(|s| (i, s))) {
            self.set_active(i);
            self.rebuild_panes(&ctx, &s);
        }
        if let Some(l) = self.relink
            && self.panes.iter().all(App::is_ready)
        {
            self.linked = l;
            self.relink = None;
            self.last_view = None;
        }
        self.record_workspace();
        self.autosave();
        // Window-level actions asked for from a pane (palette, menus).
        if let Some(a) = self.panes[self.active].shell_request.take() {
            self.shell_action(a, &ctx);
        }
        self.follow();

        if self.attach.is_some() {
            let (used, focused) = ctx.input(|i| (!i.events.is_empty(), i.viewport().focused.unwrap_or(true)));
            let settled = !ctx.input(|i| i.pointer.any_down());
            let state = self.attach_state(self.armed, settled);
            self.armed |= self.panes.iter().all(App::is_ready);
            let marks = self.bookmarks.borrow();
            for b in marks.iter().skip(self.seen_bookmarks) {
                self.pending_events.push(serde_json::json!({"kind": "bookmark", "locus": b.locus(), "note": b.note}));
            }
            self.seen_bookmarks = marks.len();
            drop(marks);
            for p in &mut self.panes {
                self.pending_events.append(&mut p.person_events);
            }
            if !self.armed {
                self.pending_events.clear();
            }
            if let Some(server) = &self.attach {
                server.publish(state, used, focused);
                for e in self.pending_events.drain(..) {
                    server.push_event(e);
                }
            }
        }
        let others_ready = self.panes.iter().skip(1).all(App::is_ready);
        self.panes[0].automate_frame(&ctx, frame, others_ready);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::{Value, json};

    /// A window whose agent commands go in through the sender (nothing listens).
    fn window(ctx: &egui::Context, args: Args) -> (Shell, std::sync::mpsc::Sender<crate::attach::Command>) {
        let (server, tx) = crate::attach::Server::detached();
        (Shell::with_pane(ctx, args, false, Some(server)).unwrap(), tx)
    }

    /// Sends a command as an agent does; its reply, and then what `state` reads.
    fn command(shell: &mut Shell, ctx: &egui::Context, tx: &std::sync::mpsc::Sender<crate::attach::Command>, method: &str, params: Value) -> (Result<Value, String>, Value) {
        let (c, rx) = crate::attach::Command::new(method, params);
        tx.send(c).unwrap();
        shell.attach_commands(ctx);
        (rx.try_recv().unwrap(), shell.attach.as_ref().unwrap().snapshot())
    }

    fn test_bam() -> Option<PathBuf> {
        let bam = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio/HG002.revio.slices.bam");
        bam.exists().then_some(bam)
    }

    /// `panes` publishes like every command: the next `state` shows the panes it made or
    /// closed. A split at a place the pane can't show is refused.
    #[test]
    fn pane_changes_are_published_and_bad_loci_refused() {
        let Some(bam) = test_bam() else { return };
        let ctx = egui::Context::default();
        let (mut shell, tx) = window(&ctx, Args { bams: vec![bam], locus: Some("chr20:31180000-31181000".into()), ..Default::default() });
        let (r, state) = command(&mut shell, &ctx, &tx, "panes", json!({"action": "split", "locus": "chr20:31185000-31186000"}));
        assert_eq!(r.unwrap()["result"]["count"], 2);
        assert_eq!((state["panes"].as_array().unwrap().len(), &state["active_pane"]), (2, &json!(2)), "published: {state}");
        let (r, state) = command(&mut shell, &ctx, &tx, "panes", json!({"action": "split", "locus": "nonsense:1-2"}));
        assert!(r.unwrap_err().contains("nonsense"));
        assert_eq!((shell.panes.len(), state["panes"].as_array().unwrap().len()), (2, 2), "nothing split");
        let (r, state) = command(&mut shell, &ctx, &tx, "panes", json!({"action": "close"}));
        assert_eq!(r.unwrap()["result"]["count"], 1);
        assert_eq!(state["panes"].as_array().unwrap().len(), 1, "published: {state}");
    }

    /// Signal tracks are files of the workspace like the others: opened or closed in one pane,
    /// in every pane, so a session saved from any pane keeps them.
    #[test]
    fn signal_tracks_are_shared_across_panes() {
        let Some(bam) = test_bam() else { return };
        let dir = std::env::temp_dir().join(format!("ggv-signal-panes-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let sig = dir.join("sig.bedgraph");
        std::fs::write(&sig, "chr20\t31180000\t31181000\t5\nchr20\t31181000\t31182000\t9\n").unwrap();
        let ctx = egui::Context::default();
        let (mut shell, tx) = window(&ctx, Args { bams: vec![bam], ..Default::default() });
        command(&mut shell, &ctx, &tx, "panes", json!({"action": "split", "locus": "chr20:31180000-31182000"})).0.unwrap();
        command(&mut shell, &ctx, &tx, "open", json!({"pane": 1, "files": [sig.display().to_string()]})).0.unwrap();
        for p in &shell.panes {
            assert_eq!(p.current_session().signals, vec![sig.clone()]);
        }
        command(&mut shell, &ctx, &tx, "tracks", json!({"pane": 2, "action": "close", "track": "sig.bedgraph"})).0.unwrap();
        for p in &shell.panes {
            assert!(p.current_session().signals.is_empty());
        }
        drop(shell);
        let _ = std::fs::remove_dir_all(dir);
    }

    fn labels(shell: &Shell) -> Vec<Option<String>> {
        shell.panes.iter().map(App::locus_label).collect()
    }

    /// Three panes at three places, with `goto` in each.
    fn three_panes(ctx: &egui::Context, bam: PathBuf) -> (Shell, std::sync::mpsc::Sender<crate::attach::Command>) {
        let (mut shell, tx) = window(ctx, Args { bams: vec![bam], ..Default::default() });
        command(&mut shell, ctx, &tx, "panes", json!({"action": "split"})).0.unwrap();
        command(&mut shell, ctx, &tx, "panes", json!({"action": "split"})).0.unwrap();
        for (pane, locus) in [(1, "chr20:31180000-31181000"), (2, "chr20:31185000-31186000"), (3, "chr20:31190000-31191000")] {
            command(&mut shell, ctx, &tx, "goto", json!({"pane": pane, "target": locus})).0.unwrap();
        }
        (shell, tx)
    }

    /// Closing the active pane of linked panes moves none of the others (the closed pane's
    /// view was taken for a move of the next one).
    #[test]
    fn closing_a_linked_pane_moves_no_other() {
        let Some(bam) = test_bam() else { return };
        let ctx = egui::Context::default();
        let (mut shell, tx) = three_panes(&ctx, bam);
        command(&mut shell, &ctx, &tx, "panes", json!({"action": "link", "linked": true})).0.unwrap();
        command(&mut shell, &ctx, &tx, "panes", json!({"action": "activate", "pane": 1})).0.unwrap();
        shell.follow();
        let before = labels(&shell);
        command(&mut shell, &ctx, &tx, "panes", json!({"action": "close", "pane": 1})).0.unwrap();
        shell.follow();
        assert_eq!(labels(&shell), before[1..].to_vec());
        // Still linked: a pan of the active pane moves the other by as much.
        let (c, x, z) = shell.panes[0].view_state().unwrap();
        shell.panes[0].follow(&c, x + 1000.0, z);
        shell.follow();
        assert_eq!(labels(&shell)[1].as_deref(), Some("chr20:31,191,000-31,192,000"));
    }

    /// Closing a pane (its Close button) reports nothing the person did: the panes after it
    /// aren't compared with its outlined base and selection.
    #[test]
    fn closing_a_pane_makes_no_person_events() {
        let Some(bam) = test_bam() else { return };
        let ctx = egui::Context::default();
        let (mut shell, tx) = window(&ctx, Args { bams: vec![bam], ..Default::default() });
        command(&mut shell, &ctx, &tx, "panes", json!({"action": "split"})).0.unwrap();
        command(&mut shell, &ctx, &tx, "goto", json!({"pane": 1, "target": "chr20:31181500"})).0.unwrap();
        command(&mut shell, &ctx, &tx, "goto", json!({"pane": 2, "target": "chr20:31185500"})).0.unwrap();
        command(&mut shell, &ctx, &tx, "select", json!({"pane": 2, "target": "chr20:31185490-31185510"})).0.unwrap();
        shell.attach_state(true, true);
        shell.pending_events.clear();
        shell.close(0);
        let state = shell.attach_state(true, true);
        assert!(shell.pending_events.is_empty(), "{:?}", shell.pending_events);
        assert_eq!(state["panes"][0]["latest"], "selection", "its own, not the closed pane's");
    }

    /// Per-file filter expressions, variant site colours, Clear filter and track order are
    /// the workspace's: made in one pane, every pane has them, and a save from any pane keeps
    /// them.
    #[test]
    fn per_file_settings_are_shared_across_panes() {
        let Some(bam) = test_bam() else { return };
        let data = bam.parent().unwrap().to_path_buf();
        let vcfs = vec![data.join("HG002.revio.deepvariant.phased.vcf.gz"), data.join("HG002.revio.pbsv.vcf.gz")];
        let ctx = egui::Context::default();
        let (mut shell, tx) = window(&ctx, Args { bams: vec![bam.clone()], vcfs: vcfs.clone(), ..Default::default() });
        command(&mut shell, &ctx, &tx, "panes", json!({"action": "split"})).0.unwrap();
        command(&mut shell, &ctx, &tx, "display", json!({"pane": 1, "file": bam.display().to_string(), "expression": "hp == 1"})).0.unwrap();
        command(&mut shell, &ctx, &tx, "display", json!({"pane": 1, "variant_color": "SVLEN", "variant_file": "HG002.revio.pbsv.vcf.gz"})).0.unwrap();
        command(&mut shell, &ctx, &tx, "tracks", json!({"pane": 1, "action": "move", "track": "HG002.revio.pbsv.vcf.gz", "direction": "up"})).0.unwrap();
        for p in &shell.panes {
            let s = p.current_session();
            assert_eq!(s.track_filters.get(&bam).map(String::as_str), Some("hp == 1"));
            assert_eq!(s.variant_colors.get(&vcfs[1]).map(String::as_str), Some("SVLEN"));
            assert_eq!(s.vcfs, vec![vcfs[1].clone(), vcfs[0].clone()]);
        }
        command(&mut shell, &ctx, &tx, "panes", json!({"action": "activate", "pane": 2})).0.unwrap();
        assert_eq!(shell.workspace().track_filters.len(), 1, "saved from the other pane too");
        command(&mut shell, &ctx, &tx, "action", json!({"pane": 2, "name": crate::app::Action::ClearFilter.id()})).0.unwrap();
        assert!(shell.panes.iter().all(|p| p.current_session().track_filters.is_empty()));
    }

    /// A sort at a base belongs to the pane it was made in, as group-by-base does: the other
    /// pane sorts the same way at its own base.
    #[test]
    fn sort_by_base_stays_in_each_pane() {
        let Some(bam) = test_bam() else { return };
        let ctx = egui::Context::default();
        let (mut shell, tx) = window(&ctx, Args { bams: vec![bam], ..Default::default() });
        command(&mut shell, &ctx, &tx, "panes", json!({"action": "split"})).0.unwrap();
        command(&mut shell, &ctx, &tx, "goto", json!({"pane": 1, "target": "chr20:31180455"})).0.unwrap();
        command(&mut shell, &ctx, &tx, "goto", json!({"pane": 2, "target": "chr20:31185500"})).0.unwrap();
        let sorts = |shell: &Shell| shell.panes.iter().map(|p| p.current_session().sort.map(|s| (s.pos, s.by))).collect::<Vec<_>>();
        let base = |pos| Some((pos, "base".to_owned()));
        command(&mut shell, &ctx, &tx, "display", json!({"pane": 2, "sort": "base", "group": "base"})).0.unwrap();
        assert_eq!(sorts(&shell), vec![base(31180455), base(31185500)], "at each pane's own outlined base");
        command(&mut shell, &ctx, &tx, "display", json!({"pane": 1, "sort": "base", "group": "base"})).0.unwrap();
        assert_eq!(sorts(&shell), vec![base(31180455), base(31185500)]);
        command(&mut shell, &ctx, &tx, "display", json!({"pane": 1, "sort": "strand"})).0.unwrap();
        assert_eq!(sorts(&shell), vec![Some((31180455, "strand".to_owned())), Some((31185500, "strand".to_owned()))], "the kind is shared");
        command(&mut shell, &ctx, &tx, "display", json!({"pane": 2, "sort": "none"})).0.unwrap();
        assert_eq!(sorts(&shell), vec![None, None]);
    }

    /// The same list reviewed in two panes: neither pane's save drops the other's decisions,
    /// and each pane shows them.
    #[test]
    fn a_queue_in_two_panes_keeps_both_panes_decisions() {
        let Some(bam) = test_bam() else { return };
        let dir = std::env::temp_dir().join(format!("ggv-queue-panes-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let list = dir.join("list.bed");
        std::fs::write(&list, "chr20\t31180000\t31180010\ta\nchr20\t31185000\t31185010\tb\nchr20\t31190000\t31190010\tc\n").unwrap();
        let ctx = egui::Context::default();
        let (mut shell, tx) = window(&ctx, Args { bams: vec![bam], ..Default::default() });
        command(&mut shell, &ctx, &tx, "panes", json!({"action": "split"})).0.unwrap();
        for pane in [1, 2] {
            command(&mut shell, &ctx, &tx, "queue", json!({"pane": pane, "action": "open", "path": list.display().to_string()})).0.unwrap();
        }
        command(&mut shell, &ctx, &tx, "queue", json!({"pane": 1, "action": "decide", "decision": "pass"})).0.unwrap();
        // Pane 2 sees item 1 decided now, and skips it.
        let (r, _) = command(&mut shell, &ctx, &tx, "queue", json!({"pane": 2, "action": "next_undecided"}));
        assert_eq!(r.unwrap()["result"]["position"], 2);
        let (r, _) = command(&mut shell, &ctx, &tx, "queue", json!({"pane": 2, "action": "decide", "decision": "fail"}));
        let r = r.unwrap();
        assert_eq!(r["result"]["counts"], json!({"pass": 1, "fail": 1, "unsure": 0}), "pane 2 kept pane 1's decision");
        assert_eq!(r["result"]["position"], 3, "and went past it");
        let (r, _) = command(&mut shell, &ctx, &tx, "queue", json!({"pane": 1}));
        assert_eq!(r.unwrap()["result"]["counts"], json!({"pass": 1, "fail": 1, "unsure": 0}), "pane 1 shows pane 2's");
        command(&mut shell, &ctx, &tx, "queue", json!({"pane": 1, "action": "go", "index": 3})).0.unwrap();
        command(&mut shell, &ctx, &tx, "queue", json!({"pane": 1, "action": "decide", "decision": "unsure"})).0.unwrap();
        let saved = std::fs::read_to_string(crate::store::sidecar(&list, "ggv-queue.tsv")).unwrap();
        assert_eq!(["\tpass\t", "\tfail\t", "\tunsure\t"].map(|d| saved.matches(d).count()), [1, 1, 1], "{saved}");
        drop(shell);
        let _ = std::fs::remove_dir_all(dir);
    }

    /// An export from one pane is refused over a file another pane has open: its queue list
    /// and that list's decisions.
    #[test]
    fn exports_refuse_another_panes_queue_files() {
        let Some(bam) = test_bam() else { return };
        let dir = std::env::temp_dir().join(format!("ggv-export-panes-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let list = dir.join("list.bed");
        let mine = "chr20\t31180000\t31180010\tmine\n";
        std::fs::write(&list, mine).unwrap();
        let ctx = egui::Context::default();
        let (mut shell, tx) = window(&ctx, Args { bams: vec![bam], ..Default::default() });
        command(&mut shell, &ctx, &tx, "panes", json!({"action": "split"})).0.unwrap();
        command(&mut shell, &ctx, &tx, "queue", json!({"pane": 2, "action": "open", "path": list.display().to_string()})).0.unwrap();
        command(&mut shell, &ctx, &tx, "queue", json!({"pane": 2, "action": "decide", "decision": "pass"})).0.unwrap();
        let decisions = crate::store::sidecar(&list, "ggv-queue.tsv");
        command(&mut shell, &ctx, &tx, "queue", json!({"pane": 1, "action": "paste", "text": "chr20:31185000\n"})).0.unwrap();
        for (method, path) in [("bookmarks", &list), ("queue", &list), ("queue", &decisions)] {
            let err = command(&mut shell, &ctx, &tx, method, json!({"pane": 1, "action": "export", "path": path.display().to_string()})).0.unwrap_err();
            assert!(err.contains("open in GGV"), "{method} {}: {err}", path.display());
        }
        assert_eq!(std::fs::read_to_string(&list).unwrap(), mine);
        assert!(std::fs::read_to_string(&decisions).unwrap().contains("\tpass\t"));
        drop(shell);
        let _ = std::fs::remove_dir_all(dir);
    }

    /// Clear filter forgets the filters remembered for the files: a split (built from the
    /// session, which applies remembered filters) doesn't bring them back.
    #[test]
    fn a_cleared_filter_stays_cleared_after_a_split() {
        let Some(data) = test_bam() else { return };
        let dir = std::env::temp_dir().join(format!("ggv-clear-filter-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        // A copy of its own: what's remembered is per file.
        let bam = dir.join("clear.bam");
        std::fs::copy(&data, &bam).unwrap();
        std::fs::copy(data.with_added_extension("bai"), bam.with_added_extension("bai")).unwrap();
        // Remembered filters in a folder of its own: tests run in parallel, and the person's
        // settings are not for tests.
        let config = dir.join("config");
        std::fs::create_dir_all(&config).unwrap();
        crate::store::TEST_CONFIG.with_borrow_mut(|c| *c = Some(config.clone()));
        let remembered = || {
            let text = std::fs::read_to_string(config.join("file-filters.json")).unwrap_or_default();
            serde_json::from_str::<std::collections::BTreeMap<PathBuf, String>>(&text).unwrap_or_default()
        };
        let ctx = egui::Context::default();
        let (mut shell, tx) = window(&ctx, Args { bams: vec![bam.clone()], ..Default::default() });
        command(&mut shell, &ctx, &tx, "display", json!({"file": bam.display().to_string(), "expression": "hp == 1"})).0.unwrap();
        assert!(remembered().contains_key(&bam), "remembered");
        command(&mut shell, &ctx, &tx, "action", json!({"name": crate::app::Action::ClearFilter.id()})).0.unwrap();
        assert!(!remembered().contains_key(&bam), "forgotten");
        command(&mut shell, &ctx, &tx, "panes", json!({"action": "split"})).0.unwrap();
        assert!(shell.panes.iter().all(|p| p.current_session().track_filters.is_empty()));
        assert!(shell.workspace().track_filters.is_empty());
        drop(shell);
        crate::store::TEST_CONFIG.with_borrow_mut(|c| *c = None);
        let _ = std::fs::remove_dir_all(dir);
    }

    /// Group-by-base is per pane, like a sort's position: making another pane active changes
    /// neither pane's grouping.
    #[test]
    fn activating_a_pane_keeps_every_panes_grouping() {
        let Some(bam) = test_bam() else { return };
        let ctx = egui::Context::default();
        let (mut shell, tx) = window(&ctx, Args { bams: vec![bam], ..Default::default() });
        command(&mut shell, &ctx, &tx, "panes", json!({"action": "split"})).0.unwrap();
        command(&mut shell, &ctx, &tx, "goto", json!({"pane": 1, "target": "chr20:31180455"})).0.unwrap();
        command(&mut shell, &ctx, &tx, "goto", json!({"pane": 2, "target": "chr20:31185500"})).0.unwrap();
        let groups = |shell: &Shell| shell.panes.iter().map(|p| p.shared_style().0.group).collect::<Vec<_>>();
        command(&mut shell, &ctx, &tx, "display", json!({"pane": 1, "group": "base"})).0.unwrap();
        let both = vec![ggv_layout::GroupBy::BaseAt(31_180_454), ggv_layout::GroupBy::BaseAt(31_185_499)];
        assert_eq!(groups(&shell), both, "each at its own outlined base");
        for pane in [2, 1, 2] {
            command(&mut shell, &ctx, &tx, "panes", json!({"action": "activate", "pane": pane})).0.unwrap();
            shell.sync(&ctx);
            assert_eq!(groups(&shell), both, "after activating pane {pane}");
        }
        command(&mut shell, &ctx, &tx, "display", json!({"pane": 2, "group": "none"})).0.unwrap();
        assert_eq!(groups(&shell), vec![ggv_layout::GroupBy::None; 2], "turning it off is shared");
    }
}
