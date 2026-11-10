//! The side column: Evidence, Region, Variant review and Review queue as stacked sections
//! instead of separate floating windows over the reads. Each section can be collapsed, pinned
//! (kept open, with a hint, when it has nothing to show), popped out into a floating window and
//! docked back, and dragged by its title to another place in the column. The layout is kept in
//! `prefs.json`.

use std::collections::BTreeMap;

use eframe::egui::{self, RichText};
use serde::{Deserialize, Serialize};

use super::App;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum Section {
    Evidence,
    Region,
    Review,
    Queue,
}

impl Section {
    pub const ALL: [Section; 4] = [Section::Evidence, Section::Region, Section::Review, Section::Queue];

    pub fn title(self) -> &'static str {
        match self {
            Section::Evidence => "Evidence",
            Section::Region => "Region",
            Section::Review => "Variant review",
            Section::Queue => "Review queue",
        }
    }

    /// What a pinned section says when it has nothing to show.
    fn idle_hint(self) -> &'static str {
        match self {
            Section::Evidence => "Counts at the outlined base or the view centre.",
            Section::Region => "No region selected: shift-drag across bases, or click a base and shift-click another.",
            Section::Review => "Not reviewing: Variants > Review.",
            Section::Queue => "No queue open: View > Review queue.",
        }
    }
}

/// How long a review or queue that needs the user is kept scrolled into view, in seconds.
const ATTENTION_SECS: f64 = 1.5;

/// A popped-out section's width.
const FLOAT_WIDTH: f32 = 380.0;

/// Where a section is and how it looks.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Place {
    /// In its own window rather than the column.
    #[serde(default)]
    pub floating: bool,
    /// Shown even when it has nothing to show.
    #[serde(default)]
    pub pinned: bool,
    #[serde(default)]
    pub collapsed: bool,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Layout {
    pub order: Vec<Section>,
    #[serde(default)]
    pub places: BTreeMap<Section, Place>,
}

impl Default for Layout {
    /// Variant review and the queue first: their buttons are what a review needs on screen, and
    /// Evidence (one block per alignment track) grows with every BAM, so below it they would
    /// start out of sight.
    fn default() -> Self {
        Self { order: vec![Section::Review, Section::Queue, Section::Evidence, Section::Region], places: BTreeMap::new() }
    }
}

/// Where a window first opens over `pane` (a pane's tracks): at its top, on the side away from
/// `avoid_x` (the outlined base, or the view centre where a jump puts what it found), held by
/// its outer corner (the returned pivot) so that however wide it turns out it hugs the pane's
/// edge rather than growing over what it avoids. `height` keeps its top within `screen`.
pub(super) fn beside(pane: egui::Rect, avoid_x: f32, height: f32, screen: egui::Rect) -> (egui::Pos2, egui::Align2) {
    const MARGIN: f32 = 8.0;
    let pane = if pane.is_positive() && pane.is_finite() { pane.intersect(screen) } else { screen };
    let pane = if pane.is_positive() { pane } else { screen };
    let y = (pane.top() + MARGIN).min(screen.bottom() - height).max(screen.top());
    if avoid_x <= pane.center().x {
        (egui::pos2(pane.right() - MARGIN, y), egui::Align2::RIGHT_TOP)
    } else {
        (egui::pos2(pane.left() + MARGIN, y), egui::Align2::LEFT_TOP)
    }
}

/// Scroll bars in `ui` that show whenever there is more than fits: solid (drawn beside the
/// content, not only while hovered as egui's default floating ones) and in the text colour (the
/// solid style's own light grey handle is hardly visible on a light panel).
pub(super) fn visible_scroll_bars(ui: &mut egui::Ui) {
    ui.spacing_mut().scroll = egui::style::ScrollStyle { foreground_color: true, bar_width: 8.0, ..egui::style::ScrollStyle::solid() };
}

/// `window`, with `id`, first opened at a spot from `beside`. It is held there (its corner
/// at the spot, by its size as drawn) until its pane has been laid out (`settled`: windows are
/// drawn before the panes, so on GGV's first frame the pane isn't known) and its own size has
/// settled; from then on it stays wherever it is put. The hold writes egui's stored position
/// itself: a window dragged by its title ignores `current_pos` (egui 0.37), and egui's own
/// pivot isn't used (a right-pivoted window jumps by its width when GGV's window is resized).
pub(super) fn placed<'a>(ctx: &egui::Context, window: egui::Window<'a>, id: egui::Id, (pos, pivot): (egui::Pos2, egui::Align2), settled: bool) -> egui::Window<'a> {
    /// Frames to hold it, for its size to settle (its first is an unseen sizing pass).
    const HOLD: u8 = 3;
    let key = id.with("placed");
    let window = window.id(id);
    let held = ctx.data(|d| d.get_temp::<u8>(key)).unwrap_or(0);
    if held >= HOLD {
        return window;
    }
    if settled {
        ctx.data_mut(|d| d.insert_temp(key, held + 1));
    }
    let size = ctx.memory(|m| m.area_rect(id)).map_or(egui::Vec2::ZERO, |r| r.size());
    let left_top = egui::pos2(pos.x - pivot.x().to_factor() * size.x, pos.y - pivot.y().to_factor() * size.y);
    if let Some(mut state) = egui::AreaState::load(ctx, id) {
        state.pivot = egui::Align2::LEFT_TOP;
        state.pivot_pos = Some(left_top);
        ctx.memory_mut(|m| m.areas_mut().set_state(egui::LayerId::new(egui::Order::Middle, id), state));
    }
    window.default_pos(left_top)
}

impl Layout {
    /// The saved layout, with any section it lacks added at the end.
    pub fn load() -> Self {
        let mut l: Layout = super::tools::load_prefs()
            .get("inspector")
            .and_then(|v| serde_json::from_value(v.clone()).ok())
            .unwrap_or_default();
        l.order.dedup();
        for s in Section::ALL {
            if !l.order.contains(&s) {
                l.order.push(s);
            }
        }
        l
    }

    fn save(&self) {
        let mut prefs = super::tools::load_prefs();
        if let Ok(v) = serde_json::to_value(self) {
            prefs.insert("inspector".into(), v);
            super::tools::save_prefs(&prefs);
        }
    }

    pub fn place(&self, s: Section) -> Place {
        self.places.get(&s).copied().unwrap_or_default()
    }

    fn set(&mut self, s: Section, p: Place) {
        self.places.insert(s, p);
        self.save();
    }

    /// Moves `s` to just before `before`.
    fn move_before(&mut self, s: Section, before: Section) {
        if s == before {
            return;
        }
        self.order.retain(|x| *x != s);
        let i = self.order.iter().position(|x| *x == before).unwrap_or(self.order.len());
        self.order.insert(i, s);
        self.save();
    }
}

impl App {
    /// This pane's tracks have been laid out (so `pane_rect` is theirs), at a real size (not
    /// in a minimized window, whose panes are a sliver).
    pub(super) fn pane_settled(&self) -> bool {
        self.canvas_rect.is_finite() && self.canvas_rect.width() >= 100.0 && self.canvas_rect.height() >= 100.0
    }

    /// This pane's tracks on screen (as last drawn), or GGV's whole window before they are.
    pub(super) fn pane_rect(&self, ctx: &egui::Context) -> egui::Rect {
        if self.pane_settled() { self.canvas_rect } else { ctx.content_rect() }
    }

    /// Where a popped-out section first opens: over this pane, clear of the base
    /// Evidence describes (the outlined base, else the view centre).
    pub(super) fn window_pos(&self, ctx: &egui::Context, height: f32) -> (egui::Pos2, egui::Align2) {
        let pane = self.pane_rect(ctx);
        beside(pane, self.base_x(ctx).unwrap_or(pane.center().x), height, ctx.content_rect())
    }

    /// Where a window first opens: over this pane, clear of the view's centre, where a
    /// jump puts what was found (a search hit, a read, a difference).
    pub(super) fn centre_window_pos(&self, ctx: &egui::Context, height: f32) -> (egui::Pos2, egui::Align2) {
        let pane = self.pane_rect(ctx);
        beside(pane, pane.center().x, height, ctx.content_rect())
    }

    /// The screen x of the base Evidence describes (the outlined base, else the view centre).
    pub(super) fn base_x(&self, ctx: &egui::Context) -> Option<f32> {
        let (_, pos) = self.evidence_base()?;
        let (vs, _) = self.visible();
        Some(self.x_of(self.pane_rect(ctx), vs, pos as f64 + 0.5))
    }

    /// Something in a review or queue section the user has to see now (it started, a site is
    /// waiting for a decision, a save failed): when this changes, the column scrolls to it.
    fn attention(&self, s: Section) -> Option<u64> {
        use std::hash::{Hash, Hasher};
        let mut h = std::collections::hash_map::DefaultHasher::new();
        match s {
            Section::Review => {
                let r = self.review.as_ref()?;
                (r.path.as_path(), self.highlight.is_some(), r.record, &r.error, &r.notice).hash(&mut h);
            }
            Section::Queue => {
                let q = self.queue.as_ref()?;
                (&q.source, q.pos, &q.error, &q.notice).hash(&mut h);
            }
            Section::Evidence | Section::Region => return None,
        }
        Some(h.finish())
    }

    /// Whether `s` has something to show now.
    fn section_active(&self, s: Section) -> bool {
        match s {
            Section::Evidence => self.show_evidence,
            Section::Region => self.selection.is_some() && !self.selecting,
            Section::Review => self.review.is_some(),
            Section::Queue => self.queue.is_some(),
        }
    }

    fn section_shown(&self, s: Section) -> bool {
        !self.contigs.is_empty() && (self.section_active(s) || self.inspector.place(s).pinned)
    }

    /// The column has something docked to show.
    pub(super) fn column_wanted(&self) -> bool {
        Section::ALL.iter().any(|&s| self.section_shown(s) && !self.inspector.place(s).floating)
    }

    fn section_body(&mut self, s: Section, ui: &mut egui::Ui) {
        if !self.section_active(s) {
            ui.label(RichText::new(s.idle_hint()).weak());
            return;
        }
        match s {
            Section::Evidence => self.evidence_panel(ui),
            Section::Region => self.region_ui(ui),
            Section::Review => self.review_ui(ui),
            Section::Queue => self.queue_ui(ui),
        }
    }

    /// Closing a section (its window's ×, or Close in the column) ends what it shows.
    fn close_section(&mut self, s: Section) {
        let mut p = self.inspector.place(s);
        p.pinned = false;
        self.inspector.set(s, p);
        match s {
            Section::Evidence => self.show_evidence = false,
            Section::Region => {
                self.selection = None;
                for t in &mut self.tracks {
                    t.summary = None;
                    t.summary_for = None;
                }
            }
            Section::Review => self.review = None,
            Section::Queue => self.queue = None,
        }
    }

    /// The side column: the docked sections in a scroll area whose bar shows whenever there is
    /// more than fits (egui's default floating bar is invisible until hovered, so nothing said the
    /// review's buttons were below the bottom).
    pub(super) fn column(&mut self, ui: &mut egui::Ui) {
        visible_scroll_bars(ui);
        egui::ScrollArea::vertical()
            .id_salt(("side-column", self.reference_slot))
            .scroll_bar_visibility(egui::scroll_area::ScrollBarVisibility::VisibleWhenNeeded)
            .show(ui, |ui| self.inspector_column(ui));
    }

    /// The docked sections, in order.
    pub(super) fn inspector_column(&mut self, ui: &mut egui::Ui) {
        let order = self.inspector.order.clone();
        let mut moved: Option<(Section, Section)> = None;
        let mut first = true;
        for s in order {
            let place = self.inspector.place(s);
            if place.floating || !self.section_shown(s) {
                // Shown again later (a new review): it is new then.
                ui.data_mut(|d| d.remove::<(Option<u64>, f64)>(egui::Id::unique(("column-attention", s, self.reference_slot))));
                continue;
            }
            if !first {
                ui.add_space(4.0);
                ui.separator();
            }
            first = false;
            // The section on the panel's own background (egui's drop zone filled it grey), and
            // a drop target only while a section title is being dragged over it.
            let area = ui
                .scope(|ui| {
                    self.section_header(s, ui, true);
                    if !self.inspector.place(s).collapsed {
                        ui.push_id(("section", s), |ui| self.section_body(s, ui));
                    }
                })
                .response
                .rect;
            // A review or queue that just started, or now needs the user (a site to decide, a
            // failed save): scrolled into view. A failed save's message is at the section's end.
            let attention = self.attention(s);
            let seen_id = egui::Id::unique(("column-attention", s, self.reference_slot));
            let now = ui.input(|i| i.time);
            if attention.is_some() && ui.data(|d| d.get_temp::<(Option<u64>, f64)>(seen_id)).map(|(a, _)| a) != Some(attention) {
                ui.data_mut(|d| d.insert_temp(seen_id, (attention, now)));
            }
            // Kept in view for a moment after: what loads meanwhile (Evidence for the new site,
            // above it) would otherwise push it back out.
            if let Some((_, since)) = ui.data(|d| d.get_temp::<(Option<u64>, f64)>(seen_id))
                && attention.is_some()
                && now - since < ATTENTION_SECS
            {
                let warn = match s {
                    Section::Review => self.review.as_ref().is_some_and(|r| r.error.is_some() || r.notice.is_some()),
                    Section::Queue => self.queue.as_ref().is_some_and(|q| q.error.is_some() || q.notice.is_some()),
                    _ => false,
                };
                let target = if warn { egui::Rect::from_min_max(egui::pos2(area.left(), (area.bottom() - 80.0).max(area.top())), area.max) } else { area };
                ui.scroll_to_rect(target, None);
                ui.ctx().request_repaint();
            }
            let zone = ui.interact(area, egui::Id::unique(("section-drop", s, self.reference_slot)), egui::Sense::hover());
            if zone.dnd_hover_payload::<Section>().is_some_and(|d| *d != s) {
                ui.painter().line_segment(
                    [area.left_top() - egui::vec2(0.0, 3.0), area.right_top() - egui::vec2(0.0, 3.0)],
                    egui::Stroke::new(2.0, ui.visuals().selection.stroke.color),
                );
            }
            if let Some(d) = zone.dnd_release_payload::<Section>() {
                moved = Some((*d, s));
            }
        }
        if let Some((s, before)) = moved {
            self.inspector.move_before(s, before);
        }
    }

    /// A section's title row: collapse, its name (drag it to reorder), pin, pop out or dock.
    fn section_header(&mut self, s: Section, ui: &mut egui::Ui, docked: bool) {
        let mut p = self.inspector.place(s);
        let mut changed = false;
        let mut close = false;
        // Too narrow for the title and the buttons on one line (a pane of a four-way split):
        // the buttons go on their own line rather than over the title.
        let narrow = ui.available_width() < 240.0;
        if docked {
            ui.horizontal(|ui| {
                if ui.small_button(if p.collapsed { "+" } else { "−" }).on_hover_text(if p.collapsed { "Expand" } else { "Collapse" }).clicked() {
                    p.collapsed = !p.collapsed;
                    changed = true;
                }
                ui.dnd_drag_source(egui::Id::unique(("section-drag", s, self.reference_slot)), s, |ui| {
                    ui.label(RichText::new(s.title()).strong());
                })
                .response
                .on_hover_text("Drag to move it up or down the column");
                if !narrow {
                    self.section_buttons(s, ui, docked, &mut p, &mut changed, &mut close);
                }
            });
        }
        if narrow || !docked {
            ui.horizontal(|ui| self.section_buttons(s, ui, docked, &mut p, &mut changed, &mut close));
        }
        if changed {
            self.inspector.set(s, p);
        }
        if close {
            self.close_section(s);
        }
    }

    /// A section's Pin, Pop out / Dock and close buttons, right-aligned.
    fn section_buttons(&mut self, s: Section, ui: &mut egui::Ui, docked: bool, p: &mut Place, changed: &mut bool, close: &mut bool) {
        let _ = s;
        {
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if docked && ui.small_button("×").on_hover_text("Close").clicked() {
                    *close = true;
                }
                let float_label = if docked { "Pop out" } else { "Dock" };
                let float_hint = if docked { "Show it in a floating window over the tracks, which you can drag anywhere in GGV's window" } else { "Put it back in the side column" };
                if ui.small_button(float_label).on_hover_text(float_hint).clicked() {
                    p.floating = docked;
                    *changed = true;
                }
                if ui
                    .selectable_label(p.pinned, "Pin")
                    .on_hover_text("Keep it open even when it has nothing to show (no selection, no review)")
                    .clicked()
                {
                    p.pinned = !p.pinned;
                    *changed = true;
                }
            });
        }
    }

    /// Popped-out sections, each in its own window.
    pub(super) fn floating_sections(&mut self, ctx: &egui::Context) {
        if !self.is_active {
            return;
        }
        for s in self.inspector.order.clone() {
            if !self.inspector.place(s).floating || !self.section_shown(s) {
                continue;
            }
            // Narrower in a narrower pane (a split), to fit beside the base at its centre.
            let width = (self.pane_rect(ctx).width() / 2.0 - 40.0).clamp(260.0, FLOAT_WIDTH);
            let mut open = true;
            placed(ctx, egui::Window::new(s.title()), egui::Id::unique((s.title(), self.reference_slot)), self.window_pos(ctx, 300.0), self.pane_settled())
                // First shown over this pane's tracks on the side away from the outlined base (or
                // the view centre), so it doesn't hide the reads it counts; then wherever it was
                // dragged.
                .default_width(width)
                .open(&mut open)
                .show(ctx, |ui| {
                    self.section_header(s, ui, false);
                    ui.separator();
                    egui::ScrollArea::vertical().max_height(ctx.content_rect().height() * 0.7).show(ui, |ui| self.section_body(s, ui));
                });
            if !open {
                self.close_section(s);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use egui::{Pos2, Rect, Vec2, pos2, vec2};

    fn data() -> Option<std::path::PathBuf> {
        let d = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio");
        d.join("HG002.revio.slices.bam").exists().then_some(d)
    }

    /// A scratch folder for one test, with GGV's settings kept in it.
    fn scratch(name: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("ggv-inspector-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(dir.join("config")).unwrap();
        crate::store::TEST_CONFIG.with_borrow_mut(|c| *c = Some(dir.join("config")));
        dir
    }

    fn frame(ctx: &egui::Context, size: Vec2, mut ui: impl FnMut(&mut egui::Ui)) -> egui::FullOutput {
        let raw = egui::RawInput { screen_rect: Some(Rect::from_min_size(Pos2::ZERO, size)), ..Default::default() };
        let mut out = ctx.run_ui(raw, |root| ui(root));
        out.textures_delta.clear();
        out
    }

    /// A window opens beside what it must not hide, inside its pane and GGV's window.
    #[test]
    fn windows_open_beside_the_base() {
        let screen = Rect::from_min_size(Pos2::ZERO, vec2(1400.0, 900.0));
        let pane = Rect::from_min_max(pos2(0.0, 60.0), pos2(1050.0, 870.0));
        let size = vec2(400.0, 300.0);
        let rect = |(pos, pivot): (Pos2, egui::Align2)| pivot.anchor_size(pos, size);
        // The base at the centre or left of the pane: the window on the right, clear of it.
        for base in [525.0, 100.0] {
            let r = rect(beside(pane, base, size.y, screen));
            assert!(r.left() > base && r.right() <= pane.right() && r.top() >= pane.top(), "{base}: {r:?}");
        }
        // Right of the centre: on the left.
        let r = rect(beside(pane, 900.0, size.y, screen));
        assert!(r.left() >= pane.left() && r.right() < 900.0, "{r:?}");
        // The last pane of a four-way split, narrower than the window: at its right edge,
        // within GGV's window.
        let pane4 = Rect::from_min_max(pos2(1050.0, 60.0), pos2(1400.0, 870.0));
        assert_eq!(rect(beside(pane4, 1225.0, size.y, screen)).right(), 1392.0);
        // Before anything is drawn (no pane yet): within GGV's window.
        let r = rect(beside(Rect::NOTHING, 0.0, size.y, screen));
        assert!(screen.contains_rect(r), "{r:?}");
    }

    /// Popped-out Evidence, Find sequence and Compare open over their own pane (the active one
    /// in a split), clear of the outlined base and of the view centre where a jump puts a hit.
    #[test]
    fn floating_windows_open_in_their_pane_clear_of_the_base() {
        let Some(data) = data() else { return };
        let dir = scratch("float");
        let ctx = egui::Context::default();
        let mut bams = Vec::new();
        for name in ["a.bam", "b.bam"] {
            std::fs::copy(data.join("HG002.revio.slices.bam"), dir.join(name)).unwrap();
            std::fs::copy(data.join("HG002.revio.slices.bam.bai"), dir.join(format!("{name}.bai"))).unwrap();
            bams.push(dir.join(name));
        }
        let mut app = App::new(&ctx, crate::app::Args { bams, locus: Some("chr20:58888838-58889038".into()), ..Default::default() }).unwrap();
        assert!(!app.contigs.is_empty(), "contigs from the BAM header");
        // This pane's tracks are the right three quarters of GGV's window (the other pane, or a
        // side column, on the left), the base outlined at their centre.
        let screen = vec2(1400.0, 900.0);
        let pane = Rect::from_min_max(pos2(350.0, 80.0), pos2(1400.0, 870.0));
        app.canvas_rect = pane;
        app.width = pane.width();
        app.highlight = Some((app.contig, app.center as u64));
        let base = app.base_x(&ctx).unwrap();
        assert!((base - pane.center().x).abs() < 20.0, "{base}");
        app.show_evidence = true;
        app.inspector.places.insert(Section::Evidence, Place { floating: true, ..Default::default() });
        app.seq = Some(super::super::tools::SeqUi::default());
        app.open_compare();
        assert!(app.compare.is_some());
        // Windows are drawn before the panes: on GGV's first frame the pane isn't known yet, and
        // the windows move to it once it is (they first opened over the middle of GGV's window).
        app.canvas_rect = Rect::NOTHING;
        frame(&ctx, screen, |ui| {
            let ctx = ui.ctx().clone();
            app.floating_sections(&ctx);
            app.seq_window(&ctx);
            app.compare_window(&ctx);
        });
        app.canvas_rect = pane;
        for _ in 0..3 {
            frame(&ctx, screen, |ui| {
                let ctx = ui.ctx().clone();
                app.floating_sections(&ctx);
                app.seq_window(&ctx);
                app.compare_window(&ctx);
            });
        }
        for id in [
            egui::Id::unique((Section::Evidence.title(), app.reference_slot)),
            egui::Id::unique(("seq-search", app.reference_slot)),
            egui::Id::unique(("compare", app.reference_slot)),
        ] {
            let r = ctx.memory(|m| m.area_rect(id)).expect("shown");
            assert!(r.left() >= pane.left() && r.right() <= pane.right() + 1.0 && r.top() >= pane.top(), "{id:?} in its pane: {r:?}");
            assert!(!(r.left()..=r.right()).contains(&base), "{id:?} clear of the base at {base}: {r:?}");
        }
        // Then they stay put when GGV's window is resized a little (a right-pivoted window
        // jumped left by its own width).
        let before: Vec<Rect> = ["seq-search", "compare"].iter().map(|n| ctx.memory(|m| m.area_rect(egui::Id::unique((*n, app.reference_slot)))).unwrap()).collect();
        for _ in 0..3 {
            frame(&ctx, screen - vec2(16.0, 39.0), |ui| {
                let ctx = ui.ctx().clone();
                app.seq_window(&ctx);
                app.compare_window(&ctx);
            });
        }
        for (n, b) in ["seq-search", "compare"].iter().zip(before) {
            let r = ctx.memory(|m| m.area_rect(egui::Id::unique((*n, app.reference_slot)))).unwrap();
            assert!((r.left() - b.left()).abs() <= 16.0, "{n} stayed: {b:?} -> {r:?}");
        }
        crate::store::TEST_CONFIG.with_borrow_mut(|c| *c = None);
        let _ = std::fs::remove_dir_all(dir);
    }

    /// A review started with Evidence above it and the column too short for both: its buttons
    /// are scrolled into view (they were below the bottom with no scroll bar), and so is a
    /// failed save's message.
    #[test]
    fn a_review_below_the_fold_scrolls_into_view() {
        let Some(data) = data() else { return };
        let dir = scratch("review");
        let vcf = dir.join("calls.vcf.gz");
        std::fs::copy(data.join("HG002.revio.deepvariant.phased.vcf.gz"), &vcf).unwrap();
        std::fs::copy(data.join("HG002.revio.deepvariant.phased.vcf.gz.tbi"), dir.join("calls.vcf.gz.tbi")).unwrap();
        let bam = data.join("HG002.revio.slices.bam");
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let mut app = App::new(&ctx, crate::app::Args { bams: vec![bam.clone(), bam.clone(), bam], vcfs: vec![vcf], locus: Some("chr20:58888938".into()), ..Default::default() }).unwrap();
        // The old order, Evidence first, as a saved layout may still have it.
        app.inspector.order = vec![Section::Evidence, Section::Region, Section::Review, Section::Queue];
        app.show_evidence = true;
        let size = vec2(330.0, 260.0);
        let show = |app: &mut App| frame(&ctx, size, |ui| app.column(ui)).platform_output.accesskit_update.expect("accesskit");
        let bounds = |tree: &egui::accesskit::TreeUpdate, prefix: &str| {
            tree.nodes.iter().find(|(_, n)| n.label().or(n.value()).is_some_and(|l| l.starts_with(prefix))).and_then(|(_, n)| n.bounds())
        };
        for _ in 0..3 {
            show(&mut app);
        }
        app.start_review(0);
        let mut tree = show(&mut app);
        for _ in 0..30 {
            tree = show(&mut app);
        }
        let next = bounds(&tree, "Next (N)").expect("the review's Next button");
        assert!(next.y0 >= 0.0 && next.y1 <= f64::from(size.y), "Next in view: {next:?}");
        // A failed save: its message is shown too.
        app.review.as_mut().unwrap().error = Some("disk full".into());
        for _ in 0..30 {
            tree = show(&mut app);
        }
        let err = bounds(&tree, "Not saved").expect("the error");
        assert!(err.y0 >= 0.0 && err.y1 <= f64::from(size.y), "error in view: {err:?}");
        crate::store::TEST_CONFIG.with_borrow_mut(|c| *c = None);
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn review_and_queue_come_first_by_default() {
        assert_eq!(Layout::default().order[..2], [Section::Review, Section::Queue]);
        for s in Section::ALL {
            assert!(Layout::default().order.contains(&s));
        }
    }

    #[test]
    fn reorder_and_round_trip() {
        let mut l = Layout::default();
        l.order.retain(|s| *s != Section::Queue);
        l.order.insert(0, Section::Queue);
        assert_eq!(l.order[0], Section::Queue);
        l.places.insert(Section::Region, Place { floating: true, pinned: true, collapsed: false });
        let back: Layout = serde_json::from_value(serde_json::to_value(&l).unwrap()).unwrap();
        assert_eq!(back, l);
        // Older prefs without some sections, or none at all, still load whole.
        let partial: Layout = serde_json::from_str(r#"{"order":["Review"]}"#).unwrap();
        assert_eq!(partial.order, vec![Section::Review]);
    }
}
