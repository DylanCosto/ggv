//! Everything the app can do, by name: one list feeds the command palette, the keyboard
//! shortcuts overlay and the key handling, so they can't disagree.

use eframe::egui::{self, Align2, Key, Modifiers, Vec2};
use ggv_layout::GroupBy;

use super::App;
use crate::{
    session::{color_name, group_name},
    tracks::{ColorBy, CoverageScale},
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Action {
    Palette,
    Shortcuts,
    GoToLocus,
    Back,
    Forward,
    ZoomIn,
    ZoomOut,
    AddBam,
    AddVcf,
    AddGenes,
    SetReference,
    SaveSession,
    OpenSession,
    RestoreSession,
    ExportPng,
    CopyImage,
    CopyLink,
    CopyCommand,
    CopyFasta,
    CopySam,
    Report,
    AddBookmark,
    Bookmarks,
    ExportBookmarks,
    QueueOpen,
    QueuePaste,
    QueueStop,
    ReviewVariants,
    NextVariant,
    PrevVariant,
    Filter,
    ClearFilter,
    FindRead,
    SeqSearch,
    Compare,
    NextDiff,
    PrevDiff,
    ToggleEvidence,
    ToggleSquish,
    ToggleSoftClips,
    ToggleDark,
    ToggleAllTranscripts,
    ToggleSharedScale,
    ToggleAutoReload,
    ReloadAll,
    Color(ColorBy),
    Group(GroupBy),
    Scale(CoverageScale),
    Split,
    ClosePane,
    LinkPanes,
    ClearSelection,
    OpenDemo,
    About,
    Assembly,
}

/// A shortcut: modifiers and key.
#[derive(Clone, Copy)]
pub struct Keys(pub Modifiers, pub Key);

const CTRL: Modifiers = Modifiers::COMMAND;
const CTRL_SHIFT: Modifiers = Modifiers { shift: true, ..Modifiers::COMMAND };
const ALT: Modifiers = Modifiers::ALT;
const NONE: Modifiers = Modifiers::NONE;

impl Action {
    pub fn all() -> Vec<Self> {
        use Action::*;
        let mut v = vec![
            Palette, Shortcuts, GoToLocus, Back, Forward, ZoomIn, ZoomOut, AddBam, AddVcf, AddGenes, SetReference, SaveSession,
            OpenSession, RestoreSession, ExportPng, CopyImage, CopyLink, CopyCommand, CopyFasta, CopySam, Report, AddBookmark,
            Bookmarks, ExportBookmarks, QueueOpen, QueuePaste, QueueStop, ReviewVariants, NextVariant, PrevVariant, Filter,
            ClearFilter, FindRead, SeqSearch, Compare, NextDiff, PrevDiff, ToggleEvidence, ToggleSquish, ToggleSoftClips,
            ToggleDark, ToggleAllTranscripts, ToggleSharedScale, ToggleAutoReload, ReloadAll,
        ];
        v.extend([ColorBy::None, ColorBy::Strand, ColorBy::Haplotype, ColorBy::InsertSize, ColorBy::BaseMods].map(Color));
        v.extend([GroupBy::None, GroupBy::Strand, GroupBy::Haplotype, GroupBy::ReadGroup, GroupBy::Alleles].map(Group));
        v.extend([CoverageScale::View, CoverageScale::Loaded, CoverageScale::Log].map(Scale));
        v.extend([Split, ClosePane, LinkPanes, ClearSelection, OpenDemo, About, Assembly]);
        v
    }

    /// The name in a menu: shorter than the palette's where the menu says the rest.
    pub fn menu_name(self) -> String {
        match self {
            Action::AddBookmark => "Bookmark this place".into(),
            a => a.name(),
        }
    }

    pub fn name(self) -> String {
        use Action::*;
        match self {
            Palette => "Command palette".into(),
            Shortcuts => "Keyboard shortcuts".into(),
            GoToLocus => "Go to a locus or gene (locus box)".into(),
            Back => "Back".into(),
            Forward => "Forward".into(),
            ZoomIn => "Zoom in".into(),
            ZoomOut => "Zoom out".into(),
            AddBam => "Open alignments (BAM/CRAM)…".into(),
            AddVcf => "Open variants (VCF/BCF)…".into(),
            AddGenes => "Open gene annotation (GTF/GFF/BED)…".into(),
            SetReference => "Open reference FASTA…".into(),
            SaveSession => "Save session as…".into(),
            OpenSession => "Open session…".into(),
            RestoreSession => "Restore last session".into(),
            ExportPng => "Export view as PNG…".into(),
            CopyImage => "Copy view image to clipboard".into(),
            CopyLink => "Copy link to this view (ggv://)".into(),
            CopyCommand => "Copy command that reopens this view".into(),
            CopyFasta => "Copy region as FASTA (selection or view)".into(),
            CopySam => "Copy reads in selection as SAM".into(),
            Report => "Evidence report (HTML)…".into(),
            AddBookmark => "Bookmark the selected region (or the outlined base, or the view)".into(),
            Bookmarks => "Bookmarks".into(),
            ExportBookmarks => "Export bookmarks as BED…".into(),
            QueueOpen => "Review queue: open a list (BED, VCF, loci)…".into(),
            QueuePaste => "Review queue: paste loci…".into(),
            QueueStop => "Review queue: stop".into(),
            ReviewVariants => "Review variants (triage the first VCF)".into(),
            NextVariant => "Next variant".into(),
            PrevVariant => "Previous variant".into(),
            Filter => "Filter reads with an expression…".into(),
            ClearFilter => "Clear filter expressions".into(),
            FindRead => "Find read by name…".into(),
            SeqSearch => "Find sequence in the reference (primer or oligo)…".into(),
            Compare => "Compare two tracks (tumor/normal)…".into(),
            NextDiff => "Next difference (compare mode)".into(),
            PrevDiff => "Previous difference (compare mode)".into(),
            ToggleEvidence => "Toggle evidence panel".into(),
            ToggleSquish => "Toggle squished reads".into(),
            ToggleSoftClips => "Toggle soft-clipped bases".into(),
            ToggleDark => "Toggle dark theme".into(),
            ToggleAllTranscripts => "Toggle all transcripts".into(),
            ToggleSharedScale => "Toggle one coverage scale for all tracks".into(),
            ToggleAutoReload => "Toggle reloading files when they change".into(),
            ReloadAll => "Reload all files".into(),
            Color(c) => format!("Colour reads by {}", color_name(c)),
            Group(g) => format!("Group reads by {}", group_name(g)),
            Scale(CoverageScale::View) => "Coverage: autoscale to view".into(),
            Scale(CoverageScale::Loaded) => "Coverage: scale to loaded data".into(),
            Scale(CoverageScale::Log) => "Coverage: log scale".into(),
            Split => "Split view (add a pane)".into(),
            ClosePane => "Close this pane".into(),
            LinkPanes => "Link panes (scroll and zoom together)".into(),
            ClearSelection => "Clear selection, outline and search marks".into(),
            About => "About GGV and copy diagnostics (for bug reports)".into(),
            Assembly => "Assembly and reference checks (which build each file is)".into(),
            OpenDemo => "Open the demo (HG002 long reads, 4.3 MB download) with a tour".into(),
        }
    }

    pub fn keys(self) -> Option<Keys> {
        use Action::*;
        Some(match self {
            Palette => Keys(CTRL, Key::K),
            Shortcuts => Keys(NONE, Key::Questionmark),
            GoToLocus => Keys(CTRL, Key::L),
            Back => Keys(ALT, Key::ArrowLeft),
            Forward => Keys(ALT, Key::ArrowRight),
            CopyImage => Keys(CTRL_SHIFT, Key::C),
            Report => Keys(CTRL_SHIFT, Key::R),
            AddBookmark => Keys(CTRL, Key::D),
            Bookmarks => Keys(CTRL_SHIFT, Key::B),
            NextVariant => Keys(CTRL, Key::F),
            PrevVariant => Keys(CTRL, Key::B),
            Filter => Keys(CTRL, Key::E),
            FindRead => Keys(CTRL_SHIFT, Key::F),
            SeqSearch => Keys(CTRL_SHIFT, Key::S),
            NextDiff => Keys(NONE, Key::CloseBracket),
            PrevDiff => Keys(NONE, Key::OpenBracket),
            ToggleEvidence => Keys(NONE, Key::E),
            ToggleSquish => Keys(NONE, Key::S),
            SaveSession => Keys(CTRL, Key::S),
            OpenSession => Keys(CTRL, Key::O),
            ClearSelection => Keys(NONE, Key::Escape),
            _ => return None,
        })
    }
}

pub fn keys_label(k: Keys) -> String {
    let Keys(m, key) = k;
    let mut s = String::new();
    if m.command || m.ctrl {
        s.push_str(if cfg!(target_os = "macos") { "Cmd+" } else { "Ctrl+" });
    }
    if m.alt {
        s.push_str("Alt+");
    }
    if m.shift {
        s.push_str("Shift+");
    }
    s.push_str(match key {
        Key::ArrowLeft => "Left",
        Key::ArrowRight => "Right",
        Key::Questionmark => "?",
        Key::OpenBracket => "[",
        Key::CloseBracket => "]",
        Key::Escape => "Esc",
        k => k.name(),
    });
    s
}

/// Keys handled elsewhere (not actions), for the overlay.
const OTHER_KEYS: &[(&str, &str)] = &[
    ("Left / Right", "Pan"),
    ("+ / -, mouse wheel, Ctrl+wheel", "Zoom (around the pointer); over rows that don't all fit, only Ctrl+wheel"),
    ("Up / Down, PgUp / PgDn", "Scroll the rows of the track under the pointer"),
    ("Wheel over rows, their scrollbar", "Scroll reads, gene or sample rows that don't all fit"),
    ("Drag", "Pan (and scroll rows)"),
    ("Shift+drag", "Select a region: statistics"),
    ("Ctrl+drag", "Measure: length, GC, depth, reads (from what's loaded)"),
    ("Double-click", "Zoom in, or fit a gene"),
    ("Right-click", "Sort, group, colour, copy, go to mate"),
    ("N / P", "Next / previous item (variant review or review queue)"),
    ("1 / 2 / 3", "Mark real/pass, artifact/fail, unsure, and move on"),
    ("U", "Next undecided item (review queue)"),
];

impl App {
    /// Runs the shortcut of any action whose keys were pressed. `typing`: a text field has
    /// focus, so only shortcuts with Ctrl (or Alt) count.
    pub(super) fn action_keys(&mut self, ctx: &egui::Context, typing: bool) {
        // The palette handles its own keys; Escape closes the shortcuts overlay first.
        if self.palette.is_some() {
            return;
        }
        if self.show_shortcuts && ctx.input_mut(|i| i.consume_key(NONE, Key::Escape)) {
            self.show_shortcuts = false;
            return;
        }
        // Escape closes the open tool windows first (they piled up otherwise); with none open
        // it clears the selection and marks as before.
        let windows = self.variant_popup.is_some()
            || self.show_bookmarks
            || self.find.is_some()
            || self.seq.is_some()
            || self.filter_ui.is_some()
            || self.show_about
            || self.assembly.show;
        if windows && ctx.input_mut(|i| i.consume_key(NONE, Key::Escape)) {
            self.variant_popup = None;
            self.show_bookmarks = false;
            self.find = None;
            if let Some(s) = self.seq.take()
                && let Some(x) = &s.search
            {
                x.cancel();
            }
            self.filter_ui = None;
            self.show_about = false;
            self.assembly.show = false;
            return;
        }
        // egui matches modifiers loosely (Ctrl+F also fires on Ctrl+Shift+F), so try the
        // shortcuts with more modifiers first.
        let mut shortcuts: Vec<(Action, Keys)> = Action::all().into_iter().filter_map(|a| a.keys().map(|k| (a, k))).collect();
        shortcuts.sort_by_key(|(_, Keys(m, _))| std::cmp::Reverse(u8::from(m.command) + u8::from(m.shift) + u8::from(m.alt)));
        let hit = ctx.input_mut(|i| {
            shortcuts
                .into_iter()
                .filter(|(_, Keys(m, _))| !typing || m.command || m.alt)
                .find(|(_, Keys(m, k))| i.consume_key(*m, *k))
                .map(|(a, _)| a)
        });
        if let Some(a) = hit {
            self.run(a, ctx);
        }
    }

    pub(super) fn run(&mut self, a: Action, ctx: &egui::Context) {
        use Action::*;
        let mut style = self.style;
        match a {
            Palette => self.palette = Some(super::palette::Palette::default()),
            Shortcuts => self.show_shortcuts = !self.show_shortcuts,
            GoToLocus => self.focus_locus = true,
            Back => self.go_back(),
            Forward => self.go_forward(),
            ZoomIn => self.zoom_at(2.0, self.width / 2.0),
            ZoomOut => self.zoom_at(0.5, self.width / 2.0),
            AddBam => {
                self.ask(rfd::FileDialog::new().add_filter("Alignments", &["bam", "cram"]), super::dialogs::Kind::PickFiles, super::dialogs::Then::AddTracks);
            }
            AddVcf => {
                self.ask(rfd::FileDialog::new().add_filter("Variants", &["gz", "bcf", "vcf"]), super::dialogs::Kind::PickFiles, super::dialogs::Then::AddVariants);
            }
            AddGenes => {
                self.ask(rfd::FileDialog::new().add_filter("Annotation", &["gtf", "gff", "gff3", "bed", "gz"]), super::dialogs::Kind::PickFiles, super::dialogs::Then::AddGenes);
            }
            SetReference => {
                self.ask(rfd::FileDialog::new().add_filter("Reference", &["fa", "fasta", "fna", "gz", "bgz", "2bit"]), super::dialogs::Kind::PickFile, super::dialogs::Then::Reference);
            }
            SaveSession => self.save_session_as(),
            OpenSession => {
                self.ask(rfd::FileDialog::new().add_filter("GGV session", &["json"]), super::dialogs::Kind::PickFile, super::dialogs::Then::OpenSession);
            }
            RestoreSession => {
                if let Some(last) = crate::session::last_session_path().filter(|p| p.exists()) {
                    let _ = self.open_session(&last, ctx);
                }
            }
            ExportPng => {
                self.ask(rfd::FileDialog::new().add_filter("PNG", &["png"]).set_file_name("ggv-view.png"), super::dialogs::Kind::Save, super::dialogs::Then::ExportPng);
            }
            CopyImage => self.capture(ctx, super::tools::Capture::Clipboard),
            CopyLink => {
                ctx.copy_text(crate::links::to_link(&self.session()));
                self.notify("Copied a ggv:// link to this view");
            }
            CopyCommand => {
                ctx.copy_text(crate::links::to_command(&self.session()));
                self.notify("Copied the command that reopens this view");
            }
            CopyFasta => self.copy_fasta(ctx),
            CopySam => self.copy_selection_sam(ctx),
            Report => self.start_report(),
            AddBookmark => self.add_bookmark(),
            Bookmarks => self.show_bookmarks = !self.show_bookmarks,
            ExportBookmarks => self.export_bookmarks(),
            QueueOpen => self.open_queue_dialog(),
            QueuePaste => {
                self.paste_queue = Some(String::new());
                ctx.memory_mut(|m| m.data.insert_temp(egui::Id::unique("paste-queue-focus"), true));
            }
            QueueStop => self.queue = None,
            ReviewVariants => {
                if !self.variants.is_empty() {
                    self.queue = None;
                    self.start_review(0);
                }
            }
            NextVariant => self.seek_variant(true, self.hover_variant),
            PrevVariant => self.seek_variant(false, self.hover_variant),
            Filter => self.open_filter_window(None),
            ClearFilter => {
                style.filter.expr = None;
                // Forgotten too, as clearing one file's filter does: else a split pane or the
                // next start applies it again.
                let open: Vec<std::path::PathBuf> = self.tracks.iter().map(|t| t.loader.bam_path.clone()).collect();
                super::tools::forget_file_filters(std::mem::take(&mut self.track_filters).into_keys().chain(open));
            }
            FindRead => self.find = Some(super::tools::FindUi::default()),
            SeqSearch => self.seq = Some(super::tools::SeqUi::default()),
            Compare => self.open_compare(),
            NextDiff => self.step_diff(true, ctx),
            PrevDiff => self.step_diff(false, ctx),
            ToggleEvidence => self.show_evidence = !self.show_evidence,
            ToggleSquish => style.squished = !style.squished,
            ToggleSoftClips => style.soft_clips = !style.soft_clips,
            ToggleDark => {
                style.dark = !style.dark;
                ctx.set_theme(if style.dark { egui::Theme::Dark } else { egui::Theme::Light });
            }
            ToggleAllTranscripts => self.all_transcripts = !self.all_transcripts,
            ToggleSharedScale => style.shared_scale = !style.shared_scale,
            ToggleAutoReload => {
                self.auto_reload = !self.auto_reload;
                self.notify(if self.auto_reload { "Files reload when they change" } else { "Files no longer reload when they change" });
            }
            ReloadAll => self.reload_all(ctx),
            Color(c) => style.color = c,
            Group(g) => style.group = g,
            Scale(s) => style.coverage_scale = s,
            Split | ClosePane | LinkPanes => self.shell_request = Some(a),
            About => self.show_about = true,
            Assembly => self.assembly.show = !self.assembly.show,
            OpenDemo => {
                // Replaces the display settings with the demo's own: don't write `style` back.
                self.open_demo(ctx);
                return;
            }
            ClearSelection => {
                self.selection = None;
                self.measure = None;
                self.selected_read = None;
                self.highlight = None;
                if let Some(s) = &mut self.seq {
                    s.search = None;
                }
                for t in &mut self.tracks {
                    t.summary = None;
                }
            }
        }
        if style != self.style {
            self.style = style;
        }
    }

    /// The `?` overlay: every shortcut, from the same list the keys use. Centred and within
    /// GGV's window whatever its size: it scrolls, can be moved and resized, and its two tables
    /// stack when there isn't room for them side by side; long names wrap rather than overlap.
    pub(super) fn shortcuts_window(&mut self, ctx: &egui::Context) {
        if !self.show_shortcuts || !self.is_active {
            return;
        }
        let screen = ctx.content_rect();
        let (default, max) = shortcuts_size(screen.size());
        let mut open = true;
        egui::Window::new("Keyboard shortcuts")
            .id(egui::Id::unique(("shortcuts", self.reference_slot)))
            .open(&mut open)
            .pivot(Align2::CENTER_CENTER)
            .default_pos(screen.center())
            .default_size(default)
            .max_size(max)
            .collapsible(false)
            .resizable(true)
            .show(ctx, |ui| {
                ui.add(egui::Label::new(egui::RichText::new("Everything else is in the command palette (Ctrl+K).").weak()).wrap());
                ui.separator();
                super::inspector::visible_scroll_bars(ui);
                egui::ScrollArea::vertical().auto_shrink(false).show(ui, |ui| {
                    let actions: Vec<(String, String)> = Action::all().into_iter().filter_map(|a| a.keys().map(|k| (keys_label(k), a.name()))).collect();
                    let other: Vec<(String, String)> = OTHER_KEYS.iter().map(|(k, what)| ((*k).to_owned(), (*what).to_owned())).collect();
                    let table = |ui: &mut egui::Ui, id: &str, rows: &[(String, String)]| {
                        // A finite column width makes the names wrap at the table's right edge.
                        egui::Grid::new(id).striped(true).num_columns(2).max_col_width(ui.available_width().max(120.0)).show(ui, |ui| {
                            for (k, what) in rows {
                                ui.add(egui::Label::new(egui::RichText::new(k).monospace()).wrap_mode(egui::TextWrapMode::Extend));
                                ui.add(egui::Label::new(what.as_str()).wrap());
                                ui.end_row();
                            }
                        });
                    };
                    if ui.available_width() >= SHORTCUTS_TWO_COLUMNS {
                        ui.columns(2, |cols| {
                            table(&mut cols[0], "shortcut-actions", &actions);
                            table(&mut cols[1], "shortcut-other", &other);
                        });
                    } else {
                        table(ui, "shortcut-actions", &actions);
                        ui.separator();
                        table(ui, "shortcut-other", &other);
                    }
                });
            });
        if !open {
            self.show_shortcuts = false;
        }
    }
}

/// The shortcuts window's tables go side by side from this width.
const SHORTCUTS_TWO_COLUMNS: f32 = 760.0;

/// The shortcuts window's first size and its largest, for GGV's window `screen` big: as wide as
/// both tables need (about 1000 px) where that fits, and never past GGV's window.
fn shortcuts_size(screen: Vec2) -> (Vec2, Vec2) {
    let max = Vec2::new((screen.x - 24.0).max(200.0), (screen.y - 48.0).max(150.0));
    (Vec2::new(1000.0f32.min(max.x), 560.0f32.min(max.y)), max)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shortcuts_are_unique_and_named() {
        let all = Action::all();
        let mut seen = Vec::new();
        for a in &all {
            assert!(!a.name().is_empty());
            if let Some(Keys(m, k)) = a.keys() {
                let key = (m.command, m.shift, m.alt, k);
                assert!(!seen.contains(&key), "{a:?} reuses {}", keys_label(Keys(m, k)));
                seen.push(key);
            }
        }
        assert!(all.len() > 50);
        assert_eq!(keys_label(Action::Palette.keys().unwrap()), if cfg!(target_os = "macos") { "Cmd+K" } else { "Ctrl+K" });
    }

    /// The shortcuts window fits GGV's window at any size, opening centred, its key column
    /// on screen (it was anchored, 1020 px wide, off both sides of an 800 px window).
    #[test]
    fn shortcuts_window_fits_small_windows() {
        for (w, h) in [(800.0, 560.0), (1280.0, 720.0), (1400.0, 900.0), (500.0, 400.0)] {
            let (default, max) = shortcuts_size(Vec2::new(w, h));
            assert!(default.x <= w - 24.0 && default.y <= h - 48.0, "{w}x{h}: {default:?}");
            assert!(max.x <= w && max.y <= h);
            let ctx = egui::Context::default();
            let mut app = App::new(&ctx, super::super::Args::default()).unwrap();
            app.show_shortcuts = true;
            let screen = egui::Rect::from_min_size(egui::Pos2::ZERO, Vec2::new(w, h));
            for _ in 0..3 {
                let raw = egui::RawInput { screen_rect: Some(screen), ..Default::default() };
                ctx.run_pass(raw, |ctx| app.shortcuts_window(ctx)).textures_delta.clear();
            }
            let rect = ctx.memory(|m| m.area_rect(egui::Id::unique(("shortcuts", app.reference_slot)))).expect("shown");
            assert!(screen.expand(1.0).contains_rect(rect), "{w}x{h}: {rect:?}");
            assert!((rect.center().x - w / 2.0).abs() < 2.0, "{w}x{h}: centred: {rect:?}");
        }
    }
}
