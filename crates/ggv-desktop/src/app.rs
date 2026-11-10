mod about;
mod actions;
mod assembly;
pub(crate) mod arrange;
mod attach;
mod coding;
mod control;
mod demo;
mod dialogs;
mod genome;
mod inspector;
mod palette;
mod tools;
mod variant_info;
mod workspace;

pub use actions::Action;
pub use control::{API_VERSION, CAPABILITIES};

use std::{collections::HashMap, path::PathBuf, sync::Arc};

use eframe::egui::{self, Align2, Color32, FontId, Pos2, Rect, Sense, Stroke, Vec2};
use ggv_layout::{AlignedRead, Coverage, ViewData, ViewOptions};
use ggv_model::{Contig, Region};

use crate::{
    genes::{self, GeneTrack},
    session::{self, Session, color_name, group_name, sort_name},
    gpu::{self, InstanceSet, Uniforms},
    track::{Span, Track},
    variants::{self, VariantTrack},
    tracks::{self, COVERAGE_H, COVERAGE_TOP, HEADER_H, LETTERS_PPB, READS_TOP, REF_H, SEPARATOR_Y, Style},
    util::{commas, nice_step},
};

/// Reads load for views up to this wide (IGV's visibility window).
/// Reads are loaded for views up to this many bases (the "visibility window"); wider views
/// show coverage only. 1 Mb by default: on 35x short reads (265k reads) panning stays above
/// 60 fps with about 1 GB per track; 2 Mb starts to stall and 5 Mb takes 3.5 GB. People choose
/// their own (toolbar, `--read-window`, MCP `display`), kept in `prefs.json`.
pub const DEFAULT_READ_WINDOW: u64 = 1_000_000;
pub const MIN_READ_WINDOW: u64 = 10_000;
pub const MAX_READ_WINDOW: u64 = 50_000_000;
static READ_WINDOW_BP: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(DEFAULT_READ_WINDOW);

pub fn read_window() -> u64 {
    READ_WINDOW_BP.load(std::sync::atomic::Ordering::Relaxed)
}

pub fn set_read_window(bp: u64) {
    READ_WINDOW_BP.store(bp.clamp(MIN_READ_WINDOW, MAX_READ_WINDOW), std::sync::atomic::Ordering::Relaxed);
}

/// The saved read window, if any (`prefs.json`).
pub fn load_read_window_pref() {
    if let Some(bp) = tools::load_prefs().get("read_window").and_then(serde_json::Value::as_u64) {
        set_read_window(bp);
    }
}

/// Sets and remembers the read window.
pub fn save_read_window(bp: u64) {
    set_read_window(bp);
    let mut prefs = tools::load_prefs();
    prefs.insert("read_window".into(), serde_json::json!(read_window()));
    tools::save_prefs(&prefs);
}

/// `30 kb`, `1 Mb`, `2.5 Mb`.
/// An SV's size, exact below 1 kb (`292 bp`), then to a tenth (`1.4 kb`, `3.0 kb`, `12.5 kb`).
fn sv_size(bp: u64) -> String {
    match bp {
        0..1000 => format!("{bp} bp"),
        1000..1_000_000 => format!("{:.1} kb", bp as f64 / 1000.0),
        _ => format!("{:.1} Mb", bp as f64 / 1e6),
    }
}

pub fn fmt_bp(bp: u64) -> String {
    if bp >= 1_000_000 {
        let mb = bp as f64 / 1e6;
        if (mb - mb.round()).abs() < 1e-9 { format!("{mb:.0} Mb") } else { format!("{mb:.1} Mb") }
    } else {
        format!("{} kb", bp / 1000)
    }
}
/// Binned depth loads for views up to this wide.
const COVERAGE_WINDOW: u64 = 5_000_000;
/// Most zoomed-in: points per base.
const MAX_PPB: f64 = 40.0;
const MAX_BINS: usize = 16_384;
/// Load at least this much when reads are shown, so panning at high zoom stays local.
const MIN_READ_SPAN: f64 = 20_000.0;
/// Bases a reads span may exceed the read window by: integer ends around fractional view edges.
const SPAN_SLACK: u64 = 2;

/// Shared rows above the tracks, in points from the canvas top.
const RULER_H: f32 = 30.0;
const REF_TOP: f32 = 32.0;
const TRACKS_TOP: f32 = REF_TOP + REF_H + 4.0;
/// A drag on the ruler or reference row marking less than this many points is taken for a slip.
const MIN_ZOOM_DRAG: f64 = 4.0;
/// Tracks split the remaining height evenly, but never get shorter than this.
const MIN_SECTION_H: f32 = 120.0;

#[derive(Default)]
pub struct Args {
    /// A session to open instead of the files below.
    pub session: Option<Session>,
    /// Open the built-in demo, with its tour.
    pub demo: bool,
    pub bams: Vec<PathBuf>,
    pub vcfs: Vec<PathBuf>,
    pub genes: Vec<PathBuf>,
    /// Signal tracks (bigWig, bedGraph, SEG).
    pub signals: Vec<PathBuf>,
    pub all_transcripts: bool,
    pub show_evidence: bool,
    pub review: bool,
    pub select: Option<String>,
    pub select_read: Option<String>,
    pub style: Style,
    /// Sort at the locus once it is resolved (it must be a single base).
    pub sort: Option<ggv_layout::SortBy>,
    pub reference: Option<PathBuf>,
    pub locus: Option<String>,
    /// A base to outline once loaded (`--highlight`).
    pub highlight: Option<String>,
    /// A built-in genome to open (`--genome`).
    pub genome: Option<String>,
    /// Built-in tracks to open once the genome's build is known (`--track`).
    pub tracks: Vec<String>,
    pub automation: Option<crate::automation::Automation>,
    /// The window's bookmarks (split panes share them).
    pub bookmarks: Option<crate::bookmarks::Shared>,
    /// The window's overlays (split panes share them).
    pub overlays: Option<crate::overlays::Shared>,
    /// Start a review queue over this list (`--queue`).
    pub queue: Option<PathBuf>,
    /// A filter expression for every track (`--filter`).
    pub filter: Option<String>,
    /// Open the sequence search with this pattern, searching the genome (`--find-seq`).
    pub find_seq: Option<String>,
    /// Find this read by name (`--find-read`).
    pub find_read: Option<String>,
    /// Compare the first two tracks (`--compare`).
    pub compare: bool,
    /// Open the command palette with this text (`--palette`, for screenshots).
    pub palette: Option<String>,
    /// Show the shortcuts overlay (`--shortcuts`, for screenshots).
    pub shortcuts: bool,
    /// Show About GGV (`--about`, for screenshots).
    pub about: bool,
    /// Write an evidence report here once loaded (`--report`).
    pub report: Option<PathBuf>,
    /// Measure this region (`--measure`, for screenshots).
    pub measure: Option<String>,
}

pub struct App {
    tracks: Vec<Track>,
    variants: Vec<VariantTrack>,
    genes: Vec<GeneTrack>,
    signals: Vec<crate::signals::SignalTrack>,
    /// Show every transcript, not just one per gene.
    all_transcripts: bool,
    /// No BAM yet: the genome (contigs and lengths) comes from an annotation.
    contigs_from_annotation: bool,
    reference_path: Option<PathBuf>,
    /// The view's genome: the first track's contigs. Other tracks map names via aliases.
    contigs: Vec<Contig>,
    contig: usize,
    /// Centre of the view, in bases (0-based, fractional).
    center: f64,
    /// Points per base.
    ppb: f64,
    width: f32,
    style: Style,
    next_request_id: u64,
    /// This pane's GPU slot for the reference row.
    reference_slot: u64,
    /// Whether keyboard input goes to this pane (the last one clicked, in split view).
    pub is_active: bool,
    generation: u64,
    /// GPU slots of closed tracks, freed on the next paint.
    retired: Vec<u64>,
    /// The reference row's data and instances, and which (track, load) they came from.
    reference: Option<(Arc<ViewData>, InstanceSet, (u64, u64))>,
    /// With no alignment track to bring the reference row, the window reads it itself.
    own_reference: Option<(PathBuf, ggv_io::ReferenceFasta)>,
    locus_text: String,
    pending_goto: Option<String>,
    /// The tag typed in the colour/group menus.
    tag_input: String,
    /// The order and heights of the tracks, as the person arranged them.
    arrangement: arrange::Arrangement,
    arranging: Option<arrange::Arranging>,
    /// The first track shown when they don't all fit (the track stack scrolls by whole tracks).
    stack_first: usize,
    /// The track whose rows' scrollbar is being dragged.
    row_bar_drag: Option<Hit>,
    /// Built-in genomes and tracks being downloaded, and a genome to open once the window is up
    /// (`--genome`).
    downloads: Vec<genome::Download>,
    pending_genome: Option<String>,
    pending_tracks: Vec<String>,
    error: Option<String>,
    /// Files whose warning (a BAM that looks truncated; see `warning`) the person dismissed.
    dismissed_warnings: Vec<PathBuf>,
    /// The view (contig index, start, end) loads were last requested for.
    requested_view: Option<(usize, f64, f64)>,
    /// The error the last failed `goto` set, so a later successful one clears only that.
    goto_error: Option<String>,
    automation: Option<crate::automation::Automation>,
    /// A base to outline (contig index, 0-based position), set by going to a single base.
    highlight: Option<(usize, u64)>,
    /// The track whose rows a vertical drag scrolls.
    drag_track: Option<Hit>,
    /// Where the context menu was opened: contig index, base, and the read under it.
    menu_at: Option<(usize, u64, Option<AlignedRead>, Option<Hit>)>,
    /// The clicked read's name: every alignment with it is outlined (mate, split parts).
    selected_read: Option<String>,
    /// The variant whose details are open (clicked).
    variant_popup: Option<variant_info::VariantPopup>,
    /// A shift-dragged region: contig index, 0-based start and end; and whether still dragging.
    selection: Option<(usize, u64, u64)>,
    selecting: bool,
    /// A drag on the ruler or reference row in progress: the base it started on and the base
    /// under the pointer. Zooms to what it marks on release.
    zoom_drag: Option<(u64, u64)>,
    /// A sort requested on the command line, applied at the first outlined base.
    pending_sort: Option<ggv_layout::SortBy>,
    /// An outlined base from a session, applied once the genome is known.
    pending_highlight: Option<String>,
    /// Where a scripted zoom cycle started (points per base, centre).
    cycle_start: Option<(f64, f64)>,
    /// Places jumped away from (contig, centre, points per base), and ones gone back from.
    back: Vec<(usize, f64, f64)>,
    forward: Vec<(usize, f64, f64)>,
    recent: Vec<PathBuf>,
    /// Files that need an index: offered, being built, or failed.
    indexing: Vec<IndexJob>,
    /// A screenshot being taken (PNG export, clipboard, report).
    capture: Option<tools::CaptureJob>,
    /// A report waiting for the evidence at its base to load.
    pending_report: Option<PathBuf>,
    last_report: Option<PathBuf>,
    canvas_rect: Rect,
    bookmarks: crate::bookmarks::Shared,
    overlays: crate::overlays::Shared,
    /// What the person did that an attached agent may be waiting for (triage decisions),
    /// taken by the window each frame.
    pub person_events: Vec<serde_json::Value>,
    show_bookmarks: bool,
    /// A bookmark just added, whose note field should take focus.
    edit_bookmark: Option<usize>,
    queue: Option<crate::queue::Queue>,
    /// The paste-a-list window's text.
    paste_queue: Option<String>,
    palette: Option<palette::Palette>,
    show_shortcuts: bool,
    focus_locus: bool,
    /// Per-file filter expressions (replacing the shared one for that track).
    track_filters: HashMap<PathBuf, crate::expr::ExprId>,
    saved_filters: Vec<String>,
    filter_ui: Option<tools::FilterUi>,
    find: Option<tools::FindUi>,
    seq: Option<tools::SeqUi>,
    compare: Option<tools::CompareUi>,
    /// A ctrl-dragged span: contig index, start, end.
    measure: Option<(usize, u64, u64)>,
    measuring: bool,
    /// The variant track under the pointer, for Ctrl+F / Ctrl+B.
    hover_variant: Option<usize>,
    jobs: (std::sync::mpsc::Sender<tools::Apply>, std::sync::mpsc::Receiver<tools::Apply>),
    notice: Option<(String, std::time::Instant)>,
    watcher: crate::watch::Watcher,
    auto_reload: bool,
    /// An action for the window rather than this pane (split, link), taken by the shell.
    pub shell_request: Option<actions::Action>,
    /// The files the window's other panes have open (set by the shell before a command), which
    /// an export from this pane must not overwrite either.
    pub other_panes_files: Vec<PathBuf>,
    /// Startup requests that wait for the genome (`--queue`, `--find-seq`…).
    startup: Option<Startup>,
    /// A FASTA handle for the odd base the loaded reference doesn't cover.
    ref_reader: std::cell::RefCell<Option<(PathBuf, ggv_io::ReferenceFasta)>>,
    /// Consequence lines for the evidence base, and what they were computed from.
    #[allow(clippy::type_complexity)]
    consequence_cache: std::cell::RefCell<Option<((usize, u64, String, usize, Option<PathBuf>), Vec<String>)>>,
    show_evidence: bool,
    review: Option<crate::review::Review>,
    /// Start reviewing once the genome and view are set (from `--review`).
    pending_review: bool,
    /// The demo tour's current step, while it's open.
    tour: Option<usize>,
    /// The side column's sections: order, docked or floating, pinned, collapsed.
    inspector: inspector::Layout,
    /// A file dialog that is open (on its own thread).
    dialog: Option<dialogs::Pending>,
    /// For waking the window from other threads.
    egui_ctx: egui::Context,
    show_about: bool,
    /// Copy diagnostics with full file paths (off: file names only).
    diagnostics_paths: bool,
    /// Recent errors, for diagnostics.
    error_log: about::ErrorLog,
    /// Open the built-in demo on the first frame (`--demo`).
    pending_demo: bool,
    /// A region to summarize once the genome is known (from `--select`).
    pending_select: Option<String>,
    /// Which assembly each open file is, and where they disagree.
    assembly: assembly::Assembly,
    /// Files the last session named that weren't there, to locate.
    missing: Option<workspace::Missing>,
    /// A session with several panes, for the window to rebuild them.
    pending_panes: Option<Session>,
    /// Alignment files asked for while the reference waits for its index, opened once it has one.
    deferred_bams: Vec<PathBuf>,
    /// A reference asked for that has no index yet: the tracks stay on the one they have until
    /// it has (or the offer is dismissed).
    pending_reference: Option<PathBuf>,
}

/// Per-file filter expressions, each variant track's site colour, and the track arrangement.
pub type FileSettings = (HashMap<PathBuf, crate::expr::ExprId>, Vec<(PathBuf, String)>, arrange::Arrangement);

/// Where each track is drawn this frame: alignment tracks, then gene tracks at the bottom.
struct Layout {
    signals: Vec<Rect>,
    variants: Vec<Rect>,
    sections: Vec<Rect>,
    genes: Vec<Rect>,
    /// Every track shown, top to bottom, with its key.
    order: Vec<(arrange::Key, Rect)>,
    /// When the tracks don't all fit: the scrollbar along the right edge.
    stack: Option<arrange::StackBar>,
    /// The height the tracks have to share.
    avail: f32,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Hit {
    Variants(usize),
    Track(usize),
    Genes(usize),
    Signal(usize),
}

impl Layout {
    fn section_at(&self, y: f32) -> Option<usize> {
        self.sections.iter().position(|r| (r.top()..r.bottom()).contains(&y))
    }

    fn hit(&self, y: f32) -> Option<Hit> {
        let within = |rs: &[Rect]| rs.iter().position(|r| (r.top()..r.bottom()).contains(&y));
        within(&self.variants)
            .map(Hit::Variants)
            .or_else(|| self.section_at(y).map(Hit::Track))
            .or_else(|| self.genes.iter().position(|r| (r.top()..r.bottom()).contains(&y)).map(Hit::Genes))
            .or_else(|| within(&self.signals).map(Hit::Signal))
    }

    fn close_button(section: Rect) -> Rect {
        Rect::from_min_size(Pos2::new(section.right() - 18.0, section.top() + 1.0), Vec2::new(14.0, HEADER_H - 2.0))
    }

    /// Where track `hit` is drawn.
    fn rect_of(&self, hit: Hit) -> Rect {
        let get = |rs: &[Rect], i: usize| rs.get(i).copied().unwrap_or(Rect::NOTHING);
        match hit {
            Hit::Variants(i) => get(&self.variants, i),
            Hit::Track(i) => get(&self.sections, i),
            Hit::Genes(i) => get(&self.genes, i),
            Hit::Signal(i) => get(&self.signals, i),
        }
    }
}

/// A track's rows that scroll (reads, gene rows, sample rows) when they don't all fit: where
/// they are drawn, how far they are scrolled and can scroll, and how many rows of what pitch.
#[derive(Clone, Copy, Debug, PartialEq)]
struct RowScroll {
    area: Rect,
    scroll: f32,
    max: f32,
    rows: usize,
    pitch: f32,
    /// What the rows are, for the header: "rows", "samples".
    noun: &'static str,
}

impl RowScroll {
    /// The scrollbar's trough, along the right edge of the rows.
    fn bar(&self) -> Rect {
        Rect::from_min_max(Pos2::new(self.area.right() - 8.0, self.area.top() + 2.0), Pos2::new(self.area.right() - 2.0, self.area.bottom() - 2.0))
    }

    /// Where the pointer takes hold of the scrollbar: a little wider than drawn.
    fn bar_hit(&self) -> Rect {
        let bar = self.bar();
        Rect::from_min_max(Pos2::new(bar.left() - 3.0, bar.top()), Pos2::new(self.area.right(), bar.bottom()))
    }

    fn thumb(&self) -> Rect {
        let bar = self.bar();
        let h = (bar.height() * self.area.height() / (self.max + self.area.height())).clamp(bar.height().min(14.0), bar.height());
        let t = bar.top() + (bar.height() - h) * (self.scroll / self.max).clamp(0.0, 1.0);
        Rect::from_min_size(Pos2::new(bar.left(), t), Vec2::new(bar.width(), h))
    }

    /// Points the rows scroll per point the thumb moves.
    fn ratio(&self) -> f32 {
        let free = self.bar().height() - self.thumb().height();
        if free > 0.0 { self.max / free } else { 0.0 }
    }

    /// The first and last rows shown, from 1.
    fn shown(&self) -> (usize, usize) {
        let scroll = self.scroll.clamp(0.0, self.max);
        let a = ((scroll / self.pitch).floor() as usize + 1).min(self.rows);
        let b = (((scroll + self.area.height()) / self.pitch).floor() as usize).clamp(a, self.rows);
        (a, b)
    }
}

impl App {
    /// A pane. The GPU pipeline must already be set up (`gpu::init`, done once by the shell).
    pub fn new(ctx: &egui::Context, args: Args) -> anyhow::Result<Self> {
        let mut app = Self {
            tracks: Vec::new(),
            variants: Vec::new(),
            genes: Vec::new(),
            signals: Vec::new(),
            all_transcripts: args.all_transcripts,
            contigs_from_annotation: false,
            reference_path: None,
            contigs: Vec::new(),
            contig: 0,
            center: 0.0,
            ppb: 1.0,
            width: 1200.0,
            style: args.style,
            next_request_id: 1,
            reference_slot: gpu::new_slot(),
            is_active: true,
            generation: 0,
            retired: Vec::new(),
            reference: None,
            own_reference: None,
            locus_text: String::new(),
            pending_goto: args.locus,
            tag_input: "CB".to_owned(),
            arrangement: arrange::Arrangement::default(),
            arranging: None,
            stack_first: 0,
            row_bar_drag: None,
            downloads: Vec::new(),
            pending_genome: args.genome.clone(),
            pending_tracks: args.tracks.clone(),
            error: None,
            dismissed_warnings: Vec::new(),
            requested_view: None,
            goto_error: None,
            automation: args.automation,
            highlight: None,
            drag_track: None,
            menu_at: None,
            selected_read: args.select_read.clone(),
            variant_popup: None,
            selection: None,
            selecting: false,
            zoom_drag: None,
            pending_sort: args.sort,
            pending_highlight: args.highlight.clone(),
            cycle_start: None,
            back: Vec::new(),
            forward: Vec::new(),
            recent: crate::files::load_recent(),
            indexing: Vec::new(),
            capture: None,
            pending_report: None,
            last_report: None,
            canvas_rect: Rect::NOTHING,
            bookmarks: args.bookmarks.clone().unwrap_or_default(),
            overlays: args.overlays.clone().unwrap_or_default(),
            person_events: Vec::new(),
            show_bookmarks: false,
            edit_bookmark: None,
            queue: None,
            paste_queue: None,
            palette: args.palette.as_ref().map(|q| palette::Palette { query: q.clone(), ..Default::default() }),
            show_shortcuts: args.shortcuts,
            focus_locus: false,
            track_filters: HashMap::new(),
            saved_filters: tools::load_saved_filters(),
            filter_ui: None,
            find: None,
            seq: None,
            compare: None,
            measure: None,
            measuring: false,
            hover_variant: None,
            jobs: std::sync::mpsc::channel(),
            notice: None,
            watcher: crate::watch::Watcher::spawn(ctx),
            auto_reload: true,
            shell_request: None,
            other_panes_files: Vec::new(),
            ref_reader: Default::default(),
            consequence_cache: Default::default(),
            startup: Some(Startup {
                queue: args.queue.clone(),
                find_seq: args.find_seq.clone(),
                find_read: args.find_read.clone(),
                compare: args.compare,
                report: args.report.clone(),
                measure: args.measure.clone(),
            }),
            show_evidence: args.show_evidence,
            review: None,
            pending_review: false,
            tour: None,
            inspector: inspector::Layout::load(),
            dialog: None,
            egui_ctx: ctx.clone(),
            show_about: args.about,
            diagnostics_paths: false,
            error_log: about::ErrorLog::default(),
            pending_demo: args.demo,
            pending_select: None,
            missing: None,
            pending_panes: None,
            deferred_bams: Vec::new(),
            pending_reference: None,
            assembly: assembly::Assembly::default(),
        };
        if let Some(r) = args.reference {
            app.set_reference(r, ctx);
        }
        match args.session {
            Some(s) => app.apply_session(s, false, ctx),
            None => {
                app.add_tracks(args.bams, ctx);
                app.add_genes(args.genes, ctx);
                app.add_variants(args.vcfs, ctx);
                app.add_signals(args.signals, ctx);
            }
        }
        if args.review && !app.variants.is_empty() {
            app.pending_review = true;
        }
        app.pending_select = args.select;
        if let Some(f) = &args.filter {
            if let Some(e) = app.apply_filter_text(f, None) {
                app.error = Some(format!("--filter: {e}"));
            }
        }
        // Filters remembered for these files.
        let remembered = tools::load_file_filters();
        let mut bad = Vec::new();
        for t in &app.tracks {
            if !app.track_filters.contains_key(&t.loader.bam_path)
                && let Some(text) = remembered.get(&t.loader.bam_path)
            {
                match crate::expr::intern(text) {
                    Ok(e) => {
                        app.track_filters.insert(t.loader.bam_path.clone(), e);
                    }
                    Err(e) => bad.push(format!("{}: remembered filter {text:?} doesn't parse ({e}); showing it unfiltered", t.title())),
                }
            }
        }
        if !bad.is_empty() {
            app.error = Some(bad.join("\n"));
        }
        Ok(app)
    }

    /// The current files, view and display options.
    fn session(&self) -> Session {
        let mut s = Session {
            version: session::VERSION,
            bams: self.tracks.iter().map(|t| t.loader.bam_path.clone()).collect(),
            vcfs: self.variants.iter().map(|v| v.path.clone()).collect(),
            genes: self.genes.iter().map(|g| g.path.clone()).collect(),
            signals: self.signals.iter().map(|s| s.path.clone()).collect(),
            reference: self.reference_path.clone(),
            all_transcripts: self.all_transcripts,
            ..Default::default()
        };
        if let Some(c) = self.contigs.get(self.contig) {
            let (vs, ve) = self.visible();
            s.locus = Some(format!("{}:{}-{}", c.name, vs.max(0.0) as u64 + 1, ve.min(self.contig_len()) as u64));
            s.highlight = self
                .highlight
                .and_then(|(ci, pos)| self.contigs.get(ci).map(|c| format!("{}:{}", c.name, pos + 1)));
        }
        s.set_style(self.style);
        s.bookmarks = self.bookmarks.borrow().clone();
        s.overlays = self.overlays.borrow().clone();
        s.track_filters = self.track_filters.iter().map(|(p, e)| (p.clone(), crate::expr::text(*e))).collect();
        s.variant_colors = self
            .variants
            .iter()
            .filter(|v| v.color != variants::SiteColor::default_for(&v.info_fields))
            .map(|v| (v.path.clone(), v.color.name().to_owned()))
            .collect();
        if self.arrangement != arrange::Arrangement::default() {
            s.track_layout = self.arrangement.saved(&self.arrangement.ordered(&self.track_keys()));
        }
        s
    }

    /// Replaces everything open with `s`. With `whole_window` (a session opened by the person or
    /// an agent) the window's panes are rebuilt to the session's, one pane when it saved one.
    fn apply_session(&mut self, mut s: Session, whole_window: bool, ctx: &egui::Context) {
        self.note_missing(&mut s);
        if s.recovered {
            self.notify("The session file was damaged; opened its last complete save instead");
        }
        if s.panes.len() > 1 {
            // The first pane is this one; the window adds the others.
            s.locus = s.panes[0].locus.clone().or(s.locus);
            s.highlight = s.panes[0].highlight.clone();
        }
        if whole_window || s.panes.len() > 1 {
            self.pending_panes = Some(s.clone());
        }
        // What was picked in the old workspace names contigs by index: in another genome the
        // same index is another contig.
        self.selection = None;
        self.selecting = false;
        self.zoom_drag = None;
        self.measure = None;
        self.measuring = false;
        self.selected_read = None;
        self.variant_popup = None;
        self.menu_at = None;
        self.drag_track = None;
        self.back.clear();
        self.forward.clear();
        self.deferred_bams.clear();
        self.retired.extend(self.tracks.iter().map(|t| t.id));
        self.retired.extend(self.genes.iter().map(|g| g.id));
        self.retired.extend(self.variants.iter().map(|v| v.id));
        self.tracks.clear();
        self.genes.clear();
        self.variants.clear();
        self.signals.clear();
        self.reference = None;
        self.contigs.clear();
        self.contigs_from_annotation = false;
        self.highlight = None;
        self.reference_path = None;
        self.pending_reference = None;
        if let Some(r) = s.reference.clone() {
            self.set_reference(r, ctx);
        }
        self.style = s.style();
        self.all_transcripts = s.all_transcripts;
        self.add_tracks(s.bams.clone(), ctx);
        self.add_genes(s.genes.clone(), ctx);
        self.add_variants(s.vcfs.clone(), ctx);
        for v in &mut self.variants {
            if let Some(c) = s.variant_colors.get(&v.path) {
                v.color = variants::SiteColor::parse(c);
            }
        }
        self.add_signals(s.signals.clone(), ctx);
        // A review of a VCF the session opens again carries on (its track is a new one).
        self.rebind_review();
        self.pending_goto = s.locus.clone().or_else(|| self.contigs.first().map(|c| c.name.clone()));
        self.pending_highlight = s.highlight;
        if !s.bookmarks.is_empty() || !self.bookmarks.borrow().is_empty() {
            *self.bookmarks.borrow_mut() = s.bookmarks;
        }
        if !s.overlays.is_empty() || !self.overlays.borrow().is_empty() {
            *self.overlays.borrow_mut() = s.overlays;
        }
        self.track_filters = s.track_filters.iter().filter_map(|(p, e)| Some((p.clone(), crate::expr::intern(e).ok()?))).collect();
        self.arrangement = arrange::Arrangement::from_saved(&s.track_layout);
    }

    fn save_session_as(&mut self) {
        self.ask(rfd::FileDialog::new().add_filter("GGV session", &["json"]).set_file_name("session.ggv.json"), dialogs::Kind::Save, dialogs::Then::SaveSession);
    }

    /// Opens the session at `path` in place of the workspace; one that can't be read leaves
    /// the workspace as it was and says why.
    fn open_session(&mut self, path: &std::path::Path, ctx: &egui::Context) -> Result<(), String> {
        match Session::load(path) {
            Ok(s) => {
                self.apply_session(s, true, ctx);
                Ok(())
            }
            Err(e) => {
                let e = format!("{e:#}");
                self.error = Some(e.clone());
                Err(e)
            }
        }
    }

    /// Moves track `hit` up or down within its kind (variants, alignments or genes).
    fn move_track(&mut self, hit: Hit, up: bool) {
        fn shift<T>(v: &mut [T], i: usize, up: bool) {
            let j = if up { i.checked_sub(1) } else { (i + 1 < v.len()).then_some(i + 1) };
            if let Some(j) = j {
                v.swap(i, j);
            }
        }
        match hit {
            Hit::Track(i) => shift(&mut self.tracks, i, up),
            Hit::Genes(i) => shift(&mut self.genes, i, up),
            Hit::Variants(i) => shift(&mut self.variants, i, up),
            Hit::Signal(i) => shift(&mut self.signals, i, up),
        }
        self.sync_compare();
    }

    fn close_hit(&mut self, hit: Hit) {
        match hit {
            Hit::Track(i) => self.close_track(i),
            Hit::Genes(i) => self.close_genes(i),
            Hit::Variants(i) => self.close_variants(i),
            Hit::Signal(i) => self.close_signal(i),
        }
    }

    /// `paths` without files already open (compared as the files they name, so a relative path
    /// or another spelling counts) or given twice: opening one again made a second, identical
    /// track that couldn't then be closed by path.
    fn not_open(&mut self, paths: Vec<PathBuf>, open: &[PathBuf]) -> Vec<PathBuf> {
        let same = |p: &std::path::Path| std::fs::canonicalize(p).unwrap_or_else(|_| p.to_owned());
        let mut seen: Vec<PathBuf> = open.iter().map(|p| same(p)).collect();
        let mut again = Vec::new();
        let mut keep = Vec::new();
        for p in paths {
            let c = same(&p);
            if seen.contains(&c) {
                again.push(crate::files::display_name(&p));
            } else {
                seen.push(c);
                keep.push(p);
            }
        }
        if !again.is_empty() {
            self.notify(format!("Already open: {}", again.join(", ")));
        }
        keep
    }

    fn add_signals(&mut self, paths: Vec<PathBuf>, ctx: &egui::Context) {
        let paths = self.existing(paths, workspace::OpenedAs::Signal);
        for path in paths {
            if !self.signals.iter().any(|s| s.path == path) {
                self.signals.push(crate::signals::SignalTrack::open(path, ctx));
            }
        }
    }

    fn close_signal(&mut self, i: usize) {
        self.signals.remove(i);
        self.drag_track = None;
    }

    fn add_variants(&mut self, paths: Vec<PathBuf>, ctx: &egui::Context) {
        let paths = self.existing(paths, workspace::OpenedAs::Variants);
        let open: Vec<PathBuf> = self.variants.iter().map(|v| v.path.clone()).collect();
        let paths = self.not_open(paths, &open);
        let paths: Vec<PathBuf> = paths.into_iter().filter(|p| !self.needs_index(p)).collect();
        for path in paths {
            let id = gpu::new_slot();
            match VariantTrack::open(id, path.clone(), ctx) {
                Ok(v) => self.variants.push(v),
                Err(e) => self.error = Some(path_error(&path, format!("{e:#}"))),
            }
        }
    }

    fn close_variants(&mut self, i: usize) {
        let v = self.variants.remove(i);
        // Its review ends with it (every decision is already saved).
        if self.review.as_ref().is_some_and(|r| r.track == v.id) {
            self.review = None;
        }
        self.retired.push(v.id);
        self.drag_track = None;
    }

    /// Asks the variant track under the pointer (or the first) for the next/previous site from
    /// the outlined base, or from the view's centre.
    fn seek_variant(&mut self, forward: bool, track: Option<usize>) {
        if self.contigs.is_empty() || self.variants.is_empty() {
            return;
        }
        let i = track.filter(|&i| i < self.variants.len()).unwrap_or(0);
        let from = match self.highlight {
            Some((c, pos)) if c == self.contig => pos,
            _ => self.center.max(0.0) as u64,
        };
        // This contig, then the following (or preceding) ones in genome order.
        let pair = |c: &Contig| (c.name.clone(), c.length);
        let contigs: Vec<(String, u64)> = if forward {
            self.contigs[self.contig..].iter().map(pair).collect()
        } else {
            self.contigs[..=self.contig].iter().rev().map(pair).collect()
        };
        self.variants[i].seek(contigs, from, forward, &mut self.next_request_id);
    }

    /// Centres a found variant and outlines it, zooming in to read level if the view is wide.
    /// The "Missing indexes" window: offers to build indexes, runs them on a thread, and opens
    /// each file once its index exists.
    fn index_window(&mut self, ctx: &egui::Context) {
        // Finished builds: open the file (or re-set the reference).
        let mut ready = Vec::new();
        for job in &mut self.indexing {
            if let IndexState::Building(rx) = &job.state
                && let Ok(result) = rx.try_recv()
            {
                job.state = match result {
                    Ok(()) => IndexState::Done,
                    Err(e) => IndexState::Failed(e),
                };
                if matches!(job.state, IndexState::Done) {
                    ready.push(job.path.clone());
                }
            }
        }
        self.indexing.retain(|j| !matches!(j.state, IndexState::Done));
        for path in ready {
            if Some(&path) == self.pending_reference.as_ref() {
                // Switches now, opening the alignments asked for with it too.
                self.set_reference(path, ctx);
            } else {
                self.open_paths(vec![path], ctx);
            }
        }
        if self.indexing.is_empty() {
            return;
        }

        let mut start = Vec::new();
        let mut dismiss = Vec::new();
        egui::Window::new("Missing indexes")
            .id(egui::Id::unique(("missing-indexes", self.reference_slot)))
            .collapsible(false)
            .resizable(false)
            .anchor(Align2::CENTER_CENTER, Vec2::ZERO)
            .show(ctx, |ui| {
            ui.label("These files have no index, which GGV needs for fast region queries.");
            ui.label("Building one reads the whole file once (large BAMs take a few minutes).");
            ui.separator();
            for (k, job) in self.indexing.iter().enumerate() {
                ui.horizontal(|ui| {
                    ui.label(crate::files::display_name(&job.path)).on_hover_text(job.path.display().to_string());
                    match &job.state {
                        IndexState::Offered => {
                            if ggv_io::index_kind(&job.path).is_some() {
                                if ui.button(format!("Build .{}", ggv_io::index_kind(&job.path).unwrap_or(""))).clicked() {
                                    start.push(k);
                                }
                            } else {
                                ui.label("can't index this here (bgzip/samtools first)");
                            }
                            if ui.button("Dismiss").clicked() {
                                dismiss.push(k);
                            }
                        }
                        IndexState::Building(_) => {
                            ui.spinner();
                            ui.label("indexing…");
                        }
                        IndexState::Failed(e) => {
                            ui.colored_label(Color32::from_rgb(200, 60, 60), e);
                            if ui.button("Dismiss").clicked() {
                                dismiss.push(k);
                            }
                        }
                        IndexState::Done => {}
                    }
                });
            }
            if self.indexing.iter().filter(|j| matches!(j.state, IndexState::Offered)).count() > 1 && ui.button("Build all").clicked() {
                start.extend(self.indexing.iter().enumerate().filter(|(_, j)| matches!(j.state, IndexState::Offered)).map(|(k, _)| k));
            }
        });
        for k in start {
            let path = self.indexing[k].path.clone();
            let (tx, rx) = std::sync::mpsc::channel();
            let c = ctx.clone();
            std::thread::spawn(move || {
                let _ = tx.send(ggv_io::build_index(&path).map(|_| ()).map_err(|e| e.to_string()));
                c.request_repaint();
            });
            self.indexing[k].state = IndexState::Building(rx);
        }
        for k in dismiss.into_iter().rev() {
            let job = self.indexing.remove(k);
            if Some(&job.path) == self.pending_reference.as_ref() {
                // The tracks stay on the reference they have; the alignments asked for with this
                // one open on it rather than being lost.
                self.pending_reference = None;
                let deferred = std::mem::take(&mut self.deferred_bams);
                self.add_tracks(deferred, ctx);
            }
        }
        ctx.request_repaint_after(std::time::Duration::from_millis(250));
    }

    /// Saves a requested PNG export once the screenshot arrives, cropped to the tracks.
    fn start_review(&mut self, track: usize) {
        let Some(v) = self.variants.get(track) else { return };
        self.queue = None; // one review at a time: they share the keys
        self.review = Some(crate::review::Review::open(v.id, &v.path));
        // Start at the outlined site if it is one; otherwise at the next site from here.
        if self.current_site().is_none() {
            self.seek_variant(true, Some(track));
        }
    }

    /// Where the reviewed VCF's track is now (tracks move; a closed one is gone).
    fn review_track(&self) -> Option<usize> {
        let id = self.review.as_ref()?.track;
        self.variants.iter().position(|v| v.id == id)
    }

    /// The reviewed track's records starting at the outlined base (file order), and which of
    /// them the review is on.
    fn review_records(&self) -> Option<((usize, u64), Vec<&ggv_io::Variant>, usize)> {
        let (c, pos) = self.highlight.filter(|(c, _)| *c == self.contig)?;
        let v = self.variants.get(self.review_track()?)?;
        let here: Vec<&ggv_io::Variant> = v.current.as_ref()?.variants.iter().filter(|x| x.start == pos).collect();
        let last = here.len().checked_sub(1)?;
        let k = self.review.as_ref()?.record.filter(|r| (r.0, r.1) == (c, pos)).map_or(0, |r| r.2.min(last));
        Some(((c, pos), here, k))
    }

    /// The reviewed track's variant at the outlined base, with its key.
    fn current_site(&self) -> Option<(crate::review::SiteKey, ggv_io::Variant)> {
        let ((c, pos), here, k) = self.review_records()?;
        let var = here[k].clone();
        let key = crate::review::SiteKey {
            contig: self.contigs[c].name.clone(),
            pos: pos + 1,
            reference: var.reference.clone(),
            alt: var.alts.join(","),
        };
        Some((key, var))
    }

    /// The review's next or previous record: the other records at the outlined base first
    /// (seeking only finds other bases), then the nearest base with any.
    fn review_step(&mut self, forward: bool) {
        if let Some(((c, pos), here, k)) = self.review_records() {
            let to = if forward { (k + 1 < here.len()).then_some(k + 1) } else { k.checked_sub(1) };
            if let (Some(j), Some(r)) = (to, &mut self.review) {
                r.record = Some((c, pos, j));
                return;
            }
        }
        if let Some(r) = &mut self.review {
            r.back = !forward;
        }
        let track = self.review_track();
        self.seek_variant(forward, track);
    }

    /// A record clicked in the reviewed track: the review goes to it.
    fn review_pick(&mut self, var: &ggv_io::Variant) {
        let Some(v) = self.review_track().and_then(|i| self.variants.get(i)) else { return };
        let k = v.current.as_ref().and_then(|c| c.variants.iter().filter(|x| x.start == var.start).position(|x| x.reference == var.reference && x.alts == var.alts));
        self.highlight = Some((self.contig, var.start));
        if let Some(r) = &mut self.review {
            r.record = Some((self.contig, var.start, k.unwrap_or(0)));
        }
    }

    /// Keeps the review on its VCF when that was reopened (a session, a reload: a new track
    /// for the same file), or ends it when the file is no longer open.
    fn rebind_review(&mut self) {
        let Some(r) = &mut self.review else { return };
        if self.variants.iter().any(|v| v.id == r.track) {
            return;
        }
        match self.variants.iter().find(|v| crate::files::same_file(&v.path, &r.source)) {
            Some(v) => r.track = v.id,
            None => {
                self.review = None;
                self.notify("The variant review ended: its VCF is no longer open (every decision was saved)");
            }
        }
    }

    /// Records a decision for the outlined site and moves to the next one.
    fn decide(&mut self, label: crate::review::Label) {
        let Some((key, _)) = self.current_site() else { return };
        let Some(review) = &mut self.review else { return };
        // This site's own note into the editor first: the panel does it when drawn, but keys
        // and agents decide while it's collapsed or the window is hidden, and the note typed for
        // the previous site would be saved onto this one.
        review.arrive(&key);
        self.person_events.push(serde_json::json!({
            "kind": "review", "locus": format!("{}:{}", key.contig, key.pos), "ref": key.reference, "alt": key.alt, "decision": label.name(),
        }));
        review.set(key, label);
        self.review_step(true);
    }

    /// The variant review: the site, its decision and note, and the controls.
    pub(super) fn review_ui(&mut self, ui: &mut egui::Ui) {
        if self.review.is_none() {
            return;
        }
        let site = self.current_site();
        let mut decide = None;
        let mut step = None;
        let mut close = false;
        let reviewed = self.review_track();
        let file = reviewed.and_then(|i| self.variants.get(i)).map(|v| v.file_name()).unwrap_or_default();
        let review = self.review.as_mut().expect("checked above");
        if let Some((key, _)) = &site {
            review.arrive(key);
        }
        let [real, artifact, unsure] = review.counts();
        ui.label(format!("{file}: {} reviewed ({real} real, {artifact} artifact, {unsure} unsure)", real + artifact + unsure));
        ui.label(egui::RichText::new(format!("saved to {}", crate::files::display_name(&review.path))).weak());
        ui.separator();
        match &site {
            Some((key, var)) => {
                ui.label(egui::RichText::new(format!("{}:{}  {} > {}", key.contig, commas(key.pos), allele(&key.reference), alleles(&key.alt))).strong());
                ui.label(format!(
                    "{}  QUAL {}  FILTER {}",
                    var.kind(),
                    var.qual.map_or(".".to_owned(), |q| format!("{q:.1}")),
                    if var.filters.is_empty() { ".".to_owned() } else { var.filters.join(",") }
                ));
                if let Some(d) = review.decisions.get(key) {
                    ui.colored_label(
                        Color32::from_rgb(d.label.rgb()[0], d.label.rgb()[1], d.label.rgb()[2]),
                        format!("marked {}{}", d.label.name(), if d.note.is_empty() { String::new() } else { format!(": {}", d.note) }),
                    );
                }
            }
            // Looking for the next site, or none outlined (Escape cleared it, or the view moved
            // off the sites): then say how to get one, rather than a spinner that never stops.
            None if reviewed.and_then(|i| self.variants.get(i)).is_some_and(|v| v.is_seeking() || v.is_loading()) => {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label("finding the site…");
                });
            }
            None => {
                ui.label(egui::RichText::new("No site outlined: Next (N) or Prev (P) goes to one, or click a variant.").weak());
            }
        }
        ui.add(egui::TextEdit::singleline(&mut review.note).hint_text("note (optional)").desired_width(260.0));
        ui.horizontal(|ui| {
            use crate::review::Label;
            for (label, key) in [(Label::Real, "1"), (Label::Artifact, "2"), (Label::Unsure, "3")] {
                let [r, g, b] = label.rgb();
                if ui.add(egui::Button::new(egui::RichText::new(format!("{} ({key})", label.name())).color(Color32::from_rgb(r, g, b)))).clicked() {
                    decide = Some(label);
                }
            }
        });
        ui.horizontal(|ui| {
            if ui.button("< Prev (P)").clicked() {
                step = Some(false);
            }
            if ui.button("Next (N) >").clicked() {
                step = Some(true);
            }
            if ui.button("Stop").clicked() {
                close = true;
            }
        });
        if let Some(e) = &review.error {
            ui.colored_label(Color32::from_rgb(200, 60, 60), format!("Not saved: {e}"));
        }
        if let Some(n) = &review.notice {
            ui.colored_label(Color32::from_rgb(190, 120, 0), n);
        }
        let mut choose_dir = false;
        // Where it saves is said above (the full path on hover there); here, choosing another folder.
        if ui.small_button("Save reviews in…").on_hover_text(format!("Saved to {}. Keep review and queue decisions in a folder you choose (for data on a read-only drive)", review.path.display())).clicked() {
            choose_dir = true;
        }
        if choose_dir {
            self.ask(rfd::FileDialog::new().set_title("Folder for review decisions"), dialogs::Kind::PickFolder, dialogs::Then::ReviewDir);
        }
        if let Some(label) = decide {
            self.decide(label);
        }
        if let Some(forward) = step {
            self.review_step(forward);
        }
        if close {
            self.review = None;
        }
    }

    /// The outlined base, while it is on screen.
    fn outlined_on_screen(&self) -> Option<(usize, u64)> {
        let (vs, ve) = self.visible();
        self.highlight.filter(|&(c, p)| c == self.contig && (p as f64) + 1.0 > vs && (p as f64) < ve)
    }

    /// The base the evidence panel describes: the outlined base while it is on screen, else the
    /// view's centre (so panning away from an outlined base follows the view again).
    fn evidence_base(&self) -> Option<(usize, u64)> {
        self.outlined_on_screen().or_else(|| (!self.contigs.is_empty()).then(|| (self.contig, self.center.max(0.0) as u64)))
    }

    fn evidence_panel(&mut self, ui: &mut egui::Ui) {
        let Some((c, pos)) = self.evidence_base() else { return };
        let contig = self.contigs[c].name.clone();
        let filters: Vec<crate::tracks::ReadFilter> = self.tracks.iter().map(|t| self.track_style(&t.loader.bam_path).filter).collect();
        for (t, filter) in self.tracks.iter_mut().zip(filters) {
            t.request_evidence(&contig, pos, filter, &mut self.next_request_id);
        }
        ui.horizontal(|ui| {
            ui.heading(format!("{contig}:{}", commas(pos + 1)));
            if ui.small_button("Copy").on_hover_text("Copy the locus").clicked() {
                ui.ctx().copy_text(format!("{contig}:{}", pos + 1));
            }
        });
        let how = match (self.outlined_on_screen().is_some(), self.highlight.is_some()) {
            (true, _) => "outlined base · click another base to move it, click it again (or Esc) to clear",
            (false, true) => "view centre (the outlined base is off screen) · click a base to outline it",
            (false, false) => "view centre · click a base to outline it and pin the counts there",
        };
        ui.label(egui::RichText::new(how).weak());
        let exprs: Vec<String> = self.tracks.iter().filter_map(|t| self.track_style(&t.loader.bam_path).filter.expr).map(crate::expr::text).collect::<std::collections::BTreeSet<_>>().into_iter().collect();
        if !exprs.is_empty() {
            ui.label(egui::RichText::new(format!("counts use the filter: {}", exprs.join(" / "))).weak());
        }

        // Genes at this base.
        for l in self.genes_at(pos) {
            ui.label(l);
        }

        // What the changes here do to the genes.
        let lines = self.consequence_lines(pos);
        if !lines.is_empty() {
            ui.separator();
            ui.label(egui::RichText::new("Consequences").strong());
            for l in lines {
                ui.monospace(l);
            }
        }

        // Variant records here, with genotypes.
        for v in &self.variants {
            let Some(cur) = &v.current else { continue };
            for var in cur.variants.iter().filter(|x| x.start <= pos && pos < x.end) {
                ui.separator();
                ui.label(egui::RichText::new(format!("{}  {} > {}", v.file_name(), allele(&var.reference), alleles(&var.alts.join(",")))).strong());
                ui.label(format!(
                    "QUAL {}  FILTER {}  AF {}",
                    var.qual.map_or(".".to_owned(), |q| format!("{q:.1}")),
                    if var.filters.is_empty() { ".".to_owned() } else { var.filters.join(",") },
                    var.af_text()
                ));
                for (s, call) in v.samples.iter().zip(&var.samples) {
                    let ad = call.ad_text().map_or_else(String::new, |ad| format!("  AD {ad}"));
                    ui.monospace(format!("{s}: {}{}{}", call.gt, call.dp.map_or(String::new(), |d| format!("  DP {d}")), ad));
                }
            }
        }

        // Exact pileup per alignment track.
        let mut copy: Option<String> = None;
        for t in &self.tracks {
            ui.separator();
            ui.horizontal(|ui| {
                let name = ui
                    .add(egui::Label::new(egui::RichText::new(t.title()).strong()).sense(egui::Sense::click()))
                    .on_hover_text(format!("{}\nsample: {}\n{}\nclick to copy", t.file_name(), t.label(), t.loader.bam_path.display()));
                if name.clicked() {
                    copy = Some(t.full_title());
                }
                if let Some(s) = t.evidence.as_ref().and_then(|e| e.stats.as_ref())
                    && ui.small_button("Copy").on_hover_text("Copy these counts as text").clicked()
                {
                    copy = Some(evidence_text(&t.full_title(), &format!("{contig}:{}", pos + 1), s));
                }
            });
            let local = t.contig_name(&contig).map(|c| c.name.clone());
            let Some(ev) = t.evidence.as_ref().filter(|e| e.pos == pos && Some(&e.contig) == local.as_ref()) else {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label("counting…");
                });
                continue;
            };
            if let Some(e) = &ev.error {
                ui.colored_label(Color32::from_rgb(200, 60, 60), e);
                continue;
            }
            let Some(s) = &ev.stats else {
                ui.label("no reads here");
                continue;
            };
            let depth = f64::from(s.depth.max(1));
            ui.label(format!("depth {}  ref {}", s.depth, s.ref_base.unwrap_or('.')));
            // The allele table keeps its columns; if the panel is narrower it scrolls sideways
            // on its own while the text around it wraps.
            egui::ScrollArea::horizontal().id_salt(("evidence-table", t.id)).show(ui, |ui| egui::Grid::new(("evidence", t.id)).striped(true).num_columns(5).show(ui, |ui| {
                for h in ["", "reads", "%", "fwd", "rev"] {
                    ui.label(egui::RichText::new(h).weak());
                }
                ui.end_row();
                let b = &s.bases;
                let rows = [("A", b.a), ("C", b.c), ("G", b.g), ("T", b.t), ("N", b.n), ("DEL", s.deletions), ("INS after", s.insertions_after)];
                for (name, c) in rows.into_iter().filter(|(_, c)| c.total > 0) {
                    let is_ref = s.ref_base.is_some_and(|r| r.to_string() == name);
                    let label = if is_ref { format!("{name} (ref)") } else { name.to_owned() };
                    ui.label(label);
                    ui.monospace(c.total.to_string());
                    ui.monospace(format!("{:.1}", 100.0 * f64::from(c.total) / depth));
                    ui.monospace(c.fwd.to_string());
                    ui.monospace(c.rev.to_string());
                    ui.end_row();
                }
            }));
            ui.label(format!(
                "MAPQ mean {:.1} ({} at 0) · base Q mean {:.1}{}",
                s.mean_mapping_quality,
                s.mapping_quality_zero,
                s.mean_base_quality,
                if s.low_base_quality > 0 { format!(" ({} below the filter)", s.low_base_quality) } else { String::new() }
            ));
            if s.soft_clipped_left + s.soft_clipped_right > 0 || s.ref_skips > 0 {
                ui.label(format!("soft clips starting {} / ending {} · ref skips {}", s.soft_clipped_left, s.soft_clipped_right, s.ref_skips));
            }
            if let Some(alt) = &s.top_alt {
                ui.label(format!("top alt {} at {:.1}% ({} fwd / {} rev){}", alt.allele, alt.fraction * 100.0, alt.count.fwd, alt.count.rev, strand_note(s)));
            }
        }
        if let Some(text) = copy {
            ui.ctx().copy_text(text);
            self.notify("Copied");
        }
    }

    /// Centres `pos` on `contig` (any naming), zooming to read level if the view is wide.
    fn jump_to(&mut self, contig: &str, pos: u64) {
        let aliases = ggv_model::contig_aliases(contig);
        let Some(c) = self.contigs.iter().position(|x| x.name == contig || aliases.contains(&x.name)) else {
            self.error = Some(format!("{contig} is not in this genome"));
            return;
        };
        self.remember_place();
        self.contig = c;
        self.center = pos as f64 + 0.5;
        if f64::from(self.width) / self.ppb > read_window() as f64 {
            self.ppb = f64::from(self.width) / 1001.0;
        }
        self.clamp_view();
    }

    /// The region-statistics window for a shift-dragged selection.
    /// The selected region: copy its sequence, each track's statistics, bookmark it.
    pub(super) fn region_ui(&mut self, ui: &mut egui::Ui) {
        let ctx = &ui.ctx().clone();
        let Some((c, start, end)) = self.selection else { return };
        if self.selecting {
            return;
        }
        let contig = self.contigs.get(c).map(|x| x.name.clone()).unwrap_or_default();
        // Statistics not yet asked for this selection (it was set by a search hit, a bookmark,
        // a link...), or a filter changed since: ask now.
        let filters: Vec<crate::tracks::ReadFilter> = self.tracks.iter().map(|t| self.track_style(&t.loader.bam_path).filter).collect();
        let key = (contig.clone(), start, end);
        for (t, filter) in self.tracks.iter_mut().zip(filters) {
            if end - start <= 10_000_000 && (t.summary_for.as_ref() != Some(&key) || (t.summary.is_some() && t.summary_filter != Some(filter))) {
                t.request_summary(&contig, start, end, filter, &mut self.next_request_id);
            }
        }
        let mut close = false;
        // Copy the reference (None) or track i's sequence (Some(i)).
        let mut copy: Option<Option<usize>> = None;
        let mut bookmark = false;
        let tracks: Vec<(String, Option<(u64, Option<Result<ggv_query::RegionSummary, String>>)>)> = self
            .tracks
            .iter()
            .map(|t| (t.title(), t.summary.clone()))
            .collect();
        ui.label(egui::RichText::new(format!("{contig}:{}-{}  ({} bp)", commas(start + 1), commas(end), commas(end - start))).strong());
        // The selected bases: the reference's, or each sample's own (from its reads).
        ui.horizontal_wrapped(|ui| {
            ui.label("Copy sequence:");
            if ui.add_enabled(self.reference_path.is_some(), egui::Button::new("Reference")).on_hover_text("The reference genome's bases here, as FASTA").on_disabled_hover_text("Open a reference first").clicked() {
                copy = Some(None);
            }
            for (i, (label, _)) in tracks.iter().enumerate() {
                if ui
                    .button(label.as_str())
                    .on_hover_text("This sample's bases here, from its reads, as FASTA: the base most reads show, IUPAC codes at heterozygous bases, N under 3 reads, deleted bases left out. Its read filter applies (filter hp == 1 for one haplotype).")
                    .clicked()
                {
                    copy = Some(Some(i));
                }
            }
        });
        for (label, summary) in &tracks {
            ui.separator();
            ui.label(egui::RichText::new(label).strong());
            match summary {
                Some((_, Some(Ok(s)))) => {
                    let r = &s.reads;
                    ui.label(format!(
                        "reads {} used of {} (dup {}, secondary {}, supplementary {}, low MAPQ {})",
                        commas(r.used),
                        commas(r.seen),
                        commas(r.duplicate),
                        commas(r.secondary),
                        commas(r.supplementary),
                        commas(r.low_mapping_quality)
                    ));
                    ui.label(format!("depth mean {:.1}  median {}  min {}  max {}", s.mean_depth, s.median_depth, s.min_depth, s.max_depth));
                    let b: Vec<String> = s.breadth.iter().map(|(k, v)| format!("{k}x {:.1}%", v * 100.0)).collect();
                    ui.label(format!("breadth {}", b.join("  ")));
                    ui.label(format!("{} variable sites (≥20% non-reference, depth ≥8)", s.variable_site_count));
                    for v in &s.variable_sites {
                        ui.monospace(format!(
                            "  {}  {}  {:.0}%  ({}/{} fwd/rev, depth {})",
                            commas(v.pos),
                            v.change(),
                            v.alt_fraction * 100.0,
                            v.alt_fwd,
                            v.alt_rev,
                            v.depth
                        ));
                    }
                    if ui.button("Copy as JSON").clicked() {
                        ui.ctx().copy_text(serde_json::to_string_pretty(s).unwrap_or_default());
                    }
                }
                Some((_, Some(Err(e)))) => {
                    ui.colored_label(Color32::from_rgb(200, 60, 60), e);
                }
                Some((_, None)) => {
                    ui.horizontal(|ui| {
                        ui.spinner();
                        ui.label("summarizing…");
                    });
                }
                None => {
                    ui.label("not in this file");
                }
            }
        }
        ui.separator();
        ui.horizontal(|ui| {
            if ui.button("Bookmark this region").on_hover_text("Ctrl+D does the same while a region is selected").clicked() {
                bookmark = true;
            }
            if ui.button("Close").clicked() {
                close = true;
            }
        });
        if bookmark {
            self.add_bookmark();
        }
        match copy {
            Some(None) => self.copy_fasta(ctx),
            Some(Some(i)) => self.copy_sample_sequence(ctx, i),
            None => {}
        }
        if close {
            self.selection = None;
            for t in &mut self.tracks {
                t.summary = None;
                t.summary_for = None;
            }
        }
    }

    /// Records the current place before a jump, for Back.
    fn remember_place(&mut self) {
        if self.contigs.is_empty() {
            return;
        }
        let here = (self.contig, self.center, self.ppb);
        let same = |p: &(usize, f64, f64)| p.0 == here.0 && (p.1 - here.1).abs() < 1.0 && (p.2 / here.2 - 1.0).abs() < 1e-6;
        if !self.back.last().is_some_and(same) {
            self.back.push(here);
            if self.back.len() > 100 {
                self.back.remove(0);
            }
        }
        self.forward.clear();
    }

    fn go_back(&mut self) {
        if let Some(prev) = self.back.pop() {
            self.forward.push((self.contig, self.center, self.ppb));
            (self.contig, self.center, self.ppb) = prev;
            self.clamp_view();
        }
    }

    fn go_forward(&mut self) {
        if let Some(next) = self.forward.pop() {
            self.back.push((self.contig, self.center, self.ppb));
            (self.contig, self.center, self.ppb) = next;
            self.clamp_view();
        }
    }

    /// Opens dropped or recent files by kind: alignments, variants, genes, a reference, or a
    /// session.
    fn open_paths(&mut self, paths: Vec<PathBuf>, ctx: &egui::Context) {
        let sorted = crate::files::sort(paths.clone());
        if let Some(session) = sorted.session {
            let _ = self.open_session(&session, ctx);
        }
        if let Some(reference) = sorted.reference {
            self.set_reference(reference, ctx);
        }
        self.add_tracks(sorted.bams, ctx);
        self.add_variants(sorted.vcfs, ctx);
        self.add_genes(sorted.genes, ctx);
        self.add_signals(sorted.signals, ctx);
        if !sorted.unknown.is_empty() {
            let names: Vec<String> = sorted.unknown.iter().map(|p| crate::files::display_name(p)).collect();
            self.error = Some(format!("Not a BAM/CRAM, VCF/BCF, GTF/GFF/BED, FASTA, bigWig/bedGraph/SEG or session: {}", names.join(", ")));
        }
        let known: Vec<PathBuf> = paths.into_iter().filter(|p| p.exists() && crate::files::kind_of(&p.to_string_lossy()).is_some()).collect();
        crate::files::remember(&mut self.recent, &known);
    }

    fn show_variant(&mut self, contig: &str, start: u64, end: u64) {
        let Some(c) = self.contigs.iter().position(|c| c.name == contig) else { return };
        self.remember_place();
        let from = self.contigs.get(self.contig).map(|c| c.name.clone());
        self.contig = c;
        self.reanchor(from.as_deref(), start);
        self.center = (start + end) as f64 / 2.0;
        if f64::from(self.width) / self.ppb > read_window() as f64 {
            self.ppb = f64::from(self.width) / 101.0;
        }
        self.highlight = Some((c, start));
        self.clamp_view();
    }

    fn add_genes(&mut self, paths: Vec<PathBuf>, ctx: &egui::Context) {
        let paths = self.existing(paths, workspace::OpenedAs::Genes);
        let open: Vec<PathBuf> = self.genes.iter().map(|g| g.path.clone()).collect();
        let paths = self.not_open(paths, &open);
        for path in &paths {
            self.check_annotation_build(path);
        }
        for path in paths {
            let id = gpu::new_slot();
            self.genes.push(GeneTrack::open(id, path, ctx));
        }
    }

    fn close_genes(&mut self, i: usize) {
        let g = self.genes.remove(i);
        self.retired.push(g.id);
        self.drag_track = None;
    }

    /// Switches the view's genome to `contigs`, staying on the same contig when it exists there.
    /// Everything else that names a contig by index (the selection, the outlined base, history…)
    /// follows its contig by name, or goes: the same index is another contig in another genome.
    fn set_genome(&mut self, contigs: Vec<Contig>) {
        let old = std::mem::replace(&mut self.contigs, contigs);
        let new = &self.contigs;
        let map = |c: usize| {
            let n = &old.get(c)?.name;
            let aliases = ggv_model::contig_aliases(n);
            new.iter().position(|x| x.name == *n || aliases.contains(&x.name))
        };
        self.contig = map(self.contig).unwrap_or(0);
        self.selection = self.selection.and_then(|(c, a, b)| Some((map(c)?, a, b)));
        self.highlight = self.highlight.and_then(|(c, p)| Some((map(c)?, p)));
        self.measure = self.measure.and_then(|(c, a, b)| Some((map(c)?, a, b)));
        self.back = self.back.iter().filter_map(|&(c, x, z)| Some((map(c)?, x, z))).collect();
        self.forward = self.forward.iter().filter_map(|&(c, x, z)| Some((map(c)?, x, z))).collect();
        if let Some(r) = &mut self.review {
            r.record = r.record.and_then(|(c, p, k)| Some((map(c)?, p, k)));
        }
        if self.selection.is_none() {
            self.selecting = false;
        }
        if self.measure.is_none() {
            self.measuring = false;
        }
        self.menu_at = None;
        self.requested_view = None;
        *self.consequence_cache.borrow_mut() = None;
    }

    /// Queues an offer to index `path` if it has no index; true if it was queued (so the caller
    /// should open it later, when the index exists).
    fn needs_index(&mut self, path: &std::path::Path) -> bool {
        if !missing_index(path) {
            return false;
        }
        if !self.indexing.iter().any(|j| j.path == path) {
            self.indexing.push(IndexJob { path: path.to_owned(), state: IndexState::Offered });
        }
        true
    }

    fn add_tracks(&mut self, paths: Vec<PathBuf>, ctx: &egui::Context) {
        let paths = self.existing(paths, workspace::OpenedAs::Alignments);
        crate::files::remember(&mut self.recent, &paths);
        let open: Vec<PathBuf> = self.tracks.iter().map(|t| t.loader.bam_path.clone()).collect();
        let paths = self.not_open(paths, &open);
        if self.pending_reference.is_some() {
            // Opened once the reference they come with is indexed.
            for p in paths {
                if !self.deferred_bams.contains(&p) {
                    self.deferred_bams.push(p);
                }
            }
            return;
        }
        let paths: Vec<PathBuf> = paths.into_iter().filter(|p| !self.needs_index(p)).collect();
        let opening = !paths.is_empty();
        let mut errors = Vec::new();
        for path in paths {
            let id = gpu::new_slot();
            match Track::open(id, path.clone(), self.reference_path.clone(), ctx) {
                // A file that opens but looks truncated is still shown, with its warning.
                Ok(track) => self.tracks.push(track),
                Err(e) => errors.push(path_error(&path, format!("{e:#}"))),
            }
        }
        // Nothing opened leaves an earlier error (a session that didn't load) in place.
        if opening {
            self.error = (!errors.is_empty()).then(|| errors.join("; "));
        }
        if self.contigs_from_annotation
            && let Some(first) = self.tracks.first()
        {
            // A BAM defines the genome better than an annotation's extents.
            let contigs = first.contigs().to_vec();
            self.set_genome(contigs);
            self.contigs_from_annotation = false;
        }
        if self.contigs.is_empty()
            && let Some(first) = self.tracks.first()
        {
            self.contigs = first.contigs().to_vec();
            if self.pending_goto.is_none() {
                // Start on the first contig, whole; type a locus to zoom in.
                self.pending_goto = self.contigs.first().map(|c| c.name.clone());
            }
        }
    }

    fn close_track(&mut self, i: usize) {
        let track = self.tracks.remove(i);
        self.retired.push(track.id);
        if self.reference.as_ref().is_some_and(|(_, _, (id, _))| *id == track.id) {
            self.reference = None;
        }
        if self.tracks.is_empty() {
            // No BAM left: fall back to an annotation's genome, the reference's or a signal
            // file's (as a window opened with only signal tracks has), else the pane is unusable.
            let contigs = self
                .genes
                .iter()
                .find_map(|g| g.annotation().map(|a| a.contigs()))
                .filter(|c| !c.is_empty())
                .or_else(|| self.reference_path.as_ref().and_then(|p| ggv_io::ReferenceFasta::open(p).ok()).map(|r| r.contigs().to_vec()).filter(|c| !c.is_empty()))
                .or_else(|| self.signals.iter().map(|s| s.contigs()).find(|c| !c.is_empty()))
                .unwrap_or_default();
            self.contigs_from_annotation = !contigs.is_empty();
            self.set_genome(contigs);
        } else if i == 0 {
            // The next track defines the genome now; keep the same place if it has the contig.
            let contigs = self.tracks.first().map(|t| t.contigs().to_vec()).unwrap_or_default();
            self.set_genome(contigs);
        }
        self.drag_track = None;
        self.sync_compare();
    }

    /// Reopens every track with a new reference (each loader keeps its own handle). One without
    /// an index waits for it (offered in Missing indexes), the tracks staying on the reference
    /// they have: switching then would close them all.
    fn set_reference(&mut self, path: PathBuf, ctx: &egui::Context) {
        if self.existing(vec![path.clone()], workspace::OpenedAs::Reference).is_empty() {
            return;
        }
        if self.needs_index(&path) {
            self.pending_reference = Some(path);
            return;
        }
        self.pending_reference = None;
        self.reference_path = Some(path);
        let paths: Vec<PathBuf> = self.tracks.iter().map(|t| t.loader.bam_path.clone()).collect();
        self.retired.extend(self.tracks.iter().map(|t| t.id));
        self.tracks.clear();
        self.reference = None;
        let contigs = std::mem::take(&mut self.contigs);
        let goto = self.pending_goto.clone();
        self.add_tracks(paths, ctx);
        if !contigs.is_empty() {
            // Same genome: the view stays where it was (add_tracks, seeing no contigs, would
            // start on the first one).
            self.contigs = contigs;
            self.pending_goto = goto;
        }
        // The alignments that waited for this reference's index (offering theirs next if they
        // lack them).
        let deferred = std::mem::take(&mut self.deferred_bams);
        self.add_tracks(deferred, ctx);
    }

    /// The warnings of the open files (a BAM that looks truncated) not dismissed: shown until
    /// dismissed, but unlike `error` they don't fail an automated run.
    fn warning(&self) -> Option<String> {
        let w: Vec<&str> = self
            .tracks
            .iter()
            .filter(|t| !self.dismissed_warnings.contains(&t.loader.bam_path))
            .filter_map(|t| t.loader.warning.as_deref())
            .collect();
        (!w.is_empty()).then(|| w.join("\n"))
    }

    fn contig_len(&self) -> f64 {
        self.contigs.get(self.contig).map_or(1.0, |c| c.length as f64)
    }

    fn visible(&self) -> (f64, f64) {
        let half = f64::from(self.width) / self.ppb / 2.0;
        (self.center - half, self.center + half)
    }

    fn clamp_view(&mut self) {
        let len = self.contig_len();
        let w = f64::from(self.width.max(1.0));
        // A contig shorter than the view even at the deepest zoom is shown at that zoom.
        self.ppb = self.ppb.clamp((w / len).min(MAX_PPB), MAX_PPB);
        let half = w / self.ppb / 2.0;
        self.center = if half * 2.0 >= len { len / 2.0 } else { self.center.clamp(half, len - half) };
    }

    /// Resolves a locus, or failing that a gene name or id from the loaded annotations (with a
    /// little padding). `Err(None)` means "not known yet": an annotation is still loading.
    fn resolve(&self, text: &str) -> Result<Region, Option<String>> {
        let locus_err = match ggv_io::resolve_locus_str(&self.contigs, text) {
            Ok(region) => return Ok(region),
            Err(e) => e.to_string(),
        };
        if let Some((_, region)) = self.gene_hits(text).into_iter().next() {
            return Ok(region);
        }
        if self.genes.iter().any(GeneTrack::is_loading) {
            return Err(None);
        }
        // A bare word that is no contig was meant as a gene: say so rather than "no such contig".
        let word = !text.contains(':') && !text.trim().is_empty();
        Err(Some(match (word, self.genes.is_empty()) {
            (true, false) => format!("`{}` is neither a contig nor a gene in the open annotations", text.trim()),
            (true, true) => format!("{locus_err} (open a gene annotation to find genes by name)"),
            (false, _) => locus_err,
        }))
    }

    /// Selects `[a, b)` on contig `c` and asks each track for its statistics there.
    fn select_region(&mut self, c: usize, a: u64, b: u64) {
        if b - a > 10_000_000 {
            self.error = Some("Select at most 10 Mb for region statistics".into());
            self.selection = None;
            return;
        }
        self.selection = Some((c, a, b));
        let contig = self.contigs[c].name.clone();
        let filters: Vec<crate::tracks::ReadFilter> = self.tracks.iter().map(|t| self.track_style(&t.loader.bam_path).filter).collect();
        for (t, filter) in self.tracks.iter_mut().zip(filters) {
            t.request_summary(&contig, a, b, filter, &mut self.next_request_id);
        }
    }

    /// Where gene `text` is in each open annotation, best first: annotations whose file name
    /// names the open genome's build (`hg19.ncbiRefSeq.gtf.gz` with hg19 data) before those
    /// naming no build, and those naming another build last. With an hg38 and an hg19 annotation
    /// both open on hg19 data, a gene went to the first one opened (GRCh38 coordinates).
    fn gene_hits(&self, text: &str) -> Vec<(usize, Region)> {
        let build = self.build();
        let mut hits: Vec<(u8, usize, Region)> = Vec::new();
        for (i, g) in self.genes.iter().enumerate() {
            let Some(span) = g.annotation().and_then(|a| a.find(text)) else { continue };
            let aliases = ggv_model::contig_aliases(&span.contig);
            let Some(contig) = self.contigs.iter().find(|c| c.name == span.contig || aliases.contains(&c.name)) else {
                continue;
            };
            let pad = ((span.end - span.start) / 20).max(50);
            let region = Region { contig: contig.name.clone(), start: span.start.saturating_sub(pad), end: (span.end + pad).min(contig.length) };
            let rank = match (build, annotation_build(&g.path)) {
                (Some(b), Some(a)) if a == b => 0,
                (_, None) => 1,
                _ => 2,
            };
            hits.push((rank, i, region));
        }
        hits.sort_by_key(|h| (h.0, h.1));
        hits.into_iter().map(|(_, i, r)| (i, r)).collect()
    }

    /// A gene track whose file names another build than the open genome's: its genes would be
    /// drawn (and found) in the wrong places. Said once, when it opens.
    fn check_annotation_build(&mut self, path: &std::path::Path) {
        if let (Some(open), Some(file)) = (self.build(), annotation_build(path))
            && open != file
        {
            self.error = Some(format!(
                "{} looks like {file}, but the open genome is {open}: its genes will be in the wrong places. Close it, or open the {open} annotation (Open > Tracks).",
                crate::files::display_name(path)
            ));
        }
    }

    /// Grouping by base and sorting name a position but not a contig: after moving to another
    /// contig (from `from`), they follow the view to `base` rather than point into it at random.
    fn reanchor(&mut self, from: Option<&str>, base: u64) {
        if from.is_none_or(|f| self.contigs.get(self.contig).is_some_and(|c| c.name == f)) {
            return;
        }
        if matches!(self.style.group, ggv_layout::GroupBy::BaseAt(p) if p != u64::MAX) {
            self.style.group = ggv_layout::GroupBy::BaseAt(base);
        }
        if let Some((_, by)) = self.style.sort {
            self.style.sort = Some((base, by));
        }
    }

    fn goto(&mut self, text: &str) {
        match self.resolve(text) {
            Err(None) => self.pending_goto = Some(text.to_owned()), // try again once genes load
            Err(Some(e)) => {
                self.error = Some(e.clone());
                self.goto_error = Some(e);
            }
            Ok(region) => {
                self.remember_place();
                let from = self.contigs.get(self.contig).map(|c| c.name.clone());
                self.contig = self.contigs.iter().position(|c| c.name == region.contig).unwrap_or(0);
                let base = if region.len() == 1 { region.start } else { (region.start + region.end) / 2 };
                self.reanchor(from.as_deref(), base);
                let w = f64::from(self.width);
                if region.len() == 1 {
                    // A single base: show it centred with ±50 bp and outlined, like `ggv render`.
                    self.center = region.start as f64 + 0.5;
                    self.ppb = w / 101.0;
                    self.highlight = Some((self.contig, region.start));
                    if let Some(by) = self.pending_sort.take() {
                        self.style.sort = Some((region.start, by));
                    }
                    if self.style.group == ggv_layout::GroupBy::BaseAt(u64::MAX) {
                        self.style.group = ggv_layout::GroupBy::BaseAt(region.start);
                    }
                } else {
                    self.center = (region.start + region.end) as f64 / 2.0;
                    self.ppb = w / region.len() as f64;
                    self.highlight = None;
                }
                for t in &mut self.tracks {
                    t.scroll_y = 0.0;
                }
                // Only the error a failed jump left: startup and loading errors (a filter that
                // didn't parse, a file that didn't open) stay until dealt with.
                if self.goto_error.take().is_some_and(|e| self.error.as_ref() == Some(&e)) {
                    self.error = None;
                }
                self.clamp_view();
            }
        }
    }

    /// Fits the view to the bases from `a` to `b`, given in either order, unless so little is
    /// marked that it was a slip.
    fn zoom_to_marked(&mut self, a: u64, b: u64) {
        let len = self.contig_len() as u64;
        let (start, last) = (a.min(b).min(len), a.max(b).min(len));
        if (last - start) as f64 * self.ppb < MIN_ZOOM_DRAG {
            return;
        }
        let end = (last + 1).min(len);
        self.remember_place();
        self.center = (start + end) as f64 / 2.0;
        self.ppb = f64::from(self.width) / (end - start) as f64;
        self.clamp_view();
    }

    fn zoom_at(&mut self, factor: f64, x_in_canvas: f32) {
        let (vs, _) = self.visible();
        let anchor = vs + f64::from(x_in_canvas) / self.ppb;
        self.ppb *= factor;
        self.clamp_view();
        let new_vs = anchor - f64::from(x_in_canvas) / self.ppb;
        self.center = new_vs + f64::from(self.width) / self.ppb / 2.0;
        self.clamp_view();
    }

    /// The span to load for the current view, in the view's contig naming.
    fn wanted(&self, pixels_per_point: f32) -> Span {
        let contig = &self.contigs[self.contig];
        let (vs, ve) = self.visible();
        let (start, end, reads) = load_span(self.center, ve - vs, contig.length);
        let region = Region { contig: contig.name.clone(), start, end };
        let bp_per_bin = (!reads).then(|| region.len() as f64 / self.bins_for(&region, pixels_per_point) as f64);
        Span { region, reads, bp_per_bin, filter: None, tag: self.style.tag(), mods: self.style.color == tracks::ColorBy::BaseMods }
    }

    fn bins_for(&self, region: &Region, pixels_per_point: f32) -> usize {
        ((region.len() as f64 * self.ppb * f64::from(pixels_per_point)) as usize).clamp(1, MAX_BINS)
    }

    /// Picks up parsed annotations; with no BAM open, the first one defines the genome.
    /// Opened signal files; with nothing else open, the first one supplies the contigs, so a
    /// window of signal tracks alone can be navigated (and become ready).
    fn poll_signals(&mut self) {
        let mut opened = None;
        for s in &mut self.signals {
            if s.poll() && opened.is_none() {
                opened = Some(s.contigs());
            }
        }
        if let Some(contigs) = opened.filter(|c| !c.is_empty())
            && self.contigs.is_empty()
        {
            self.contigs = contigs;
            self.contigs_from_annotation = true;
            if self.pending_goto.is_none() {
                self.pending_goto = self.contigs.first().map(|c| c.name.clone());
            }
        }
    }

    fn poll_genes(&mut self) {
        let mut ready = None;
        for g in &mut self.genes {
            if g.poll() && ready.is_none() {
                ready = g.annotation().map(|a| a.contigs());
            }
        }
        if let Some(contigs) = ready
            && self.contigs.is_empty()
        {
            self.contigs = contigs;
            self.contigs_from_annotation = true;
            if self.pending_goto.is_none() {
                self.pending_goto = self.contigs.first().map(|c| c.name.clone());
            }
        }
    }

    /// Loading without drawing, while the window is hidden or minimised (no `ui` runs then): an
    /// agent that moved the view can still wait for it to be ready.
    pub fn load_hidden(&mut self, pixels_per_point: f32) {
        if self.contigs.is_empty() {
            return;
        }
        self.clamp_view();
        self.update_tracks(pixels_per_point);
    }

    /// Receives finished loads and requests new ones, for every track.
    fn update_tracks(&mut self, pixels_per_point: f32) {
        if self.contigs.is_empty() {
            return;
        }
        let (vs, ve) = self.visible();
        let contig = self.contigs[self.contig].clone();
        let mut found = None;
        for v in &mut self.variants {
            v.dark = self.style.dark;
            if let Some(hit) = v.receive(&mut self.generation) {
                found = Some(hit);
            }
            v.ensure_loaded(&contig.name, contig.length, vs, ve, &mut self.next_request_id);
        }
        match found {
            Some(Some((c, s, e))) => {
                self.show_variant(&c, s, e);
                // Going back, the review arrives at the last record at that base.
                if let (Some(r), Some((ci, pos))) = (&mut self.review, self.highlight) {
                    r.record = Some((ci, pos, if std::mem::take(&mut r.back) { usize::MAX } else { 0 }));
                }
            }
            Some(None) => {
                self.error = Some("No more variants in that direction".to_owned());
                if let Some(r) = &mut self.review {
                    r.back = false;
                }
            }
            None => {}
        }
        let (vs, ve) = self.visible();
        for g in &mut self.genes {
            g.update(&contig.name, contig.length, vs, ve, self.ppb, self.all_transcripts, self.style.dark, &mut self.generation);
        }
        let want = self.wanted(pixels_per_point);
        let (vs, ve) = self.visible();
        let view_contig = self.contigs[self.contig].name.clone();
        let opts = ViewOptions {
            // Read-level loads keep everything mapped and filter in memory; wide loads filter
            // here, since no reads are kept.
            filters: if want.reads {
                ggv_query::PileupOptions { exclude_flags: 0x4, ..Default::default() }
            } else {
                self.style.filter.pileup_options()
            },
            bins: self.bins_for(&want.region, pixels_per_point),
            tag: self.style.tag(),
            mods: self.style.color == tracks::ColorBy::BaseMods,
            max_read_window: read_window() + SPAN_SLACK,
            max_coverage_window: COVERAGE_WINDOW,
            // The app packs reads itself (sorting, grouping) on the loader thread.
            max_rows: 0,
            ..Default::default()
        };
        let styles: Vec<Style> = self.tracks.iter().map(|t| self.track_style(&t.loader.bam_path)).collect();
        for (t, style) in self.tracks.iter_mut().zip(styles) {
            let local = t.contig_name(&view_contig).map(|c| c.name.clone());
            t.receive(local.as_deref(), &mut self.generation);
            t.ensure_loaded(&want, vs, ve, opts.clone(), style, &mut self.next_request_id);
            t.apply_style(style, &mut self.generation);
        }
        self.update_reference();
        let width = self.width;
        for sg in &mut self.signals {
            sg.update(&contig.name, vs, ve, width);
        }
        self.requested_view = Some((self.contig, vs, ve));
    }

    /// The reference row comes from the first track that loaded per-base data with a reference;
    /// with no alignment track, the window reads it itself (zoomed in to 100 kb or less).
    fn update_reference(&mut self) {
        if self.tracks.is_empty() {
            self.own_reference_row();
            return;
        }
        let source = self.tracks.iter().find_map(|t| {
            let cur = t.current.as_ref()?;
            cur.reference.as_ref().map(|r| (t.id, cur, r.clone()))
        });
        let Some((track_id, cur, instances)) = source else { return };
        if self.reference.as_ref().is_some_and(|(_, _, src)| *src == (track_id, cur.id)) {
            return;
        }
        self.generation += 1;
        let set = InstanceSet { generation: self.generation, instances };
        self.reference = Some((cur.view.clone(), set, (track_id, cur.id)));
    }

    fn own_reference_row(&mut self) {
        const OWN: u64 = u64::MAX; // the row's source id when the window read it
        let Some(path) = self.reference_path.clone() else { return };
        let Some(contig) = self.contigs.get(self.contig).cloned() else { return };
        let (vs, ve) = self.visible();
        if ve - vs > read_window() as f64 {
            return;
        }
        let (vs, ve) = (vs.max(0.0) as u64, (ve.max(0.0) as u64).min(contig.length));
        if let Some((view, _, (OWN, _))) = &self.reference
            && view.region.contig == contig.name
            && view.region.start <= vs
            && view.region.end >= ve
        {
            return;
        }
        if self.own_reference.as_ref().is_none_or(|(p, _)| *p != path) {
            match ggv_io::ReferenceFasta::open(&path) {
                Ok(r) => self.own_reference = Some((path, r)),
                Err(e) => {
                    self.error = Some(format!("Reference: {e}"));
                    self.reference_path = None;
                    return;
                }
            }
        }
        let width = ve.saturating_sub(vs).max(1);
        let region = ggv_model::Region { contig: contig.name.clone(), start: vs.saturating_sub(width), end: (ve + width).min(contig.length) };
        let Some((_, reader)) = &mut self.own_reference else { return };
        let Ok(seq) = reader.fetch(&region) else { return };
        let view = ViewData { region, reference: Some(seq), coverage: ggv_layout::Coverage::Unavailable, reads: None, read_counts: Default::default(), label: String::new() };
        self.generation += 1;
        let set = InstanceSet { generation: self.generation, instances: Arc::new(crate::tracks::build_reference(&view)) };
        self.reference = Some((Arc::new(view), set, (OWN, self.generation)));
    }

    fn layout(&self, rect: Rect) -> Layout {
        // In the order and at the sizes the person chose (see `arrange`).
        self.arrange(rect)
    }
}

/// Records the graphics adapter for About and diagnostics.
pub fn set_gpu_info(info: &eframe::wgpu::AdapterInfo) {
    about::set_gpu(info);
}

/// An allele for display: a long one (a structural variant's REF can be thousands of bases)
/// becomes its first bases and its length, e.g. `GCCTT… (6,802 bp)`.
pub(crate) fn allele(a: &str) -> String {
    const SHOWN: usize = 12;
    if a.len() <= 2 * SHOWN {
        a.to_owned()
    } else {
        format!("{}… ({} bp)", &a[..SHOWN], commas(a.len() as u64))
    }
}

/// Comma-separated alleles, each shortened as [`allele`] does.
pub(crate) fn alleles(list: &str) -> String {
    list.split(',').map(allele).collect::<Vec<_>>().join(",")
}

/// Pane entry points for the shell (see `shell.rs`).
impl App {
    pub fn is_scripted(&self) -> bool {
        self.automation.is_some()
    }

    pub fn current_session(&self) -> Session {
        self.session()
    }

    pub fn open_files(&mut self, paths: Vec<PathBuf>, ctx: &egui::Context) {
        self.open_paths(paths, ctx);
    }

    pub fn toolbar_ui(&mut self, ui: &mut egui::Ui, ctx: &egui::Context) {
        self.toolbar(ui, ctx);
    }

    pub fn status_ui(&mut self, ui: &mut egui::Ui, frame: &eframe::Frame) {
        self.status(ui, frame);
    }

    /// The display options panes share (the toolbar edits the active pane's).
    pub fn shared_style(&self) -> (Style, bool) {
        (self.style, self.all_transcripts)
    }

    pub fn set_shared_style(&mut self, (style, all_transcripts): (Style, bool)) {
        // Group-by-base and a sort at a base name a position in this pane's view: the kind is
        // shared, at this pane's own position, else its outlined base, else its centre. So
        // making another pane active (its style pushed to the others) changes no pane.
        let (group, sort) = (self.style.group, self.style.sort);
        let here = self.highlight.filter(|&(c, _)| c == self.contig).map_or(self.center.max(0.0) as u64, |(_, p)| p);
        self.style = style;
        if matches!(style.group, ggv_layout::GroupBy::BaseAt(_)) {
            self.style.group = if matches!(group, ggv_layout::GroupBy::BaseAt(_)) { group } else { ggv_layout::GroupBy::BaseAt(here) };
        }
        if let Some((_, by)) = style.sort {
            self.style.sort = Some((sort.map_or(here, |(p, _)| p), by));
        }
        self.all_transcripts = all_transcripts;
    }

    /// Takes in review and queue decisions another pane saved for the same VCF or list.
    pub fn refresh_decisions(&mut self) {
        if let Some(r) = &mut self.review {
            r.refresh();
        }
        if let Some(q) = &mut self.queue {
            // Read again as the genome now has it (a reference or BAM opened or closed).
            q.settle(&self.contigs);
            q.refresh();
        }
    }

    /// The per-file settings panes share: filter expressions, variant site colours and the
    /// track arrangement (sizes and order).
    pub fn file_settings(&self) -> FileSettings {
        let colors = self.variants.iter().map(|v| (v.path.clone(), v.color.name().to_owned())).collect();
        (self.track_filters.clone(), colors, self.arrangement.clone())
    }

    pub fn set_file_settings(&mut self, (filters, colors, arrangement): &FileSettings) {
        self.track_filters.clone_from(filters);
        for (path, name) in colors {
            if let Some(v) = self.variants.iter_mut().find(|v| v.path == *path)
                && v.color.name() != name
            {
                v.color = variants::SiteColor::parse(name);
            }
        }
        self.arrangement.clone_from(arrangement);
    }

    /// Every GPU slot this pane draws into, for freeing when the pane closes.
    pub fn take_slots(&mut self) -> Vec<u64> {
        let mut ids: Vec<u64> = std::mem::take(&mut self.retired);
        ids.extend(self.tracks.iter().map(|t| t.id));
        ids.extend(self.genes.iter().map(|g| g.id));
        ids.extend(self.variants.iter().map(|v| v.id));
        ids.push(self.reference_slot);
        ids
    }

    /// Slots to free at this pane's next paint.
    pub fn retire(&mut self, ids: Vec<u64>) {
        self.retired.extend(ids);
    }

    /// Where the view is: contig, centre, points per base.
    pub fn view_state(&self) -> Option<(String, f64, f64)> {
        self.contigs.get(self.contig).map(|c| (c.name.clone(), self.center, self.ppb))
    }

    /// Moves the view (for linked panes); a contig this genome lacks is ignored.
    pub fn follow(&mut self, contig: &str, center: f64, ppb: f64) {
        let aliases = ggv_model::contig_aliases(contig);
        let Some(c) = self.contigs.iter().position(|x| x.name == contig || aliases.contains(&x.name)) else { return };
        self.contig = c;
        self.center = center;
        self.ppb = ppb;
        self.clamp_view();
    }

    /// The view as `contig:start-end`, for the pane header.
    pub fn locus_label(&self) -> Option<String> {
        let c = self.contigs.get(self.contig)?;
        let (vs, ve) = self.visible();
        Some(format!("{}:{}-{}", c.name, commas(vs.max(0.0) as u64 + 1), commas(ve.min(self.contig_len()) as u64)))
    }


    /// Makes this pane's files match `src` (same tracks, reference, genes and signals), keeping
    /// its view.
    pub fn sync_files(&mut self, src: &Session, ctx: &egui::Context) {
        let mine = self.session();
        if mine.reference != src.reference
            && let Some(r) = src.reference.clone()
        {
            self.set_reference(r, ctx);
        }
        let mine = self.session();
        for (i, p) in mine.bams.iter().enumerate().rev() {
            if !src.bams.contains(p) {
                self.close_track(i);
            }
        }
        for (i, p) in mine.vcfs.iter().enumerate().rev() {
            if !src.vcfs.contains(p) {
                self.close_variants(i);
            }
        }
        for (i, p) in mine.genes.iter().enumerate().rev() {
            if !src.genes.contains(p) {
                self.close_genes(i);
            }
        }
        for (i, p) in mine.signals.iter().enumerate().rev() {
            if !src.signals.contains(p) {
                self.close_signal(i);
            }
        }
        let add = |want: &[PathBuf], have: &[PathBuf]| want.iter().filter(|p| !have.contains(p)).cloned().collect::<Vec<_>>();
        self.add_tracks(add(&src.bams, &mine.bams), ctx);
        self.add_variants(add(&src.vcfs, &mine.vcfs), ctx);
        self.add_genes(add(&src.genes, &mine.genes), ctx);
        self.add_signals(add(&src.signals, &mine.signals), ctx);
        // And in its order (a track moved there moves here).
        let rank = |want: &[PathBuf], p: &PathBuf| want.iter().position(|w| w == p).unwrap_or(usize::MAX);
        let before = self.track_keys();
        self.tracks.sort_by_key(|t| rank(&src.bams, &t.loader.bam_path));
        self.variants.sort_by_key(|v| rank(&src.vcfs, &v.path));
        self.genes.sort_by_key(|g| rank(&src.genes, &g.path));
        self.signals.sort_by_key(|s| rank(&src.signals, &s.path));
        if self.track_keys() != before {
            // Track indices name other tracks now.
            self.drag_track = None;
            self.menu_at = None;
            self.hover_variant = None;
            self.sync_compare();
        }
    }

    /// Whether `text` names a place in this pane (a locus, region or gene), as `goto` checks it.
    pub fn check_locus(&self, text: &str) -> Result<(), String> {
        if self.contigs.is_empty() {
            return Err("no files are open in this pane".into());
        }
        match self.resolve(text) {
            // Genes still loading: going there waits for them.
            Ok(_) | Err(None) => Ok(()),
            Err(Some(e)) => Err(e),
        }
    }

    /// Per-frame windows (index offers, export, review, region statistics).
    pub fn windows(&mut self, ctx: &egui::Context) {
        self.poll_jobs(ctx);
        self.watch_files(ctx);
        if self.is_active {
            let typing = ctx.memory(|m| m.focused().is_some());
            self.action_keys(ctx, typing);
        }
        self.run_startup(ctx);
        self.index_window(ctx);
        self.finish_capture(ctx);
        self.poll_report(ctx);
        if self.pending_review && !self.contigs.is_empty() && self.pending_goto.is_none() {
            self.pending_review = false;
            self.start_review(0);
        }
        if !self.contigs.is_empty()
            && self.pending_goto.is_none()
            && let Some(text) = self.pending_select.take()
            && let Ok(r) = ggv_io::resolve_locus_str(&self.contigs, &text)
            && let Some(c) = self.contigs.iter().position(|x| x.name == r.contig)
        {
            self.selection = Some((c, r.start, r.end));
            let filters: Vec<crate::tracks::ReadFilter> = self.tracks.iter().map(|t| self.track_style(&t.loader.bam_path).filter).collect();
            for (t, filter) in self.tracks.iter_mut().zip(filters) {
                t.request_summary(&r.contig, r.start, r.end, filter, &mut self.next_request_id);
            }
        }
        // Region, review and queue: docked in the side column or floating (see `inspector`).
        self.floating_sections(ctx);
        self.variant_window(ctx);
        self.bookmarks_window(ctx);
        self.filter_window(ctx);
        self.find_window(ctx);
        self.seq_window(ctx);
        self.compare_window(ctx);
        self.paste_queue_window(ctx);
        self.missing_window(ctx);
        self.update_assembly();
        self.assembly_window(ctx);
        self.shortcuts_window(ctx);
        self.palette_window(ctx);
        self.tour_window(ctx);
        self.poll_dialog(ctx);
        self.error_log.note(self.error.as_ref());
        self.about_window(ctx);
        if std::mem::take(&mut self.pending_demo) {
            self.open_demo(ctx);
        }
    }

    /// The pane itself: the side column (Evidence, Region, review; when any is docked and has
    /// something to show) and the track canvas.
    pub fn body(&mut self, ui: &mut egui::Ui) {
        if self.column_wanted() {
            // Beside the tracks when the pane has room for both (at most half of it); in a narrow
            // pane (three or four split) below them at full width instead, where the allele table
            // and statistics stay readable. Either way it scrolls both ways rather
            // than cutting off what doesn't fit.
            const MIN_SIDE: f32 = 300.0;
            let width = ui.available_width();
            if width * 0.5 >= MIN_SIDE {
                // Text wraps to the panel (it scrolls only up and down); it can't be dragged
                // narrower than its content can wrap to, so nothing is cut off.
                egui::Panel::right(("evidence", self.reference_slot))
                    .default_size(330.0f32.min(width * 0.5))
                    .min_size(240.0)
                    .max_size(width * 0.5)
                    .show(ui, |ui| {
                        self.column(ui);
                    });
            } else {
                let height = ui.available_height();
                egui::Panel::bottom(("evidence-below", self.reference_slot))
                    .resizable(true)
                    .default_size((height * 0.4).clamp(220.0, 360.0))
                    .max_size(height * 0.6)
                    .show(ui, |ui| {
                        self.column(ui);
                    });
            }
        }
        let width = self.width;
        egui::Frame::NONE.fill(self.pal().bg).show(ui, |ui| self.canvas(ui));
        if self.width != width {
            // The toolbar, drawn before the tracks, showed the locus for the old width: a resize
            // with no input after it left the locus box wrong until the next event.
            ui.ctx().request_repaint();
        }
    }

    /// Scripted runs: `others_ready` says the other panes have finished loading too.
    pub fn automate_frame(&mut self, ctx: &egui::Context, frame: &eframe::Frame, others_ready: bool) {
        self.automate(ctx, frame, others_ready);
    }

    /// Everything shown has finished loading (for scripted screenshots and benchmarks).
    pub fn is_ready(&self) -> bool {
        // A track lacking the view's contig has nothing to load, so it counts as ready.
        // An open index offer is a stable state to screenshot, even with nothing else loaded.
        // Loads are requested for the view each frame; until that has happened for this view (just
        // after a jump), what's loaded is the old place's, so not ready.
        let (vs, ve) = self.visible();
        let requested = self.contigs.is_empty()
            || self.requested_view.is_some_and(|(c, s, e)| c == self.contig && (s - vs).abs() < 1e-6 && (e - ve).abs() < 1e-6);
        // Nothing open and nothing on its way (the welcome screen) is a stable state too.
        let empty = self.tracks.is_empty() && self.genes.is_empty() && self.variants.is_empty() && self.signals.is_empty();
        let coming = self.pending_genome.is_some() || self.pending_demo || !self.pending_tracks.is_empty() || !self.downloads.is_empty();
        !self.indexing.is_empty()
            || (empty && !coming && self.automation.is_some())
            || (!self.tracks.is_empty() || !self.genes.is_empty() || !self.variants.is_empty() || !self.signals.is_empty())
            && requested
            && self.signals.iter().all(|s| !s.is_loading())
            // A track that failed to load (a CRAM without its reference, a damaged file) has
            // settled too: waiting for it hung scripted screenshots and agents (it shows its error).
            && self.tracks.iter().all(|t| (t.current.is_some() || t.missing_contig || t.error.is_some()) && !t.is_loading())
            && self.genes.iter().all(|g| !g.is_loading() && (g.annotation().is_none() || g.built.is_some()))
            && self.variants.iter().all(|v| !v.is_loading())
            && (!self.show_evidence || self.tracks.iter().all(|t| t.evidence.is_some() || t.missing_contig))
            && self.review.as_ref().is_none_or(|_| self.current_site().is_some() || self.review_track().and_then(|i| self.variants.get(i)).is_none_or(|v| !v.is_seeking()))
            && self.pending_select.is_none()
            && self.tracks.iter().all(|t| !matches!(t.summary, Some((_, None))))
            && self.pending_goto.is_none()
            && self.capture.is_none()
            && self.pending_report.is_none()
            && self.startup.is_none()
            && self.seq.as_ref().and_then(|s| s.search.as_ref()).is_none_or(|s| s.done)
            && self.find.as_ref().and_then(|f| f.search.as_ref()).is_none_or(|f| f.done)
            && self.compare.as_ref().is_none_or(|c| c.scan.is_none())
    }
}

impl App {
    fn automate(&mut self, ctx: &egui::Context, frame: &eframe::Frame, others_ready: bool) {
        use crate::automation::Step;
        let ready = self.is_ready() && others_ready;
        let loads = (self.tracks.iter().map(|t| t.loads).sum(), self.tracks.iter().map(|t| t.cache_hits).sum());
        let Some(auto) = &mut self.automation else { return };
        let step = auto.step(ctx, ready, frame.info().cpu_usage, loads);
        if matches!(step, Step::Done) && auto.bench && !auto.report_loads_done {
            auto.report_loads_done = true;
            for t in &self.tracks {
                let mut ms = t.load_ms.clone();
                ms.sort_by(f64::total_cmp);
                if let (Some(p50), Some(max)) = (ms.get(ms.len() / 2), ms.last()) {
                    println!("loads {}: {} · p50 {p50:.0} ms · max {max:.0} ms", t.title(), ms.len());
                }
            }
        }
        if auto.report {
            auto.report = false;
            let (vs, ve) = self.visible();
            println!("view {}:{}-{} ({:.0} bp, {:.4} points/base)", self.contigs[self.contig].name, vs as u64 + 1, ve as u64, ve - vs, self.ppb);
            for t in &self.tracks {
                println!("  {}", t.debug_state());
            }
        }
        match step {
            Step::Wait => {}
            Step::Cycle(n) => {
                let (start_ppb, start_center) = *self.cycle_start.get_or_insert((self.ppb, self.center));
                if n <= 60 {
                    self.zoom_at(0.8, self.width / 2.0); // all the way out (clamped at the contig)
                } else if n > 90 {
                    self.zoom_at(1.25, self.width / 2.0);
                    if self.ppb >= start_ppb {
                        self.ppb = start_ppb;
                        self.center = start_center;
                        self.clamp_view();
                        if let Some(a) = &mut self.automation {
                            a.cycle_back = Some(std::time::Instant::now());
                        }
                    }
                }
            }
            Step::Pan(fraction) => self.center += f64::from(self.width * fraction) / self.ppb,
            Step::Zoom(factor) => self.zoom_at(factor, self.width / 2.0),
            Step::Done => {
                // An error still showing (a file that didn't open, an export that failed) fails
                // the run, rather than only appearing in the screenshot.
                if let Some(e) = &self.error
                    && !crate::automation::failed()
                {
                    crate::automation::fail(e);
                }
                ctx.send_viewport_cmd(egui::ViewportCommand::Close)
            }
        }
    }

    fn toolbar(&mut self, ui: &mut egui::Ui, ctx: &egui::Context) {
        if let Some(name) = self.pending_genome.take() {
            self.use_genome(&name, ctx);
        }
        self.open_pending_tracks(ctx);
        self.poll_downloads(ctx);
        // One row: where you are (back, forward, the locus box, zoom), then four menus, chips for
        // settings changed from their defaults (× resets one), then split, Evidence and commands.
        // The locus box is always wide enough for the whole locus and takes the room left over;
        // when the row is short the chips fold into one menu, then Split, Evidence and Commands
        // into More. A pane too narrow even for that wraps.
        let chips = self.chips();
        let fit = self.toolbar_fit(ui, &chips.iter().map(|(l, _)| l.as_str()).collect::<Vec<_>>());
        ui.horizontal_wrapped(|ui| {
            if ui.add_enabled(!self.back.is_empty(), egui::Button::new("<")).on_hover_text("Back (Alt+Left)").clicked() {
                self.go_back();
            }
            if ui.add_enabled(!self.forward.is_empty(), egui::Button::new(">")).on_hover_text("Forward (Alt+Right)").clicked() {
                self.go_forward();
            }
            let locus_w = fit.locus_w;
            let edit = ui.add(
                egui::TextEdit::singleline(&mut self.locus_text).id_salt("locus")
                    .desired_width(locus_w)
                    .hint_text("chr20:1,895,892 or BRCA1")
                    .font(egui::TextStyle::Monospace),
            );
            let edit = if fit.go { edit } else { edit.on_hover_text("Enter to go there") };
            if std::mem::take(&mut self.focus_locus) {
                edit.request_focus();
                if let Some(mut state) = egui::TextEdit::load_state(ui.ctx(), edit.id) {
                    state.cursor.set_char_range(Some(egui::text::CCursorRange::two(
                        egui::text::CCursor::new(0),
                        egui::text::CCursor::new(self.locus_text.chars().count()),
                    )));
                    state.store(ui.ctx(), edit.id);
                }
            }
            let submitted = edit.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
            if (fit.go && ui.button("Go").clicked()) || submitted {
                self.pending_goto = Some(self.locus_text.clone());
            }
            if !edit.has_focus()
                && self.pending_goto.is_none()
                && let Some(text) = self.locus_now()
            {
                self.locus_text = text;
            }

            if ui.button("−").on_hover_text("Zoom out (−)").clicked() {
                self.zoom_at(0.5, self.width / 2.0);
            }
            if ui.button("+").on_hover_text("Zoom in (+)").clicked() {
                self.zoom_at(2.0, self.width / 2.0);
            }
            ui.separator();

            ui.menu_button("Open", |ui| {
                for (label, action) in [
                    ("Alignments (BAM, CRAM)…", actions::Action::AddBam),
                    ("Variants (VCF, BCF)…", actions::Action::AddVcf),
                    ("Genes (GTF, GFF3, BED)…", actions::Action::AddGenes),
                ] {
                    if ui.button(label).clicked() {
                        ui.close();
                        self.run(action, ctx);
                    }
                }
                self.reference_menu(ui, ctx);
                self.tracks_menu(ui, ctx);
                ui.separator();
                if ui.button("Save session as…").clicked() {
                    self.save_session_as();
                    ui.close();
                }
                if ui.button("Open session…").clicked() {
                    ui.close();
                    self.run(actions::Action::OpenSession, ctx);
                }
                if ui.button("Export view as PNG…").clicked() {
                    self.run(actions::Action::ExportPng, ctx);
                    ui.close();
                }
                for a in [actions::Action::CopyImage, actions::Action::CopyLink, actions::Action::CopyCommand] {
                    if ui.button(a.name()).clicked() {
                        self.run(a, ctx);
                        ui.close();
                    }
                }
                ui.menu_button("Recent files", |ui| {
                    if self.recent.is_empty() {
                        ui.label("Nothing yet");
                    }
                    for p in self.recent.clone() {
                        if ui.button(crate::files::display_name(&p)).on_hover_text(p.display().to_string()).clicked() {
                            self.open_paths(vec![p], ctx);
                            ui.close();
                        }
                    }
                });
                if let Some(last) = session::last_session_path().filter(|p| p.exists())
                    && ui.button("Restore last session").clicked()
                {
                    let _ = self.open_session(&last, ctx);
                    ui.close();
                }
                ui.separator();
                if ui.button(actions::Action::OpenDemo.name()).clicked() {
                    self.run(actions::Action::OpenDemo, ctx);
                    ui.close();
                }
            });

            let mut style = self.style;
            let here = match self.highlight {
                Some((c, p)) if c == self.contig => p,
                _ => self.center.max(0.0) as u64,
            };
            // The tag being typed, lent to the colour and group menus while Reads is open.
            let mut tag_buf = std::mem::take(&mut self.tag_input);
            let tag_input = &mut tag_buf;
            egui::containers::menu::MenuButton::new("Reads").config(keep_open()).ui(ui, |ui| {
                ui.menu_button(format!("Colour: {}", crate::session::color_text(style.color).replace("tag:", "tag ")), |ui| color_options(ui, &mut style.color, tag_input));
                ui.menu_button(format!("Group: {}", group_label(style.group)), |ui| group_options(ui, &mut style.group, here, tag_input));
                let sort_label = style.sort.map_or_else(|| "Sort".to_owned(), |(p, by)| format!("Sort: {} @ {}", sort_name(by), commas(p + 1)));
                ui.menu_button(sort_label, |ui| {
                    ui.label(format!("At {}", commas(here + 1)));
                    sort_options(ui, &mut style.sort, here);
                });
                ui.menu_button("Filters", |ui| {
                    let f = &mut style.filter;
                    ui.checkbox(&mut f.hide_duplicates, "Hide duplicates");
                    ui.checkbox(&mut f.hide_secondary, "Hide secondary");
                    ui.checkbox(&mut f.hide_supplementary, "Hide supplementary");
                    ui.checkbox(&mut f.hide_qc_fail, "Hide QC-fail");
                    ui.separator();
                    mapq_options(ui, &mut f.min_mapq);
                    ui.separator();
                    if ui.button("Reset").clicked() {
                        *f = tracks::ReadFilter::default();
                    }
                });
                ui.separator();
                ui.checkbox(&mut style.squished, "Squish rows");
                ui.checkbox(&mut style.soft_clips, "Soft clips");
                ui.checkbox(&mut style.link_splits, "Link split reads").on_hover_text("Put each split read's parts (primary and supplementary) on one row, joined, in its primary's group");
                ui.checkbox(&mut style.junctions, "Junctions").on_hover_text("Splice junctions (RNA-seq) as arcs over the coverage, with their read counts: red + strand, blue −, grey unknown");
                ui.separator();
                ui.menu_button("Coverage", |ui| {
                    ui.radio_value(&mut style.coverage_scale, tracks::CoverageScale::View, "Autoscale to view");
                    ui.radio_value(&mut style.coverage_scale, tracks::CoverageScale::Loaded, "Scale to loaded data");
                    ui.radio_value(&mut style.coverage_scale, tracks::CoverageScale::Log, "Log scale");
                    ui.separator();
                    ui.checkbox(&mut style.shared_scale, "Same scale for all tracks");
                });
                self.read_window_ui(ui);
            });
            self.tag_input = tag_buf;
            if !self.variants.is_empty() {
                ui.menu_button(if self.review.is_some() { "Variants (reviewing)" } else { "Variants" }, |ui| {
                    if ui.button("< Previous variant").on_hover_text("Previous variant (Ctrl+B)").clicked() {
                        self.seek_variant(false, None);
                    }
                    if ui.button("Next variant >").on_hover_text("Next variant (Ctrl+F)").clicked() {
                        self.seek_variant(true, None);
                    }
                    let label = if self.review.is_some() { "Reviewing" } else { "Review…" };
                    ui.menu_button(label, |ui| {
                        let names: Vec<String> = self.variants.iter().map(|v| v.file_name()).collect();
                        for (i, name) in names.iter().enumerate() {
                            if ui.button(format!("Review {name}")).clicked() {
                                self.start_review(i);
                                ui.close();
                            }
                        }
                        if self.review.is_some() && ui.button("Stop reviewing").clicked() {
                            self.review = None;
                            ui.close();
                        }
                    });
                            });
            }
            egui::containers::menu::MenuButton::new("View").config(keep_open()).ui(ui, |ui| {
                if !self.genes.is_empty() {
                    ui.checkbox(&mut self.all_transcripts, "All transcripts");
                }
                if ui.checkbox(&mut style.dark, "Dark").changed() {
                    ctx.set_theme(if style.dark { egui::Theme::Dark } else { egui::Theme::Light });
                }
                let n = self.overlays.borrow().len();
                if n > 0 && ui.button(format!("Clear overlays ({n})")).on_hover_text("Remove the regions an agent marked on the view").clicked() {
                    self.overlays.borrow_mut().clear();
                }
                use actions::Action as A;
                // In groups: find and compare, bookmarks and review, export, files, help.
                let groups: [&[A]; 5] = [
                    &[A::FindRead, A::SeqSearch, A::Filter, A::Compare],
                    &[A::Bookmarks, A::AddBookmark, A::QueueOpen, A::QueuePaste],
                    &[A::Report, A::CopyFasta, A::CopySam],
                    &[A::ReloadAll, A::Assembly],
                    &[A::Shortcuts, A::About],
                ];
                for (g, actions) in groups.iter().enumerate() {
                    ui.separator();
                    for &a in *actions {
                        let label = match a.keys() {
                            Some(k) => format!("{}    {}", a.menu_name(), actions::keys_label(k)),
                            None => a.menu_name(),
                        };
                        if ui.button(label).clicked() {
                            self.run(a, ctx);
                            ui.close();
                        }
                    }
                    if g == 3 {
                        ui.checkbox(&mut self.auto_reload, "Reload files when they change");
                    }
                }
            });
            self.style = style;

            // Changed settings, each with × to reset it; folded into one menu when they would push
            // Split, Evidence and Commands onto a second row.
            if !chips.is_empty() && fit.chips != ChipFit::Each {
                let text = if fit.chips == ChipFit::Short { format!("{} changed", chips.len()) } else { format!("{} settings changed", chips.len()) };
                let label = egui::RichText::new(text).color(ui.visuals().strong_text_color());
                egui::containers::menu::MenuButton::from_button(
                    egui::Button::new(label).fill(ui.visuals().selection.bg_fill.gamma_multiply(0.35)).corner_radius(10.0),
                )
                .config(keep_open())
                .ui(ui, |ui| {
                    ui.label(egui::RichText::new("Changed from the defaults (click to reset)").weak());
                    for (label, reset) in &chips {
                        if ui.button(format!("{label}  ×")).clicked() {
                            reset(self);
                        }
                    }
                })
                .0
                .on_hover_text(chips.iter().map(|c| c.0.as_str()).collect::<Vec<_>>().join("
"));
            } else {
                for (label, reset) in chips {
                    // Full-size text in the theme's strong text colour on a light accent fill: the
                    // small size blurred at 1x scaling.
                    let chip = egui::Button::new(egui::RichText::new(format!("{label}  ×")).color(ui.visuals().strong_text_color()))
                        .fill(ui.visuals().selection.bg_fill.gamma_multiply(0.35))
                        .stroke(egui::Stroke::new(1.0, ui.visuals().selection.bg_fill.gamma_multiply(0.8)))
                        .corner_radius(10.0);
                    if ui.add(chip).on_hover_text("Reset to the default").clicked() {
                        reset(self);
                    }
                }
            }

            ui.separator();
            let palette_key = actions::keys_label(actions::Action::Palette.keys().expect("has keys"));
            if fit.more {
                // A narrow window: these three in one menu, keeping the toolbar on one row.
                egui::containers::menu::MenuButton::new("More").config(keep_open()).ui(ui, |ui| {
                    if ui.button("Split").on_hover_text("Open another pane beside this one (same files, its own position)").clicked() {
                        ui.close();
                        self.run(actions::Action::Split, ctx);
                    }
                    ui.checkbox(&mut self.show_evidence, "Evidence").on_hover_text("Exact counts at the outlined base (or the centre)");
                    if ui.button(format!("Commands    {palette_key}")).clicked() {
                        ui.close();
                        self.run(actions::Action::Palette, ctx);
                    }
                });
            } else {
                if ui.button("Split").on_hover_text("Open another pane beside this one (same files, its own position)").clicked() {
                    self.run(actions::Action::Split, ctx);
                }
                ui.toggle_value(&mut self.show_evidence, "Evidence").on_hover_text("Exact counts at the outlined base (or the centre)");
                if ui.button("Commands").on_hover_text(format!("Command palette ({palette_key}): genes, places, every action")).clicked() {
                    self.run(actions::Action::Palette, ctx);
                }
            }
        });
    }

    /// The locus the view shows, as the locus box writes it.
    fn locus_now(&self) -> Option<String> {
        let c = self.contigs.get(self.contig)?;
        let (vs, ve) = self.visible();
        Some(format!("{}:{}-{}", c.name, commas(vs.max(0.0) as u64 + 1), commas(ve.min(self.contig_len()) as u64)))
    }

    /// How the toolbar fits the row, from the widths of its parts as `ui` draws them.
    fn toolbar_fit(&self, ui: &egui::Ui, chips: &[&str]) -> ToolbarFit {
        let gap = ui.spacing().item_spacing.x;
        let pad = 2.0 * ui.spacing().button_padding.x;
        let font = egui::TextStyle::Button.resolve(ui.style());
        let mono = egui::TextStyle::Monospace.resolve(ui.style());
        let text_w = |t: String, font: &FontId| ui.fonts_mut(|f| f.layout_no_wrap(t, font.clone(), Color32::WHITE).size().x);
        let button = |t: &str| text_w(t.to_owned(), &font) + pad + gap;
        let separator = 6.0 + gap;
        let variants = if self.variants.is_empty() {
            0.0
        } else {
            button(if self.review.is_some() { "Variants (reviewing)" } else { "Variants" })
        };
        let fixed = ["<", ">", "−", "+", "Open", "Reads", "View"].iter().map(|t| button(t)).sum::<f32>() + variants + separator + gap;
        // The longer of what's typed and the view's locus, with the box's margins and cursor.
        let longest = self.locus_now().filter(|l| l.len() > self.locus_text.len()).unwrap_or_else(|| self.locus_text.clone());
        let locus_need = (text_w(longest, &mono) + 24.0).max(160.0);
        let n = chips.len();
        let chips_w = if n == 0 {
            [0.0; 3]
        } else {
            [
                chips.iter().map(|l| text_w(format!("{l}  ×"), &font) + 24.0 + gap).sum(),
                button(&format!("{n} settings changed")),
                button(&format!("{n} changed")),
            ]
        };
        let tail = [separator + button("Split") + button("Evidence") + button("Commands"), separator + button("More")];
        fit_toolbar(ui.available_width(), fixed, button("Go"), locus_need, chips_w, tail)
    }

    /// Display settings changed from their defaults, as short labels with a reset each.
    fn chips(&self) -> Vec<(String, fn(&mut Self))> {
        let d = Style::default();
        let st = self.style;
        let mut v: Vec<(String, fn(&mut Self))> = Vec::new();
        if st.color != d.color {
            v.push((format!("colour: {}", crate::session::color_text(st.color).replace("tag:", "tag ")), |a| a.style.color = Style::default().color));
        }
        if st.group != d.group {
            v.push((format!("group: {}", group_label(st.group)), |a| a.style.group = Style::default().group));
        }
        if let Some((_, by)) = st.sort {
            v.push((format!("sort: {}", sort_name(by)), |a| a.style.sort = None));
        }
        if st.filter.min_mapq != d.filter.min_mapq {
            v.push((format!("min MAPQ {}", st.filter.min_mapq), |a| a.style.filter.min_mapq = Style::default().filter.min_mapq));
        }
        let flags = |f: &tracks::ReadFilter| (f.hide_duplicates, f.hide_secondary, f.hide_supplementary, f.hide_qc_fail, f.expr);
        if flags(&st.filter) != flags(&d.filter) {
            v.push(("filtered".into(), |a| {
                let mapq = a.style.filter.min_mapq;
                a.style.filter = Style::default().filter;
                a.style.filter.min_mapq = mapq;
            }));
        }
        if st.squished {
            v.push(("squished".into(), |a| a.style.squished = false));
        }
        if st.soft_clips {
            v.push(("soft clips".into(), |a| a.style.soft_clips = false));
        }
        if st.link_splits {
            v.push(("linked splits".into(), |a| a.style.link_splits = false));
        }
        if st.junctions {
            v.push(("junctions".into(), |a| a.style.junctions = false));
        }
        v
    }

    fn status(&mut self, ui: &mut egui::Ui, frame: &eframe::Frame) {
        self.status_bar(ui, frame.info().cpu_usage);
    }

    /// The status bar: downloads, the assembly, a notice, then the error and the warning (each
    /// with its × first, cut to fit with the full text on hover, click to copy), the counts, and
    /// the frame time pinned at the right.
    fn status_bar(&mut self, ui: &mut egui::Ui, cpu: Option<f32>) {
        ui.horizontal(|ui| {
            self.download_status(ui);
            self.assembly_status(ui);
            self.notice_ui(ui);
            let dark = ui.visuals().dark_mode;
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if let Some(cpu) = cpu {
                    ui.label(format!("frame {:.2} ms", cpu * 1000.0));
                }
                ui.with_layout(egui::Layout::left_to_right(egui::Align::Center), |ui| {
                    let warning = self.warning();
                    // The messages share the row, leaving the counts some room.
                    let shown = usize::from(self.error.is_some()) + usize::from(warning.is_some());
                    let share = ((ui.available_width() - 160.0).max(120.0) / shown.max(1) as f32).max(80.0);
                    if let Some(e) = self.error.clone() {
                        match status_message(ui, &e, error_rgb(dark), share) {
                            Some(MessageClick::Dismiss) => self.error = None,
                            Some(MessageClick::Copy) => self.notify("Copied the error"),
                            None => {}
                        }
                    }
                    if let Some(w) = warning {
                        match status_message(ui, &w, warning_rgb(dark), share) {
                            Some(MessageClick::Dismiss) => {
                                let shown = self.tracks.iter().filter(|t| t.loader.warning.is_some()).map(|t| t.loader.bam_path.clone());
                                self.dismissed_warnings.extend(shown.collect::<Vec<_>>());
                            }
                            Some(MessageClick::Copy) => self.notify("Copied the warning"),
                            None => {}
                        }
                    }
                    if self.tracks.is_empty() && self.genes.is_empty() && self.variants.is_empty() && self.signals.is_empty() {
                        ui.add(egui::Label::new("Add a BAM or a gene annotation to start (or pass them on the command line).").truncate());
                    } else {
                        ui.add(egui::Label::new(self.counts_text()).truncate());
                        if self.tracks.iter().any(Track::is_loading) {
                            ui.separator();
                            ui.spinner();
                        }
                    }
                });
            });
        });
    }

    /// The status bar's counts: the tracks on screen by kind, then the reads and GPU work.
    fn counts_text(&self) -> String {
        let mut reads = 0;
        let mut instances = 0;
        let mut slowest = 0.0f64;
        for t in &self.tracks {
            if let Some(c) = &t.current {
                reads += c.view.reads.as_ref().map_or(0, |(r, _)| r.len());
                instances += c.instances.instances.len();
                slowest = slowest.max(c.elapsed.as_secs_f64() * 1000.0);
            }
        }
        instances += self.variants.iter().filter_map(|v| v.current.as_ref()).map(|c| c.instances.instances.len()).sum::<usize>();
        instances += self.genes.iter().filter_map(|g| g.built.as_ref()).map(|b| b.instances.instances.len()).sum::<usize>();
        let mut parts = vec![tracks_text(self.tracks.len(), self.variants.len(), self.genes.len(), self.signals.len())];
        if !self.tracks.is_empty() {
            parts.push(format!("{} reads", commas(reads as u64)));
        }
        parts.push(format!("{} GPU instances", commas(instances as u64)));
        if !self.tracks.is_empty() {
            parts.push(format!("slowest load {slowest:.0} ms"));
        }
        parts.join(" · ")
    }

    fn canvas(&mut self, ui: &mut egui::Ui) {
        let (rect, response) = ui.allocate_exact_size(ui.available_size(), Sense::click_and_drag());
        self.width = rect.width().max(1.0);
        self.canvas_rect = rect;
        let ppp = ui.ctx().pixels_per_point();
        self.poll_genes();
        self.poll_signals();

        // Check for a genome before taking the request: an annotation may still be loading.
        if !self.contigs.is_empty()
            && let Some(text) = self.pending_goto.take()
        {
            self.goto(&text);
        }
        if !self.contigs.is_empty()
            && self.pending_goto.is_none()
            && let Some(text) = self.pending_highlight.take()
            && let Ok(r) = ggv_io::resolve_locus_str(&self.contigs, &text)
            && let Some(c) = self.contigs.iter().position(|c| c.name == r.contig)
        {
            self.highlight = Some((c, r.start));
        }
        if self.contigs.is_empty() {
            self.welcome(ui, rect);
            return;
        }

        let layout = self.layout(rect);
        self.handle_input(ui, &response, rect, &layout);
        if self.contigs.is_empty() {
            return; // the last track was closed
        }
        let layout = self.layout(rect);
        self.clamp_view();
        for (t, section) in self.tracks.iter_mut().zip(&layout.sections) {
            t.scroll_y = t.scroll_y.clamp(0.0, t.max_scroll(self.style, section.height()));
            t.view_h = section.height();
        }
        for (g, section) in self.genes.iter_mut().zip(&layout.genes) {
            g.scroll_y = g.scroll_y.clamp(0.0, g.max_scroll(section.height()));
        }
        for (v, section) in self.variants.iter_mut().zip(&layout.variants) {
            v.scroll_y = v.scroll_y.clamp(0.0, v.max_scroll(section.height()));
        }
        self.update_tracks(ppp);

        let painter = ui.painter_at(rect);
        let (vs, ve) = self.visible();
        let mut retire = std::mem::take(&mut self.retired);

        // Shared reference row.
        if let Some((view, set, _)) = &self.reference
            && view.region.contig == self.contigs[self.contig].name
        {
            let ref_rect = Rect::from_min_size(Pos2::new(rect.left(), rect.top() + REF_TOP), Vec2::new(rect.width(), REF_H));
            let uniforms = self.uniforms(ref_rect, vs, view.region.start, 0.0, ppp);
            gpu::paint(&painter.with_clip_rect(ref_rect), ref_rect, self.reference_slot, uniforms, set.clone(), std::mem::take(&mut retire));
        }

        // One instanced draw per track, each in its own section.
        let scales = self.coverage_scales(vs, ve);
        for ((t, section), &cov_scale) in self.tracks.iter().zip(&layout.sections).zip(&scales) {
            let Some(cur) = &t.current else { continue };
            let section = section.intersect(rect);
            if section.height() <= 0.0 {
                continue;
            }
            let mut uniforms = self.uniforms(section, vs, cur.view.region.start, t.scroll_y, ppp);
            uniforms.cov_scale = cov_scale;
            gpu::paint(&painter.with_clip_rect(section), section, t.id, uniforms, cur.instances.clone(), std::mem::take(&mut retire));
        }
        for (v, section) in self.variants.iter().zip(&layout.variants) {
            let Some(cur) = &v.current else { continue };
            if cur.region.contig != self.contigs[self.contig].name || ve - vs > variants::VARIANT_WINDOW as f64 {
                continue;
            }
            let section = section.intersect(rect);
            if section.height() <= 0.0 {
                continue;
            }
            let mut uniforms = self.uniforms(section, vs, cur.region.start, v.scroll_y, ppp);
            uniforms.scroll_top = variants::gt_top(cur.site_rows);
            gpu::paint(&painter.with_clip_rect(section), section, v.id, uniforms, cur.instances.clone(), std::mem::take(&mut retire));
        }
        for (g, section) in self.genes.iter().zip(&layout.genes) {
            let Some(built) = &g.built else { continue };
            let section = section.intersect(rect);
            if section.height() <= 0.0 {
                continue;
            }
            let mut uniforms = self.uniforms(section, vs, built.start, g.scroll_y, ppp);
            uniforms.scroll_top = genes::ROWS_TOP;
            gpu::paint(&painter.with_clip_rect(section), section, g.id, uniforms, built.instances.clone(), std::mem::take(&mut retire));
        }
        self.retired.extend(retire); // nothing painted this frame: free them next time

        self.draw_shared_overlays(&painter, rect, vs, ve);
        for (i, section) in layout.sections.iter().enumerate() {
            if !section.is_positive() {
                continue; // scrolled out of the track stack
            }
            self.draw_track_overlays(&painter.with_clip_rect(section.intersect(rect)), i, *section, vs, ve, scales.get(i).copied().unwrap_or(1.0));
        }
        for (i, section) in layout.variants.iter().enumerate() {
            if !section.is_positive() {
                continue; // scrolled out of the track stack
            }
            self.draw_variant_overlays(&painter.with_clip_rect(section.intersect(rect)), i, *section, vs, ve);
        }
        for (i, section) in layout.genes.iter().enumerate() {
            if !section.is_positive() {
                continue; // scrolled out of the track stack
            }
            self.draw_gene_overlays(&painter.with_clip_rect(section.intersect(rect)), i, *section, vs, ve);
        }
        for (i, section) in layout.signals.iter().enumerate() {
            if !section.is_positive() {
                continue; // scrolled out of the track stack
            }
            self.draw_signal_overlays(&painter.with_clip_rect(section.intersect(rect)), i, *section, vs, ve);
        }
        self.draw_row_bars(&painter, &layout, response.hover_pos());
        // On top of every track's header and labels.
        self.draw_tool_overlays(&painter, rect, &layout, vs, ve);
        self.draw_affordances(&painter, &layout, response.hover_pos());
        self.draw_arranging(&painter, &layout);
        self.draw_stack_bar(&painter, &layout);

        if response.secondary_clicked()
            && let Some(p) = response.interact_pointer_pos()
        {
            let base = (vs + f64::from(p.x - rect.left()) / self.ppb).floor().max(0.0) as u64;
            let read = match layout.hit(p.y) {
                Some(Hit::Track(i)) => self.read_at(i, layout.sections[i], p, vs).cloned(),
                _ => None,
            };
            self.menu_at = Some((self.contig, base, read, layout.hit(p.y)));
        }
        response.context_menu(|ui| self.context_menu(ui));

        // A track header shows its full name (the title may be shortened).
        let header = response.hover_pos().and_then(|pos| {
            let hit = layout.hit(pos.y)?;
            let section = match hit {
                Hit::Track(i) => layout.sections[i],
                Hit::Genes(i) => layout.genes[i],
                Hit::Variants(i) => layout.variants[i],
                Hit::Signal(i) => layout.signals[i],
            };
            (pos.y < section.top() + HEADER_H && !Layout::close_button(section).contains(pos)).then(|| match hit {
                Hit::Track(i) => self.tracks[i].full_title(),
                Hit::Genes(i) => self.genes[i].path.display().to_string(),
                Hit::Variants(i) => self.variants[i].path.display().to_string(),
                Hit::Signal(i) => self.signals[i].path.display().to_string(),
            })
        });
        let on_bar = response.hover_pos().and_then(|p| self.rows_under(&layout, p).filter(|(_, rs)| rs.bar_hit().contains(p)));
        if let Some(name) = header {
            show_tooltip(&response, vec![name]);
        } else if let Some((_, rs)) = on_bar {
            let (a, b) = rs.shown();
            show_tooltip(
                &response,
                vec![
                    format!("{} {}–{} of {}", rs.noun, commas(a as u64), commas(b as u64), commas(rs.rows as u64)),
                    "Wheel or drag here to scroll them · Ctrl+wheel zooms".to_owned(),
                ],
            );
        } else if let Some(pos) = response.hover_pos() {
            match layout.hit(pos.y) {
                Some(Hit::Track(i)) => self.hover(&response, &painter, i, layout.sections[i], pos, vs),
                Some(Hit::Genes(i)) => self.hover_genes(&response, &painter, i, layout.genes[i], pos, vs),
                Some(Hit::Variants(i)) => self.hover_variants(&response, &painter, i, layout.variants[i], pos, vs),
                Some(Hit::Signal(i)) => self.hover_signal(&response, i, layout.signals[i], pos, vs),
                None if pos.y < rect.top() + REF_TOP => self.hover_ruler(&response, rect, pos, vs),
                None => {}
            }
        }
    }

    /// A track header's title, shortened with an ellipsis to end before the close button (a long
    /// name ran underneath it). The full name is on the header's hover.
    /// It ends before `right` (what [`Self::header_rows`] leaves).
    fn header_title(&self, painter: &egui::Painter, section: Rect, text: String, right: f32) -> Rect {
        let font = FontId::proportional(11.5);
        let color = self.pal().text;
        let room = (right - (section.left() + 4.0)).max(0.0);
        let mut galley = painter.layout_no_wrap(text.clone(), font.clone(), color);
        if galley.size().x > room {
            let mut chars: Vec<char> = text.chars().collect();
            while galley.size().x > room && !chars.is_empty() {
                chars.pop();
                let s: String = chars.iter().collect();
                galley = painter.layout_no_wrap(format!("{}\u{2026}", s.trim_end()), font.clone(), color);
            }
        }
        let at = Pos2::new(section.left() + 4.0, section.top() + HEADER_H / 2.0 - galley.size().y / 2.0);
        let rect = Rect::from_min_size(at, galley.size());
        painter.galley(at, galley, color);
        rect
    }

    /// When a track's rows don't all fit: "rows 1–41 of 377" at the right of its header, before
    /// the grip. Returns where the title must end.
    fn header_rows(&self, painter: &egui::Painter, section: Rect, rows: Option<RowScroll>) -> f32 {
        let end = Layout::grip(section).left() - 6.0;
        let Some(rs) = rows else { return end };
        let (a, b) = rs.shown();
        let text = format!("{} {}–{} of {}", rs.noun, commas(a as u64), commas(b as u64), commas(rs.rows as u64));
        let galley = painter.layout_no_wrap(text, FontId::proportional(10.5), self.pal().muted);
        let w = galley.size().x;
        // Only with room left for a title.
        if end - w < section.left() + 120.0 {
            return end;
        }
        painter.galley(Pos2::new(end - w, section.top() + HEADER_H / 2.0 - galley.size().y / 2.0), galley, self.pal().muted);
        end - w - 12.0
    }

    /// The rows of a track that scroll, when they don't all fit: reads below the coverage, gene
    /// rows, or sample rows (genotypes).
    fn row_scroll(&self, hit: Hit, section: Rect) -> Option<RowScroll> {
        let below = |y: f32| Rect::from_min_max(Pos2::new(section.left(), section.top() + y), section.max);
        let rs = match hit {
            Hit::Track(i) => {
                let t = self.tracks.get(i)?;
                let cur = t.current.as_ref()?;
                let (reads, _) = cur.view.reads.as_ref()?;
                // No rows to scroll when none of the loaded reads is on screen.
                let (vs, ve) = self.visible();
                if !reads_in_view(reads, &cur.packed, vs, ve) {
                    return None;
                }
                let rl = tracks::row_layout(&cur.packed, cur.style);
                let area = below(READS_TOP);
                RowScroll { area, scroll: t.scroll_y, max: (rl.content_height() - area.height() + 12.0).max(0.0), rows: rl.rows, pitch: rl.pitch, noun: "rows" }
            }
            Hit::Genes(i) => {
                let g = self.genes.get(i)?;
                let rows = g.built.as_ref()?.rows.len();
                RowScroll { area: below(genes::ROWS_TOP), scroll: g.scroll_y, max: g.max_scroll(section.height()), rows, pitch: genes::PITCH, noun: "rows" }
            }
            Hit::Variants(i) => {
                let v = self.variants.get(i)?;
                let area = below(v.gt_top());
                RowScroll { area, scroll: v.scroll_y, max: v.max_scroll(section.height()), rows: v.samples.len(), pitch: variants::GT_PITCH, noun: "samples" }
            }
            Hit::Signal(_) => return None,
        };
        (rs.max > 0.5 && rs.area.height() >= 12.0 && rs.rows > 0).then_some(rs)
    }

    /// Scrolls the rows of a track by `dy` points (clamped when the frame lays them out).
    fn scroll_rows(&mut self, hit: Hit, dy: f32) {
        match hit {
            Hit::Track(i) => {
                if let Some(t) = self.tracks.get_mut(i) {
                    t.scroll_y = (t.scroll_y + dy).max(0.0);
                }
            }
            Hit::Genes(i) => {
                if let Some(g) = self.genes.get_mut(i) {
                    g.scroll_y = (g.scroll_y + dy).max(0.0);
                }
            }
            Hit::Variants(i) => {
                if let Some(v) = self.variants.get_mut(i) {
                    v.scroll_y = (v.scroll_y + dy).max(0.0);
                }
            }
            Hit::Signal(_) => {}
        }
    }

    /// A track's rows scrollbar: drag the thumb, or click the trough to page towards the
    /// pointer. True when it took the pointer (so the view doesn't pan or select as well).
    fn row_bar_input(&mut self, ui: &egui::Ui, response: &egui::Response, layout: &Layout) -> bool {
        let on_bar = |s: &Self, p: Pos2| s.rows_under(layout, p).filter(|(_, rs)| rs.bar_hit().contains(p));
        if response.drag_started()
            && let Some(p) = ui.input(|i| i.pointer.press_origin()).or(response.interact_pointer_pos())
            && let Some((hit, _)) = on_bar(self, p)
        {
            self.row_bar_drag = Some(hit);
        }
        if let Some(hit) = self.row_bar_drag {
            if let Some(rs) = self.row_scroll(hit, layout.rect_of(hit)) {
                self.scroll_rows(hit, response.drag_delta().y * rs.ratio());
            }
            if response.drag_stopped() || !response.dragged() {
                self.row_bar_drag = None;
            }
            return true;
        }
        if response.clicked()
            && let Some(p) = response.interact_pointer_pos()
            && let Some((hit, rs)) = on_bar(self, p)
        {
            let thumb = rs.thumb();
            if !thumb.y_range().contains(p.y) {
                let page = rs.area.height() - rs.pitch;
                self.scroll_rows(hit, if p.y < thumb.top() { -page } else { page });
            }
            return true;
        }
        false
    }

    /// The track whose overflowing rows are under `p`: there the wheel scrolls them (and
    /// Ctrl+wheel zooms), elsewhere the wheel zooms.
    fn rows_under(&self, layout: &Layout, p: Pos2) -> Option<(Hit, RowScroll)> {
        let hit = layout.hit(p.y)?;
        let rs = self.row_scroll(hit, layout.rect_of(hit))?;
        rs.area.contains(p).then_some((hit, rs))
    }

    /// A thin scrollbar at the right of every track whose rows don't all fit.
    fn draw_row_bars(&self, painter: &egui::Painter, layout: &Layout, hover: Option<Pos2>) {
        let pal = self.pal();
        let tracks = layout.sections.iter().enumerate().map(|(i, r)| (Hit::Track(i), *r));
        let genes = layout.genes.iter().enumerate().map(|(i, r)| (Hit::Genes(i), *r));
        let variants = layout.variants.iter().enumerate().map(|(i, r)| (Hit::Variants(i), *r));
        for (hit, r) in tracks.chain(genes).chain(variants).filter(|(_, r)| r.is_positive()) {
            let Some(rs) = self.row_scroll(hit, r) else { continue };
            let hot = self.row_bar_drag == Some(hit) || hover.is_some_and(|p| rs.bar_hit().contains(p));
            let bar = rs.bar();
            // Backed by the page colour, so the trough shows over reads as well as on white.
            painter.rect_filled(bar.expand(1.0), 3.0, pal.bg.gamma_multiply(0.85));
            painter.rect_filled(bar, 3.0, if hot { pal.rule_strong.gamma_multiply(0.6) } else { pal.rule });
            painter.rect_filled(rs.thumb(), 3.0, if hot { pal.icon } else { pal.muted });
        }
    }

    /// Over the ruler: every overlay label and bookmark note at the base under the pointer, so
    /// those the ruler has no room for can still be read.
    fn hover_ruler(&self, response: &egui::Response, rect: Rect, pos: Pos2, vs: f64) {
        let Some(c) = self.contigs.get(self.contig) else { return };
        let at = vs + f64::from(pos.x - rect.left()) / self.ppb;
        // Within 3 px, so a one-base mark can be hovered.
        let slack = 3.0 / self.ppb;
        let near = |s: u64, e: u64| (s as f64) - slack <= at && at < (e as f64) + slack;
        let mut lines: Vec<String> = Vec::new();
        for o in self.overlays.borrow().iter().filter(|o| o.contig == c.name && near(o.start, o.end) && !o.label.is_empty()) {
            lines.push(format!("{}  ({}:{}-{})", o.label, o.contig, commas(o.start + 1), commas(o.end)));
        }
        for b in self.bookmarks.borrow().iter().filter(|b| b.contig == c.name && near(b.start, b.end)) {
            let place = format!("{}:{}-{}", b.contig, commas(b.start + 1), commas(b.end));
            lines.push(if b.note.is_empty() { format!("bookmark {place}") } else { format!("bookmark: {}  ({place})", b.note) });
        }
        if !lines.is_empty() {
            show_tooltip(response, lines);
        }
    }

    /// The read under `pos` in alignment track `i`.
    fn read_at(&self, i: usize, section: Rect, pos: Pos2, vs: f64) -> Option<&AlignedRead> {
        let t = self.tracks.get(i)?;
        let cur = t.current.as_ref()?;
        let (reads, _) = cur.view.reads.as_ref()?;
        let y = pos.y - section.top() - READS_TOP + t.scroll_y;
        if y < 0.0 {
            return None;
        }
        let rl = tracks::row_layout(&cur.packed, cur.style);
        let row = cur.packed.rows.get((y / rl.pitch) as usize)?;
        let base = (vs + f64::from(pos.x - section.left()) / self.ppb).floor().max(0.0) as u64;
        let k = row.partition_point(|&i| reads[i].end <= base);
        let r = &reads[*row.get(k)?];
        (r.start <= base).then_some(r)
    }

    /// Right-click menu (IGV-style): sort at the clicked base, group, colour, copy.
    fn context_menu(&mut self, ui: &mut egui::Ui) {
        let Some((contig, base, read, hit)) = self.menu_at.clone() else {
            ui.close();
            return;
        };
        let name = self.contigs.get(contig).map_or("", |c| c.name.as_str()).to_owned();
        ui.label(egui::RichText::new(format!("{name}:{}", commas(base + 1))).strong());
        ui.separator();
        if let Some(Hit::Variants(i)) = hit
            && let Some(v) = self.variants.get_mut(i)
        {
            ui.menu_button(format!("Colour sites: {}", v.color.name()), |ui| {
                let mut pick = None;
                if ui.selectable_label(v.color == variants::SiteColor::AlleleFrequency, "Allele frequency (AF)").clicked() {
                    pick = Some(variants::SiteColor::AlleleFrequency);
                }
                ui.separator();
                ui.label(egui::RichText::new("INFO field").weak());
                egui::ScrollArea::vertical().max_height(360.0).show(ui, |ui| {
                    for f in v.info_fields.iter() {
                        let on = v.color == variants::SiteColor::Info(f.id.clone());
                        if ui.selectable_label(on, &f.id).on_hover_text(&f.description).clicked() {
                            pick = Some(variants::SiteColor::Info(f.id.clone()));
                        }
                    }
                });
                if let Some(c) = pick {
                    v.color = c;
                    ui.close();
                }
            });
            ui.separator();
        }
        let mut style = self.style;
        sort_options(ui, &mut style.sort, base);
        ui.separator();
        if ui.button("Group by base here").clicked() {
            style.group = ggv_layout::GroupBy::BaseAt(base);
            ui.close();
        }
        let tag_input = &mut self.tag_input;
        ui.menu_button(format!("Group: {}", group_label(style.group)), |ui| group_options(ui, &mut style.group, base, tag_input));
        let tag_input = &mut self.tag_input;
        ui.menu_button(format!("Colour: {}", crate::session::color_text(style.color).replace("tag:", "tag ")), |ui| color_options(ui, &mut style.color, tag_input));
        ui.checkbox(&mut style.soft_clips, "Show soft-clipped bases");
        ui.checkbox(&mut style.link_splits, "Link split reads");
        ui.checkbox(&mut style.squished, "Squish");
        ui.checkbox(&mut style.junctions, "Splice junctions");
        self.style = style;
        ui.separator();
        if let Some(read) = &read {
            if ui.button("Copy read name").clicked() {
                ui.ctx().copy_text(read.name.clone());
                ui.close();
            }
            if let Some(Hit::Track(i)) = hit {
                if ui.button("Copy read sequence").clicked() {
                    self.copy_read(&ui.ctx().clone(), i, read, true);
                    ui.close();
                }
                if ui.button("Copy read as SAM").clicked() {
                    self.copy_read(&ui.ctx().clone(), i, read, false);
                    ui.close();
                }
            }
            // Mate and other alignments of this read, through history (Alt+Left comes back).
            let track = match hit {
                Some(Hit::Track(i)) => self.tracks.get(i),
                _ => None,
            };
            if let Some(m) = read.mate
                && let Some(mate_contig) = track.and_then(|t| t.contigs().get(m.contig_id.max(0) as usize)).map(|c| c.name.clone())
                && ui.button(format!("Go to mate  {mate_contig}:{}", commas(m.pos.max(0) as u64 + 1))).clicked()
            {
                self.selected_read = Some(read.name.clone());
                self.jump_to(&mate_contig, m.pos.max(0) as u64);
                ui.close();
            }
            for sa in &read.supplementary {
                let label = format!("Go to split part  {}:{} ({}, MAPQ {})", sa.contig, commas(sa.pos + 1), if sa.reverse { "−" } else { "+" }, sa.mapq);
                if ui.button(label).clicked() {
                    self.selected_read = Some(read.name.clone());
                    self.jump_to(&sa.contig, sa.pos);
                    ui.close();
                }
            }
        }
        if ui.button("Copy position").clicked() {
            ui.ctx().copy_text(format!("{name}:{}", base + 1));
            ui.close();
        }
        let ctx = ui.ctx().clone();
        if ui.button(if self.selection.is_some() { "Copy selection as FASTA" } else { "Copy view as FASTA" }).clicked() {
            self.copy_fasta(&ctx);
            ui.close();
        }
        if self.selection.is_some() && ui.button("Copy selected reads as SAM").clicked() {
            self.copy_selection_sam(&ctx);
            ui.close();
        }
        if ui.button("Copy view image").clicked() {
            self.run(actions::Action::CopyImage, &ctx);
            ui.close();
        }
        if ui.button("Bookmark here").clicked() {
            let i = crate::bookmarks::add(
                &mut self.bookmarks.borrow_mut(),
                crate::bookmarks::Bookmark { contig: name.clone(), start: base, end: base + 1, note: String::new() },
            );
            self.show_bookmarks = true;
            self.edit_bookmark = Some(i);
            ui.close();
        }
        if let Some(hit) = hit {
            ui.separator();
            if let Hit::Track(i) = hit
                && ui.button("Filter this track…").clicked()
            {
                self.open_filter_window(Some(i));
                ui.close();
            }
            if ui.button("Move track up").clicked() {
                self.move_track(hit, true);
                ui.close();
            }
            if ui.button("Move track down").clicked() {
                self.move_track(hit, false);
                ui.close();
            }
            if ui.button("Close track").clicked() {
                self.close_hit(hit);
                self.menu_at = None;
                ui.close();
            }
        }
    }

    /// Points per coverage height unit for each alignment track: scaled to the deepest visible
    /// position (or loaded, or log), per track or shared across tracks.
    fn coverage_scales(&self, vs: f64, ve: f64) -> Vec<f32> {
        let units: Vec<f32> = self
            .tracks
            .iter()
            .map(|t| {
                t.current.as_ref().map_or(1.0, |c| {
                    tracks::max_units(&c.coverage, &c.depth_index, c.view.region.start, c.view.region.end, vs, ve, self.style.coverage_scale)
                })
            })
            .collect();
        let shared = units.iter().copied().fold(1.0f32, f32::max);
        let usable = COVERAGE_H - 12.0; // room for the label
        units.iter().map(|&u| usable / if self.style.shared_scale { shared } else { u }).collect()
    }

    fn pal(&self) -> Pal {
        Pal::new(self.style.dark)
    }

    fn uniforms(&self, rect: Rect, vs: f64, origin: u64, scroll_y: f32, ppp: f32) -> Uniforms {
        Uniforms {
            size: [rect.width(), rect.height()],
            offset_bp: split_offset(vs, origin).0,
            ppb: self.ppb as f32,
            scroll_y,
            scroll_top: READS_TOP,
            min_px: 1.0 / ppp,
            linear_out: 0,
            cov_scale: 1.0,
            cov_bottom: COVERAGE_TOP + COVERAGE_H,
            offset_frac: split_offset(vs, origin).1,
            _pad: 0.0,
        }
    }

    fn handle_input(&mut self, ui: &egui::Ui, response: &egui::Response, rect: Rect, layout: &Layout) {
        if self.stack_input(ui, response, layout) {
            return;
        }
        // Close buttons in track headers.
        if response.clicked()
            && let Some(p) = response.interact_pointer_pos()
        {
            if let Some(i) = layout.sections.iter().position(|s| Layout::close_button(*s).contains(p)) {
                self.close_track(i);
                return;
            }
            if let Some(i) = layout.genes.iter().position(|s| Layout::close_button(*s).contains(p)) {
                self.close_genes(i);
                return;
            }
            if let Some(i) = layout.variants.iter().position(|s| Layout::close_button(*s).contains(p)) {
                self.close_variants(i);
                return;
            }
        }
        // Moving and resizing tracks (headers and bottom edges).
        if self.arrange_input(ui, response, layout) {
            return;
        }
        if self.row_bar_input(ui, response, layout) {
            return;
        }

        // Drag pans: horizontally through the genome, vertically through the rows of the track
        // where the drag started.
        let base_at = |p: Pos2, s: &Self| {
            let (vs, _) = s.visible();
            (vs + f64::from(p.x - rect.left()) / s.ppb).floor().max(0.0) as u64
        };
        // Ctrl-drag measures, from what is loaded (no disk).
        let command = ui.input(|i| i.modifiers.command);
        if response.drag_started()
            && command
            && let Some(p) = response.interact_pointer_pos()
        {
            let b = base_at(p, self);
            self.measure = Some((self.contig, b, b + 1));
            self.measuring = true;
        }
        if self.measuring {
            if let (Some(p), Some((c, a, _))) = (response.interact_pointer_pos(), self.measure) {
                let b = base_at(p, self);
                self.measure = Some((c, a.min(b), a.max(b) + 1));
            }
            if response.drag_stopped() {
                self.measuring = false;
            }
            return;
        }
        // Shift-drag selects a region for statistics instead of panning.
        let shift = ui.input(|i| i.modifiers.shift);
        if response.drag_started()
            && shift
            && let Some(p) = response.interact_pointer_pos()
        {
            let b = base_at(p, self);
            self.selection = Some((self.contig, b, b + 1));
            self.selecting = true;
        }
        if self.selecting {
            if let (Some(p), Some((c, a, _))) = (response.interact_pointer_pos(), self.selection) {
                let b = base_at(p, self);
                self.selection = Some((c, a.min(b), a.max(b) + 1));
            }
            if response.drag_stopped() {
                self.selecting = false;
                if let Some((c, a, b)) = self.selection {
                    self.select_region(c, a, b);
                }
            }
            return;
        }
        // A plain drag from the ruler or reference row marks a region to zoom to instead of
        // panning. (Where the button went down: the drag is reported a few points later.)
        if response.drag_started()
            && let Some(p) = ui.input(|i| i.pointer.press_origin())
            && p.y < rect.top() + TRACKS_TOP
        {
            let b = base_at(p, self);
            self.zoom_drag = Some((b, b));
        }
        if let Some((from, _)) = self.zoom_drag {
            ui.ctx().set_cursor_icon(egui::CursorIcon::ResizeHorizontal);
            if let Some(p) = response.interact_pointer_pos() {
                let to = base_at(Pos2::new(p.x.clamp(rect.left(), rect.right()), p.y), self);
                self.zoom_drag = Some((from, to));
            }
            if response.drag_stopped() {
                let marked = self.zoom_drag.take().unwrap_or((from, from));
                self.zoom_to_marked(marked.0, marked.1);
            }
            return;
        }
        if response.clicked()
            && let Some(p) = response.interact_pointer_pos()
        {
            let (vs, _) = self.visible();
            let hit = layout.hit(p.y);
            self.selected_read = match hit {
                Some(Hit::Track(i)) => self.read_at(i, layout.sections[i], p, vs).map(|r| r.name.clone()),
                _ => None,
            };
            // Reads, coverage and the reference outline the base clicked; genes and variants have
            // their own clicks.
            let on_header = layout.order.iter().any(|(_, r)| (r.top()..r.top() + crate::tracks::HEADER_H).contains(&p.y));
            if let Some(Hit::Variants(i)) = hit
                && !on_header
                && let Some(var) = self.variant_at(i, layout.variants[i], p, vs).cloned()
            {
                if self.review_track() == Some(i) {
                    self.review_pick(&var);
                }
                self.open_variant_popup(i, var);
            }
            if matches!(hit, Some(Hit::Track(_)) | None) && !on_header && !self.contigs.is_empty() {
                let base = (base_at(p, self)).min(self.contig_len().max(1.0) as u64 - 1);
                match self.highlight {
                    // Shift-click: select from the outlined base to this one (as shift-drag does).
                    Some((c, from)) if shift && c == self.contig && from != base => {
                        self.select_region(c, from.min(base), from.max(base) + 1);
                    }
                    _ => self.highlight = if self.highlight == Some((self.contig, base)) { None } else { Some((self.contig, base)) },
                }
            }
        }
        if response.drag_started() {
            self.drag_track = response.interact_pointer_pos().and_then(|p| layout.hit(p.y));
        }
        let drag = response.drag_delta();
        if drag != Vec2::ZERO {
            self.center -= f64::from(drag.x) / self.ppb;
            match self.drag_track {
                Some(Hit::Track(i)) => {
                    if let Some(t) = self.tracks.get_mut(i) {
                        t.scroll_y -= drag.y;
                    }
                }
                Some(Hit::Genes(i)) => {
                    if let Some(g) = self.genes.get_mut(i) {
                        g.scroll_y -= drag.y;
                    }
                }
                Some(Hit::Variants(i)) => {
                    if let Some(v) = self.variants.get_mut(i) {
                        v.scroll_y -= drag.y;
                    }
                }
                Some(Hit::Signal(_)) | None => {}
            }
        }
        if response.drag_stopped() {
            self.drag_track = None;
        }
        if response.double_clicked()
            && let Some(p) = response.interact_pointer_pos()
        {
            // On a gene: fit the transcript (with a little room). Elsewhere: zoom in 2x.
            let gene = match layout.hit(p.y) {
                Some(Hit::Genes(i)) => self.gene_at(i, layout.genes[i], p).map(|(_, t)| (t.start, t.end)),
                _ => None,
            };
            match gene {
                Some((start, end)) => {
                    self.remember_place();
                    let pad = ((end - start) / 20).max(20) as f64;
                    self.center = (start + end) as f64 / 2.0;
                    self.ppb = f64::from(self.width) / ((end - start) as f64 + 2.0 * pad);
                    self.clamp_view();
                }
                None => self.zoom_at(2.0, p.x - rect.left()),
            }
        }
        let pointer = ui.input(|i| i.pointer.hover_pos());
        if response.hovered() {
            let (mut scroll, zoom) = ui.input(|i| (i.smooth_scroll_delta(), i.zoom_delta()));
            let x = pointer.map_or(self.width / 2.0, |p| p.x - rect.left());
            // Over rows that don't all fit (reads, gene or sample rows), the wheel scrolls them;
            // elsewhere (ruler, headers, coverage, tracks that fit) it zooms around the pointer,
            // as Ctrl+wheel does everywhere (egui turns it into `zoom`). Horizontal scroll
            // (trackpads, Shift+wheel) pans.
            if scroll.y != 0.0
                && let Some((hit, _)) = pointer.and_then(|p| self.rows_under(layout, p))
            {
                self.scroll_rows(hit, -scroll.y);
                scroll.y = 0.0;
            }
            let factor = f64::from(zoom) * (f64::from(scroll.y) * 0.004).exp();
            if (factor - 1.0).abs() > 1e-6 {
                self.zoom_at(factor, x);
            }
            if scroll.x != 0.0 {
                self.center -= f64::from(scroll.x) / self.ppb;
            }
        }

        self.hover_variant = pointer.and_then(|p| match layout.hit(p.y) {
            Some(Hit::Variants(i)) => Some(i),
            _ => None,
        });
        let typing = ui.ctx().memory(|m| m.focused().is_some());
        if typing || !self.is_active {
            return;
        }
        // Review queue: 1 / 2 / 3 pass / fail / unsure and move on; N / P step; U next undecided.
        if self.queue.is_some() {
            let keys = ui.input(|i| {
                [egui::Key::Num1, egui::Key::Num2, egui::Key::Num3, egui::Key::N, egui::Key::P, egui::Key::U].map(|k| !i.modifiers.any() && i.key_pressed(k))
            });
            for (k, v) in crate::queue::Verdict::ALL.into_iter().enumerate() {
                if keys[k] {
                    self.queue_decide(v);
                }
            }
            if keys[3] || keys[4] {
                self.queue_step(keys[3]);
            }
            if keys[5] {
                self.queue_next_undecided();
            }
        }
        // Review: 1 / 2 / 3 mark the outlined site real / artifact / unsure and move on; N / P step.
        if self.review.is_some() {
            let keys = ui.input(|i| {
                [egui::Key::Num1, egui::Key::Num2, egui::Key::Num3, egui::Key::N, egui::Key::P].map(|k| !i.modifiers.any() && i.key_pressed(k))
            });
            use crate::review::Label;
            for (k, label) in [Label::Real, Label::Artifact, Label::Unsure].into_iter().enumerate() {
                if keys[k] {
                    self.decide(label);
                }
            }
            if keys[3] || keys[4] {
                self.review_step(keys[3]);
            }
        }
        // Row-scrolling keys act on the track under the pointer (or the first).
        let target = pointer.and_then(|p| layout.section_at(p.y)).unwrap_or(0);
        let page = layout.sections.get(target).map_or(200.0, |s| s.height() - READS_TOP);
        let step = f64::from(self.width) * 0.2 / self.ppb;
        let (left, right, up, down, page_up, page_down, zoom_in, zoom_out) = ui.input(|i| {
            (
                !i.modifiers.alt && !i.modifiers.command && i.key_pressed(egui::Key::ArrowLeft),
                !i.modifiers.alt && !i.modifiers.command && i.key_pressed(egui::Key::ArrowRight),
                i.key_pressed(egui::Key::ArrowUp),
                i.key_pressed(egui::Key::ArrowDown),
                i.key_pressed(egui::Key::PageUp),
                i.key_pressed(egui::Key::PageDown),
                i.key_pressed(egui::Key::Plus) || i.key_pressed(egui::Key::Equals),
                i.key_pressed(egui::Key::Minus),
            )
        });
        if left {
            self.center -= step;
        }
        if right {
            self.center += step;
        }
        if let Some(t) = self.tracks.get_mut(target) {
            t.scroll_y += 36.0 * (f32::from(down) - f32::from(up)) + page * (f32::from(page_down) - f32::from(page_up));
        }
        if zoom_in {
            self.zoom_at(2.0, self.width / 2.0);
        }
        if zoom_out {
            self.zoom_at(0.5, self.width / 2.0);
        }
    }

    fn x_of(&self, rect: Rect, vs: f64, pos: f64) -> f32 {
        rect.left() + ((pos - vs) * self.ppb) as f32
    }

    fn draw_shared_overlays(&self, painter: &egui::Painter, rect: Rect, vs: f64, ve: f64) {
        let text = self.pal().text;
        let muted = self.pal().muted;
        let small = FontId::monospace(10.0);

        // Ruler.
        let len = (ve - vs).max(1.0);
        let label_w = 7.0 * commas(ve as u64).len() as f64 + 20.0;
        let mut step = nice_step(len / 10.0);
        while step as f64 * self.ppb < label_w {
            step = nice_step(step as f64 * 1.5);
        }
        let mut p = ((vs.max(0.0) as u64 + 1).div_ceil(step)) * step;
        while (p as f64) <= ve {
            let x = self.x_of(rect, vs, p as f64 - 0.5);
            painter.line_segment([Pos2::new(x, rect.top() + 16.0), Pos2::new(x, rect.top() + 24.0)], Stroke::new(1.0, muted));
            let label = commas(p);
            let half = 3.0 * label.len() as f32;
            if x - half >= rect.left() && x + half <= rect.right() {
                painter.text(Pos2::new(x, rect.top() + 9.0), Align2::CENTER_CENTER, label, small.clone(), text);
            }
            p += step;
        }
        painter.hline(rect.x_range(), rect.top() + RULER_H - 5.5, Stroke::new(1.0, muted));

        // Reference letters.
        if self.ppb >= f64::from(LETTERS_PPB)
            && let Some((view, _, _)) = &self.reference
            && view.region.contig == self.contigs[self.contig].name
            && let Some(reference) = &view.reference
        {
            // Proportional and growing with the zoom, as IGV's sequence track.
            let font = FontId::proportional((self.ppb as f32 * 0.85).clamp(11.0, 15.0));
            for p in vs.floor().max(view.region.start as f64) as u64..ve.ceil() as u64 {
                if let Some(&b) = reference.get((p - view.region.start) as usize) {
                    let [r, g, bl] = tracks::base_text_rgb(b, self.style.dark);
                    painter.text(
                        Pos2::new(self.x_of(rect, vs, p as f64 + 0.5), rect.top() + REF_TOP + REF_H / 2.0),
                        Align2::CENTER_CENTER,
                        (b as char).to_string(),
                        font.clone(),
                        Color32::from_rgb(r, g, bl),
                    );
                }
            }
        }

        // A shift-dragged selection, or the region being marked to zoom to, across every track.
        if let Some((c, a, b)) = self.zoom_drag.map(|(a, b)| (self.contig, a.min(b), a.max(b) + 1)).or(self.selection)
            && c == self.contig
        {
            let x0 = self.x_of(rect, vs, a as f64).max(rect.left());
            let x1 = self.x_of(rect, vs, b as f64).min(rect.right());
            if x1 > x0 {
                painter.rect_filled(
                    Rect::from_min_max(Pos2::new(x0, rect.top() + RULER_H - 6.0), Pos2::new(x1, rect.bottom())),
                    0.0,
                    Color32::from_rgba_unmultiplied(70, 130, 230, 40),
                );
            }
        }

        // Overlays: tinted bands across every track, each with a label tag in the ruler band, in
        // at most two lanes above the sequence, so tags never reach it or the tracks; labels
        // that don't fit are counted in a "+N" tag and shown on hovering the ruler.
        let contig = self.contigs.get(self.contig).map(|c| c.name.clone()).unwrap_or_default();
        let mut tags: Vec<(f32, f32, f32)> = Vec::new();
        let mut unshown = 0;
        const TAG_H: f32 = 15.0;
        let lanes = ((REF_TOP - 1.0) / TAG_H).floor().max(1.0) as usize;
        // Only regions in view: one wholly off to the side must not appear as an edge stripe
        // (clamping first made it look visible).
        for o in self.overlays.borrow().iter().filter(|o| o.contig == contig && (o.end as f64) > vs && (o.start as f64) < ve) {
            let x0 = self.x_of(rect, vs, o.start as f64).max(rect.left());
            let x1 = self.x_of(rect, vs, o.end as f64).min(rect.right()).max(x0 + 2.0);
            let [r, g, b] = crate::overlays::rgb(&o.color).unwrap_or([235, 190, 30]);
            let top = rect.top() + RULER_H;
            painter.rect_filled(Rect::from_min_max(Pos2::new(x0, top), Pos2::new(x1, rect.bottom())), 0.0, Color32::from_rgba_unmultiplied(r, g, b, 38));
            let edge = Stroke::new(1.0, Color32::from_rgba_unmultiplied(r, g, b, 170));
            painter.vline(x0, top..=rect.bottom(), edge);
            painter.vline(x1, top..=rect.bottom(), edge);
            if !o.label.is_empty() {
                let galley = painter.layout_no_wrap(o.label.clone(), FontId::proportional(11.0), Color32::WHITE);
                let w = galley.size().x + 8.0;
                let x = x0.clamp(rect.left() + 2.0, (rect.right() - w - 2.0).max(rect.left() + 2.0));
                let lane = (0..lanes).find(|&k| {
                    let y = rect.top() + 1.0 + k as f32 * TAG_H;
                    !tags.iter().any(|&(a, b, ty)| x < b && x + w > a && (ty - y).abs() < TAG_H)
                });
                let Some(lane) = lane else {
                    unshown += 1;
                    continue;
                };
                let y = rect.top() + 1.0 + lane as f32 * TAG_H;
                tags.push((x, x + w, y));
                painter.rect_filled(Rect::from_min_size(Pos2::new(x, y), Vec2::new(w, TAG_H)), 3.0, Color32::from_rgba_unmultiplied(r, g, b, 225));
                painter.galley(Pos2::new(x + 4.0, y + 1.5), galley, Color32::WHITE);
            }
        }
        if unshown > 0 {
            let galley = painter.layout_no_wrap(format!("+{unshown} more (hover the ruler)"), FontId::proportional(10.0), self.pal().text);
            let at = Pos2::new(rect.right() - galley.size().x - 8.0, rect.top() + 1.0);
            painter.rect_filled(Rect::from_min_size(at, galley.size() + Vec2::new(6.0, 3.0)), 3.0, self.pal().label_bg);
            painter.galley(at + Vec2::new(3.0, 1.5), galley, self.pal().text);
        }

        // The outlined base runs through every track.
        if let Some((contig, pos)) = self.highlight
            && contig == self.contig
        {
            let x0 = self.x_of(rect, vs, pos as f64);
            let x1 = self.x_of(rect, vs, pos as f64 + 1.0);
            let (l, r) = outline_columns(x0, x1, painter.ctx().pixels_per_point());
            let w = 1.0 / painter.ctx().pixels_per_point();
            // Two pixels clear of the reference letters, which start at the top of their row.
            let (top, bottom) = (rect.top() + REF_TOP - 2.0 - w, rect.bottom());
            let color = self.pal().outline;
            for x in [l, r] {
                painter.rect_filled(Rect::from_min_max(Pos2::new(x, top), Pos2::new(x + w, bottom)), 0.0, color);
            }
            painter.rect_filled(Rect::from_min_max(Pos2::new(l, top), Pos2::new(r + w, top + w)), 0.0, color);
        }
    }

    /// Splice junctions over a track's coverage: an arc per intron its loaded reads skip,
    /// rising from the coverage baseline, taller and heavier with more reads, labelled with the
    /// count. The colour is the strand from the splice-site motif (red +, blue −) when both
    /// ends are loaded, else grey. Only reads that pass the track's filters count (the same
    /// reads drawn), and junctions with a single read are left out when there are more than a
    /// few dozen.
    #[allow(clippy::too_many_arguments)]
    fn draw_junctions(&self, painter: &egui::Painter, section: Rect, vs: f64, ve: f64, t: &Track, view: &ggv_layout::ViewData, style: Style) {
        // Counted once per load and filter (every read's skips), not every frame.
        let all = t.junctions(style);
        let mut counts: std::collections::BTreeMap<(u64, u64), u32> =
            all.iter().filter(|((a, b), _)| (*a as f64) < ve && (*b as f64) > vs).map(|(k, n)| (*k, *n)).collect();
        if counts.len() > 40 {
            counts.retain(|_, n| *n > 1);
        }
        let Some(max) = counts.values().copied().max() else { return };
        let base = |p: u64| view.reference.as_ref().and_then(|r| r.get(p.checked_sub(view.region.start)? as usize).copied());
        let floor = section.top() + COVERAGE_TOP + COVERAGE_H;
        let clip = painter.with_clip_rect(Rect::from_min_max(Pos2::new(section.left(), section.top() + COVERAGE_TOP + 10.0), Pos2::new(section.right(), floor)).intersect(painter.clip_rect()));
        let font = FontId::monospace(10.0);
        let mut labels: Vec<Rect> = Vec::new();
        for (&(s, e), &n) in &counts {
            let strand = match (base(s), base(s + 1), base(e.saturating_sub(2)), base(e.saturating_sub(1))) {
                (Some(a), Some(b), Some(x), Some(z)) => ggv_model::splice_motif([a, b], [x, z]).1,
                _ => '.',
            };
            let color = match strand {
                '+' => Color32::from_rgb(210, 60, 60),
                '-' => Color32::from_rgb(60, 100, 220),
                _ => Color32::from_rgb(140, 140, 140),
            };
            let weight = (f64::from(n).ln_1p() / f64::from(max).ln_1p()) as f32;
            // The tallest arc and its label stay under the coverage label.
            let h = 6.0 + (COVERAGE_H - 32.0) * weight;
            let (x1, x2) = (self.x_of(section, vs, s as f64), self.x_of(section, vs, e as f64));
            const STEPS: usize = 32;
            let points: Vec<Pos2> = (0..=STEPS)
                .map(|i| {
                    let t = i as f32 / STEPS as f32;
                    Pos2::new(x1 + (x2 - x1) * t, floor - h * (std::f32::consts::PI * t).sin())
                })
                .collect();
            clip.add(egui::Shape::line(points, Stroke::new(1.0 + 2.0 * weight, color)));
            let mid = Pos2::new(((x1 + x2) / 2.0).clamp(section.left() + 10.0, section.right() - 10.0), floor - h - 1.0);
            // A count only where it doesn't overprint another: nearby junctions put their labels
            // in the same place.
            let galley = clip.layout_no_wrap(n.to_string(), font.clone(), color);
            let at = Rect::from_min_size(Pos2::new(mid.x - galley.size().x / 2.0, mid.y - galley.size().y), galley.size()).expand(1.0);
            if x2 - x1 > 14.0 && !labels.iter().any(|r: &Rect| r.intersects(at)) {
                clip.galley(at.min + Vec2::splat(1.0), galley, color);
                labels.push(at);
            }
        }
    }

    fn draw_track_overlays(&self, painter: &egui::Painter, i: usize, section: Rect, vs: f64, ve: f64, cov_scale: f32) {
        let t = &self.tracks[i];
        let muted = self.pal().muted;
        let small = FontId::monospace(10.0);
        let top = section.top();

        // Header: name, state, close button.
        painter.hline(section.x_range(), top + 0.5, Stroke::new(1.0, self.pal().rule_strong));
        painter.rect_filled(
            Rect::from_min_size(Pos2::new(section.left(), top + 1.0), Vec2::new(section.width(), HEADER_H - 1.0)),
            0.0,
            self.pal().header,
        );
        let mut header = t.title();
        if let Some(e) = t.current.as_ref().and_then(|c| c.style.filter.expr) {
            header.push_str(&format!("  ·  filter: {}", crate::expr::text(e)));
        }
        if t.missing_contig {
            header.push_str(&format!("  ·  no {} in this file", self.contigs[self.contig].name));
        } else if let Some(e) = &t.error {
            header.push_str(&format!("  ·  error: {e}"));
        }
        let right = self.header_rows(painter, section, self.row_scroll(Hit::Track(i), section));
        self.header_title(painter, section, header, right);
        let close = Layout::close_button(section);
        painter.text(close.center(), Align2::CENTER_CENTER, "×", FontId::proportional(14.0), self.pal().icon);

        let Some(cur) = &t.current else { return };
        let view = &cur.view;

        let cov_label = match &*cur.coverage {
            Coverage::Unavailable => "Coverage: zoom in (view wider than 5 Mb)".to_owned(),
            _ => {
                // The scale's top, as a depth: what the tallest bar means right now.
                let scale = self.style.coverage_scale;
                let usable = COVERAGE_H - 12.0;
                let units = usable / cov_scale;
                let depth = match scale {
                    tracks::CoverageScale::Log => units.exp2() - 1.0,
                    _ => units,
                };
                format!(
                    "Coverage [0-{}]{}{}",
                    commas(depth.round() as u64),
                    if scale == tracks::CoverageScale::Log { " log" } else { "" },
                    if matches!(&*cur.coverage, Coverage::Estimated(_)) { " · estimated from the index (zoom in under 5 Mb to count)" } else { "" }
                )
            }
        };
        painter.text(Pos2::new(section.left() + 4.0, top + COVERAGE_TOP + 1.0), Align2::LEFT_TOP, cov_label, small.clone(), muted);
        painter.hline(section.x_range(), top + SEPARATOR_Y, Stroke::new(1.0, self.pal().rule));
        if self.style.junctions && view.reads.is_some() {
            self.draw_junctions(painter, section, vs, ve, t, view, cur.style);
        }

        let packed = &cur.packed;
        let Some((reads, _)) = &view.reads else {
            painter.text(
                Pos2::new(section.left() + 4.0, top + READS_TOP + 4.0),
                Align2::LEFT_TOP,
                format!("Zoom in to see reads (view wider than {})", fmt_bp(read_window())),
                small,
                muted,
            );
            return;
        };

        let rl = tracks::row_layout(packed, cur.style);
        let reads_rect = Rect::from_min_max(Pos2::new(section.left(), top + READS_TOP), section.max).intersect(painter.clip_rect());
        let rp = painter.with_clip_rect(reads_rect);
        // Nothing to draw: say why rather than leave the track blank. Reads loaded around the view
        // (for panning) but none on screen count as none: no group labels or counts for them.
        let in_view = reads_in_view(reads, packed, vs, ve);
        if !in_view {
            let why = if reads.is_empty() {
                "No reads here".to_owned()
            } else if packed.rows.is_empty() {
                format!("No reads pass the read filters here ({} hidden)", commas(reads.len() as u64))
            } else {
                "No reads in view".to_owned()
            };
            rp.text(Pos2::new(section.left() + 4.0, top + READS_TOP + 4.0), Align2::LEFT_TOP, why, small.clone(), muted);
        }

        // Group labels, in the gap row above each group.
        for g in packed.groups.iter().filter(|_| in_view) {
            let y = top + READS_TOP + g.first_row as f32 * rl.pitch - t.scroll_y;
            if g.first_row > 0 {
                rp.hline(section.x_range(), y - rl.pitch / 2.0, Stroke::new(1.0, self.pal().rule_strong));
            }
            let label = format!("{}  ({} reads)", g.label, commas(g.reads as u64));
            let galley = rp.layout_no_wrap(label, FontId::proportional(10.5), self.pal().text);
            let at = Pos2::new(section.left() + 4.0, y + 1.0);
            rp.rect_filled(Rect::from_min_size(at - Vec2::new(2.0, 1.0), galley.size() + Vec2::new(4.0, 2.0)), 2.0, self.pal().label_bg);
            rp.galley(at, galley, self.pal().text);
        }

        // The selected read's alignments (mate, split parts), outlined wherever they are loaded.
        if let Some(name) = &self.selected_read {
            let first = (t.scroll_y / rl.pitch).floor() as usize;
            let last = (((t.scroll_y + section.height() - READS_TOP) / rl.pitch).ceil() as usize).min(packed.rows.len());
            for r in first..last {
                let y = top + READS_TOP + r as f32 * rl.pitch - t.scroll_y;
                for read in visible_in_row(reads, &packed.rows[r], vs, ve).filter(|x| &x.name == name) {
                    let x0 = self.x_of(section, vs, read.start as f64).max(section.left() - 2.0);
                    let x1 = self.x_of(section, vs, read.end as f64).min(section.right() + 2.0);
                    rp.rect_stroke(
                        Rect::from_min_max(Pos2::new(x0, y - 1.0), Pos2::new(x1, y + rl.row_h + 1.0)),
                        1.0,
                        Stroke::new(2.0, Color32::from_rgb(255, 140, 0)),
                        egui::StrokeKind::Outside,
                    );
                }
            }
        }

        // Letters and deletion labels for the visible rows only.
        if rl.row_h < 9.0 {
            return;
        }
        let first_row = (t.scroll_y / rl.pitch).floor() as usize;
        let last_row = (((t.scroll_y + section.height() - READS_TOP) / rl.pitch).ceil() as usize).min(packed.rows.len());
        let letters = self.ppb >= f64::from(LETTERS_PPB);
        let font = FontId::monospace(10.0);

        for r in first_row..last_row {
            let y = top + READS_TOP + r as f32 * rl.pitch - t.scroll_y;
            for read in visible_in_row(reads, &packed.rows[r], vs, ve) {
                if letters {
                    let clipped = if cur.style.soft_clips { &read.clip_mismatches[..] } else { &[] };
                    let from = read.mismatches.partition_point(|m| (m.pos as f64) < vs.floor());
                    let visible = read.mismatches[from..].iter().take_while(|m| (m.pos as f64) < ve);
                    for m in visible.chain(clipped.iter().filter(|m| (m.pos as f64) >= vs.floor() && (m.pos as f64) < ve)) {
                        let [cr, cg, cb] = tracks::base_rgb(m.base);
                        rp.text(
                            Pos2::new(self.x_of(section, vs, m.pos as f64 + 0.5), y + rl.row_h / 2.0),
                            Align2::CENTER_CENTER,
                            (m.base as char).to_string(),
                            font.clone(),
                            Color32::from_rgba_unmultiplied(cr, cg, cb, tracks::quality_alpha(m.quality)),
                        );
                    }
                }
                for d in read.deletions.iter().filter(|d| d.len >= 2) {
                    let x0 = self.x_of(section, vs, d.pos as f64);
                    let x1 = self.x_of(section, vs, (d.pos + u64::from(d.len)) as f64);
                    let label = d.len.to_string();
                    let lw = 6.0 * label.len() as f32;
                    if x1 - x0 > lw + 8.0 && x1 > section.left() && x0 < section.right() {
                        let c = Pos2::new((x0 + x1) / 2.0, y + rl.row_h / 2.0);
                        rp.rect_filled(Rect::from_center_size(c, Vec2::new(lw + 4.0, rl.row_h - 2.0)), 0.0, self.pal().bg);
                        rp.text(c, Align2::CENTER_CENTER, label, font.clone(), self.pal().text);
                    }
                }
            }
        }
    }

    fn hover(&self, response: &egui::Response, painter: &egui::Painter, i: usize, section: Rect, pos: Pos2, vs: f64) {
        let t = &self.tracks[i];
        let Some(cur) = &t.current else { return };
        let view = &cur.view;
        let base = (vs + f64::from(pos.x - section.left()) / self.ppb).floor();
        if base < view.region.start as f64 || base >= view.region.end as f64 {
            return;
        }
        let base = base as u64;
        let y = pos.y - section.top();
        let contig = &view.region.contig;

        // Coverage column.
        if (COVERAGE_TOP..COVERAGE_TOP + COVERAGE_H).contains(&y)
            && let Coverage::PerBase(cols) = &*cur.coverage
        {
            let c = cols[(base - view.region.start) as usize];
            let mut lines = vec![t.short_label(), format!("{contig}:{}", commas(base + 1)), format!("depth {}", c.depth)];
            if let Some(r) = c.ref_base {
                lines.push(format!("ref {}", r as char));
            }
            if c.counts.iter().any(|n| *n > 0) {
                let total = f64::from(c.depth.max(1));
                for (k, b) in ["A", "C", "G", "T", "N"].iter().enumerate() {
                    if c.counts[k] > 0 {
                        lines.push(format!("{b}: {} ({:.0}%)", c.counts[k], 100.0 * f64::from(c.counts[k]) / total));
                    }
                }
            }
            if c.deletions > 0 {
                lines.push(format!("DEL: {}", c.deletions));
            }
            show_tooltip(response, lines);
            return;
        }

        // Read under the pointer.
        let packed = &cur.packed;
        let Some((reads, _)) = &view.reads else { return };
        if y < READS_TOP {
            return;
        }
        let rl = tracks::row_layout(packed, cur.style);
        let row = ((y - READS_TOP + t.scroll_y) / rl.pitch).floor() as usize;
        let Some(row_reads) = packed.rows.get(row) else { return };
        let k = row_reads.partition_point(|&i| reads[i].end <= base);
        let Some(&idx) = row_reads.get(k) else { return };
        let r = &reads[idx];
        if r.start > base {
            return;
        }

        let ry = section.top() + READS_TOP + row as f32 * rl.pitch - t.scroll_y;
        let x0 = self.x_of(section, vs, r.start as f64).max(section.left() - 2.0);
        let x1 = self.x_of(section, vs, r.end as f64).min(section.right() + 2.0);
        painter.with_clip_rect(section).rect_stroke(
            Rect::from_min_max(Pos2::new(x0, ry), Pos2::new(x1, ry + rl.row_h)),
            0.0,
            Stroke::new(1.0, self.pal().text),
            egui::StrokeKind::Outside,
        );

        let at_base = if let Some(m) = r.mismatches.iter().find(|m| m.pos == base) {
            format!("{} (Q{}) vs ref", m.base as char, m.quality)
        } else if r.deletions.iter().any(|d| (d.pos..d.pos + u64::from(d.len)).contains(&base)) {
            "deleted".to_owned()
        } else if r.ref_skips.iter().any(|d| (d.pos..d.pos + u64::from(d.len)).contains(&base)) {
            "skipped (N)".to_owned()
        } else if let Some((_, l)) = tracks::mod_likelihoods(r).find(|(p, _)| *p == base) {
            format!("modified {:.0}%", f64::from(l) / 2.55)
        } else if view.reference.is_some() {
            "matches ref".to_owned()
        } else {
            "aligned".to_owned()
        };
        show_tooltip(
            response,
            vec![
                format!("{}  [{}]", r.name, t.short_label()),
                format!("{contig}:{}-{}  {}", commas(r.start + 1), commas(r.end), if r.reverse { "(−)" } else { "(+)" }),
                format!("MAPQ {}  flags 0x{:x}", r.mapq, r.flags),
                format!(
                    "{} mismatches · {} del · {} ins · clips {}/{}",
                    r.mismatches.len(),
                    r.deletions.len(),
                    r.insertions.len(),
                    r.soft_clip_left,
                    r.soft_clip_right
                ),
                format!("at {}: {at_base}", commas(base + 1)),
                read_extras(r),
            ],
        );
    }
}

impl App {
    fn draw_gene_overlays(&self, painter: &egui::Painter, i: usize, section: Rect, vs: f64, ve: f64) {
        let g = &self.genes[i];
        let top = section.top();
        painter.hline(section.x_range(), top + 0.5, Stroke::new(1.0, self.pal().rule_strong));
        painter.rect_filled(
            Rect::from_min_size(Pos2::new(section.left(), top + 1.0), Vec2::new(section.width(), HEADER_H - 1.0)),
            0.0,
            self.pal().header,
        );
        let right = self.header_rows(painter, section, self.row_scroll(Hit::Genes(i), section));
        self.header_title(painter, section, self.gene_header(i), right);
        let close = Layout::close_button(section);
        painter.text(close.center(), Align2::CENTER_CENTER, "×", FontId::proportional(14.0), self.pal().icon);

        let Some(built) = &g.built else { return };
        let rows_rect = Rect::from_min_max(Pos2::new(section.left(), top + genes::ROWS_TOP), section.max).intersect(painter.clip_rect());
        let rp = painter.with_clip_rect(rows_rect);
        // Nothing here: say so rather than leave the track blank (the header's file-wide counts
        // looked like a failed load).
        if self.genes_in_view(i) == Some(0) {
            let what = if built.is_all() { "No transcripts here" } else { "No genes here" };
            rp.text(Pos2::new(section.left() + 4.0, top + genes::ROWS_TOP + 2.0), Align2::LEFT_TOP, what, FontId::monospace(10.0), self.pal().muted);
        }
        let first = (g.scroll_y / genes::PITCH).floor() as usize;
        let last = (((g.scroll_y + rows_rect.height()) / genes::PITCH).ceil() as usize).min(built.rows.len());
        let font = FontId::proportional(11.0);

        for r in first..last {
            let y = top + genes::ROWS_TOP + r as f32 * genes::PITCH - g.scroll_y;
            let mid = y + genes::FEATURE_H / 2.0;
            let row = &built.rows[r];
            let from = row.partition_point(|&k| (built.items[k].end as f64) <= vs);
            let mut label_end = f32::NEG_INFINITY;
            for &k in row[from..].iter().take_while(|&&k| (built.items[k].start as f64) < ve) {
                let t = &built.items[k];
                let [cr, cg, cb] = genes::feature_rgb(t.cds.is_some(), self.style.dark);
                let color = Color32::from_rgb(cr, cg, cb);
                self.draw_codons(&rp, section, t, y, genes::FEATURE_H, vs, ve);

                // Strand chevrons along introns wide enough to hold them.
                let dir = match t.strand {
                    ggv_annot::Strand::Forward => 1.0,
                    ggv_annot::Strand::Reverse => -1.0,
                    ggv_annot::Strand::Unknown => 0.0,
                };
                if dir != 0.0 {
                    for w in t.exons.windows(2) {
                        let x0 = self.x_of(section, vs, w[0].1 as f64).max(section.left());
                        let x1 = self.x_of(section, vs, w[1].0 as f64).min(section.right());
                        let mut x = x0 + 8.0;
                        while x + 4.0 < x1 {
                            rp.line_segment([Pos2::new(x - 2.0 * dir, mid - 3.0), Pos2::new(x + dir, mid)], Stroke::new(1.0, color));
                            rp.line_segment([Pos2::new(x + dir, mid), Pos2::new(x - 2.0 * dir, mid + 3.0)], Stroke::new(1.0, color));
                            x += 24.0;
                        }
                    }
                }

                // Name under the feature, centred on its visible part, unless it would collide.
                let name = if built_all(g) { format!("{} ({})", t.gene_name, ggv_annot::without_copy_suffix(&t.id)) } else { t.gene_name.to_string() };
                let x0 = self.x_of(section, vs, t.start as f64).max(section.left());
                let x1 = self.x_of(section, vs, t.end as f64).min(section.right());
                let galley = rp.layout_no_wrap(name, font.clone(), self.pal().text);
                let w = galley.size().x;
                // In a pane narrower than the label, it starts at the left edge.
                let lx = ((x0 + x1) / 2.0 - w / 2.0).clamp(section.left() + 2.0, (section.right() - w - 2.0).max(section.left() + 2.0));
                if lx > label_end + 6.0 {
                    rp.galley(Pos2::new(lx, y + genes::FEATURE_H), galley, self.pal().text);
                    label_end = lx + w;
                }
            }
        }
    }

    /// How many transcripts gene track `i` shows in the view, once its rows are built for this
    /// contig.
    fn genes_in_view(&self, i: usize) -> Option<usize> {
        let built = self.genes.get(i)?.built.as_ref()?;
        if self.contigs.get(self.contig).is_none_or(|c| c.name != built.contig) {
            return None;
        }
        let (vs, ve) = self.visible();
        Some(built.items.iter().filter(|t| (t.end as f64) > vs && (t.start as f64) < ve).count())
    }

    /// Gene track `i`'s header: what is in view (as variant tracks say), then the file's counts.
    fn gene_header(&self, i: usize) -> String {
        let g = &self.genes[i];
        match self.genes_in_view(i) {
            Some(n) => format!("Genes  ·  {}  ·  {} in view  ·  {}", g.file_name(), commas(n as u64), g.status()),
            None => format!("Genes  ·  {}  ·  {}", g.file_name(), g.status()),
        }
    }

    /// The transcript under `pos` in gene track `i`, with its row.
    fn gene_at(&self, i: usize, section: Rect, pos: Pos2) -> Option<(usize, &ggv_annot::Transcript)> {
        let g = self.genes.get(i)?;
        let built = g.built.as_ref()?;
        let y = pos.y - section.top() - genes::ROWS_TOP + g.scroll_y;
        if y < 0.0 {
            return None;
        }
        let r = (y / genes::PITCH).floor() as usize;
        let row = built.rows.get(r)?;
        let (vs, _) = self.visible();
        let base = (vs + f64::from(pos.x - section.left()) / self.ppb).floor().max(0.0) as u64;
        let k = row.partition_point(|&k| built.items[k].end <= base);
        let t = &built.items[*row.get(k)?];
        (t.start <= base).then_some((r, t))
    }

    fn hover_genes(&self, response: &egui::Response, painter: &egui::Painter, i: usize, section: Rect, pos: Pos2, vs: f64) {
        let g = &self.genes[i];
        let Some(built) = &g.built else { return };
        let Some((r, t)) = self.gene_at(i, section, pos) else { return };
        let base = (vs + f64::from(pos.x - section.left()) / self.ppb).floor() as u64;

        let ry = section.top() + genes::ROWS_TOP + r as f32 * genes::PITCH - g.scroll_y;
        let x0 = self.x_of(section, vs, t.start as f64).max(section.left() - 2.0);
        let x1 = self.x_of(section, vs, t.end as f64).min(section.right() + 2.0);
        painter.with_clip_rect(section).rect_stroke(
            Rect::from_min_max(Pos2::new(x0, ry - 1.0), Pos2::new(x1, ry + genes::FEATURE_H + 1.0)),
            1.0,
            Stroke::new(1.0, self.pal().text),
            egui::StrokeKind::Outside,
        );

        let strand = match t.strand {
            ggv_annot::Strand::Forward => "+",
            ggv_annot::Strand::Reverse => "−",
            ggv_annot::Strand::Unknown => "?",
        };
        let mut lines = vec![
            format!("{}  {}", t.gene_name, ggv_annot::without_copy_suffix(&t.id)),
            format!("{}:{}-{} ({strand})", built.contig, commas(t.start + 1), commas(t.end)),
            format!("{} · {} exons · {} bp coding", if t.biotype.is_empty() { "-" } else { &t.biotype }, t.exons.len(), commas(t.cds_len())),
        ];
        let mut tags = Vec::new();
        if t.tags & ggv_annot::TAG_MANE_SELECT != 0 {
            tags.push("MANE Select");
        }
        if t.tags & ggv_annot::TAG_CANONICAL != 0 {
            tags.push("Ensembl canonical");
        }
        if !tags.is_empty() {
            lines.push(tags.join(" · "));
        }
        // Exons are numbered in the direction of transcription.
        if let Some(loc) = ggv_annot::coding::locate(t, base) {
            lines.push(format!("at {}: {}", commas(base + 1), loc.feature));
            if let Some(c) = self.base_in_transcript(t, base) {
                lines.push(format!("  {c}"));
            }
        }
        show_tooltip(response, lines);
    }
}

impl App {
    /// "Reads up to 1 Mb": how wide a view still shows reads (presets, or any size).
    fn read_window_ui(&mut self, ui: &mut egui::Ui) {
        let rw = read_window();
        let mut new = None;
        ui.menu_button(format!("Reads up to {}", fmt_bp(rw)), |ui| {
                for bp in [30_000, 100_000, 250_000, 500_000, 1_000_000, 2_000_000, 5_000_000] {
                    let hint = match bp {
                        30_000 => " (IGV's default)",
                        1_000_000 => " (default)",
                        5_000_000 => " (slow on deep data)",
                        _ => "",
                    };
                    if ui.selectable_label(rw == bp, format!("{}{hint}", fmt_bp(bp))).clicked() {
                        new = Some(bp);
                    }
                }
                ui.separator();
                ui.horizontal(|ui| {
                    ui.label("Any size:");
                    let mut kb = rw / 1000;
                    if ui
                        .add(egui::DragValue::new(&mut kb).range(MIN_READ_WINDOW / 1000..=MAX_READ_WINDOW / 1000).speed(10).suffix(" kb"))
                        .changed()
                    {
                        new = Some(kb * 1000);
                    }
                });
            })
            .response
            .on_hover_text(
                "Views up to this wide show reads; wider views show coverage only.\n\
                 Cost grows with the reads in view. On 35x short reads (about 265,000 reads per Mb):\n\
                 1 Mb: about 1 GB per track, panning above 60 fps\n\
                 2 Mb: about 1.5 GB, occasional pauses when zooming\n\
                 5 Mb: about 3.5 GB, slower panning\n\
                 Sparse data (long reads, exomes, low-pass) goes much further at the same cost.",
            );
        if let Some(bp) = new {
            save_read_window(bp);
            self.notify(format!("Reads are shown for views up to {}", fmt_bp(read_window())));
        }
    }

    fn draw_signal_overlays(&self, painter: &egui::Painter, i: usize, section: Rect, vs: f64, ve: f64) {
        let s = &self.signals[i];
        let top = section.top();
        painter.hline(section.x_range(), top + 0.5, Stroke::new(1.0, self.pal().rule_strong));
        painter.rect_filled(Rect::from_min_size(Pos2::new(section.left(), top + 1.0), Vec2::new(section.width(), HEADER_H - 1.0)), 0.0, self.pal().header);
        let n = section.width().max(1.0) as usize;
        let body = Rect::from_min_max(Pos2::new(section.left(), top + HEADER_H + 2.0), Pos2::new(section.right(), section.bottom() - 2.0));
        let mut state = if let Some(e) = &s.error {
            format!("error: {e}")
        } else if s.is_loading() {
            "loading…".to_owned()
        } else {
            String::new()
        };
        if !s.is_loading() && s.error.is_none() {
            if s.is_seg {
                for (k, row) in s.rows.iter().enumerate() {
                    let y0 = body.top() + k as f32 * crate::signals::SEG_ROW_H;
                    for (c, b) in s.columns(k, vs, ve, n).iter().enumerate() {
                        if let Some(b) = b {
                            painter.rect_filled(Rect::from_min_size(Pos2::new(section.left() + c as f32, y0), Vec2::new(1.0, crate::signals::SEG_ROW_H - 1.0)), 0.0, crate::signals::seg_color(b.mean));
                        }
                    }
                    painter.text(Pos2::new(section.left() + 3.0, y0 + crate::signals::SEG_ROW_H / 2.0), Align2::LEFT_CENTER, row, FontId::monospace(8.0), self.pal().muted);
                }
                state = format!("{} samples · blue loss, red gain (±1.5)", s.rows.len());
            } else {
                let cols = s.columns(0, vs, ve, n);
                let (lo, hi) = cols.iter().flatten().fold((0.0f64, 0.0f64), |(lo, hi), b| (lo.min(b.min), hi.max(b.max)));
                let span = (hi - lo).max(1e-9);
                let y_of = |v: f64| body.bottom() - ((v - lo) / span) as f32 * body.height();
                let base = y_of(0.0);
                let color = egui::Color32::from_rgb(70, 110, 200);
                for (c, b) in cols.iter().enumerate() {
                    if let Some(b) = b {
                        let (y_mean, x) = (y_of(b.mean), section.left() + c as f32);
                        painter.rect_filled(Rect::from_min_max(Pos2::new(x, y_mean.min(base)), Pos2::new(x + 1.0, y_mean.max(base).max(y_mean.min(base) + 1.0))), 0.0, color);
                    }
                }
                state = format!("range {} – {}", fmt_value(lo), fmt_value(hi));
            }
        }
        self.header_title(painter, section, format!("Signal  ·  {}  ·  {state}", s.file_name()), Layout::grip(section).left() - 6.0);
        let close = Layout::close_button(section);
        painter.text(close.center(), Align2::CENTER_CENTER, "×", FontId::proportional(14.0), self.pal().icon);
    }

    fn hover_signal(&self, response: &egui::Response, i: usize, section: Rect, pos: Pos2, vs: f64) {
        let s = &self.signals[i];
        let base = vs + f64::from(pos.x - section.left()) / self.ppb;
        let row = if s.is_seg { ((pos.y - section.top() - HEADER_H - 2.0) / crate::signals::SEG_ROW_H).floor().max(0.0) as usize } else { 0 };
        let Some(b) = s.at(row, base) else { return };
        let contig = self.contigs.get(self.contig).map(|c| c.name.clone()).unwrap_or_default();
        let mut lines = vec![s.file_name(), format!("{contig}:{}", commas(base.max(0.0) as u64 + 1))];
        if s.is_seg && let Some(name) = s.rows.get(row) {
            lines.push(format!("sample {name}"));
        }
        lines.push(if b.min == b.max { format!("value {}", fmt_value(b.mean)) } else { format!("mean {} (min {}, max {})", fmt_value(b.mean), fmt_value(b.min), fmt_value(b.max)) });
        show_tooltip(response, lines);
    }

    fn draw_variant_overlays(&self, painter: &egui::Painter, i: usize, section: Rect, vs: f64, ve: f64) {
        let v = &self.variants[i];
        let top = section.top();
        painter.hline(section.x_range(), top + 0.5, Stroke::new(1.0, self.pal().rule_strong));
        painter.rect_filled(
            Rect::from_min_size(Pos2::new(section.left(), top + 1.0), Vec2::new(section.width(), HEADER_H - 1.0)),
            0.0,
            self.pal().header,
        );
        let state = if let Some(e) = &v.error {
            format!("error: {e}")
        } else if ve - vs > variants::VARIANT_WINDOW as f64 {
            "zoom in to see variants (view wider than 10 Mb)".to_owned()
        } else if let Some(cur) = &v.current {
            let shown = cur.variants.iter().filter(|x| (x.end as f64) > vs && (x.start as f64) < ve).count();
            let more = if cur.overflow > 0 { format!(" · {} overlapping in the last row", commas(cur.overflow as u64)) } else { String::new() };
            format!("{} in view · loaded {} in {:.0} ms{more}", commas(shown as u64), commas(cur.variants.len() as u64), cur.elapsed.as_secs_f64() * 1000.0)
        } else {
            "loading…".to_owned()
        };
        let right = self.header_rows(painter, section, self.row_scroll(Hit::Variants(i), section));
        let title = self.header_title(
            painter,
            section,
            format!("Variants  ·  {}  ·  {} sample{}  ·  {state}", v.file_name(), v.samples.len(), if v.samples.len() == 1 { "" } else { "s" }),
            right,
        );
        let close = Layout::close_button(section);
        // What the site colours mean, after the title while it fits; nothing (not a bare
        // "CLNSIG:") when no loaded record has the field.
        let legend = v.legend(6);
        if let variants::SiteColor::Info(key) = &v.color
            && !legend.is_empty()
        {
            let font = FontId::proportional(10.5);
            let cy = top + HEADER_H / 2.0;
            let mut x = title.right() + 14.0;
            let galley = painter.layout_no_wrap(format!("{key}:"), font.clone(), self.pal().text);
            if x + galley.size().x < right - 2.0 {
                let w = galley.size().x;
                painter.galley(Pos2::new(x, cy - galley.size().y / 2.0), galley, self.pal().text);
                x += w + 6.0;
                for (label, [r, g, b]) in legend {
                    let galley = painter.layout_no_wrap(label, font.clone(), self.pal().text);
                    if x + 12.0 + galley.size().x > right - 2.0 {
                        break;
                    }
                    painter.rect_filled(Rect::from_center_size(Pos2::new(x + 4.0, cy), Vec2::splat(8.0)), 1.0, Color32::from_rgb(r, g, b));
                    let w = galley.size().x;
                    painter.galley(Pos2::new(x + 11.0, cy - galley.size().y / 2.0), galley, self.pal().text);
                    x += 11.0 + w + 10.0;
                }
            }
        }
        painter.text(close.center(), Align2::CENTER_CENTER, "×", FontId::proportional(14.0), self.pal().icon);

        // Review decisions: a coloured tick above each reviewed site.
        if let Some(review) = self.review.as_ref().filter(|r| r.track == v.id)
            && let Some(cur) = self.shown_variants(i)
        {
            let contig = &self.contigs[self.contig].name;
            for (var, &row) in cur.variants.iter().zip(cur.rows.iter()).filter(|(x, _)| (x.end as f64) > vs && (x.start as f64) < ve) {
                let key = crate::review::SiteKey {
                    contig: contig.clone(),
                    pos: var.start + 1,
                    reference: var.reference.clone(),
                    alt: var.alts.join(","),
                };
                if let Some(d) = review.decisions.get(&key) {
                    let (s, e) = var.shown();
                    let x = self.x_of(section, vs, (s + e) as f64 / 2.0);
                    let [r, g, b] = d.label.rgb();
                    painter.rect_filled(
                        Rect::from_center_size(Pos2::new(x, top + variants::site_y(row) - 1.0), Vec2::new(7.0, 4.0)),
                        1.0,
                        Color32::from_rgb(r, g, b),
                    );
                }
            }
        }

        // Structural variants: a type-and-size tag in the site bar where it fits (deletions,
        // duplications, inversions span their bases already), and for insertions, which take
        // no reference bases, a purple caret at the insertion point with the size beside it.
        if let Some(cur) = self.shown_variants(i) {
            let small = FontId::proportional(10.0);
            let all = Rect::from_min_max(
                Pos2::new(section.left(), top + variants::SITE_TOP),
                Pos2::new(section.right(), top + variants::site_y(cur.site_rows as u16 - 1) + variants::SITE_H),
            );
            let sp = painter.with_clip_rect(all.expand2(Vec2::new(0.0, 2.0)).intersect(painter.clip_rect()));
            let mut last_label_end = vec![f32::NEG_INFINITY; cur.site_rows];
            for ((var, &row), sv) in cur.variants.iter().zip(cur.rows.iter()).zip(cur.sv.iter()).filter(|((x, _), _)| (x.end as f64) > vs && (x.start as f64) < ve) {
                let Some((kind, len)) = *sv else { continue };
                let site = Rect::from_min_size(Pos2::new(section.left(), top + variants::site_y(row)), Vec2::new(section.width(), variants::SITE_H));
                let last_label_end = &mut last_label_end[usize::from(row)];
                let label = match len {
                    Some(n) => format!("{kind} {}", sv_size(n)),
                    None => kind.to_owned(),
                };
                let (s, e) = var.shown();
                let (x0, x1) = (self.x_of(section, vs, s as f64), self.x_of(section, vs, e as f64));
                let ins = kind == "INS";
                if ins {
                    let purple = Color32::from_rgb(140, 60, 200);
                    let cx = x0;
                    sp.add(egui::Shape::convex_polygon(
                        vec![Pos2::new(cx - 4.0, site.top()), Pos2::new(cx + 4.0, site.top()), Pos2::new(cx, site.top() + 6.0)],
                        purple,
                        Stroke::NONE,
                    ));
                    sp.vline(cx, site.y_range(), Stroke::new(1.5, purple));
                }
                let galley = sp.layout_no_wrap(label, small.clone(), Color32::WHITE);
                let w = galley.size().x + 6.0;
                // Inside a span wide enough, else just right of the site; never over another tag.
                // Centred on the part of the span on screen: a long deletion's middle can be off it.
                let (v0, v1) = (x0.max(site.left()), x1.min(site.right()));
                let lx = if !ins && v1 - v0 >= w + 4.0 { ((v0 + v1) / 2.0 - w / 2.0).max(site.left()) } else { x1.max(x0) + 4.0 };
                if lx < *last_label_end + 4.0 || lx + w > site.right() {
                    continue;
                }
                *last_label_end = lx + w;
                // Coloured by an INFO field (ClinVar's significance), the tag keeps the site's colour.
                let bg = match &v.color {
                    variants::SiteColor::Info(key) => {
                        let [r, g, b] = variants::info_rgb(key, var.info_value(key));
                        Color32::from_rgb(r, g, b)
                    }
                    _ if ins => Color32::from_rgb(140, 60, 200),
                    _ => Color32::from_rgba_unmultiplied(40, 40, 40, 200),
                };
                sp.rect_filled(Rect::from_min_size(Pos2::new(lx, site.top() + 0.5), Vec2::new(w, variants::SITE_H - 1.0)), 2.0, bg);
                sp.galley(Pos2::new(lx + 3.0, site.top() + (variants::SITE_H - galley.size().y) / 2.0), galley, Color32::WHITE);
            }
        }

        // Sample names at the left edge of their genotype rows, never over a genotype: when a
        // site sits under the label it moves to the right edge, and with sites under both it
        // is left out (hovering a genotype names its sample).
        let gt_top = v.gt_top();
        let rows = Rect::from_min_max(Pos2::new(section.left(), top + gt_top), section.max).intersect(painter.clip_rect());
        let rp = painter.with_clip_rect(rows);
        let font = FontId::proportional(10.0);
        // Screen x-extents of the sites in view, as drawn (at least 3 px, centred).
        let sites: Vec<(f32, f32)> = self.shown_variants(i).map_or_else(Vec::new, |cur| {
            cur.variants
                .iter()
                .filter_map(|var| {
                    let (s, e) = var.shown();
                    ((e as f64) >= vs && (s as f64) <= ve).then(|| {
                        let (x0, x1) = (self.x_of(section, vs, s as f64), self.x_of(section, vs, e as f64));
                        let mid = (x0 + x1) / 2.0;
                        (x0.min(mid - 1.5), x1.max(mid + 1.5))
                    })
                })
                .collect()
        });
        let clear = |l: f32, r: f32| !sites.iter().any(|&(a, b)| a < r && b > l);
        for (s, name) in v.samples.iter().enumerate() {
            let y = top + gt_top + s as f32 * variants::GT_PITCH - v.scroll_y;
            let galley = rp.layout_no_wrap(name.clone(), font.clone(), self.pal().text);
            let w = galley.size().x + 6.0;
            let left = if clear(section.left(), section.left() + w + 1.0) {
                section.left()
            } else if clear(section.right() - w - 1.0, section.right()) {
                section.right() - w
            } else {
                continue;
            };
            rp.rect_filled(Rect::from_min_size(Pos2::new(left, y), Vec2::new(w, variants::GT_H)), 0.0, self.pal().label_bg);
            rp.galley(Pos2::new(left + 3.0, y + (variants::GT_H - galley.size().y) / 2.0), galley, self.pal().text);
        }
    }

    /// What variant track `i` draws: what it loaded, when that is on the contig in view and the
    /// view is narrow enough to show variants (as the GPU pass decides). Tags, review ticks and
    /// clicks use only this, so nothing loaded for an earlier view shows elsewhere.
    fn shown_variants(&self, i: usize) -> Option<&variants::Loaded> {
        let cur = self.variants.get(i)?.current.as_ref()?;
        let (vs, ve) = self.visible();
        (self.contigs.get(self.contig).is_some_and(|c| c.name == cur.region.contig) && ve - vs <= variants::VARIANT_WINDOW as f64).then_some(cur)
    }

    /// The variant under `pos` in variant track `i` (drawn at least 3 px wide, so within 2 px):
    /// in the site row under the pointer, or over the genotype rows, the nearest in any row
    /// (the top row first).
    fn variant_at(&self, i: usize, section: Rect, pos: Pos2, vs: f64) -> Option<&ggv_io::Variant> {
        let cur = self.shown_variants(i)?;
        let slack = 2.0 / self.ppb;
        let at = vs + f64::from(pos.x - section.left()) / self.ppb;
        let y = pos.y - section.top() - variants::SITE_TOP;
        let row = (y >= 0.0 && y < cur.site_rows as f32 * variants::SITE_PITCH).then(|| (y / variants::SITE_PITCH) as u16);
        cur.variants
            .iter()
            .zip(cur.rows.iter())
            .filter(|(x, r)| {
                let (s, e) = x.shown();
                (s as f64) - slack <= at && at < (e as f64) + slack && row.is_none_or(|row| **r == row)
            })
            .min_by(|(a, ra), (b, rb)| ra.cmp(rb).then((a.shown().0 as f64 - at).abs().total_cmp(&(b.shown().0 as f64 - at).abs())))
            .map(|(v, _)| v)
    }

    fn hover_variants(&self, response: &egui::Response, painter: &egui::Painter, i: usize, section: Rect, pos: Pos2, vs: f64) {
        let v = &self.variants[i];
        let Some(cur) = &v.current else { return };
        let Some(var) = self.variant_at(i, section, pos, vs) else { return };

        let (s, e) = var.shown();
        let mid = self.x_of(section, vs, (s + e) as f64 / 2.0);
        let x0 = self.x_of(section, vs, s as f64).min(mid - 1.5);
        let x1 = self.x_of(section, vs, e as f64).max(mid + 1.5);
        painter.with_clip_rect(section).rect_stroke(
            Rect::from_min_max(Pos2::new(x0 - 1.0, section.top() + variants::SITE_TOP - 1.0), Pos2::new(x1 + 1.0, section.bottom() - 1.0)),
            0.0,
            Stroke::new(1.0, self.pal().text),
            egui::StrokeKind::Outside,
        );

        let contig = &cur.region.contig;
        let mut lines = vec![
            format!(
                "{contig}:{}  {}  {}",
                commas(var.start + 1),
                match var.structural(50) {
                    Some((k, Some(n))) => format!("{k} · {} (structural)", sv_size(n)),
                    Some((k, None)) => format!("{k} (structural)"),
                    None => var.kind(),
                },
                if var.id == "." { "" } else { &var.id }
            ),
            format!("{} > {}", ellipsize(&var.reference, 40), ellipsize(&var.alts.join(","), 40)),
            format!(
                "QUAL {}  FILTER {}  AF {}",
                var.qual.map_or_else(|| ".".to_owned(), |q| format!("{q:.1}")),
                if var.filters.is_empty() { ".".to_owned() } else { var.filters.join(",") },
                var.af_text(),
            ),
        ];
        let y = pos.y - section.top() - v.gt_top() + v.scroll_y;
        if y >= 0.0
            && let Some(call) = var.samples.get((y / variants::GT_PITCH) as usize)
        {
            let s = (y / variants::GT_PITCH) as usize;
            let mut l = format!("{}: {}", v.samples[s], call.gt);
            if let Some(dp) = call.dp {
                l.push_str(&format!("  DP {dp}"));
            }
            if let Some(gq) = call.gq {
                l.push_str(&format!("  GQ {gq}"));
            }
            if let Some(ad) = call.ad_text() {
                l.push_str(&format!("  AD {ad}"));
            }
            lines.push(l);
        } else {
            let count = |c| var.samples.iter().filter(|s| s.class == c).count();
            lines.push(format!(
                "het {} · hom-alt {} · hom-ref {} · no-call {}",
                count(ggv_io::GenotypeClass::Het),
                count(ggv_io::GenotypeClass::HomAlt),
                count(ggv_io::GenotypeClass::HomRef),
                count(ggv_io::GenotypeClass::NoCall)
            ));
        }
        // What it does to the genes here.
        let contig = self.contigs.get(self.contig).map(|c| c.name.clone()).unwrap_or_default();
        for a in var.alts.iter().filter(|a| !a.starts_with('<') && *a != "*").take(2) {
            for c in self.consequences(&contig, var.start, &var.reference, a).into_iter().take(3) {
                lines.push(format!("{} {}  {}{}  {}", c.gene, c.transcript, c.hgvs_c, c.hgvs_p.map(|p| format!("  {p}")).unwrap_or_default(), c.effect));
            }
        }
        // ClinVar's summary, or the first INFO fields; everything is a click away.
        let clinvar = variant_info::clinvar_lines(var);
        if clinvar.is_empty() {
            for (k, val) in var.info_pairs().take(6) {
                lines.push(ellipsize(&if val.is_empty() { k.to_owned() } else { format!("{k}={val}") }, 80));
            }
        } else {
            lines.extend(clinvar.iter().take(5).map(|(k, val)| ellipsize(&format!("{k}: {val}"), 80)));
        }
        lines.push("click for every INFO field".to_owned());
        show_tooltip(response, lines);
    }
}

fn ellipsize(s: &str, max: usize) -> String {
    if s.chars().count() <= max { s.to_owned() } else { format!("{}…", s.chars().take(max).collect::<String>()) }
}

fn built_all(g: &GeneTrack) -> bool {
    g.built.as_ref().is_some_and(genes::Built::is_all)
}

/// Haplotype, read group, insert size and mate, when present.
fn read_extras(r: &AlignedRead) -> String {
    let mut parts = Vec::new();
    if let Some(hp) = r.haplotype {
        parts.push(format!("HP {hp}"));
    }
    if let Some(rg) = &r.read_group {
        parts.push(format!("RG {rg}"));
    }
    if let Some(m) = r.mate {
        if m.same_contig {
            parts.push(format!("insert {}  mate {} {}", r.template_length, commas(m.pos.max(0) as u64 + 1), if m.reverse { "(−)" } else { "(+)" }));
        } else {
            parts.push("mate on another contig".to_owned());
        }
    }
    if !r.mods.is_empty() {
        // Calls per code, with how many are likely modified (> 50%).
        let mut codes: Vec<(u8, usize, usize)> = Vec::new();
        for m in &r.mods {
            match codes.iter_mut().find(|c| c.0 == m.code) {
                Some(c) => {
                    c.1 += 1;
                    c.2 += usize::from(m.prob >= 128);
                }
                None => codes.push((m.code, 1, usize::from(m.prob >= 128))),
            }
        }
        let list: Vec<String> = codes.iter().map(|(c, n, hi)| format!("{}: {hi}/{n} called", *c as char)).collect();
        parts.push(format!("mods {}", list.join(", ")));
    }
    if parts.is_empty() { "unpaired, no HP/RG".to_owned() } else { parts.join(" · ") }
}




fn color_options(ui: &mut egui::Ui, value: &mut tracks::ColorBy, tag: &mut String) {
    for c in [tracks::ColorBy::None, tracks::ColorBy::Strand, tracks::ColorBy::Haplotype, tracks::ColorBy::InsertSize, tracks::ColorBy::BaseMods] {
        ui.selectable_value(value, c, color_name(c));
    }
    if let Some(t) = tag_row(ui, matches!(value, tracks::ColorBy::Tag(_)), tag) {
        *value = tracks::ColorBy::Tag(t);
    }
}

/// "tag [CB]": the tag to colour or group by, chosen when the name is two characters.
fn tag_row(ui: &mut egui::Ui, selected: bool, tag: &mut String) -> Option<[u8; 2]> {
    let mut chosen = None;
    ui.horizontal(|ui| {
        let parsed = ggv_layout::parse_tag_name(tag);
        if ui.add_enabled(parsed.is_some(), egui::Button::selectable(selected, "tag")).on_hover_text("Any read tag: CB (cell barcode), PS (phase set), UB, MI, an aligner's own...").clicked() {
            chosen = parsed;
        }
        let edit = ui.add(egui::TextEdit::singleline(tag).desired_width(28.0).char_limit(2).hint_text("CB"));
        if edit.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
            chosen = ggv_layout::parse_tag_name(tag);
        }
    });
    chosen
}

/// Grouping choices; "base" groups at `pos` (the outlined or clicked base).
fn group_options(ui: &mut egui::Ui, value: &mut ggv_layout::GroupBy, pos: u64, tag: &mut String) {
    use ggv_layout::GroupBy;
    for g in [GroupBy::None, GroupBy::Strand, GroupBy::Haplotype, GroupBy::ReadGroup, GroupBy::Alleles] {
        ui.selectable_value(value, g, group_name(g));
    }
    ui.selectable_value(value, GroupBy::BaseAt(pos), format!("base at {}", commas(pos + 1)));
    if let Some(t) = tag_row(ui, matches!(value, GroupBy::Tag(_)), tag) {
        *value = GroupBy::Tag(t);
    }
}

/// Menus whose settings you change several at a time stay open until you click elsewhere
/// (or press Escape); an action that opens a window still closes them itself.
fn keep_open() -> egui::containers::menu::MenuConfig {
    egui::containers::menu::MenuConfig::new().close_behavior(egui::PopupCloseBehavior::CloseOnClickOutside)
}

/// The build an annotation file's name declares (`hg38.ncbiRefSeq.gtf.gz`, `gencode.v44.GRCh38...`,
/// `Homo_sapiens.GRCh37.87.gtf`), as GGV's genome id; `None` when the name doesn't say.
fn annotation_build(path: &std::path::Path) -> Option<&'static str> {
    build_named(&crate::files::display_name(path))
}

use ggv_io::assembly::build_named;

/// Minimum mapping quality: one click on a common value, or type any (0-255). A drag box in a
/// menu drew with no frame and was hard to change, so this uses buttons and a text field.
fn mapq_options(ui: &mut egui::Ui, value: &mut u8) {
    ui.label("Hide reads with mapping quality below");
    ui.horizontal(|ui| {
        for q in [0u8, 1, 5, 10, 20, 30, 60] {
            if ui.selectable_label(*value == q, q.to_string()).clicked() {
                *value = q;
            }
        }
    });
    ui.horizontal(|ui| {
        ui.label("Other");
        let id = ui.make_persistent_id("min-mapq-text");
        let mut text = ui.data_mut(|d| d.get_temp::<String>(id)).unwrap_or_else(|| value.to_string());
        let edit = ui.add(egui::TextEdit::singleline(&mut text).desired_width(44.0).hint_text("0-255"));
        if edit.changed()
            && let Ok(v) = text.trim().parse::<u8>()
        {
            *value = v;
        }
        if !edit.has_focus() {
            text = value.to_string();
        }
        ui.data_mut(|d| d.insert_temp(id, text));
    });
    ui.label(egui::RichText::new("MAPQ 0: the read fits another place equally well (repeats, duplicated genes).").weak().small());
}

fn group_label(g: ggv_layout::GroupBy) -> String {
    match g {
        ggv_layout::GroupBy::BaseAt(p) => format!("base @ {}", commas(p + 1)),
        ggv_layout::GroupBy::Tag(t) => format!("tag {}", String::from_utf8_lossy(&t)),
        g => group_name(g).to_owned(),
    }
}

fn sort_options(ui: &mut egui::Ui, value: &mut Option<(u64, ggv_layout::SortBy)>, pos: u64) {
    use ggv_layout::SortBy;
    for by in [SortBy::Base, SortBy::Strand, SortBy::MappingQuality, SortBy::Start, SortBy::InsertSize] {
        if ui.button(format!("Sort by {}", sort_name(by))).clicked() {
            *value = Some((pos, by));
            ui.close();
        }
    }
    if value.is_some() && ui.button("Clear sort").clicked() {
        *value = None;
        ui.close();
    }
}

fn show_tooltip(response: &egui::Response, lines: Vec<String>) {
    response.clone().on_hover_ui_at_pointer(|ui| {
        for l in lines {
            ui.monospace(l);
        }
    });
}

/// The span to load around a view of `vis` bases centred on `center`, and whether it is a
/// reads span: about 3 view widths (at least 20 kb for reads), capped by the read or coverage
/// window. Any view up to the read window shows reads (a typed 1 Mb locus with "Reads up to
/// 1 Mb"); its span may run `SPAN_SLACK` bases over, which the loader allows, and always covers
/// the view, fractional edges included.
/// Bases from `origin` to `vs`, as whole bases and the fraction, for the shader. As one f32,
/// 2.5 Mb out (a reads span loaded zoomed out, kept when zooming in) resolves only a quarter
/// base: tracks drew up to an eighth of a base off the ruler and outline, and the error
/// changed with every pan.
fn split_offset(vs: f64, origin: u64) -> (f32, f32) {
    let d = vs - origin as f64;
    (d.floor() as f32, d.rem_euclid(1.0) as f32)
}

fn load_span(center: f64, vis: f64, contig_len: u64) -> (u64, u64, bool) {
    let vis = vis.max(1.0);
    let reads = vis <= read_window() as f64;
    let cap = if reads { read_window() + SPAN_SLACK } else { COVERAGE_WINDOW } as f64;
    let floor = if reads { MIN_READ_SPAN } else { 0.0 };
    let span = (vis * 3.0).max(floor).max(vis + 2.0).min(cap.max(vis + 2.0)).floor() as u64;
    let start = (center - span as f64 / 2.0).max(0.0).floor() as u64;
    let end = (start + span).min(contig_len).max(start + 1);
    (start, end, reads)
}

/// Reads in one packed row that overlap `[vs, ve)`. Rows hold non-overlapping reads sorted by
/// start, so their ends are sorted too and a binary search finds the first.
/// Whether any packed (filtered) read overlaps `vs..ve`.
fn reads_in_view(reads: &[AlignedRead], packed: &ggv_layout::Packed, vs: f64, ve: f64) -> bool {
    packed.rows.iter().any(|row| visible_in_row(reads, row, vs, ve).next().is_some())
}

fn visible_in_row<'a>(reads: &'a [AlignedRead], row: &'a [usize], vs: f64, ve: f64) -> impl Iterator<Item = &'a AlignedRead> {
    let first = row.partition_point(|&i| (reads[i].end as f64) <= vs);
    row[first..].iter().map(move |&i| &reads[i]).take_while(move |r| (r.start as f64) < ve)
}

/// Where to draw the outline around a base spanning `x0..x1` (points): the physical pixel just
/// left of the coverage bar and the one just right of it, as left edges in points. The bar is
/// what the GPU draws (tracks.wgsl): a pixel is filled when its centre falls inside the quad, and
/// bases 4 points or wider lose 1 point on the right as a gap. Snapping to those pixels keeps the
/// outline equally close on both sides at any zoom and pixel density.
fn outline_columns(x0: f32, x1: f32, pixels_per_point: f32) -> (f32, f32) {
    let ppp = pixels_per_point.max(0.1);
    let right = if x1 - x0 >= 4.0 { x1 - 1.0 } else { x1.max(x0 + 1.0 / ppp) };
    let first = (x0 * ppp - 0.5).ceil();
    let last = ((right * ppp - 0.5).ceil() - 1.0).max(first);
    ((first - 1.0) / ppp, (last + 1.0) / ppp)
}

/// A note on the top alt allele's strand support, judged against the other reads: one-strand
/// support is a classic artifact only where the other reads come from the strand it lacks (at
/// a capture edge every read can be on one strand).
fn strand_note(s: &ggv_query::PositionStats) -> String {
    let Some(alt) = &s.top_alt else { return String::new() };
    if alt.count.total < 4 || (alt.count.fwd > 0 && alt.count.rev > 0) {
        return String::new();
    }
    let b = &s.bases;
    let (all_f, all_r) = [b.a, b.c, b.g, b.t, b.n, s.deletions].iter().fold((0u32, 0u32), |(f, r), c| (f + c.fwd, r + c.rev));
    let (other_f, other_r) = (all_f.saturating_sub(alt.count.fwd), all_r.saturating_sub(alt.count.rev));
    let lacking = if alt.count.fwd == 0 { other_f } else { other_r };
    let other = other_f + other_r;
    if other >= 4 && f64::from(lacking) >= 0.2 * f64::from(other) {
        format!("  (one strand only while other reads are on both: suspicious)")
    } else if lacking == 0 {
        "  (every read here is on one strand: strand bias can't be judged)".to_owned()
    } else {
        String::new()
    }
}

/// The evidence for one file at one base, as text to paste.
fn evidence_text(name: &str, locus: &str, s: &ggv_query::PositionStats) -> String {
    let depth = f64::from(s.depth.max(1));
    let mut out = format!("{name}\n{locus}  depth {}  ref {}\nallele\treads\t%\tfwd\trev\n", s.depth, s.ref_base.unwrap_or('.'));
    let b = &s.bases;
    for (n, c) in [("A", b.a), ("C", b.c), ("G", b.g), ("T", b.t), ("N", b.n), ("DEL", s.deletions), ("INS after", s.insertions_after)] {
        if c.total > 0 {
            out.push_str(&format!("{n}\t{}\t{:.1}\t{}\t{}\n", c.total, 100.0 * f64::from(c.total) / depth, c.fwd, c.rev));
        }
    }
    out.push_str(&format!("MAPQ mean {:.1} ({} at 0), base Q mean {:.1}\n", s.mean_mapping_quality, s.mapping_quality_zero, s.mean_base_quality));
    if let Some(alt) = &s.top_alt {
        out.push_str(&format!("top alt {} at {:.1}% ({} fwd / {} rev){}\n", alt.allele, alt.fraction * 100.0, alt.count.fwd, alt.count.rev, strand_note(s)));
    }
    out
}

#[cfg(test)]
mod tests {
    /// Reads loaded around the view for panning, but none on screen: the track says "No reads
    /// in view" and shows no row count or scrollbar for them.
    #[test]
    fn reads_loaded_beside_the_view_are_not_in_view() {
        let data = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio");
        let Ok(mut bam) = ggv_io::BamFile::open(data.join("HG002.revio.slices.bam")) else { return };
        let region = ggv_model::Region { contig: "chr20".into(), start: 31_165_000, end: 31_184_000 };
        let view = ggv_layout::load_view(&mut bam, &region, None, String::new(), &ggv_layout::ViewOptions { max_rows: 0, ..Default::default() }).unwrap();
        let (reads, _) = view.reads.as_ref().unwrap();
        // Packed as the app packs them (the loader leaves that to it).
        let packed = &ggv_layout::pack(reads, 10_000, 2);
        assert!(!reads.is_empty());
        assert!(super::reads_in_view(reads, packed, 31_170_000.0, 31_170_100.0));
        // Past the last read's end, still inside what a pan might load.
        let end = reads.iter().map(|r| r.end).max().unwrap() as f64;
        assert!(!super::reads_in_view(reads, packed, end + 10.0, end + 110.0));
    }

    #[test]
    fn annotation_builds_from_file_names() {
        use std::path::Path;
        assert_eq!(super::annotation_build(Path::new("C:/x/hg38.ncbiRefSeq.gtf.gz")), Some("hg38"));
        assert_eq!(super::annotation_build(Path::new("/g/hg19/tracks/hg19.ncbiRefSeq.gtf.gz")), Some("hg19"));
        assert_eq!(super::annotation_build(Path::new("Homo_sapiens.GRCh37.87.gtf")), Some("hg19"));
        assert_eq!(super::annotation_build(Path::new("gencode.v44.basic.annotation.GRCh38.gff3")), Some("hg38"));
        assert_eq!(super::annotation_build(Path::new("genes.gtf")), None);
        assert_eq!(super::annotation_build(Path::new("hg19ToHg38.chain.bed")), None);
        assert_eq!(super::annotation_build(Path::new("lift.hg19.hg38.bed")), None);
    }

    use super::*;

    #[test]
    fn sv_sizes() {
        assert_eq!([sv_size(292), sv_size(1411), sv_size(3039), sv_size(71_507), sv_size(2_400_000)], ["292 bp", "1.4 kb", "3.0 kb", "71.5 kb", "2.4 Mb"]);
    }

    #[test]
    fn outline_hugs_the_coverage_bar_on_both_sides() {
        // 14 points per base starting at 692.6 (as in the reproduction): the bar fills pixels
        // 693..=705 (gap at 706); the outline must sit at 692 and 706.
        assert_eq!(outline_columns(692.6, 706.6, 1.0), (692.0, 706.0));
        // Pixel-aligned base.
        assert_eq!(outline_columns(100.0, 110.0, 1.0), (99.0, 109.0));
        // Narrow base (no gap): the lines straddle it.
        assert_eq!(outline_columns(100.0, 102.0, 1.0), (99.0, 102.0));
        // 1.5x display: the bar spans 150..163.5 physical pixels, so it fills 150..=162; lines at
        // 149 and 163.
        let (l, r) = outline_columns(100.0, 110.0, 1.5);
        assert_eq!((l * 1.5, r * 1.5), (149.0, 163.0));
    }

    #[test]
    fn reads_spans_fit_the_read_window_and_cover_the_view() {
        let contig_len = 248_956_422;
        let mut checked = 0;
        let w = read_window() as f64;
        for vis in [101.0, 20_000.0, w / 3.0, w / 3.0 + 1.0, w / 2.0 + 0.5, w - 2.1, w - 2.0, w - 1.9, w, w * 1.5] {
            for center in [43_092_919.5, 43_092_919.123, 1_000_000.0] {
                let (start, end, reads) = load_span(center, vis, contig_len);
                let (vs, ve) = (center - vis / 2.0, center + vis / 2.0);
                if reads {
                    assert!(end - start <= read_window() + SPAN_SLACK, "vis {vis}: reads span {} > window", end - start);
                }
                assert!(start as f64 <= vs && end as f64 >= ve, "vis {vis} centre {center}: [{start}, {end}) misses [{vs}, {ve})");
                checked += 1;
            }
        }
        assert_eq!(checked, 30);
        // The boundary itself: the widest view that still loads reads.
        assert!(load_span(5e7, w, contig_len).2);
        assert!(!load_span(5e7, w + 0.1, contig_len).2);
    }

    #[test]
    fn far_origins_stay_on_the_ruler() {
        // Where the shader puts base `pos`, in bases from the view's left edge, against f64.
        for origin in [0u64, 56_354_246, 58_844_246] {
            for vs in [58_854_221.0, 58_854_221.37, 58_854_221.5, 58_854_221.937, 58_854_222.0625] {
                let (whole, frac) = split_offset(vs, origin);
                for pos in [58_854_221u64, 58_854_245, 58_854_260] {
                    let x = (pos - origin) as f32;
                    let shader = f64::from((x - whole) - frac);
                    let exact = pos as f64 - vs;
                    // Positions past 2^24 bases from the origin aren't exact in f32 at all; the
                    // app keeps spans far shorter than that.
                    if pos - origin < 1 << 24 {
                        assert!((shader - exact).abs() < 1e-3, "origin {origin} vs {vs} pos {pos}: {shader} vs {exact}");
                    }
                }
            }
        }
        // The single-f32 offset this replaced: an eighth of a base off at 2.5 Mb.
        let old = f64::from((58_854_221.37f64 - 56_354_246.0) as f32);
        assert!((old - 2_499_975.37).abs() > 0.1);
    }

    #[test]
    fn spans_clip_to_the_contig() {
        let (start, end, _) = load_span(100.0, 1000.0, 16_569);
        assert_eq!(start, 0);
        assert!(end <= 16_569);
        let (_, end, _) = load_span(16_500.0, 1000.0, 16_569);
        assert_eq!(end, 16_569);
    }
}

/// Requests from the command line that wait until the genome is known.
struct Startup {
    queue: Option<PathBuf>,
    find_seq: Option<String>,
    find_read: Option<String>,
    compare: bool,
    report: Option<PathBuf>,
    measure: Option<String>,
}

impl App {
    fn run_startup(&mut self, ctx: &egui::Context) {
        if self.contigs.is_empty() || self.pending_goto.is_some() || self.startup.is_none() {
            return;
        }
        // Wait for the first loads, so on-screen searches and reports see data.
        if self.tracks.iter().any(|t| t.current.is_none() && !t.missing_contig && t.error.is_none()) {
            return;
        }
        let st = self.startup.take().expect("checked");
        if let Some(q) = st.queue {
            self.open_queue_file(&q);
        }
        if let Some(text) = st.find_seq {
            let mut s = tools::SeqUi { text, ..Default::default() };
            if let Some(m) = s.text.split_once('/').map(|(a, b)| (a.to_owned(), b.parse::<u8>().unwrap_or(0))) {
                s.text = m.0;
                s.mismatches = m.1;
            }
            self.seq = Some(s);
            self.run_seq_now(ctx);
        }
        if let Some(name) = st.find_read {
            self.find = Some(tools::FindUi { text: name, ..Default::default() });
            self.run_find_now(ctx);
        }
        if st.compare {
            self.open_compare();
        }
        if let Some(m) = st.measure
            && let Ok(r) = ggv_io::resolve_locus_str(&self.contigs, &m)
            && let Some(c) = self.contigs.iter().position(|x| x.name == r.contig)
        {
            self.measure = Some((c, r.start, r.end));
        }
        if let Some(p) = st.report {
            self.report_to(p);
        }
    }
}

struct IndexJob {
    path: PathBuf,
    state: IndexState,
}

enum IndexState {
    Offered,
    Building(std::sync::mpsc::Receiver<Result<(), String>>),
    Failed(String),
    Done,
}

/// How the changed-settings chips show in the toolbar.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ChipFit {
    /// One chip each, with its ×.
    Each = 0,
    /// One menu, "3 settings changed".
    Folded = 1,
    /// One menu, "3 changed".
    Short = 2,
}

/// The toolbar's layout for a row: the locus box's width and what folds away.
#[derive(Clone, Copy, Debug, PartialEq)]
struct ToolbarFit {
    locus_w: f32,
    chips: ChipFit,
    /// Split, Evidence and Commands in a More menu.
    more: bool,
    /// The Go button (Enter in the locus box does the same).
    go: bool,
}

/// Fits the toolbar into a `row` wide row: `fixed` for the buttons and menus that always show,
/// `go` for the Go button, `chips` for the chips each / folded / short, `tail` for Split,
/// Evidence and Commands and for More. The locus box gets at least `locus_need` (the whole
/// locus) and the room left, up to 520; the chips fold, then the tail, then Go goes (Enter
/// still goes), until the row fits. A row that can't fit even then wraps, the locus still whole.
fn fit_toolbar(row: f32, fixed: f32, go: f32, locus_need: f32, chips: [f32; 3], tail: [f32; 2]) -> ToolbarFit {
    const SLACK: f32 = 8.0;
    const LOCUS_MAX: f32 = 520.0;
    let options = [
        (ChipFit::Each, false, true),
        (ChipFit::Folded, false, true),
        (ChipFit::Folded, true, true),
        (ChipFit::Short, true, true),
        (ChipFit::Short, true, false),
    ];
    let rest = |(c, more, with_go): (ChipFit, bool, bool)| fixed + if with_go { go } else { 0.0 } + chips[c as usize] + tail[usize::from(more)] + SLACK;
    let last = options[options.len() - 1];
    let o = options.into_iter().find(|&o| rest(o) + locus_need <= row).unwrap_or(last);
    let locus_w = (row - rest(o)).clamp(locus_need, LOCUS_MAX.max(locus_need));
    ToolbarFit { locus_w, chips: o.0, more: o.1, go: o.2 }
}

/// What a click on a status-bar message did.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum MessageClick {
    Dismiss,
    Copy,
}

/// A message as the status bar shows it: on one line, each line's leading file path cut to the
/// file's name, so the reason after it fits a narrow bar. Hover and copy have it whole.
fn status_line(text: &str) -> String {
    let short = |line: &str| -> String {
        // "C:\data\x.bam: why": the path ends at the first ": " (a drive's colon has no space).
        if let Some((head, rest)) = line.split_once(": ")
            && (head.contains('/') || head.contains('\\'))
            && let Some(name) = std::path::Path::new(head).file_name()
        {
            return format!("{}: {rest}", name.to_string_lossy());
        }
        line.to_owned()
    };
    text.lines().map(short).collect::<Vec<_>>().join("  ·  ")
}

/// One status-bar message: its × first (so it stays on screen), then the text on one line, cut
/// with "…" to `max_w`; hovering shows the whole text and clicking it copies it.
fn status_message(ui: &mut egui::Ui, text: &str, color: Color32, max_w: f32) -> Option<MessageClick> {
    let mut click = None;
    if ui.small_button("×").on_hover_text("Dismiss").clicked() {
        click = Some(MessageClick::Dismiss);
    }
    let line = status_line(text);
    let label = egui::Label::new(egui::RichText::new(line).color(color)).truncate().show_tooltip_when_elided(false).sense(Sense::click());
    let r = ui
        .scope(|ui| {
            ui.set_max_width(max_w);
            ui.add(label)
        })
        .inner
        .on_hover_cursor(egui::CursorIcon::PointingHand)
        .on_hover_ui(|ui| {
            ui.set_max_width(640.0);
            ui.label(text);
            ui.label(egui::RichText::new("Click to copy").weak());
        });
    if r.clicked() {
        ui.ctx().copy_text(text.to_owned());
        click = Some(MessageClick::Copy);
    }
    ui.separator();
    click
}

/// The error text colour: a darker red on light panels, a light red on dark ones (both at least
/// 4.5:1 against the panel).
fn error_rgb(dark: bool) -> Color32 {
    if dark { Color32::from_rgb(255, 120, 110) } else { Color32::from_rgb(180, 30, 30) }
}

/// The warning text colour, amber in both themes.
fn warning_rgb(dark: bool) -> Color32 {
    if dark { Color32::from_rgb(240, 175, 60) } else { Color32::from_rgb(160, 95, 0) }
}

/// The tracks on screen by kind: "1 alignment track", or "4 tracks (1 alignment, 2 variant,
/// 1 gene)". Counting alignments alone said "1 track" with four on screen.
fn tracks_text(alignments: usize, variants: usize, genes: usize, signals: usize) -> String {
    let total = alignments + variants + genes + signals;
    let kinds: Vec<String> = [(alignments, "alignment"), (variants, "variant"), (genes, "gene"), (signals, "signal")]
        .into_iter()
        .filter(|(n, _)| *n > 0)
        .map(|(n, kind)| format!("{n} {kind}"))
        .collect();
    let plural = if total == 1 { "" } else { "s" };
    match kinds.len() {
        0 => "no tracks".to_owned(),
        1 => format!("{} track{plural}", kinds[0]),
        _ => format!("{total} tracks ({})", kinds.join(", ")),
    }
}

/// `e` about the file at `path`, naming the file once: ggv-io's errors already start with the
/// path, and prefixing it again pushed the reason off the status bar.
pub(crate) fn path_error(path: &std::path::Path, e: impl std::fmt::Display) -> String {
    let e = e.to_string();
    let shown = path.display().to_string();
    if e.contains(&shown) { e } else { format!("{shown}: {e}") }
}

/// Whether a file GGV reads through an index has none next to it. A file that isn't there has
/// nothing to index: opening it says it wasn't found instead.
fn missing_index(path: &std::path::Path) -> bool {
    if !path.exists() {
        return false;
    }
    let exists = |ext: &str| path.with_added_extension(ext).exists() || path.with_extension(ext).exists();
    let name = path.to_string_lossy().to_ascii_lowercase();
    if name.ends_with(".bam") {
        !(exists("bai") || exists("csi"))
    } else if name.ends_with(".cram") {
        !exists("crai")
    } else if name.ends_with(".vcf.gz") || name.ends_with(".vcf.bgz") {
        !(exists("tbi") || exists("csi"))
    } else if name.ends_with(".bcf") {
        !exists("csi")
    } else if ggv_io::index_kind(path) == Some("fai") {
        // A bgzipped FASTA is read through its .gzi too; Build .fai writes both.
        !exists("fai") || ((name.ends_with(".gz") || name.ends_with(".bgz")) && !exists("gzi"))
    } else {
        false
    }
}

/// Overlay colours for the light and dark themes.
#[derive(Clone, Copy)]
struct Pal {
    bg: Color32,
    text: Color32,
    muted: Color32,
    icon: Color32,
    header: Color32,
    rule: Color32,
    rule_strong: Color32,
    label_bg: Color32,
    outline: Color32,
}

impl Pal {
    fn new(dark: bool) -> Self {
        if dark {
            Self {
                bg: Color32::from_gray(24),
                text: Color32::from_gray(225),
                muted: Color32::from_gray(150),
                icon: Color32::from_gray(170),
                header: Color32::from_gray(42),
                rule: Color32::from_gray(60),
                rule_strong: Color32::from_gray(110),
                label_bg: Color32::from_black_alpha(200),
                outline: Color32::from_white_alpha(190),
            }
        } else {
            Self {
                bg: Color32::WHITE,
                text: Color32::from_gray(30),
                muted: Color32::from_gray(120),
                icon: Color32::from_gray(90),
                header: Color32::from_gray(246),
                rule: Color32::from_gray(225),
                rule_strong: Color32::from_gray(160),
                label_bg: Color32::from_white_alpha(220),
                outline: Color32::from_black_alpha(170),
            }
        }
    }
}

/// A signal value for display: up to 4 significant digits.
fn fmt_value(v: f64) -> String {
    if v == 0.0 || !v.is_finite() {
        return format!("{v}");
    }
    let digits = (3 - v.abs().log10().floor() as i32).clamp(0, 6) as usize;
    format!("{v:.digits$}")
}

