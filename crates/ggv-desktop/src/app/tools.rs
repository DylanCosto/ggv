//! The tools around the view: review queue, bookmarks, filter expressions, read and sequence
//! search, track comparison, the measure tool, copy-as, evidence reports and file reloading.
//! Anything that touches the disk runs on a background thread and reports back through
//! [`App::spawn_job`], so the UI never waits.

use std::{
    collections::BTreeMap,
    path::{Path, PathBuf},
    sync::{Arc, Mutex},
    time::{Duration, Instant},
};

use eframe::egui::{self, Align2, Color32, FontId, Pos2, Rect, Stroke, Vec2};
use ggv_layout::{AlignedRead, Coverage};
use ggv_model::Region;

use super::{App, Layout};
use crate::{
    bookmarks::{self, Bookmark},
    diff, expr, finder, queue, seqsearch,
    tracks::{COVERAGE_TOP, Style},
    util::commas,
};

/// Work finished on another thread, applied to the app on the UI thread.
pub type Apply = Box<dyn FnOnce(&mut App, &egui::Context) + Send>;

pub enum Capture {
    Png(PathBuf),
    Clipboard,
    Report(PathBuf),
}

pub struct CaptureJob {
    what: Capture,
    rect: Rect,
    slot: Arc<Mutex<Option<Arc<egui::ColorImage>>>>,
}

#[derive(Default)]
pub struct FindUi {
    pub text: String,
    pub search: Option<finder::Search>,
    /// Found among the loaded reads, without a scan.
    pub loaded_hit: Option<String>,
    pub(super) focused: bool,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Scope {
    View,
    Contig,
    Genome,
}

pub struct SeqUi {
    pub text: String,
    pub mismatches: u8,
    pub scope: Scope,
    pub search: Option<seqsearch::Search>,
    pub error: Option<String>,
    pub(super) focused: bool,
}

impl Default for SeqUi {
    fn default() -> Self {
        Self { text: String::new(), mismatches: 0, scope: Scope::Genome, search: None, error: None, focused: false }
    }
}

pub struct FilterUi {
    pub text: String,
    pub error: Option<String>,
    /// A track's file (its own filter), or `None` for every track.
    pub target: Option<PathBuf>,
    /// The field has had focus once (it takes it when the window opens).
    pub focused: bool,
}

pub struct CompareUi {
    pub a: usize,
    pub b: usize,
    /// The compared files: `a` and `b` follow them when tracks move or close.
    pub files: (PathBuf, PathBuf),
    pub settings: diff::Settings,
    /// Differences in what both tracks have loaded, and what they were computed from.
    cache: Option<((u64, u64, u64, u64, u64), Arc<Vec<diff::Site>>)>,
    pub scan: Option<diff::Scan>,
    /// Where the running scan started and with what, to restart it when that changes.
    job: Option<ScanJob>,
    pub status: Option<String>,
}

/// What a comparison scan reads: both files, their filters, the reference and the settings.
type ScanKey = (PathBuf, PathBuf, crate::tracks::ReadFilter, crate::tracks::ReadFilter, Option<PathBuf>, diff::Settings);

struct ScanJob {
    forward: bool,
    contig: String,
    from: u64,
    key: ScanKey,
}

/// Measure-tool numbers for a span, from what is loaded.
struct Measured {
    len: u64,
    gc: Option<f64>,
    tracks: Vec<(String, Option<f64>, usize)>,
}

fn config_file(name: &str) -> Option<PathBuf> {
    // Through store's folder, so tests that set one keep off the person's settings.
    crate::store::config_dir().map(|d| d.join(name))
}

/// App-wide preferences (`prefs.json`): a JSON object, read and written whole.
pub fn load_prefs() -> serde_json::Map<String, serde_json::Value> {
    config_file("prefs.json").and_then(|p| std::fs::read_to_string(p).ok()).and_then(|t| serde_json::from_str(&t).ok()).unwrap_or_default()
}

pub fn save_prefs(prefs: &serde_json::Map<String, serde_json::Value>) {
    save_json("prefs.json", prefs);
}

pub fn load_saved_filters() -> Vec<String> {
    config_file("filters.json").and_then(|p| std::fs::read_to_string(p).ok()).and_then(|t| serde_json::from_str(&t).ok()).unwrap_or_default()
}

pub fn load_file_filters() -> BTreeMap<PathBuf, String> {
    config_file("file-filters.json").and_then(|p| std::fs::read_to_string(p).ok()).and_then(|t| serde_json::from_str(&t).ok()).unwrap_or_default()
}

/// Forgets the filters remembered for `paths`.
pub fn forget_file_filters(paths: impl IntoIterator<Item = PathBuf>) {
    let mut saved = load_file_filters();
    let before = saved.len();
    for p in paths {
        saved.remove(&p);
    }
    if saved.len() != before {
        save_json("file-filters.json", &saved);
    }
}

fn save_json<T: serde::Serialize>(name: &str, value: &T) {
    if let Some(p) = config_file(name) {
        let _ = crate::store::write_atomic(&p, serde_json::to_string_pretty(value).unwrap_or_default().as_bytes());
    }
}

fn open_in_system(path: &Path) {
    let p = path.as_os_str();
    let _ = if cfg!(windows) {
        std::process::Command::new("explorer").arg(p).spawn()
    } else if cfg!(target_os = "macos") {
        std::process::Command::new("open").arg(p).spawn()
    } else {
        std::process::Command::new("xdg-open").arg(p).spawn()
    };
}

fn rgb([r, g, b]: [u8; 3]) -> Color32 {
    Color32::from_rgb(r, g, b)
}

impl App {
    pub(super) fn notify(&mut self, msg: impl Into<String>) {
        self.notice = Some((msg.into(), Instant::now()));
    }

    /// Runs `work` on a new thread; what it returns is applied on the UI thread.
    pub(super) fn spawn_job(&self, ctx: &egui::Context, work: impl FnOnce() -> Apply + Send + 'static) {
        let tx = self.jobs.0.clone();
        let ctx = ctx.clone();
        let _ = std::thread::Builder::new().name("ggv-job".into()).spawn(move || {
            let _ = tx.send(work());
            ctx.request_repaint();
        });
    }

    pub(super) fn poll_jobs(&mut self, ctx: &egui::Context) {
        let done: Vec<Apply> = self.jobs.1.try_iter().collect();
        for f in done {
            f(self, ctx);
        }
    }

    /// The selection, else the outlined base, else the view: contig index, start, end.
    pub(super) fn place(&self) -> Option<(usize, u64, u64)> {
        if self.contigs.is_empty() {
            return None;
        }
        if let Some(s) = self.selection {
            return Some(s);
        }
        if let Some((c, p)) = self.highlight.filter(|(c, _)| *c == self.contig) {
            return Some((c, p, p + 1));
        }
        let (vs, ve) = self.visible();
        Some((self.contig, vs.max(0.0) as u64, (ve.min(self.contig_len())) as u64))
    }

    // ---- Capture: PNG, clipboard, report ----

    pub(super) fn capture(&mut self, ctx: &egui::Context, what: Capture) {
        let slot: Arc<Mutex<Option<Arc<egui::ColorImage>>>> = Arc::default();
        let (s2, c2) = (slot.clone(), ctx.clone());
        ctx.request_screenshot(move |image| {
            *s2.lock().unwrap() = Some(image);
            c2.request_repaint();
        });
        self.capture = Some(CaptureJob { what, rect: self.canvas_rect, slot });
    }

    pub(super) fn finish_capture(&mut self, ctx: &egui::Context) {
        let Some(job) = &self.capture else { return };
        let Some(image) = job.slot.lock().unwrap().take() else { return };
        let job = self.capture.take().expect("checked");
        let cropped = crate::automation::crop(&image, &job.rect, ctx.pixels_per_point());
        match job.what {
            Capture::Png(path) => match crate::automation::encode_png(&cropped).and_then(|png| std::fs::write(&path, png).map_err(|e| e.to_string())) {
                Ok(()) => self.notify(format!("Saved {}", path.display())),
                Err(e) => self.error = Some(format!("export: {e}")),
            },
            Capture::Clipboard => {
                ctx.copy_image(cropped);
                self.notify("Copied the view to the clipboard");
            }
            Capture::Report(path) => {
                let png = crate::automation::encode_png(&cropped).ok();
                let html = self.build_report(png).html();
                match std::fs::write(&path, html) {
                    Ok(()) => {
                        self.last_report = Some(path.clone());
                        self.notify(format!("Report saved to {}", crate::files::display_name(&path)));
                    }
                    Err(e) => {
                        self.error = Some(format!("report: {e}"));
                        if self.is_scripted() {
                            crate::automation::fail(&format!("report {}: {e}", path.display()));
                        }
                    }
                }
            }
        }
    }

    // ---- Evidence report ----

    pub(super) fn start_report(&mut self) {
        let Some((c, pos)) = self.evidence_base() else { return };
        let name = format!("ggv-{}-{}.html", self.contigs[c].name, pos + 1);
        self.ask(rfd::FileDialog::new().add_filter("HTML", &["html"]).set_file_name(name), super::dialogs::Kind::Save, super::dialogs::Then::Report);
    }

    /// Writes a report to `path` once the evidence at the base has loaded.
    pub(super) fn report_to(&mut self, path: PathBuf) {
        let Some((c, pos)) = self.evidence_base() else { return };
        let contig = self.contigs[c].name.clone();
        let filters: Vec<crate::tracks::ReadFilter> = self.tracks.iter().map(|t| self.track_style(&t.loader.bam_path).filter).collect();
        for (t, filter) in self.tracks.iter_mut().zip(filters) {
            t.request_evidence(&contig, pos, filter, &mut self.next_request_id);
        }
        self.pending_report = Some(path);
    }

    pub(super) fn poll_report(&mut self, ctx: &egui::Context) {
        if self.pending_report.is_none() || self.capture.is_some() {
            return;
        }
        let Some((c, pos)) = self.evidence_base() else { return };
        let contig = &self.contigs[c].name;
        let ready = self.tracks.iter().all(|t| {
            t.missing_contig || t.evidence.as_ref().is_some_and(|e| e.pos == pos && t.contig_name(contig).is_some_and(|x| x.name == e.contig))
        }) && !self.tracks.iter().any(crate::track::Track::is_loading);
        if ready {
            let path = self.pending_report.take().expect("checked");
            self.capture(ctx, Capture::Report(path));
        }
    }

    fn build_report(&self, png: Option<Vec<u8>>) -> ggv_render::report::Report {
        use ggv_render::report::{Report, Section};
        let (c, pos) = self.evidence_base().unwrap_or((self.contig, 0));
        let contig = self.contigs.get(c).map(|x| x.name.clone()).unwrap_or_default();
        let mut title = format!("{contig}:{}", commas(pos + 1));
        let mut sections = Vec::new();

        // Variant records at the base.
        for v in &self.variants {
            let Some(cur) = &v.current else { continue };
            for var in cur.variants.iter().filter(|x| x.start <= pos && pos < x.end) {
                if !title.contains('>') {
                    title = format!("{title}  {} > {}", super::allele(&var.reference), super::alleles(&var.alts.join(",")));
                }
                let mut s = Section { heading: format!("Variant record · {}", v.file_name()), ..Default::default() };
                s.lines.push(format!(
                    "{}  {} > {}  {}  QUAL {}  FILTER {}  AF {}",
                    if var.id.is_empty() || var.id == "." { "-" } else { &var.id },
                    super::allele(&var.reference),
                    super::alleles(&var.alts.join(",")),
                    var.kind(),
                    var.qual.map_or(".".to_owned(), |q| format!("{q:.1}")),
                    if var.filters.is_empty() { ".".to_owned() } else { var.filters.join(",") },
                    var.af_text()
                ));
                if !var.samples.is_empty() {
                    s.table.push(["sample", "GT", "DP", "GQ", "AD"].map(String::from).to_vec());
                    for (name, call) in v.samples.iter().zip(&var.samples) {
                        s.table.push(vec![
                            name.clone(),
                            call.gt.clone(),
                            call.dp.map_or(".".into(), |d| d.to_string()),
                            call.gq.map_or(".".into(), |d| d.to_string()),
                            call.ad_text().unwrap_or_else(|| ".".into()),
                        ]);
                    }
                }
                sections.push(s);
            }
        }

        let lines = self.consequence_lines(pos);
        if !lines.is_empty() {
            sections.push(Section { heading: "Consequences".into(), lines, ..Default::default() });
        }

        // Pileup at the base, per track.
        for t in &self.tracks {
            let mut s = Section { heading: t.title(), ..Default::default() };
            match t.evidence.as_ref().filter(|e| e.pos == pos) {
                Some(ev) if ev.stats.is_some() => {
                    let st = ev.stats.as_ref().expect("checked");
                    let depth = f64::from(st.depth.max(1));
                    s.lines.push(format!("depth {}  ·  reference {}", st.depth, st.ref_base.unwrap_or('.')));
                    s.table.push(["allele", "reads", "%", "forward", "reverse"].map(String::from).to_vec());
                    let b = &st.bases;
                    for (name, n) in [("A", b.a), ("C", b.c), ("G", b.g), ("T", b.t), ("N", b.n), ("DEL", st.deletions), ("INS after", st.insertions_after)] {
                        if n.total > 0 {
                            let is_ref = st.ref_base.is_some_and(|r| r.to_string() == name);
                            s.table.push(vec![
                                if is_ref { format!("{name} (ref)") } else { name.to_owned() },
                                n.total.to_string(),
                                format!("{:.1}", 100.0 * f64::from(n.total) / depth),
                                n.fwd.to_string(),
                                n.rev.to_string(),
                            ]);
                        }
                    }
                    s.lines.push(format!(
                        "MAPQ mean {:.1} ({} at 0)  ·  base quality mean {:.1}  ·  soft clips starting {} / ending {}",
                        st.mean_mapping_quality, st.mapping_quality_zero, st.mean_base_quality, st.soft_clipped_left, st.soft_clipped_right
                    ));
                    if let Some(alt) = &st.top_alt {
                        s.lines.push(format!("top alternate {} at {:.1}% ({} forward / {} reverse)", alt.allele, alt.fraction * 100.0, alt.count.fwd, alt.count.rev));
                        if alt.count.total >= 4 && (alt.count.fwd == 0 || alt.count.rev == 0) {
                            s.warning = Some("The alternate allele is supported by one strand only — a common artifact.".into());
                        }
                    }
                }
                Some(ev) => s.lines.push(ev.error.clone().unwrap_or_else(|| "no reads here".into())),
                None => s.lines.push("not in this file".into()),
            }
            if let Some(Ok(sum)) = t.summary.as_ref().and_then(|(_, r)| r.as_ref())
                && let Some((sc, a, b)) = self.selection
            {
                s.lines.push(format!(
                    "selected region {}:{}-{}: mean depth {:.1}, median {}, {} variable sites",
                    self.contigs[sc].name,
                    commas(a + 1),
                    commas(b),
                    sum.mean_depth,
                    sum.median_depth,
                    sum.variable_site_count
                ));
            }
            sections.push(s);
        }

        // Notes from review, the queue and bookmarks.
        let mut notes = Section { heading: "Review notes".into(), ..Default::default() };
        if let Some(r) = &self.review
            && let Some((key, _)) = self.current_site()
            && let Some(d) = r.decisions.get(&key)
        {
            notes.lines.push(format!("variant review: {}{}", d.label.name(), if d.note.is_empty() { String::new() } else { format!(" — {}", d.note) }));
        }
        if let Some(q) = &self.queue
            && let Some(item) = q.current()
            && let Some(d) = q.decision(q.pos)
        {
            notes.lines.push(format!("queue ({}): {} {}{}", q.source, item.key(), d.verdict.name(), if d.note.is_empty() { String::new() } else { format!(" — {}", d.note) }));
        }
        for b in self.bookmarks.borrow().iter().filter(|b| b.contig == contig && b.start <= pos && pos < b.end) {
            notes.lines.push(format!("bookmark {}{}", b.locus(), if b.note.is_empty() { String::new() } else { format!(": {}", b.note) }));
        }
        if notes.lines.is_empty() {
            notes.lines.push("none".into());
        }
        sections.push(notes);

        // How this was made, so a collaborator can tell which inputs and settings produced it.
        let mut how = Section { heading: "Provenance".into(), ..Default::default() };
        how.lines.push(super::about::build_line());
        if let Some(r) = &self.reference_path {
            how.lines.push(format!("reference: {}", crate::files::display_name(r)));
        }
        for t in &self.tracks {
            let f = self.track_style(&t.loader.bam_path).filter;
            how.lines.push(format!(
                "{}: exclude flags 0x{:x}, min MAPQ {}{}",
                t.full_title(),
                f.exclude_flags(),
                f.min_mapq,
                f.expr.map_or(String::new(), |e| format!(", filter {}", crate::expr::text(e)))
            ));
        }
        for v in &self.variants {
            how.lines.push(format!("variants: {}", v.file_name()));
        }
        how.lines.push("counts: exact (every read at the base is counted, whatever the view draws)".into());
        sections.push(how);

        let files: Vec<String> = self.tracks.iter().map(|t| crate::files::display_name(&t.loader.bam_path)).collect();
        let now = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(0, |d| d.as_secs());
        let mut subtitle = format!("{} · view {}", ggv_render::report::utc_stamp(now), self.locus_label().unwrap_or_default());
        let exprs: std::collections::BTreeSet<String> = self.tracks.iter().filter_map(|t| self.track_style(&t.loader.bam_path).filter.expr).map(crate::expr::text).collect();
        if !exprs.is_empty() {
            subtitle.push_str(&format!(" · reads filtered by {}", exprs.into_iter().collect::<Vec<_>>().join(" / ")));
        }
        if !files.is_empty() {
            subtitle.push_str(&format!(" · {}", files.join(", ")));
        }
        Report { title, subtitle, png, sections }
    }

    // ---- Copy as ----

    pub(super) fn copy_fasta(&mut self, ctx: &egui::Context) {
        let Some((c, start, end)) = self.place() else { return };
        let Some(reference) = self.reference_path.clone() else {
            self.error = Some("Copy as FASTA needs a reference (Reference… or drop a FASTA)".into());
            return;
        };
        if end - start > 10_000_000 {
            self.error = Some("Select at most 10 Mb to copy as FASTA".into());
            return;
        }
        let region = Region { contig: self.contigs[c].name.clone(), start, end };
        self.spawn_job(ctx, move || {
            let result = ggv_io::ReferenceFasta::open(&reference).and_then(|mut f| f.fetch(&region));
            Box::new(move |app: &mut App, ctx: &egui::Context| match result {
                Ok(seq) => {
                    let mut out = format!(">{}:{}-{}\n", region.contig, region.start + 1, region.end);
                    for line in seq.chunks(60) {
                        out.push_str(&String::from_utf8_lossy(line));
                        out.push('\n');
                    }
                    ctx.copy_text(out);
                    app.notify(format!("Copied {} bp as FASTA", commas(seq.len() as u64)));
                }
                Err(e) => app.error = Some(format!("FASTA: {e}")),
            })
        });
    }

    /// Copies track `track`'s own sequence over the selection (or the view) as FASTA: at each
    /// base what most of its reads show, IUPAC codes at heterozygous bases, `N` where under 3
    /// reads cover it, deleted bases left out (see `consensus`). The track's read filter applies,
    /// so filtering by `hp == 1` gives one haplotype's sequence.
    pub(super) fn copy_sample_sequence(&mut self, ctx: &egui::Context, track: usize) {
        let Some((c, start, end)) = self.place() else { return };
        let Some(t) = self.tracks.get(track) else { return };
        if end - start > 1_000_000 {
            self.error = Some("Select at most 1 Mb to copy a sample's sequence".into());
            return;
        }
        let view_contig = self.contigs[c].name.clone();
        let Some(local) = t.contig_name(&view_contig).map(|x| x.name.clone()) else {
            self.error = Some(format!("{} has no reads on {view_contig}", t.title()));
            return;
        };
        let filter = self.track_style(&t.loader.bam_path).filter;
        let (path, reference, title) = (t.loader.bam_path.clone(), self.reference_path.clone(), t.title());
        let header_filter = filter.expr.map(|e| format!(" filter={}", crate::expr::text(e))).unwrap_or_default();
        let region = Region { contig: local, start, end };
        self.spawn_job(ctx, move || {
            let result = (|| -> Result<crate::consensus::Made, String> {
                let mut bam = ggv_io::BamFile::open_with_reference(&path, reference.as_ref()).map_err(|e| e.to_string())?;
                let seq = match &reference {
                    Some(r) => Some(ggv_io::ReferenceFasta::open(r).and_then(|mut f| f.fetch(&region)).map_err(|e| e.to_string())?),
                    None => None,
                };
                let mut stats = Vec::new();
                ggv_query::pileup(&mut bam, &region, seq.as_deref(), &filter.pileup_options_with(reference.as_deref()), |s| {
                    stats.push(s);
                    Ok(())
                })
                .map_err(|e| e.to_string())?;
                Ok(crate::consensus::call(&stats, region.start, region.end))
            })();
            Box::new(move |app: &mut App, ctx: &egui::Context| match result {
                Ok(m) => {
                    let header = format!("{title} {}:{}-{} consensus{header_filter}", region.contig, region.start + 1, region.end);
                    ctx.copy_text(crate::consensus::fasta(&header.replace(' ', "|"), &m.sequence));
                    let mut note = format!("Copied {title}'s sequence: {} bp", commas(m.sequence.len() as u64));
                    let extras: Vec<String> = [(m.ambiguous, "heterozygous (IUPAC)"), (m.low_depth, "N (under 3 reads)"), (m.deleted, "deleted, left out")]
                        .into_iter()
                        .filter(|(n, _)| *n > 0)
                        .map(|(n, what)| format!("{} {what}", commas(n as u64)))
                        .collect();
                    if !extras.is_empty() {
                        note.push_str(&format!(" ({})", extras.join(", ")));
                    }
                    app.notify(note);
                }
                Err(e) => app.error = Some(format!("{title}: copying its sequence: {e}")),
            })
        });
    }

    /// Copies one read as its SAM line (`sequence_only`: just the bases).
    pub(super) fn copy_read(&mut self, ctx: &egui::Context, track: usize, read: &AlignedRead, sequence_only: bool) {
        let Some(t) = self.tracks.get(track) else { return };
        let Some(cur) = &t.current else { return };
        let (path, reference) = (t.loader.bam_path.clone(), self.reference_path.clone());
        let region = Region { contig: cur.view.region.contig.clone(), start: read.start, end: read.end.max(read.start + 1) };
        let (name, flags, start) = (read.name.clone(), read.flags, read.start);
        self.spawn_job(ctx, move || {
            let result = ggv_io::BamFile::open_with_reference(&path, reference.as_ref()).and_then(|mut bam| {
                bam.sam_lines(&region, 1, |r| {
                    r.name().is_some_and(|n| <[u8]>::eq(n.as_ref(), name.as_bytes()))
                        && u16::from(r.flags()) == flags
                        && r.alignment_start().and_then(|p| p.ok()).map(|p| usize::from(p) as u64 - 1) == Some(start)
                })
            });
            Box::new(move |app: &mut App, ctx: &egui::Context| match result {
                Ok(lines) if !lines.is_empty() => {
                    let line = &lines[0];
                    if sequence_only {
                        let seq = line.split('\t').nth(9).unwrap_or("").to_owned();
                        app.notify(format!("Copied the read's sequence ({} bases)", seq.len()));
                        ctx.copy_text(seq);
                    } else {
                        ctx.copy_text(line.clone());
                        app.notify("Copied the read as SAM");
                    }
                }
                Ok(_) => app.error = Some("That read is no longer in the file".into()),
                Err(e) => app.error = Some(format!("copy read: {e}")),
            })
        });
    }

    pub(super) fn copy_selection_sam(&mut self, ctx: &egui::Context) {
        let Some((c, start, end)) = self.selection else {
            self.error = Some("Shift-drag to select a region first".into());
            return;
        };
        let contig = self.contigs[c].name.clone();
        // Each file with the filter its reads are drawn with (its own expression, if it has one).
        let files: Vec<(PathBuf, Option<String>, ggv_query::PileupOptions)> = self
            .tracks
            .iter()
            .map(|t| (t.loader.bam_path.clone(), t.contig_name(&contig).map(|x| x.name.clone()), self.track_style(&t.loader.bam_path).filter.pileup_options_with(self.reference_path.as_deref())))
            .collect();
        let reference = self.reference_path.clone();
        self.spawn_job(ctx, move || {
            const LIMIT: usize = 20_000;
            let mut out = String::new();
            let mut n = 0;
            let mut errors = Vec::new();
            for (path, local, opts) in files {
                let Some(local) = local else { continue };
                let region = Region { contig: local, start, end };
                // The reference under the selection, for expressions about mismatches.
                let bases = opts
                    .read_filter
                    .as_ref()
                    .and(reference.as_ref())
                    .and_then(|r| ggv_io::ReferenceFasta::open(r).ok()?.fetch(&region).ok());
                let keep = |r: &ggv_io::pasta::bam::Record| keeps_record(&opts, r, &region.contig, bases.as_deref().map(|b| (b, start)));
                match ggv_io::BamFile::open_with_reference(&path, reference.as_ref()).and_then(|mut b| b.sam_lines(&region, LIMIT, keep)) {
                    Ok(lines) => {
                        out.push_str(&format!("@CO\tggv: {} {}:{}-{}\n", path.display(), region.contig, start + 1, end));
                        n += lines.len();
                        for l in lines {
                            out.push_str(&l);
                            out.push('\n');
                        }
                    }
                    Err(e) => errors.push(e.to_string()),
                }
            }
            Box::new(move |app: &mut App, ctx: &egui::Context| {
                if !errors.is_empty() {
                    app.error = Some(errors.join("; "));
                }
                ctx.copy_text(out);
                app.notify(format!("Copied {} reads as SAM{}", commas(n as u64), if n >= LIMIT { " (first 20,000 per file)" } else { "" }));
            })
        });
    }
}

/// Whether record `r` passes `opts` as the reads drawn do (flags, MAPQ, then the filter
/// expression, given the reference over the span when there is one).
fn keeps_record(opts: &ggv_query::PileupOptions, r: &ggv_io::pasta::bam::Record, contig: &str, reference: Option<(&[u8], u64)>) -> bool {
    let mut counts = ggv_query::ReadCounts::default();
    if opts.read_filter.is_none() {
        return counts.admit_fields(u16::from(r.flags()), r.mapping_quality().map_or(255, u8::from), opts);
    }
    ggv_io::pasta::pileup::Read::new(r.as_ref().to_vec()).is_ok_and(|read| counts.admit_in(&read, opts, contig, reference))
}

impl App {
    // ---- Bookmarks ----

    pub(super) fn add_bookmark(&mut self) {
        let Some((c, start, end)) = self.place() else { return };
        let i = bookmarks::add(&mut self.bookmarks.borrow_mut(), Bookmark { contig: self.contigs[c].name.clone(), start, end, note: String::new() });
        self.show_bookmarks = true;
        self.edit_bookmark = Some(i);
        let what = if end - start > 1 { format!("{}:{}-{}", self.contigs[c].name, commas(start + 1), commas(end)) } else { format!("{}:{}", self.contigs[c].name, commas(start + 1)) };
        self.notify(format!("Bookmarked {what}; add a note in the Bookmarks window"));
    }

    pub(super) fn go_to_bookmark(&mut self, i: usize) {
        let Some(b) = self.bookmarks.borrow().get(i).cloned() else { return };
        if b.end == b.start + 1 {
            self.pending_goto = Some(format!("{}:{}", b.contig, b.start + 1));
        } else {
            self.pending_goto = Some(format!("{}:{}-{}", b.contig, b.start + 1, b.end));
        }
    }

    pub(super) fn export_bookmarks(&mut self) {
        if self.bookmarks.borrow().is_empty() {
            self.error = Some("No bookmarks yet (Ctrl+D adds one)".into());
            return;
        }
        let bed = bookmarks::to_bed(&self.bookmarks.borrow());
        self.ask(rfd::FileDialog::new().add_filter("BED", &["bed"]).set_file_name("bookmarks.bed"), super::dialogs::Kind::Save, super::dialogs::Then::SaveText(bed));
    }

    pub(super) fn bookmarks_window(&mut self, ctx: &egui::Context) {
        if !self.show_bookmarks || !self.is_active {
            return;
        }
        let mut open = true;
        let mut go = None;
        let mut remove = None;
        let mut add = false;
        let mut export = false;
        let edit = self.edit_bookmark.take();
        egui::Window::new("Bookmarks")
            .id(egui::Id::unique(("bookmarks", self.reference_slot)))
            .open(&mut open)
            .default_pos(Pos2::new(ctx.content_rect().right() - 380.0, 90.0))
            .default_width(340.0)
            .show(ctx, |ui| {
                let mut list = self.bookmarks.borrow_mut();
                if list.is_empty() {
                    ui.label("No bookmarks yet. Ctrl+D bookmarks the selection, the outlined base, or the view.");
                }
                egui::ScrollArea::vertical().max_height(320.0).show(ui, |ui| {
                    for (i, b) in list.iter_mut().enumerate() {
                        ui.horizontal(|ui| {
                            if ui.link(b.locus()).on_hover_text("Go there").clicked() {
                                go = Some(i);
                            }
                            let r = ui.add(egui::TextEdit::singleline(&mut b.note).hint_text("note").desired_width(150.0));
                            if edit == Some(i) {
                                r.request_focus();
                            }
                            if ui.small_button("Delete").clicked() {
                                remove = Some(i);
                            }
                        });
                    }
                });
                ui.separator();
                ui.horizontal(|ui| {
                    add = ui.button("Add this place (Ctrl+D)").clicked();
                    export = ui.add_enabled(!list.is_empty(), egui::Button::new("Export BED…")).clicked();
                });
            });
        if let Some(i) = remove {
            self.bookmarks.borrow_mut().remove(i);
        }
        if let Some(i) = go {
            self.go_to_bookmark(i);
        }
        if add {
            self.add_bookmark();
        }
        if export {
            self.export_bookmarks();
        }
        if !open {
            self.show_bookmarks = false;
        }
    }

    // ---- Review queue ----

    pub(super) fn open_queue_dialog(&mut self) {
        self.ask(rfd::FileDialog::new().add_filter("Places", &["bed", "vcf", "gz", "txt", "tsv", "csv", "list"]), super::dialogs::Kind::PickFile, super::dialogs::Then::OpenQueue);
    }

    pub(super) fn open_queue_file(&mut self, path: &Path) {
        match queue::read_file(path) {
            Ok(text) => self.start_queue(queue::Queue::for_contigs(crate::files::display_name(path), &text, Some(path), &self.contigs)),
            Err(e) => self.error = Some(super::path_error(path, e)),
        }
    }

    pub(super) fn start_queue(&mut self, q: queue::Queue) {
        if q.items.is_empty() {
            self.error = Some(format!("{}: no places found ({} lines skipped)", q.source, q.skipped));
            return;
        }
        self.review = None; // one review at a time: they share the keys
        self.queue = Some(q);
        self.show_queue_item();
    }

    /// Shows the queue's current item: a single base centred with ±50 bp and outlined, a
    /// region with a little room either side.
    pub(super) fn show_queue_item(&mut self) {
        // The place as it is meant in this genome (a BED row on a contig with `:` in its name).
        if let Some(q) = &mut self.queue {
            q.settle(&self.contigs);
        }
        let Some(item) = self.queue.as_ref().and_then(|q| q.current()).cloned() else { return };
        let text = match item.end {
            Some(e) if e == item.start + 1 => format!("{}:{}", item.contig, item.start + 1),
            Some(e) => {
                let pad = ((e - item.start) / 10).max(20);
                format!("{}:{}-{}", item.contig, item.start.saturating_sub(pad) + 1, e + pad)
            }
            None => item.contig.clone(),
        };
        self.pending_goto = Some(text);
    }

    pub(super) fn queue_decide(&mut self, v: queue::Verdict) {
        if let Some(q) = &mut self.queue {
            q.decide(v);
        }
        self.show_queue_item();
    }

    pub(super) fn queue_step(&mut self, forward: bool) {
        if let Some(q) = &mut self.queue {
            q.step(forward);
        }
        self.show_queue_item();
    }

    pub(super) fn queue_next_undecided(&mut self) {
        if let Some(q) = &mut self.queue
            && let Some(i) = q.next_undecided(q.pos + 1)
        {
            q.go(i);
        }
        self.show_queue_item();
    }

    pub(super) fn paste_queue_window(&mut self, ctx: &egui::Context) {
        let Some(mut text) = self.paste_queue.take() else { return };
        if !self.is_active {
            return;
        }
        let mut open = true;
        let mut start = false;
        egui::Window::new("Review queue: paste places")
            .id(egui::Id::unique(("paste-queue", self.reference_slot)))
            .open(&mut open)
            .anchor(Align2::CENTER_CENTER, Vec2::ZERO)
            .collapsible(false)
            .show(ctx, |ui| {
                ui.label("One per line: chr1:100-200, chr1:12345, BED (chr1 99 200 name) or VCF lines. Text after a locus is its label.");
                let edit = ui.add(egui::TextEdit::multiline(&mut text).desired_rows(12).desired_width(460.0).font(egui::TextStyle::Monospace));
                // Typing goes into the box when the window opens, not to the tracks as shortcuts.
                let first = ui.memory_mut(|m| m.data.remove_temp::<bool>(egui::Id::unique("paste-queue-focus")).unwrap_or(false));
                if first {
                    edit.request_focus();
                }
                let (items, skipped) = queue::parse(&text);
                ui.horizontal(|ui| {
                    start = ui.add_enabled(!items.is_empty(), egui::Button::new(format!("Review {} places", items.len()))).clicked();
                    if skipped > 0 {
                        ui.label(egui::RichText::new(format!("{skipped} lines not understood")).weak());
                    }
                });
            });
        if start {
            self.start_queue(queue::Queue::for_contigs("pasted list".into(), &text, None, &self.contigs));
        } else if open {
            self.paste_queue = Some(text);
        }
    }

    /// The review queue: the current place, its decision and note, the controls and the list.
    pub(super) fn queue_ui(&mut self, ui: &mut egui::Ui) {
        if self.queue.is_none() {
            return;
        }
        let mut decide = None;
        let mut step = None;
        let mut undecided = false;
        let mut stop = false;
        let mut export: Option<bool> = None;
        let mut go = None;
        let q = self.queue.as_mut().expect("checked");
        let [pass, fail, unsure] = q.counts();
        let done = pass + fail + unsure;
        ui.horizontal(|ui| {
            ui.label(format!("{} of {} · {done} decided ({pass} pass, {fail} fail, {unsure} unsure)", q.pos + 1, q.items.len()));
        });
        ui.add(egui::ProgressBar::new(done as f32 / q.items.len().max(1) as f32).desired_width(300.0));
        if let Some(p) = &q.save_path {
            ui.label(egui::RichText::new(format!("saved to {}", crate::files::display_name(p))).weak());
        }
        ui.separator();
        if let Some(item) = q.current() {
            ui.label(egui::RichText::new(item.key()).strong().monospace());
            if !item.label.is_empty() {
                ui.label(&item.label);
            }
            if let Some(d) = q.decision(q.pos) {
                ui.colored_label(rgb(d.verdict.rgb()), format!("marked {}{}", d.verdict.name(), if d.note.is_empty() { String::new() } else { format!(": {}", d.note) }));
            }
        }
        ui.add(egui::TextEdit::singleline(&mut q.note).hint_text("note (optional)").desired_width(300.0));
        ui.horizontal(|ui| {
            for (v, key) in queue::Verdict::ALL.into_iter().zip(["1", "2", "3"]) {
                if ui.button(egui::RichText::new(format!("{} ({key})", v.name())).color(rgb(v.rgb()))).clicked() {
                    decide = Some(v);
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
            undecided = ui.button("Next undecided (U)").clicked();
        });
        ui.horizontal(|ui| {
            if ui.button("Export TSV…").clicked() {
                export = Some(false);
            }
            if ui.button("Export BED…").clicked() {
                export = Some(true);
            }
            stop = ui.button("Stop").clicked();
        });
        egui::CollapsingHeader::new("All places").id_salt(("queue-list", self.reference_slot)).show(ui, |ui| {
            egui::ScrollArea::vertical().max_height(200.0).show_rows(ui, 18.0, q.items.len(), |ui, rows| {
                for i in rows {
                    let it = &q.items[i];
                    let mark = q.decision(i).map_or(Color32::GRAY, |d| rgb(d.verdict.rgb()));
                    ui.horizontal(|ui| {
                        ui.colored_label(mark, "●");
                        if ui.selectable_label(i == q.pos, format!("{}  {}", it.key(), it.label)).clicked() {
                            go = Some(i);
                        }
                    });
                }
            });
        });
        if let Some(e) = &q.error {
            ui.colored_label(Color32::from_rgb(200, 60, 60), format!("Not saved: {e}"));
        }
        if let Some(n) = &q.notice {
            ui.colored_label(Color32::from_rgb(190, 120, 0), n);
        }
        if let Some(v) = decide {
            self.queue_decide(v);
        }
        if let Some(f) = step {
            self.queue_step(f);
        }
        if undecided {
            self.queue_next_undecided();
        }
        if let Some(i) = go {
            if let Some(q) = &mut self.queue {
                q.go(i);
            }
            self.show_queue_item();
        }
        if let Some(bed) = export
            && let Some(q) = &self.queue
        {
            let (ext, body) = if bed { ("bed", q.to_bed()) } else { ("tsv", q.to_tsv(true)) };
            self.ask(rfd::FileDialog::new().add_filter(ext, &[ext]).set_file_name(format!("review.{ext}")), super::dialogs::Kind::Save, super::dialogs::Then::SaveText(body));
        }
        if stop {
            self.queue = None;
        }
    }

    // ---- Filter expressions ----

    /// The style a track draws with: its own filter expression, if it has one.
    pub(super) fn track_style(&self, path: &Path) -> Style {
        let mut s = self.style;
        s.detail = crate::tracks::Detail::for_ppb(self.ppb);
        if let Some(e) = self.track_filters.get(path) {
            s.filter.expr = Some(*e);
        }
        s
    }

    pub(super) fn open_filter_window(&mut self, track: Option<usize>) {
        let target = track.and_then(|i| self.tracks.get(i)).map(|t| t.loader.bam_path.clone());
        let current = match &target {
            Some(p) => self.track_filters.get(p).copied(),
            None => self.style.filter.expr,
        };
        self.filter_ui = Some(FilterUi { text: current.map(expr::text).unwrap_or_default(), error: None, target, focused: false });
    }

    /// Applies `text` to every track (`target: None`) or one file; empty text clears it.
    pub(super) fn apply_filter_text(&mut self, text: &str, target: Option<PathBuf>) -> Option<String> {
        let id = if text.trim().is_empty() {
            None
        } else {
            match expr::intern(text) {
                Ok(id) => Some(id),
                Err(e) => return Some(e),
            }
        };
        match target {
            None => self.style.filter.expr = id,
            Some(p) => {
                match id {
                    Some(id) => self.track_filters.insert(p.clone(), id),
                    None => self.track_filters.remove(&p),
                };
                // Remembered per file, for the next time it is opened.
                let mut saved = load_file_filters();
                match id {
                    Some(id) => saved.insert(p, expr::text(id)),
                    None => saved.remove(&p),
                };
                save_json("file-filters.json", &saved);
            }
        }
        None
    }

    pub(super) fn filter_window(&mut self, ctx: &egui::Context) {
        let Some(mut f) = self.filter_ui.take() else { return };
        if !self.is_active {
            return;
        }
        let mut open = true;
        let mut apply = false;
        let mut save = false;
        let mut load = None;
        let mut forget = None;
        let title = match &f.target {
            Some(p) => format!("Filter {}", crate::files::display_name(p)),
            None => "Filter reads (all tracks)".to_owned(),
        };
        let saved = self.saved_filters.clone();
        egui::Window::new(title)
            .id(egui::Id::unique(("filter", self.reference_slot)))
            .open(&mut open)
            .default_pos(Pos2::new(ctx.content_rect().center().x - 260.0, 90.0))
            .default_width(520.0)
            .show(ctx, |ui| {
                let edit = ui.add(
                    egui::TextEdit::singleline(&mut f.text)
                        .hint_text("mapq >= 20 && hp == 1 && !dup")
                        .desired_width(f32::INFINITY)
                        .font(egui::TextStyle::Monospace),
                );
                // Typing goes into the field when the window opens (it went to the tracks as
                // shortcuts: `-` zoomed out, `s` squished the rows).
                if !f.focused {
                    edit.request_focus();
                    f.focused = true;
                }
                if edit.changed() {
                    f.error = (!f.text.trim().is_empty()).then(|| expr::Expr::parse(&f.text).err()).flatten();
                }
                apply |= edit.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
                match &f.error {
                    Some(e) => {
                        ui.colored_label(Color32::from_rgb(200, 60, 60), e);
                    }
                    None => {
                        ui.label(egui::RichText::new("Applies where reads are loaded (views under 100 kb), on top of the flag and MAPQ filters.").weak());
                    }
                }
                ui.horizontal(|ui| {
                    apply |= ui.add_enabled(f.error.is_none(), egui::Button::new("Apply")).clicked();
                    if ui.button("Clear").clicked() {
                        f.text.clear();
                        f.error = None;
                        apply = true;
                    }
                    save = ui.add_enabled(f.error.is_none() && !f.text.trim().is_empty(), egui::Button::new("Save for reuse")).clicked();
                });
                if !saved.is_empty() {
                    ui.separator();
                    ui.label("Saved filters");
                    for s in &saved {
                        ui.horizontal(|ui| {
                            if ui.link(egui::RichText::new(s).monospace()).clicked() {
                                load = Some(s.clone());
                            }
                            if ui.small_button("Forget").clicked() {
                                forget = Some(s.clone());
                            }
                        });
                    }
                }
                egui::CollapsingHeader::new("Fields and operators").id_salt(("filter-help", self.reference_slot)).show(ui, |ui| {
                    egui::Grid::new("filter-fields").striped(true).show(ui, |ui| {
                        for (name, what) in expr::FIELDS {
                            ui.monospace(*name);
                            ui.label(*what);
                            ui.end_row();
                        }
                    });
                    ui.label("Operators: && || ! (or and, or, not), == != < <= > >=, ~ (contains), & (bits), parentheses.");
                });
            });
        if let Some(s) = load {
            f.text = s;
            f.error = None;
            apply = true;
        }
        if let Some(s) = forget {
            self.saved_filters.retain(|x| *x != s);
            save_json("filters.json", &self.saved_filters);
        }
        if save {
            let t = f.text.trim().to_owned();
            if !self.saved_filters.contains(&t) {
                self.saved_filters.insert(0, t);
                save_json("filters.json", &self.saved_filters);
            }
        }
        if apply {
            f.error = self.apply_filter_text(&f.text, f.target.clone());
            if f.error.is_none() && !f.text.trim().is_empty() {
                self.notify(format!("Filtering reads: {}", f.text.trim()));
            }
        }
        if open {
            self.filter_ui = Some(f);
        }
    }

    // ---- Find a read ----

    pub(super) fn find_window(&mut self, ctx: &egui::Context) {
        let Some(mut f) = self.find.take() else { return };
        if !self.is_active {
            return;
        }
        if let Some(s) = &mut f.search
            && s.poll()
            && let Some(hit) = s.found.clone()
        {
            self.selected_read = Some(s.name.clone());
            self.jump_to(&hit.contig, hit.pos);
        }
        let mut open = true;
        let mut go: Option<(String, u64)> = None;
        let mut search = false;
        let mut copy = None;
        // Off to the side of the pane: what is found is centred in the view.
        let spot = self.centre_window_pos(ctx, 200.0);
        super::inspector::placed(ctx, egui::Window::new("Find read by name"), egui::Id::unique(("find-read", self.reference_slot)), spot, self.pane_settled())
            .open(&mut open)
            .default_width(460.0)
            .show(ctx, |ui| {
                let edit = ui.add(egui::TextEdit::singleline(&mut f.text).hint_text("read name").desired_width(f32::INFINITY).font(egui::TextStyle::Monospace));
                if !f.focused {
                    edit.request_focus();
                    f.focused = true;
                }
                search |= edit.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
                ui.horizontal(|ui| {
                    search |= ui.add_enabled(!f.text.trim().is_empty() && !self.tracks.is_empty(), egui::Button::new("Find")).clicked();
                    ui.label(egui::RichText::new("Loaded reads first, then every file, all cores.").weak());
                });
                if let Some(name) = &f.loaded_hit {
                    ui.label(format!("{name}: on screen (outlined)"));
                }
                if let Some(s) = &f.search {
                    match (&s.found, s.done) {
                        (Some(hit), _) => {
                            ui.separator();
                            ui.label(egui::RichText::new(format!("{}:{}  ·  {}", hit.contig, commas(hit.pos + 1), crate::files::display_name(&hit.file))).strong());
                            ui.label(format!("flags {} (0x{:x})  MAPQ {}", hit.flags, hit.flags, hit.mapq));
                            ui.horizontal(|ui| {
                                if ui.button("Go to it").clicked() {
                                    go = Some((hit.contig.clone(), hit.pos));
                                }
                                if let Some((mc, mp)) = &hit.mate
                                    && ui.button(format!("Go to mate {mc}:{}", commas(mp + 1))).clicked()
                                {
                                    go = Some((mc.clone(), *mp));
                                }
                            });
                            for (sc, sp) in &hit.splits {
                                if ui.button(format!("Go to split part {sc}:{}", commas(sp + 1))).clicked() {
                                    go = Some((sc.clone(), *sp));
                                }
                            }
                            if ui.button("Copy as SAM").clicked() {
                                copy = Some(hit.sam.clone());
                            }
                        }
                        (None, false) => {
                            let (done, total) = s.progress();
                            ui.horizontal(|ui| {
                                ui.spinner();
                                ui.label(format!("searching {} ({done} of {total} contigs)", s.at.lock().unwrap()));
                            });
                        }
                        (None, true) => {
                            ui.label(format!("{} not found in any open file (unplaced unmapped reads aren't searched)", s.name));
                        }
                    }
                    if let Some(e) = &s.error {
                        ui.colored_label(Color32::from_rgb(200, 60, 60), e);
                    }
                }
            });
        if search {
            self.start_find(&mut f, ctx);
        }
        if let Some((c, p)) = go {
            self.jump_to(&c, p);
        }
        if let Some(sam) = copy {
            ctx.copy_text(sam);
            self.notify("Copied the read as SAM");
        }
        if open {
            self.find = Some(f);
        }
    }

    fn start_find(&mut self, f: &mut FindUi, ctx: &egui::Context) {
        let name = f.text.trim().to_owned();
        f.loaded_hit = None;
        f.search = None;
        let view_contig = self.contigs.get(self.contig).map(|c| c.name.clone()).unwrap_or_default();
        let (vs, ve) = self.visible();
        // Loaded already? No disk needed. Reads load well beyond the view, so only one in the
        // visible span is on screen; one loaded off screen is gone to.
        let (mut on_screen, mut off_screen) = (false, None);
        'tracks: for t in &self.tracks {
            let Some(c) = &t.current else { continue };
            let Some((reads, _)) = &c.view.reads else { continue };
            let here = t.contig_name(&view_contig).is_some_and(|x| x.name == c.view.region.contig);
            for r in reads.iter().filter(|r| r.name == name) {
                if here && (r.start as f64) < ve && (r.end as f64) > vs {
                    on_screen = true;
                    break 'tracks;
                }
                off_screen.get_or_insert((c.view.region.contig.clone(), r.start));
            }
        }
        if on_screen || off_screen.is_some() {
            if !on_screen && let Some((contig, pos)) = off_screen {
                self.jump_to(&contig, pos);
            }
            self.selected_read = Some(name.clone());
            f.loaded_hit = Some(name);
            return;
        }
        let files = self
            .tracks
            .iter()
            .map(|t| {
                let local = t.contig_name(&view_contig).map(|c| c.name.clone());
                let mut regions = Vec::new();
                if let Some(l) = &local {
                    let len = t.contig_name(l).map_or(0, |c| c.length);
                    regions.push(Region { contig: l.clone(), start: (vs.max(0.0) as u64).min(len.saturating_sub(1)), end: (ve.max(1.0) as u64).min(len).max(1) });
                    regions.push(Region { contig: l.clone(), start: 0, end: len });
                }
                for c in t.contigs() {
                    if Some(&c.name) != local.as_ref() && c.length > 0 {
                        regions.push(Region { contig: c.name.clone(), start: 0, end: c.length });
                    }
                }
                (t.loader.bam_path.clone(), regions)
            })
            .collect();
        f.search = Some(finder::Search::start(name, files, self.reference_path.clone(), ctx));
    }

    pub(super) fn run_seq_now(&mut self, ctx: &egui::Context) {
        if let Some(mut s) = self.seq.take() {
            self.start_seq_search(&mut s, ctx);
            self.seq = Some(s);
        }
    }

    pub(super) fn run_find_now(&mut self, ctx: &egui::Context) {
        if let Some(mut f) = self.find.take() {
            self.start_find(&mut f, ctx);
            self.find = Some(f);
        }
    }

    // ---- Sequence search ----

    pub(super) fn seq_window(&mut self, ctx: &egui::Context) {
        let Some(mut s) = self.seq.take() else { return };
        if !self.is_active {
            return;
        }
        if let Some(search) = &mut s.search {
            search.poll();
        }
        let mut open = true;
        let mut start = false;
        let mut go = None;
        // Off to the side of the pane: what is found is centred in the view.
        let spot = self.centre_window_pos(ctx, 300.0);
        super::inspector::placed(ctx, egui::Window::new("Find sequence in the reference"), egui::Id::unique(("seq-search", self.reference_slot)), spot, self.pane_settled())
            .open(&mut open)
            .default_width(480.0)
            .show(ctx, |ui| {
                if self.reference_path.is_none() {
                    ui.colored_label(Color32::from_rgb(200, 60, 60), "Open a reference first (Open > Reference).");
                }
                // It searches the genome sequence only; a match is selected across every track
                // so you can see the reads there, but the reads themselves aren't searched.
                ui.label(
                    egui::RichText::new("Searches the reference genome on both strands, not the reads. Choosing a match selects it across every track.").weak(),
                );
                let edit = ui.add(
                    egui::TextEdit::singleline(&mut s.text)
                        .hint_text("primer or oligo, 4–64 bases (IUPAC codes allowed)")
                        .desired_width(f32::INFINITY)
                        .font(egui::TextStyle::Monospace),
                );
                if !s.focused {
                    edit.request_focus();
                    s.focused = true;
                }
                start |= edit.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
                ui.horizontal(|ui| {
                    ui.label("Mismatches");
                    ui.add(egui::DragValue::new(&mut s.mismatches).range(0..=seqsearch::MAX_MISMATCHES));
                    ui.separator();
                    ui.radio_value(&mut s.scope, Scope::View, "View");
                    ui.radio_value(&mut s.scope, Scope::Contig, "Contig");
                    ui.radio_value(&mut s.scope, Scope::Genome, "Genome");
                });
                ui.horizontal(|ui| {
                    start |= ui.add_enabled(self.reference_path.is_some() && !s.text.trim().is_empty(), egui::Button::new("Search both strands")).clicked();
                    if let Some(sr) = &s.search
                        && !sr.done
                        && ui.button("Stop").clicked()
                    {
                        sr.cancel();
                    }
                });
                if let Some(e) = &s.error {
                    ui.colored_label(Color32::from_rgb(200, 60, 60), e);
                }
                if let Some(sr) = &s.search {
                    ui.separator();
                    if sr.done {
                        ui.label(format!(
                            "{} matches for {} (≤{} mismatches) in {:.1} s{}",
                            commas(sr.hits.len() as u64),
                            sr.pattern,
                            sr.mismatches,
                            sr.elapsed.map_or(0.0, |d| d.as_secs_f64()),
                            if sr.truncated { " — stopped at the first 5,000" } else { "" }
                        ));
                    } else {
                        ui.add(egui::ProgressBar::new(sr.progress()).show_percentage().text(format!("{} matches so far", sr.hits.len())));
                    }
                    if let Some(e) = &sr.error {
                        ui.colored_label(Color32::from_rgb(200, 60, 60), e);
                    }
                    // Short enough that the window fits GGV's window (560 px high and up).
                    let list = (ctx.content_rect().height() - 360.0).clamp(90.0, 260.0);
                    egui::ScrollArea::vertical().max_height(list).show_rows(ui, 18.0, sr.hits.len(), |ui, rows| {
                        for h in &sr.hits[rows] {
                            let label = format!(
                                "{}:{}  {}  {}",
                                h.contig,
                                commas(h.start + 1),
                                if h.reverse { "−" } else { "+" },
                                if h.mismatches == 0 { "exact".to_owned() } else { format!("{} mismatch{}", h.mismatches, if h.mismatches == 1 { "" } else { "es" }) }
                            );
                            if ui.selectable_label(false, egui::RichText::new(label).monospace()).clicked() {
                                go = Some(h.clone());
                            }
                        }
                    });
                }
            });
        if start {
            self.start_seq_search(&mut s, ctx);
        }
        if let Some(h) = go {
            let len = s.search.as_ref().map_or(20, |x| x.len);
            self.jump_to(&h.contig, h.start + len / 2);
            self.zoom_to_show(len);
            if let Some(c) = self.contigs.iter().position(|x| x.name == h.contig) {
                self.selection = Some((c, h.start, h.start + len));
            }
        }
        if open {
            self.seq = Some(s);
        }
    }

    /// After a jump to a `len`-bp match: zoomed in, if need be, so the match (centred) spans a
    /// twentieth of the view, at least 150 bp across, rather than a pixel of a wide view.
    pub(super) fn zoom_to_show(&mut self, len: u64) {
        let span = (len as f64 * 20.0).max(150.0);
        if f64::from(self.width) / self.ppb > span {
            self.ppb = f64::from(self.width) / span;
            self.clamp_view();
        }
    }

    pub(super) fn start_seq_search(&mut self, s: &mut SeqUi, ctx: &egui::Context) {
        let Some(reference) = self.reference_path.clone() else { return };
        let Ok(fasta) = ggv_io::ReferenceFasta::open(&reference) else {
            s.error = Some("can't open the reference".into());
            return;
        };
        let local = |name: &str| ggv_io::resolve_locus_str(fasta.contigs(), name).ok();
        let regions: Vec<Region> = match s.scope {
            Scope::Genome => fasta.contigs().iter().map(|c| Region { contig: c.name.clone(), start: 0, end: c.length }).collect(),
            Scope::Contig => self.contigs.get(self.contig).and_then(|c| local(&c.name)).into_iter().collect(),
            Scope::View => {
                let (vs, ve) = self.visible();
                self.contigs
                    .get(self.contig)
                    .and_then(|c| local(&c.name))
                    .map(|r| Region { start: (vs.max(0.0) as u64).min(r.end), end: (ve.max(0.0) as u64).min(r.end), ..r })
                    .into_iter()
                    .collect()
            }
        };
        match seqsearch::Search::start(reference, &s.text, s.mismatches, regions, ctx) {
            Ok(search) => {
                s.error = None;
                s.search = Some(search);
            }
            Err(e) => s.error = Some(e),
        }
    }

    // ---- Compare two tracks ----

    pub(super) fn open_compare(&mut self) {
        if self.tracks.len() < 2 {
            self.error = Some("Open two alignment tracks to compare them".into());
            return;
        }
        if self.compare.is_none() {
            let files = (self.tracks[0].loader.bam_path.clone(), self.tracks[1].loader.bam_path.clone());
            self.compare = Some(CompareUi { a: 0, b: 1, files, settings: diff::Settings::default(), cache: None, scan: None, job: None, status: None });
        }
    }

    /// Compares alignment tracks `a` and `b` (indexes).
    pub(super) fn set_compare_tracks(&mut self, a: usize, b: usize) {
        let (Some(ta), Some(tb)) = (self.tracks.get(a), self.tracks.get(b)) else { return };
        let files = (ta.loader.bam_path.clone(), tb.loader.bam_path.clone());
        if let Some(c) = &mut self.compare {
            (c.a, c.b, c.files) = (a, b, files);
            c.cache = None;
        }
    }

    fn scan_key(&self) -> Option<ScanKey> {
        let c = self.compare.as_ref()?;
        let filter = |p: &Path| self.track_style(p).filter;
        Some((c.files.0.clone(), c.files.1.clone(), filter(&c.files.0), filter(&c.files.1), self.reference_path.clone(), c.settings))
    }

    /// Keeps the comparison on the files picked when tracks move or close, and restarts a
    /// running scan whose files, filters, reference or settings changed (its answer would be
    /// for the old ones).
    pub(super) fn sync_compare(&mut self) {
        let Some(c) = &self.compare else { return };
        let find = |p: &Path| self.tracks.iter().position(|t| t.loader.bam_path == p);
        let (a, b) = match (find(&c.files.0), find(&c.files.1)) {
            (Some(a), Some(b)) if a != b => (a, b),
            _ => {
                self.compare = None;
                self.notify("A compared track was closed, so the comparison was too");
                return;
            }
        };
        let key = self.scan_key();
        let c = self.compare.as_mut().expect("checked");
        if (c.a, c.b) != (a, b) {
            (c.a, c.b) = (a, b);
            c.cache = None;
        }
        let Some(job) = c.job.take_if(|j| c.scan.is_some() && Some(&j.key) != key.as_ref()) else { return };
        c.scan = None;
        if let Some(ci) = self.contigs.iter().position(|x| x.name == job.contig) {
            self.start_diff_scan(job.forward, ci, job.from, None);
        }
    }

    /// The span both compared tracks have checked base by base, `[lo, hi)`: read-level loads
    /// on the view's contig with a reference (without one, allele counts aren't filled in, and
    /// binned coverage has none).
    fn compare_checked(&self) -> Option<(u64, u64)> {
        let c = self.compare.as_ref()?;
        self.reference_path.as_ref()?;
        let contig = &self.contigs.get(self.contig)?.name;
        let span = |i: usize| {
            let t = self.tracks.get(i)?;
            let cur = t.current.as_ref().filter(|_| !t.missing_contig)?;
            let here = t.contig_name(contig)?;
            (cur.view.region.contig == here.name && matches!(&*cur.coverage, Coverage::PerBase(_))).then_some((cur.view.region.start, cur.view.region.end))
        };
        let ((a0, a1), (b0, b1)) = (span(c.a)?, span(c.b)?);
        let (lo, hi) = (a0.max(b0), a1.min(b1));
        (lo < hi).then_some((lo, hi))
    }

    /// Differences in the loaded data of both compared tracks (recomputed only when either
    /// track loads new data or the settings change).
    pub(super) fn diff_sites(&mut self) -> Option<Arc<Vec<diff::Site>>> {
        // Coverage another contig left loaded (a track lacking this one) isn't compared.
        self.compare_checked()?;
        let c = self.compare.as_mut()?;
        let (ta, tb) = (self.tracks.get(c.a)?, self.tracks.get(c.b)?);
        let (ca, cb) = (ta.current.as_ref()?, tb.current.as_ref()?);
        let key = (ca.id, cb.id, c.settings.threshold.to_bits(), u64::from(c.settings.min_depth), ta.id);
        if let Some((k, sites)) = &c.cache
            && *k == key
        {
            return Some(sites.clone());
        }
        let sites = match (&*ca.coverage, &*cb.coverage) {
            (Coverage::PerBase(a), Coverage::PerBase(b)) => {
                diff::sites(a, ca.view.region.start, b, cb.view.region.start, 0, u64::MAX, c.settings)
            }
            _ => Vec::new(),
        };
        let sites = Arc::new(sites);
        c.cache = Some((key, sites.clone()));
        Some(sites)
    }

    pub(super) fn step_diff(&mut self, forward: bool, _ctx: &egui::Context) {
        if self.compare.is_none() {
            self.open_compare();
        }
        self.sync_compare();
        if self.compare.is_none() || self.contigs.is_empty() {
            return;
        }
        let at = match self.highlight {
            Some((hc, p)) if hc == self.contig => p,
            _ => self.center.max(0.0) as u64,
        };
        // In what both tracks checked first; otherwise scan the files from the first base not
        // checked.
        let checked = self.compare_checked();
        let sites = checked.and_then(|_| self.diff_sites());
        let checked = checked.filter(|_| sites.is_some());
        match diff::step(sites.as_deref().map_or(&[], Vec::as_slice), checked, at, forward) {
            diff::Step::Site(pos) => {
                let contig = self.contigs[self.contig].name.clone();
                // Says what was found, as a scan does: "No more differences" from an earlier
                // step stayed beside the site.
                let site = sites.as_deref().and_then(|s| s.iter().find(|x| x.pos == pos));
                if let Some(c) = &mut self.compare {
                    c.status = site.map(|s| site_status(&contig, s));
                }
                self.show_variant(&contig, pos, pos + 1);
            }
            diff::Step::Scan(from) => self.start_diff_scan(forward, self.contig, from, None),
        }
    }

    /// Scans the compared files from `from` on contig `ci` (see [`diff::Scan::start`]), each with
    /// its own read filter. `status` is shown while it runs.
    fn start_diff_scan(&mut self, forward: bool, ci: usize, from: u64, status: Option<String>) {
        let Some(key) = self.scan_key() else { return };
        let Some(c) = &self.compare else { return };
        let (Some(a), Some(b)) = (self.tracks.get(c.a), self.tracks.get(c.b)) else { return };
        let Some(contig) = self.contigs.get(ci).map(|x| x.name.clone()) else { return };
        let opts = [key.2.pileup_options_with(self.reference_path.as_deref()), key.3.pileup_options_with(self.reference_path.as_deref())];
        let files = diff::Files { a: a.loader.bam_path.clone(), b: b.loader.bam_path.clone(), reference: self.reference_path.clone(), opts };
        let pair = |x: &ggv_model::Contig| (x.name.clone(), x.length);
        let contigs: Vec<(String, u64)> = if forward {
            self.contigs[ci..].iter().map(pair).collect()
        } else {
            self.contigs[..=ci].iter().rev().map(pair).collect()
        };
        let scan = diff::Scan::start(files, contigs, from, forward, c.settings, &self.egui_ctx);
        if let Some(c) = &mut self.compare {
            c.scan = Some(scan);
            c.job = Some(ScanJob { forward, contig, from, key });
            c.status = status;
        }
    }

    pub(super) fn compare_window(&mut self, ctx: &egui::Context) {
        self.sync_compare();
        // A finished scan: go there.
        if let Some(c) = &mut self.compare
            && let Some(found) = c.scan.as_ref().and_then(diff::Scan::poll)
        {
            c.scan = None;
            c.job = None;
            match found {
                diff::Found::Site(contig, s) => {
                    c.status = Some(site_status(&contig, &s));
                    self.show_variant(&contig, s.pos, s.pos + 1);
                }
                diff::Found::None => c.status = Some("No more differences in that direction".into()),
                diff::Found::Error(e) => c.status = Some(format!("error: {e}")),
            }
        }
        if self.compare.is_none() || !self.is_active {
            return;
        }
        let sites = self.diff_sites();
        let (vs, ve) = self.visible();
        let in_view = sites.as_ref().map_or(0, |s| s.iter().filter(|x| (x.pos as f64) >= vs && (x.pos as f64) < ve).count());
        let names: Vec<String> = self.tracks.iter().map(|t| t.title()).collect();
        let has_reference = self.reference_path.is_some();
        let reads_level = ve - vs <= super::read_window() as f64;
        let mut open = true;
        let mut step = None;
        let (spot, settled) = (self.centre_window_pos(ctx, 220.0), self.pane_settled());
        let c = self.compare.as_mut().expect("checked");
        let picked = (c.a, c.b);
        super::inspector::placed(ctx, egui::Window::new("Compare tracks"), egui::Id::unique(("compare", self.reference_slot)), spot, settled)
            .open(&mut open)
            // Over its own pane, clear of the view centre (placed above); it can be dragged.
            .resizable(false)
            .show(ctx, |ui| {
                egui::Grid::new(("compare-grid", c.a)).num_columns(2).show(ui, |ui| {
                    for (label, which) in [("A", &mut c.a), ("B", &mut c.b)] {
                        ui.label(label);
                        egui::ComboBox::from_id_salt(("compare", label)).width(280.0).selected_text(names.get(*which).cloned().unwrap_or_default()).show_ui(ui, |ui| {
                            for (i, n) in names.iter().enumerate() {
                                ui.selectable_value(which, i, n);
                            }
                        });
                        ui.end_row();
                    }
                    ui.label("Difference");
                    ui.add(egui::Slider::new(&mut c.settings.threshold, diff::THRESHOLDS).custom_formatter(|v, _| format!("{:.0}%", v * 100.0)));
                    ui.end_row();
                    ui.label("Min depth");
                    ui.add(egui::DragValue::new(&mut c.settings.min_depth).range(diff::MIN_DEPTHS));
                    ui.end_row();
                });
                if !has_reference {
                    ui.label(egui::RichText::new("Without a reference, differences are only found by stepping ([ and ]).").weak());
                } else if !reads_level {
                    ui.label(egui::RichText::new("Zoom in under 100 kb to mark differences in the view.").weak());
                } else {
                    ui.label(format!("{in_view} differences in view (orange)"));
                }
                ui.horizontal(|ui| {
                    if ui.button("< Prev ([)").clicked() {
                        step = Some(false);
                    }
                    if ui.button("Next (]) >").clicked() {
                        step = Some(true);
                    }
                    if let Some(scan) = &c.scan {
                        ui.spinner();
                        ui.label(format!("scanning {}", scan.at.lock().unwrap()));
                    }
                });
                if let Some(s) = &c.status {
                    ui.label(s);
                }
            });
        if c.a == c.b {
            c.b = (c.a + 1) % names.len().max(1);
        }
        // Tracks picked here: follow those files from now on (a running scan restarts).
        let (a, b) = (c.a, c.b);
        if (a, b) != picked {
            self.set_compare_tracks(a, b);
        }
        self.sync_compare();
        if let Some(f) = step {
            self.step_diff(f, ctx);
        }
        if !open {
            self.compare = None;
        }
    }

    // ---- Measure ----

    fn measured(&self, a: u64, b: u64) -> Measured {
        let gc = self.reference.as_ref().and_then(|(view, _, _)| {
            let r = &view.region;
            let seq = view.reference.as_ref()?;
            if r.contig != self.contigs[self.contig].name || a < r.start || b > r.end {
                return None;
            }
            let s = &seq[(a - r.start) as usize..(b - r.start) as usize];
            let acgt = s.iter().filter(|c| matches!(c, b'A' | b'C' | b'G' | b'T')).count();
            (acgt > 0).then(|| s.iter().filter(|c| matches!(c, b'G' | b'C')).count() as f64 / acgt as f64)
        });
        let tracks = self
            .tracks
            .iter()
            .map(|t| {
                let Some(cur) = &t.current else { return (t.short_label(), None, 0) };
                let r = &cur.view.region;
                let depth = match &*cur.coverage {
                    Coverage::PerBase(cols) if a >= r.start && b <= r.end && b > a => {
                        Some(cols[(a - r.start) as usize..(b - r.start) as usize].iter().map(|c| f64::from(c.depth)).sum::<f64>() / (b - a) as f64)
                    }
                    _ => None,
                };
                let reads = cur.view.reads.as_ref().map_or(0, |(reads, _)| {
                    let rows = &cur.packed.rows;
                    rows.iter().flatten().filter(|&&i| reads[i].start < b && reads[i].end > a).count()
                });
                (t.short_label(), depth, reads)
            })
            .collect();
        Measured { len: b - a, gc, tracks }
    }

    // ---- Overlays ----

    /// Bookmarks and the queue item on the ruler, sequence-search hits on the reference row,
    /// compare-mode differences through the alignment tracks, and the measure readout.
    pub(super) fn draw_tool_overlays(&mut self, painter: &egui::Painter, rect: Rect, layout: &Layout, vs: f64, ve: f64) {
        let contig = self.contigs[self.contig].name.clone();
        let diff = if self.compare.is_some() { self.diff_sites() } else { None };
        let x_of = |p: f64| rect.left() + ((p - vs) * self.ppb) as f32;
        let span = |a: u64, b: u64| (x_of(a as f64).max(rect.left()), x_of(b as f64).min(rect.right()).max(x_of(a as f64) + 2.0));

        // Notes share one line: each takes the room up to the next note (shortened with an
        // ellipsis), and one with no room is left to the ruler's hover rather than printed over
        // another.
        let mut marks: Vec<(f32, String)> = Vec::new();
        for b in self.bookmarks.borrow().iter().filter(|b| b.contig == contig && (b.end as f64) > vs && (b.start as f64) < ve) {
            let (x0, x1) = span(b.start, b.end);
            marks.push((x0, b.note.clone()));
            // Between the ruler's baseline and the reference row.
            let band = Rect::from_min_max(Pos2::new(x0, rect.top() + super::RULER_H - 4.5), Pos2::new(x1.max(x0 + 3.0), rect.top() + super::REF_TOP - 0.5));
            painter.rect_filled(band, 1.0, Color32::from_rgb(150, 90, 200));
        }
        marks.retain(|m| !m.1.is_empty());
        marks.sort_by(|a, b| a.0.total_cmp(&b.0));
        let note_color = Color32::from_rgb(130, 70, 180);
        let font = FontId::proportional(10.0);
        let mut free_from = rect.left();
        for k in 0..marks.len() {
            let (x0, note) = (marks[k].0, &marks[k].1);
            if x0 + 2.0 < free_from {
                continue;
            }
            let limit = marks.get(k + 1).map_or(rect.right(), |m| m.0) - 4.0;
            let room = limit - (x0 + 2.0);
            if room < 24.0 {
                continue;
            }
            let mut text = note.clone();
            let mut galley = painter.layout_no_wrap(text.clone(), font.clone(), note_color);
            while galley.size().x > room && text.chars().count() > 1 {
                text.pop();
                galley = painter.layout_no_wrap(format!("{}\u{2026}", text.trim_end()), font.clone(), note_color);
            }
            let w = galley.size().x;
            painter.galley(Pos2::new(x0 + 2.0, rect.top() + 24.0 - galley.size().y), galley, note_color);
            free_from = x0 + 2.0 + w + 6.0;
        }

        if let Some(item) = self.queue.as_ref().and_then(|q| q.current())
            && item.contig == contig
            && let Some(end) = item.end.filter(|&e| e > item.start + 1)
        {
            let (x0, x1) = span(item.start, end);
            painter.rect_stroke(
                Rect::from_min_max(Pos2::new(x0, rect.top() + super::REF_TOP - 1.0), Pos2::new(x1, rect.bottom())),
                0.0,
                Stroke::new(1.5, Color32::from_rgb(230, 140, 20)),
                egui::StrokeKind::Inside,
            );
        }

        if let Some(search) = self.seq.as_ref().and_then(|s| s.search.as_ref()) {
            let len = search.len;
            let from = search.hits.partition_point(|h| h.contig < contig);
            for h in search.hits[from..].iter().take_while(|h| h.contig == contig) {
                if (h.start + len) as f64 <= vs || h.start as f64 >= ve {
                    continue;
                }
                let (x0, x1) = span(h.start, h.start + len);
                let y = rect.top() + super::REF_TOP - 3.0;
                let color = if h.mismatches == 0 { Color32::from_rgb(20, 150, 90) } else { Color32::from_rgb(200, 150, 20) };
                painter.rect_filled(Rect::from_min_max(Pos2::new(x0, y - 3.0), Pos2::new(x1, y)), 0.0, color);
                let tip = if h.reverse { Pos2::new(x0, y - 1.5) } else { Pos2::new(x1, y - 1.5) };
                painter.circle_filled(tip, 2.0, color);
            }
        }

        if let Some(sites) = diff
            && let Some((ia, ib)) = self.compare.as_ref().map(|c| (c.a, c.b))
        {
            let orange = Color32::from_rgba_unmultiplied(255, 140, 0, 70);
            let from = sites.partition_point(|s| (s.pos as f64) < vs.floor());
            for s in sites[from..].iter().take_while(|s| (s.pos as f64) < ve) {
                let (x0, x1) = span(s.pos, s.pos + 1);
                for i in [ia, ib] {
                    if let Some(sec) = layout.sections.get(i) {
                        painter.rect_filled(Rect::from_min_max(Pos2::new(x0, sec.top() + COVERAGE_TOP), Pos2::new(x1, sec.bottom())), 0.0, orange);
                    }
                }
                painter.rect_filled(Rect::from_min_max(Pos2::new(x0, rect.top() + 16.0), Pos2::new(x1, rect.top() + 24.0)), 0.0, Color32::from_rgb(255, 140, 0));
            }
        }

        if let Some((c, a, b)) = self.measure
            && c == self.contig
            && b > a
        {
            let (x0, x1) = span(a, b);
            let top = rect.top() + super::REF_TOP - 1.0;
            painter.rect_filled(Rect::from_min_max(Pos2::new(x0, top), Pos2::new(x1, rect.bottom())), 0.0, Color32::from_rgba_unmultiplied(20, 160, 140, 28));
            for x in [x0, x1] {
                painter.vline(x, top..=rect.bottom(), Stroke::new(1.0, Color32::from_rgb(20, 140, 120)));
            }
            let m = self.measured(a, b);
            let mut lines = vec![format!("{} bp", commas(m.len))];
            if let Some(gc) = m.gc {
                lines.push(format!("GC {:.1}%", gc * 100.0));
            }
            for (name, depth, reads) in &m.tracks {
                lines.push(match depth {
                    Some(d) => format!("{name}: depth {d:.1}, {} reads", commas(*reads as u64)),
                    None => format!("{name}: zoom in to measure"),
                });
            }
            // Wrapped to the pane, so every line stays readable in a narrow pane.
            let font = FontId::monospace(11.0);
            let wrap = (rect.width() - 16.0).clamp(80.0, 420.0);
            let galley = painter.layout(lines.join("\n"), font, self.pal().text, wrap);
            let at = Pos2::new((x1 + 6.0).min(rect.right() - galley.size().x - 6.0).max(rect.left() + 4.0), top + 4.0);
            painter.rect_filled(Rect::from_min_size(at - Vec2::splat(4.0), galley.size() + Vec2::splat(8.0)), 3.0, self.pal().label_bg);
            painter.galley(at, galley, self.pal().text);
        }
    }

    // ---- Reloading changed files ----

    pub(super) fn watch_files(&mut self, ctx: &egui::Context) {
        self.watcher.enabled.store(self.auto_reload, std::sync::atomic::Ordering::Relaxed);
        if !self.auto_reload {
            return;
        }
        let mut files: Vec<PathBuf> = self.tracks.iter().map(|t| t.loader.bam_path.clone()).collect();
        files.extend(self.variants.iter().map(|v| v.path.clone()));
        files.extend(self.genes.iter().map(|g| g.path.clone()));
        files.extend(self.signals.iter().map(|s| s.path.clone()));
        files.extend(self.reference_path.clone());
        self.watcher.set(files);
        let changed: Vec<PathBuf> = self.watcher.rx.try_iter().collect();
        for p in changed {
            self.reload_path(&p, ctx);
        }
    }

    pub(super) fn reload_all(&mut self, ctx: &egui::Context) {
        let mut files: Vec<PathBuf> = self.tracks.iter().map(|t| t.loader.bam_path.clone()).collect();
        files.extend(self.variants.iter().map(|v| v.path.clone()));
        files.extend(self.genes.iter().map(|g| g.path.clone()));
        files.extend(self.signals.iter().map(|s| s.path.clone()));
        for p in files {
            self.reload_path(&p, ctx);
        }
    }

    /// Reopens `path` in place (same position in the track list, same row scroll), keeping the
    /// view where it is.
    pub(super) fn reload_path(&mut self, path: &Path, ctx: &egui::Context) {
        let name = crate::files::display_name(path);
        // An earlier failed reload of this file is settled by this one (a new failure says so
        // again): it stayed on screen after the file was fixed.
        if self.error.as_deref().is_some_and(|e| e.starts_with(&format!("reloading {name}:"))) {
            self.error = None;
        }
        if self.reference_path.as_deref() == Some(path) {
            self.set_reference(path.to_owned(), ctx);
            self.notify(format!("Reloaded {name}"));
            return;
        }
        if let Some(i) = self.tracks.iter().position(|t| t.loader.bam_path == path) {
            match crate::track::Track::open(crate::gpu::new_slot(), path.to_owned(), self.reference_path.clone(), ctx) {
                Ok(mut t) => {
                    t.scroll_y = self.tracks[i].scroll_y;
                    let old = std::mem::replace(&mut self.tracks[i], t);
                    self.retired.push(old.id);
                    if self.reference.as_ref().is_some_and(|(_, _, (id, _))| *id == old.id) {
                        self.reference = None;
                    }
                    // The new file's own warning, if any, shows again even if the old one's was dismissed.
                    self.dismissed_warnings.retain(|p| p != path);
                    self.notify(format!("Reloaded {name}"));
                }
                Err(e) => self.error = Some(format!("reloading {name}: {e:#} (will retry when it changes again)")),
            }
        }
        if let Some(i) = self.variants.iter().position(|v| v.path == path) {
            match crate::variants::VariantTrack::open(crate::gpu::new_slot(), path.to_owned(), ctx) {
                Ok(mut v) => {
                    // As the person left it: rows scrolled and sites coloured.
                    v.scroll_y = self.variants[i].scroll_y;
                    v.color = self.variants[i].color.clone();
                    let old = std::mem::replace(&mut self.variants[i], v);
                    self.retired.push(old.id);
                    // A review of it carries on with the new track.
                    self.rebind_review();
                    self.notify(format!("Reloaded {name}"));
                }
                Err(e) => self.error = Some(format!("reloading {name}: {e:#}")),
            }
        }
        if let Some(i) = self.genes.iter().position(|g| g.path == path) {
            let mut g = crate::genes::GeneTrack::open(crate::gpu::new_slot(), path.to_owned(), ctx);
            g.scroll_y = self.genes[i].scroll_y;
            let old = std::mem::replace(&mut self.genes[i], g);
            self.retired.push(old.id);
            self.notify(format!("Reloading {name}"));
        }
        // Signal files are read when opened (a bedGraph or SEG whole): open it again.
        if let Some(i) = self.signals.iter().position(|s| s.path == path) {
            self.signals[i] = crate::signals::SignalTrack::open(path.to_owned(), ctx);
            self.notify(format!("Reloading {name}"));
        }
    }

    /// The transient notice and the latest report, for the status bar.
    pub(super) fn notice_ui(&mut self, ui: &mut egui::Ui) {
        if let Some((msg, at)) = &self.notice {
            let age = at.elapsed();
            if age < Duration::from_secs(5) {
                ui.label(egui::RichText::new(msg).color(Color32::from_rgb(40, 140, 80)));
                ui.ctx().request_repaint_after(Duration::from_secs(5) - age);
                if let Some(p) = self.last_report.clone()
                    && msg.starts_with("Report saved")
                    && ui.small_button("Open").clicked()
                {
                    open_in_system(&p);
                }
                ui.separator();
            } else {
                self.notice = None;
            }
        }
    }
}

/// A difference found by Compare, for its status line.
fn site_status(contig: &str, s: &diff::Site) -> String {
    format!("{contig}:{}  {:.0}% difference in {}", commas(s.pos + 1), s.diff * 100.0, diff::ALLELES[s.allele])
}

#[cfg(test)]
mod tests {
    use super::keeps_record;

    /// A short search hit chosen in a wide view is zoomed to (it was a 1-px line at 14 kb across);
    /// one already wide enough on screen keeps the zoom.
    #[test]
    fn a_short_hit_is_zoomed_to() {
        let bam = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio/HG002.revio.slices.bam");
        if !bam.exists() {
            return;
        }
        let ctx = eframe::egui::Context::default();
        let mut app = super::App::new(&ctx, super::super::Args { bams: vec![bam], locus: Some("chr20:31150000-31164500".into()), ..Default::default() }).unwrap();
        app.width = 1050.0;
        app.ppb = 1050.0 / 14_500.0;
        app.jump_to("chr20", 31_157_048);
        app.zoom_to_show(12);
        let (vs, ve) = app.visible();
        assert!(ve - vs <= 241.0 && (vs..ve).contains(&31_157_048.0), "{vs}-{ve}");
        let ppb = app.ppb;
        app.zoom_to_show(12);
        assert_eq!(app.ppb, ppb, "already shown: no change");
    }

    #[test]
    fn copied_reads_pass_the_filter_expression() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio");
        let mut bam = ggv_io::BamFile::open(dir.join("HG002.revio.slices.bam")).unwrap();
        let region = ggv_model::Region { contig: "chr20".into(), start: 31_166_000, end: 31_166_050 };
        let with = |text: Option<&str>| crate::tracks::ReadFilter { expr: text.map(|t| crate::expr::intern(t).unwrap()), ..Default::default() }.pileup_options();
        let mut count = |opts: &ggv_query::PileupOptions| bam.sam_lines(&region, 10_000, |r| keeps_record(opts, r, "chr20", None)).unwrap().len();
        let all = count(&with(None));
        assert!(all > 0);
        // No read here has MAPQ over 60: nothing is copied.
        assert_eq!(count(&with(Some("mapq > 60"))), 0);
        // One haplotype: some, not all.
        let hp1 = count(&with(Some("hp == 1")));
        assert!(hp1 > 0 && hp1 < all, "{hp1} of {all}");
        // Flags still apply: supplementary records only when not hidden.
        let mut f = crate::tracks::ReadFilter { hide_supplementary: true, ..Default::default() };
        let no_supp = count(&f.pileup_options());
        f.hide_supplementary = false;
        assert!(no_supp <= count(&f.pileup_options()));
    }
}
