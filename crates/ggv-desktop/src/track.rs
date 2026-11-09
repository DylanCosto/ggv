//! One alignment track: its loader thread, what it has on screen, what it has asked for, and
//! its own vertical scroll. Tracks share the view (locus and zoom) but load independently, in
//! parallel.

use std::{path::PathBuf, sync::Arc, time::Duration};

use eframe::egui;
use ggv_layout::{ViewData, ViewOptions};
use ggv_model::{Contig, Region, contig_aliases};

use crate::{
    gpu::InstanceSet,
    loader::{Loaded, Loader, Request},
    tracks::{self, Style},
};

/// A read-level span more than this many times longer than the one wanted is reloaded rather
/// than kept: zoomed in on it, every read in megabases off screen would be rebuilt at full
/// detail (modification colours per base, letter backings).
pub const MAX_OVERSIZE: u64 = 8;

/// Reads per splice junction, by (start, end).
pub type JunctionCounts = std::sync::Arc<std::collections::BTreeMap<(u64, u64), u32>>;

/// What a load covers, so we can tell whether it satisfies the current view. Regions use the
/// track's own contig names.
#[derive(Clone, Debug)]
pub struct Span {
    pub region: Region,
    pub reads: bool,
    /// Bases per coverage bin, for binned loads.
    pub bp_per_bin: Option<f64>,
    /// The read filter applied at load (binned loads only; read-level loads filter in memory).
    pub filter: Option<tracks::ReadFilter>,
    /// The tag reads were loaded with (for colouring or grouping by it).
    pub tag: Option<[u8; 2]>,
    /// Whether reads were loaded with their base-modification calls (to colour by them).
    pub mods: bool,
}

impl Span {
    /// Whether this span can draw `want` (same contig and mode) over `[vs, ve)` without
    /// reloading.
    pub fn satisfies(&self, want: &Span, vs: f64, ve: f64) -> bool {
        self.region.contig == want.region.contig
            && self.reads == want.reads
            && self.filter == want.filter
            && self.tag == want.tag
            && (self.mods || !want.mods)
            && (!self.reads || self.region.len() <= want.region.len().saturating_mul(MAX_OVERSIZE))
            && self.region.start as f64 <= vs.max(want.region.start as f64)
            && self.region.end as f64 >= ve.min(want.region.end as f64)
            && match (self.bp_per_bin, want.bp_per_bin) {
                // Binned depth is fine until it's twice as coarse as the screen.
                (Some(have), Some(w)) => have <= w * 2.0,
                _ => true,
            }
    }
}

pub struct Current {
    pub view: Arc<ViewData>,
    /// Reads packed for `style` (filtering, sorting and grouping).
    pub packed: Arc<ggv_layout::Packed>,
    /// Coverage for `style`'s read filter.
    pub coverage: Arc<ggv_layout::Coverage>,
    pub depth_index: Arc<tracks::DepthIndex>,
    pub instances: InstanceSet,
    pub style: Style,
    pub span: Span,
    pub elapsed: Duration,
    /// The request id this data answered.
    pub id: u64,
    /// Reference-row instances for this span (built by the loader).
    pub reference: Option<Arc<Vec<crate::gpu::Instance>>>,
}

pub struct Track {
    /// GPU slot id; unique for the app's lifetime.
    pub id: u64,
    pub loader: Loader,
    pub current: Option<Current>,
    requested: Option<(u64, Span)>,
    /// Requests sent and not yet answered (the worker may skip some), oldest first.
    inflight: Vec<(u64, Span)>,
    pub scroll_y: f32,
    /// The height of the track's section when last drawn (0 before), for reporting which read
    /// rows are on screen.
    pub view_h: f32,
    /// How long each load took (for benchmarks).
    pub load_ms: Vec<f64>,
    /// Splice junctions in what's loaded, for the filter they were counted with.
    junction_cache: std::cell::RefCell<Option<(u64, tracks::ReadFilter, JunctionCounts)>>,
    pub loads: u64,
    /// Loads served from the loader's memory (cache or prefetch) rather than the disk.
    pub cache_hits: u64,
    /// The span whose load failed last (not asked for again while it would serve the view).
    failed: Option<Span>,
    pub error: Option<String>,
    /// The view's contig is not in this file under any alias.
    pub missing_contig: bool,
    /// The latest pileup at the evidence base, and what it was asked for.
    pub evidence: Option<crate::loader::Evidence>,
    evidence_key: Option<(String, u64, tracks::ReadFilter)>,
    evidence_id: u64,
    /// The summary of a selected region: the request id and, once back, its result.
    pub summary: Option<(u64, Option<Result<ggv_query::RegionSummary, String>>)>,
    /// The filter the summary was asked for with (it is redone when the filter changes).
    pub summary_filter: Option<tracks::ReadFilter>,
    /// The region the summary was asked for (view contig name, start, end), also when the
    /// contig isn't in this file, so a missing summary for it means "not in this file".
    pub summary_for: Option<(String, u64, u64)>,
}

impl Track {
    pub fn open(id: u64, bam: PathBuf, reference: Option<PathBuf>, ctx: &egui::Context) -> anyhow::Result<Self> {
        Ok(Self {
            id,
            loader: Loader::spawn(bam, reference, ctx.clone())?,
            current: None,
            requested: None,
            inflight: Vec::new(),
            scroll_y: 0.0,
            view_h: 0.0,
            load_ms: Vec::new(),
            junction_cache: Default::default(),
            loads: 0,
            cache_hits: 0,
            error: None,
            failed: None,
            missing_contig: false,
            evidence: None,
            evidence_key: None,
            evidence_id: 0,
            summary: None,
            summary_filter: None,
            summary_for: None,
        })
    }

    /// The sample name from the file's read groups (`SM`), else the file name.
    pub fn label(&self) -> &str {
        &self.loader.label
    }

    pub fn file_name(&self) -> String {
        crate::files::display_name(&self.loader.bam_path)
    }

    /// A short name for tight places: the sample, unless it is an opaque ID (then the file).
    pub fn short_label(&self) -> String {
        if opaque_id(self.label()) { self.file_name() } else { self.label().to_owned() }
    }

    /// `file · sample` for headers, with an ID-like sample shortened (`d27993fa…`), and just the
    /// file when the sample adds nothing.
    pub fn title(&self) -> String {
        let (file, sample) = (self.file_name(), self.label());
        let stem = file.split('.').next().unwrap_or(&file);
        if sample == file || sample == stem || sample.is_empty() {
            file
        } else if opaque_id(sample) {
            format!("{file}  ·  {}…", sample.chars().take(8).collect::<String>())
        } else {
            format!("{file}  ·  {sample}")
        }
    }

    /// `file · full sample name`, for copying.
    pub fn full_title(&self) -> String {
        let (file, sample) = (self.file_name(), self.label());
        if sample == file { file } else { format!("{file}  ·  {sample}") }
    }

    pub fn contigs(&self) -> &[Contig] {
        &self.loader.contigs
    }

    pub fn is_loading(&self) -> bool {
        self.requested.is_some()
    }

    /// This file's name for `contig` (exact, or via `chr` aliases).
    pub fn contig_name(&self, contig: &str) -> Option<&Contig> {
        let contigs = self.contigs();
        contigs.iter().find(|c| c.name == contig).or_else(|| {
            contig_aliases(contig).into_iter().find_map(|a| contigs.iter().find(|c| c.name == a))
        })
    }

    /// Translates a span in the view's naming into this file's naming.
    pub fn localize(&self, want: &Span) -> Option<Span> {
        let contig = self.contig_name(&want.region.contig)?;
        let end = want.region.end.min(contig.length);
        let start = want.region.start.min(end.saturating_sub(1));
        Some(Span { region: Region { contig: contig.name.clone(), start, end }, ..want.clone() })
    }

    /// Requests data for `want` (view naming) unless what's shown or pending already covers
    /// `[vs, ve)`.
    pub fn ensure_loaded(&mut self, want: &Span, vs: f64, ve: f64, opts: ViewOptions, style: Style, next_id: &mut u64) {
        let Some(mut want) = self.localize(want) else {
            self.missing_contig = true;
            return;
        };
        self.missing_contig = false;
        want.filter = (!want.reads).then_some(style.filter);
        want.tag = style.tag();
        if self.current.as_ref().is_some_and(|c| c.span.satisfies(&want, vs, ve))
            || self.requested.as_ref().is_some_and(|(_, s)| s.satisfies(&want, vs, ve))
            || self.failed.as_ref().is_some_and(|s| s.satisfies(&want, vs, ve))
        {
            return;
        }
        let id = *next_id;
        *next_id += 1;
        let v0 = (vs.max(0.0).floor() as u64).clamp(want.region.start, want.region.end.saturating_sub(1));
        let v1 = (ve.max(0.0).ceil() as u64).clamp(v0 + 1, want.region.end.max(v0 + 1));
        self.loader.request(Request { id, region: want.region.clone(), view: (v0, v1), opts, style });
        if self.inflight.len() >= 64 {
            self.inflight.remove(0);
        }
        self.inflight.push((id, want.clone()));
        self.requested = Some((id, want));
    }

    /// Asks for the pileup at `pos` (view naming) unless it is already shown or pending.
    pub fn request_evidence(&mut self, contig: &str, pos: u64, filter: tracks::ReadFilter, next_id: &mut u64) {
        let Some(local) = self.contig_name(contig).map(|c| c.name.clone()) else {
            self.evidence = None;
            return;
        };
        let key = (local.clone(), pos, filter);
        if self.evidence_key.as_ref() == Some(&key) {
            return;
        }
        self.evidence_id = *next_id;
        *next_id += 1;
        self.loader.request_pileup(self.evidence_id, local, pos, filter.pileup_options());
        self.evidence_key = Some(key);
    }

    /// Asks for a summary of `[start, end)` on `contig` (view naming).
    pub fn request_summary(&mut self, contig: &str, start: u64, end: u64, filter: tracks::ReadFilter, next_id: &mut u64) {
        self.summary_for = Some((contig.to_owned(), start, end));
        let Some(c) = self.contig_name(contig) else {
            self.summary = None;
            return;
        };
        let region = Region { contig: c.name.clone(), start, end: end.min(c.length) };
        let id = *next_id;
        *next_id += 1;
        let opts = ggv_query::SummaryOptions { pileup: filter.pileup_options(), max_variable_sites: 10, ..Default::default() };
        self.loader.request_summary(id, region, opts);
        self.summary = Some((id, None));
        self.summary_filter = Some(filter);
    }

    /// Takes finished loads. `view_contig` is the view's contig in this file's naming. Call
    /// [`Track::apply_style`] after [`Track::ensure_loaded`].
    pub fn receive(&mut self, view_contig: Option<&str>, generation: &mut u64) {
        for s in self.loader.summary_rx.try_iter().collect::<Vec<_>>() {
            if let Some((id, slot)) = &mut self.summary
                && *id == s.id
            {
                *slot = Some(s.result);
            }
        }
        for ev in self.loader.evidence_rx.try_iter().collect::<Vec<_>>() {
            if ev.id == self.evidence_id {
                self.evidence = Some(ev);
            }
        }
        let results: Vec<_> = self.loader.rx.try_iter().collect();
        for result in results {
            match result {
                Ok(loaded) => self.accept(loaded, view_contig, generation),
                // A failed load is settled (not loading), and isn't asked for again until the
                // view needs another span: retrying every frame would spin. Only the request it
                // belongs to: an older load failing mustn't cancel the one still in flight.
                Err((id, e)) => {
                    self.inflight.retain(|(i, _)| *i != id);
                    if self.requested.as_ref().is_some_and(|(r, _)| *r == id) {
                        self.failed = self.requested.take().map(|(_, s)| s);
                    }
                    self.error = Some(e);
                }
            }
        }
    }

    /// Rebuilds what's on screen for a new style, unless a newer load is on its way: that one
    /// comes built for the style, and the span on screen may be one being replaced for being far
    /// larger than the view (rebuilding megabases at full detail as the view zooms in).
    pub fn apply_style(&mut self, style: Style, generation: &mut u64) {
        let replacing = matches!((&self.current, &self.requested), (Some(c), Some((id, _))) if *id > c.id);
        if !replacing && self.current.as_ref().is_some_and(|c| c.style != style) {
            self.restyle(style, generation);
        }
    }

    fn accept(&mut self, loaded: Loaded, view_contig: Option<&str>, generation: &mut u64) {
        // Take any result newer than what's on screen, not just the latest request: while the
        // view keeps moving (a fling), every frame supersedes the last request, and waiting for
        // the newest would show nothing until the motion stops.
        if self.current.as_ref().is_some_and(|c| loaded.id <= c.id) {
            return;
        }
        let Some(pos) = self.inflight.iter().position(|(id, _)| *id == loaded.id) else { return };
        let mut span = self.inflight[pos].1.clone();
        // Record what actually came back, not what was asked for: if the loader returned
        // coverage only, the span must not claim reads, or the view would never reload.
        span.reads = loaded.view.reads.is_some();
        if span.reads {
            span.filter = None;
            // A cached span may be larger than, or shifted from, what was asked for.
            span.region = loaded.view.region.clone();
        }
        self.inflight.drain(..=pos);
        let is_latest = self.requested.as_ref().is_some_and(|(id, _)| *id == loaded.id);
        if !is_latest && view_contig != Some(span.region.contig.as_str()) {
            return; // an older load for a contig we've since left
        }
        if is_latest {
            self.requested = None;
        }
        *generation += 1;
        self.loads += 1;
        self.load_ms.push(loaded.elapsed.as_secs_f64() * 1000.0);
        self.cache_hits += u64::from(loaded.cached);
        self.error = None;
        self.failed = None;
        self.current = Some(Current {
            view: loaded.view,
            packed: loaded.packed,
            coverage: loaded.coverage,
            depth_index: loaded.depth_index,
            instances: InstanceSet { generation: *generation, instances: loaded.instances },
            style: loaded.style,
            span,
            elapsed: loaded.elapsed,
            id: loaded.id,
            reference: loaded.reference,
        });
    }

    /// Rebuilds instances from the loaded data after a style change; no disk access. Re-packs
    /// only when sorting or grouping changed (and only then resets the scroll).
    fn restyle(&mut self, style: Style, generation: &mut u64) {
        if let Some(cur) = &mut self.current {
            let refilter = cur.style.filter != style.filter;
            if refilter {
                cur.coverage = Arc::new(tracks::coverage_for(&cur.view, style));
                cur.depth_index = Arc::new(tracks::DepthIndex::of(&cur.coverage));
            }
            if refilter || cur.style.pack_options() != style.pack_options() {
                cur.packed = Arc::new(tracks::pack(&cur.view, style));
                self.scroll_y = 0.0;
            }
            *generation += 1;
            cur.instances = InstanceSet {
                generation: *generation,
                instances: Arc::new(tracks::build(&cur.view, &cur.coverage, &cur.packed, style)),
            };
            cur.style = style;
        }
    }

    /// Reads per splice junction (start, end) in what's loaded, for `style`'s filter: counted
    /// once per load and filter rather than every frame.
    pub fn junctions(&self, style: Style) -> JunctionCounts {
        let Some(cur) = &self.current else { return JunctionCounts::default() };
        let mut cache = self.junction_cache.borrow_mut();
        if let Some((id, f, c)) = cache.as_ref()
            && *id == cur.id
            && *f == style.filter
        {
            return c.clone();
        }
        let mut counts = std::collections::BTreeMap::new();
        if let Some((reads, _)) = &cur.view.reads {
            let expr = style.filter.expr.map(crate::expr::get);
            for r in reads.iter().filter(|r| style.filter.passes(r, expr.as_deref())) {
                for s in &r.ref_skips {
                    *counts.entry((s.pos, s.pos + u64::from(s.len))).or_default() += 1;
                }
            }
        }
        let counts = std::sync::Arc::new(counts);
        *cache = Some((cur.id, style.filter, counts.clone()));
        counts
    }

    /// What is loaded and pending, for debugging load behaviour.
    pub fn debug_state(&self) -> String {
        let span = |s: &Span| format!("{} reads={} bp/bin={:?}", s.region, s.reads, s.bp_per_bin.map(|b| b.round()));
        format!(
            "{}: current [{}] requested [{}] inflight {} loads {} error {:?}",
            self.label(),
            self.current.as_ref().map_or("none".to_owned(), |c| format!(
                "id {} {} reads loaded={} coverage={}",
                c.id,
                span(&c.span),
                c.view.reads.as_ref().map_or(0, |(r, _)| r.len()),
                match &c.view.coverage {
                    ggv_layout::Coverage::PerBase(_) => "per-base",
                    ggv_layout::Coverage::Binned(_) => "binned",
                    ggv_layout::Coverage::Estimated(_) => "estimated",
                    ggv_layout::Coverage::Unavailable => "unavailable",
                }
            )),
            self.requested.as_ref().map_or("none".to_owned(), |(id, s)| format!("id {id} {}", span(s))),
            self.inflight.len(),
            self.loads,
            self.error,
        )
    }

    pub fn row_layout(&self, style: Style) -> Option<tracks::RowLayout> {
        self.current.as_ref().map(|c| tracks::row_layout(&c.packed, style))
    }

    /// How far the read rows can scroll in a section of height `section_h`.
    pub fn max_scroll(&self, style: Style, section_h: f32) -> f32 {
        self.row_layout(style)
            .map_or(0.0, |rl| (rl.content_height() - (section_h - tracks::READS_TOP) + 12.0).max(0.0))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn span(contig: &str, start: u64, end: u64, reads: bool, bp_per_bin: Option<f64>) -> Span {
        Span { region: Region { contig: contig.into(), start, end }, reads, bp_per_bin, filter: None, tag: None, mods: false }
    }

    #[test]
    fn satisfies_needs_cover_contig_and_mode() {
        let have = span("chr1", 1000, 5000, true, None);
        let want = span("chr1", 0, 10_000, true, None);
        assert!(have.satisfies(&want, 2000.0, 3000.0));
        assert!(!have.satisfies(&want, 500.0, 3000.0), "view starts before the loaded span");
        assert!(!have.satisfies(&span("chr2", 0, 10_000, true, None), 2000.0, 3000.0));
        assert!(!have.satisfies(&span("chr1", 0, 10_000, false, Some(10.0)), 2000.0, 3000.0));
        let mut filtered = have.clone();
        filtered.filter = Some(tracks::ReadFilter::default());
        assert!(!filtered.satisfies(&have, 2000.0, 3000.0), "a load with another filter doesn't count");
    }

    #[test]
    fn satisfies_clips_to_contig_end() {
        // A view running past the contig end is covered by a span that reaches the end.
        let have = span("chrM", 0, 16569, true, None);
        let want = span("chrM", 0, 16569, true, None);
        assert!(have.satisfies(&want, 16000.0, 17000.0));
    }

    #[test]
    fn oversized_and_mod_spans() {
        // Zoomed far into a 5 Mb span: reload a small one rather than rebuild it all.
        let have = span("chr1", 0, 5_000_000, true, None);
        assert!(have.satisfies(&span("chr1", 2_000_000, 2_700_000, true, None), 2_200_000.0, 2_400_000.0));
        assert!(!have.satisfies(&span("chr1", 2_490_000, 2_510_000, true, None), 2_499_000.0, 2_501_000.0));
        // Colouring by modifications needs a span loaded with them; one with them serves either.
        let mut mods = span("chr1", 0, 10_000, true, None);
        assert!(!mods.satisfies(&Span { mods: true, ..mods.clone() }, 2000.0, 3000.0));
        mods.mods = true;
        assert!(mods.satisfies(&span("chr1", 0, 10_000, true, None), 2000.0, 3000.0));
    }

    #[test]
    fn binned_resolution() {
        let have = span("chr1", 0, 3_000_000, false, Some(1000.0));
        assert!(have.satisfies(&span("chr1", 0, 3_000_000, false, Some(600.0)), 0.0, 1e6));
        assert!(!have.satisfies(&span("chr1", 0, 3_000_000, false, Some(400.0)), 0.0, 1e6));
    }
}

/// A sample name that is an identifier rather than a name: a UUID, or a long run of hex digits
/// and dashes (as some pipelines write in `SM`).
pub fn opaque_id(s: &str) -> bool {
    let hexish = s.chars().all(|c| c.is_ascii_hexdigit() || c == '-');
    (s.len() == 36 && s.matches('-').count() == 4 && hexish) || (s.len() >= 20 && hexish)
}

#[cfg(test)]
mod name_tests {
    use super::opaque_id;

    #[test]
    fn opaque_ids() {
        assert!(opaque_id("d27993fa-916f-4305-ac18-22a29c1c4b0a"));
        assert!(opaque_id("0123456789abcdef0123456789"));
        assert!(!opaque_id("HG002"));
        assert!(!opaque_id("NA12878-tumor"));
    }
}
