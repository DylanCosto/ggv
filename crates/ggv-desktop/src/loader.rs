//! Loads regions on a background thread so the UI never waits on disk.
//!
//! The UI sends a request whenever the view leaves the loaded span; the worker always serves
//! the newest request, dropping any that were superseded while it was busy, and also builds
//! the GPU instances so the UI thread only uploads them.
//!
//! Read-level loads are kept in a small cache, and while the worker has nothing to do it
//! decodes the spans either side of the last load (the direction of travel first). Panning
//! into them, or going back to a place just left, then needs no disk at all.

use std::{
    path::PathBuf,
    sync::{Arc, mpsc},
    thread,
    time::{Duration, Instant},
};

use eframe::egui;
use ggv_io::{BamFile, ReferenceFasta};
use ggv_layout::{ViewData, ViewOptions, load_view};
use ggv_model::{Contig, Region};

use crate::{gpu::Instance, tracks};

pub struct Request {
    pub id: u64,
    pub region: Region,
    /// The visible part of `region`: a cached span that covers it (with room to pan) will do.
    pub view: (u64, u64),
    pub opts: ViewOptions,
    pub style: tracks::Style,
}

pub struct Loaded {
    pub id: u64,
    pub view: Arc<ViewData>,
    pub instances: Arc<Vec<Instance>>,
    /// Reads packed for the request's style (filtering, sorting and grouping).
    pub packed: Arc<ggv_layout::Packed>,
    /// Coverage for the request's read filter.
    pub coverage: Arc<ggv_layout::Coverage>,
    /// Block maxima of `coverage`, for the view's deepest column.
    pub depth_index: Arc<tracks::DepthIndex>,
    /// The shared reference row for this span, built here rather than on the UI thread.
    pub reference: Option<Arc<Vec<Instance>>>,
    pub style: tracks::Style,
    pub elapsed: Duration,
    /// Served from memory (the cache, or a span decoded ahead), not the disk.
    pub cached: bool,
}

/// The exact pileup at one base, for the evidence panel.
pub struct Evidence {
    pub id: u64,
    pub contig: String,
    pub pos: u64,
    pub stats: Option<ggv_query::PositionStats>,
    pub error: Option<String>,
}

enum Job {
    Load(Request),
    Pileup { id: u64, contig: String, pos: u64, opts: ggv_query::PileupOptions },
    Summary { id: u64, region: ggv_model::Region, opts: ggv_query::SummaryOptions },
}

/// `ggv query summary` over a selected region.
pub struct Summary {
    pub id: u64,
    pub result: Result<ggv_query::RegionSummary, String>,
}

pub struct Loader {
    tx: mpsc::Sender<Job>,
    /// Finished loads; a failure carries its request's id.
    pub rx: mpsc::Receiver<Result<Loaded, (u64, String)>>,
    pub evidence_rx: mpsc::Receiver<Evidence>,
    pub summary_rx: mpsc::Receiver<Summary>,
    pub contigs: Vec<Contig>,
    pub label: String,
    pub bam_path: PathBuf,
    /// A problem with the file that doesn't stop it opening (no BGZF end-of-file marker).
    pub warning: Option<String>,
}

impl Loader {
    /// Opens the files on the calling thread (so errors surface immediately), then hands them
    /// to a worker thread.
    pub fn spawn(bam_path: PathBuf, reference_path: Option<PathBuf>, ctx: egui::Context) -> anyhow::Result<Self> {
        let bam = BamFile::open_with_reference(&bam_path, reference_path.as_ref())?;
        let reference = reference_path.as_ref().map(ReferenceFasta::open).transpose()?;
        // Reference over whole reads, for reads that run past a query's own bases (a
        // `mismatches` filter, and mismatches drawn past a loaded span's edge).
        let whole = reference_path.map(ggv_query::ReadReference::from_path);
        let contigs = bam.contigs().to_vec();
        let warning = bam.warning();
        let info = bam.info();
        let label = info.read_groups.iter().find_map(|rg| rg.sample.clone()).unwrap_or_else(|| {
            bam_path.file_name().unwrap_or_default().to_string_lossy().into_owned()
        });

        let (tx, worker_rx) = mpsc::channel::<Job>();
        let (worker_tx, rx) = mpsc::channel();
        let (evidence_tx, evidence_rx) = mpsc::channel();
        let (summary_tx, summary_rx) = mpsc::channel();
        let worker_label = label.clone();

        thread::Builder::new().name("ggv-loader".into()).spawn(move || {
            let mut bam = bam;
            let mut reference = reference;
            let mut cache = Cache::default();
            // Spans to decode ahead of the view while idle, nearest last.
            let mut ahead: Vec<(Region, ViewOptions)> = Vec::new();
            let mut last_center: Option<f64> = None;
            // The last read-level span sent and how it was packed, so the next packing keeps
            // reads in their rows while panning.
            let mut last_layout: Option<(Arc<ViewData>, Arc<ggv_layout::Packed>, tracks::Style)> = None;
            loop {
                let first = match worker_rx.try_recv() {
                    Ok(job) => job,
                    Err(mpsc::TryRecvError::Disconnected) => break,
                    Err(mpsc::TryRecvError::Empty) => {
                        if let Some((region, opts)) = ahead.pop() {
                            if !cache.covers(&region) {
                                let req = Request { id: 0, view: (region.start, region.end), region, opts, style: tracks::Style::default() };
                                if let Ok(view) = load(&mut bam, reference.as_mut(), whole.as_ref(), &req, &worker_label)
                                    && view.reads.is_some()
                                {
                                    cache.put(Arc::new(view));
                                }
                            }
                            continue;
                        }
                        match worker_rx.recv() {
                            Ok(job) => job,
                            Err(_) => break,
                        }
                    }
                };
                // Take everything queued: pileups all run (they're tiny), loads collapse to the
                // newest (older ones are already out of date).
                let mut latest = None;
                for job in std::iter::once(first).chain(worker_rx.try_iter()) {
                    match job {
                        Job::Load(req) => latest = Some(req),
                        Job::Summary { id, region, mut opts } => {
                            if opts.pileup.reference.is_none() {
                                opts.pileup.reference = whole.clone();
                            }
                            let seq = reference.as_mut().and_then(|r| r.fetch(&region).ok());
                            let result = ggv_query::summarize(&mut bam, &region, seq.as_deref(), &opts).map_err(|e| e.to_string());
                            if summary_tx.send(Summary { id, result }).is_err() {
                                return;
                            }
                            ctx.request_repaint();
                        }
                        Job::Pileup { id, contig, pos, mut opts } => {
                            if opts.reference.is_none() {
                                opts.reference = whole.clone();
                            }
                            let ev = pileup_at(&mut bam, reference.as_mut(), id, contig, pos, &opts);
                            if evidence_tx.send(ev).is_err() {
                                return;
                            }
                            ctx.request_repaint();
                        }
                    }
                }
                let Some(req) = latest else { continue };
                if (req.opts.tag, req.opts.mods) != (cache.tag, cache.mods) {
                    cache.views.clear();
                    (cache.tag, cache.mods) = (req.opts.tag, req.opts.mods);
                    ahead.clear();
                }
                let t = Instant::now();
                let wants_reads = req.region.len() <= req.opts.max_read_window;
                let cached = if wants_reads { cache.get(&req.region, req.view) } else { None };
                let from_memory = cached.is_some();
                let view = match cached {
                    Some(v) => Ok(v),
                    None => load(&mut bam, reference.as_mut(), whole.as_ref(), &req, &worker_label).map(|v| {
                        let v = Arc::new(v);
                        if v.reads.is_some() {
                            cache.put(v.clone());
                        }
                        v
                    }),
                };
                // Plan the spans either side for when the worker is idle, the way the view is
                // going first (popped from the end).
                ahead.clear();
                if let Ok(v) = &view
                    && let Some((reads, _)) = &v.reads
                    && !reads.is_empty()
                    && approx_bytes(v) <= PREFETCH_BYTES
                {
                    let center = (req.view.0 + req.view.1) as f64 / 2.0;
                    let right_first = last_center.is_none_or(|c| center >= c);
                    last_center = Some(center);
                    let contig_len = contig_length(&bam, &v.region.contig);
                    let (near, far) = (neighbour(&v.region, right_first, contig_len), neighbour(&v.region, !right_first, contig_len));
                    ahead.extend(far.map(|r| (r, req.opts.clone())));
                    ahead.extend(near.map(|r| (r, req.opts.clone())));
                }
                let result = view.map(|view| {
                    let hints = match &last_layout {
                        Some((v, p, s)) if v.region.contig == view.region.contig && s.pack_options() == req.style.pack_options() && s.filter == req.style.filter => {
                            v.reads.as_ref().map(|(reads, _)| ggv_layout::row_hints(reads, p)).unwrap_or_default()
                        }
                        _ => ggv_layout::RowHints::default(),
                    };
                    let packed = Arc::new(tracks::pack_sticky(&view, req.style, &hints));
                    if view.reads.is_some() {
                        last_layout = Some((view.clone(), packed.clone(), req.style));
                    }
                    let coverage = tracks::coverage_for(&view, req.style);
                    let depth_index = Arc::new(tracks::DepthIndex::of(&coverage));
                    let instances = tracks::build(&view, &coverage, &packed, req.style);
                    let reference = (view.reference.is_some() && view.reads.is_some())
                        .then(|| Arc::new(tracks::build_reference(&view)));
                    Loaded {
                        id: req.id,
                        view,
                        instances: Arc::new(instances),
                        packed,
                        coverage: Arc::new(coverage),
                        depth_index,
                        reference,
                        style: req.style,
                        elapsed: t.elapsed(),
                        cached: from_memory,
                    }
                });
                let result = result.map_err(|e| (req.id, e));
                if worker_tx.send(result).is_err() {
                    break; // the UI is gone
                }
                ctx.request_repaint();
            }
        })?;

        Ok(Self { tx, rx, evidence_rx, summary_rx, contigs, label, bam_path, warning })
    }

    pub fn request(&self, req: Request) {
        let _ = self.tx.send(Job::Load(req));
    }

    pub fn request_summary(&self, id: u64, region: ggv_model::Region, opts: ggv_query::SummaryOptions) {
        let _ = self.tx.send(Job::Summary { id, region, opts });
    }

    pub fn request_pileup(&self, id: u64, contig: String, pos: u64, opts: ggv_query::PileupOptions) {
        let _ = self.tx.send(Job::Pileup { id, contig, pos, opts });
    }
}

/// The cache keeps at most this many spans, and about this many bytes in all (per track: four
/// deep tracks then hold at most two gigabytes between them).
const CACHE_SPANS: usize = 5;
const CACHE_BYTES: usize = 512 << 20;
/// Spans larger than this are not decoded ahead: the cache must keep the span on view and both
/// neighbours, or the neighbours are decoded only to be dropped.
const PREFETCH_BYTES: usize = CACHE_BYTES / 3;

/// Roughly how much memory a loaded span holds: the reads with everything they own, the
/// per-base coverage and the reference.
pub fn approx_bytes(v: &ViewData) -> usize {
    use ggv_layout::AlignedRead;
    let read = |r: &AlignedRead| {
        std::mem::size_of::<AlignedRead>()
            + r.name.capacity()
            + 16 * (r.blocks.len() + r.deletions.len() + r.ref_skips.len() + r.insertions.len())
            + 16 * (r.mismatches.len() + r.clip_mismatches.len())
            + 8 * (r.mods.len() + r.low_quality.len())
            + 64 * r.supplementary.len()
            + 8 // its place in a packed row
    };
    let reads = v.reads.as_ref().map_or(0, |(rs, _)| rs.iter().map(read).sum());
    let coverage = match &v.coverage {
        ggv_layout::Coverage::PerBase(c) => c.len() * std::mem::size_of::<ggv_layout::CoverageColumn>(),
        ggv_layout::Coverage::Binned(b) | ggv_layout::Coverage::Estimated(b) => b.len() * 16,
        ggv_layout::Coverage::Unavailable => 0,
    };
    reads + coverage + v.reference.as_ref().map_or(0, Vec::len)
}

fn contig_length(bam: &BamFile, contig: &str) -> u64 {
    bam.contigs().iter().find(|c| c.name == contig).map_or(0, |c| c.length)
}

/// The span next to `region`, half its length over: the view leaves a span from its middle
/// third, so this one covers it there with room to keep going.
fn neighbour(region: &Region, right: bool, contig_len: u64) -> Option<Region> {
    let len = region.len();
    let (start, end) = if right {
        if region.end >= contig_len {
            return None;
        }
        let s = region.start + len / 2;
        (s, (s + len).min(contig_len))
    } else {
        if region.start == 0 {
            return None;
        }
        let s = region.start.saturating_sub(len / 2);
        (s, s + len)
    };
    (end > start).then(|| Region { contig: region.contig.clone(), start, end })
}

/// Recent read-level loads, most recently used last.
#[derive(Default)]
struct Cache {
    /// Spans with their approximate size.
    views: Vec<(Arc<ViewData>, usize)>,
    /// The tag the cached reads carry, and whether they keep their modification calls; spans
    /// loaded for another tag, or without calls the view now draws, are no use.
    tag: Option<[u8; 2]>,
    mods: bool,
}

impl Cache {
    /// A cached span covering `view` with room to pan (half the view's width on each side, or
    /// up to the contig's start), the roomiest.
    fn get(&mut self, region: &Region, view: (u64, u64)) -> Option<Arc<ViewData>> {
        // Not one the track would reload at once for being far larger than asked for.
        let max_len = region.len().saturating_mul(crate::track::MAX_OVERSIZE);
        let regions: Vec<&Region> = self.views.iter().map(|v| &v.0.region).filter(|r| r.len() <= max_len).collect();
        let i = pick(&regions, &region.contig, view)?;
        let i = self.views.iter().position(|v| &v.0.region == regions[i]).expect("picked from views");
        let v = self.views.remove(i);
        let out = v.0.clone();
        self.views.push(v);
        Some(out)
    }

    fn covers(&self, region: &Region) -> bool {
        self.views.iter().any(|(v, _)| v.region.contig == region.contig && v.region.start <= region.start && v.region.end >= region.end)
    }

    fn put(&mut self, v: Arc<ViewData>) {
        let bytes = approx_bytes(&v);
        if bytes > CACHE_BYTES {
            return; // larger than the whole budget: not worth evicting everything for
        }
        self.views.retain(|(x, _)| x.region != v.region);
        self.views.push((v, bytes));
        while self.views.len() > CACHE_SPANS || self.views.iter().map(|x| x.1).sum::<usize>() > CACHE_BYTES {
            self.views.remove(0);
        }
    }
}

/// Which of `regions` serves `view`: it must cover the view with half a view width spare on
/// each side (or reach the contig start), and the one with the most room wins.
fn pick(regions: &[&Region], contig: &str, view: (u64, u64)) -> Option<usize> {
    let w = view.1.saturating_sub(view.0);
    regions
        .iter()
        .enumerate()
        .filter(|(_, r)| r.contig == contig && r.start <= view.0 && r.end >= view.1)
        .filter(|(_, r)| (view.0 - r.start >= w / 2 || r.start == 0) && r.end - view.1 >= w / 2)
        .max_by_key(|(_, r)| (view.0 - r.start).min(r.end - view.1))
        .map(|(i, _)| i)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(start: u64, end: u64) -> Region {
        Region { contig: "chr1".into(), start, end }
    }

    #[test]
    fn neighbours_keep_a_panning_view_covered() {
        // A 30 kb span around a 10 kb view; the view leaves it after 10 kb of panning right.
        let span = r(100_000, 130_000);
        let next = neighbour(&span, true, 1_000_000).unwrap();
        assert_eq!(next, r(115_000, 145_000));
        // The view just past the old span's end is served by the neighbour, with room.
        assert_eq!(pick(&[&span, &next], "chr1", (120_001, 130_001)), Some(1));
        let prev = neighbour(&span, false, 1_000_000).unwrap();
        assert_eq!(prev, r(85_000, 115_000));
        assert_eq!(pick(&[&span, &prev], "chr1", (99_999, 109_999)), Some(1));
        // Contig edges.
        assert_eq!(neighbour(&r(0, 30_000), false, 1_000_000), None);
        assert_eq!(neighbour(&r(970_000, 1_000_000), true, 1_000_000), None);
        assert_eq!(neighbour(&r(960_000, 990_000), true, 1_000_000), Some(r(975_000, 1_000_000)));
    }

    #[test]
    fn cached_spans_need_room_around_the_view() {
        let span = r(100_000, 130_000);
        assert_eq!(pick(&[&span], "chr1", (110_000, 120_000)), Some(0));
        assert_eq!(pick(&[&span], "chr1", (102_000, 112_000)), None, "too close to the start");
        assert_eq!(pick(&[&span], "chr2", (110_000, 120_000)), None);
        assert_eq!(pick(&[&r(0, 30_000)], "chr1", (0, 10_000)), Some(0), "the contig start counts as room");
    }

    /// Reads running past a loaded span are decoded against the reference over all of them, so
    /// their mismatches (and a `mismatches` filter) don't stop at the span's edge.
    #[test]
    fn reads_past_the_span_get_their_whole_reference() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio");
        let (bam, fa) = (dir.join("HG002.revio.slices.bam"), dir.join("GRCh38.chr20.slices.fa.gz"));
        if !bam.exists() {
            return;
        }
        let region = Region { contig: "chr20".into(), start: 58_880_000, end: 58_880_200 };
        let req = Request { id: 1, view: (region.start, region.end), region: region.clone(), opts: ViewOptions { max_rows: 0, ..Default::default() }, style: tracks::Style::default() };
        let run = |whole: Option<&ggv_query::ReadReference>| {
            let mut b = BamFile::open_with_reference(&bam, Some(&fa)).unwrap();
            let mut f = ReferenceFasta::open(&fa).unwrap();
            load(&mut b, Some(&mut f), whole, &req, "x").unwrap().reads.unwrap().0
        };
        let (edge, whole) = (run(None), run(Some(&ggv_query::ReadReference::from_path(&fa))));
        assert_eq!(edge.len(), whole.len());
        let inside = |p: u64| (region.start..region.end).contains(&p);
        assert!(edge.iter().all(|r| r.mismatches.iter().all(|m| inside(m.pos))));
        assert!(whole.iter().any(|r| r.mismatches.iter().any(|m| !inside(m.pos))), "no read has a mismatch past the span");
        // Inside the span, the same mismatches either way.
        for (a, b) in edge.iter().zip(&whole) {
            let within: Vec<u64> = b.mismatches.iter().map(|m| m.pos).filter(|&p| inside(p)).collect();
            assert_eq!(a.mismatches.iter().map(|m| m.pos).collect::<Vec<_>>(), within, "{}", a.name);
        }
    }
}

/// `ggv query pileup` at one base (the same engine and numbers, samtools-compatible).
fn pileup_at(
    bam: &mut BamFile,
    reference: Option<&mut ReferenceFasta>,
    id: u64,
    contig: String,
    pos: u64,
    opts: &ggv_query::PileupOptions,
) -> Evidence {
    let region = ggv_model::Region { contig: contig.clone(), start: pos, end: pos + 1 };
    let seq = reference.and_then(|r| r.fetch(&region).ok());
    let mut stats = None;
    let result = ggv_query::pileup(bam, &region, seq.as_deref(), opts, |s| {
        stats = Some(s);
        Ok(())
    });
    Evidence { id, contig, pos, stats, error: result.err().map(|e| e.to_string()) }
}

fn load(bam: &mut BamFile, reference: Option<&mut ReferenceFasta>, whole: Option<&ggv_query::ReadReference>, req: &Request, label: &str) -> Result<ViewData, String> {
    // The reference is only drawn per base, so skip fetching it for wide views.
    let seq = match reference {
        Some(r) if req.region.len() <= req.opts.max_read_window => {
            Some(r.fetch(&req.region).map_err(|e| e.to_string())?)
        }
        _ => None,
    };
    let mut opts = req.opts.clone();
    if opts.filters.reference.is_none() {
        opts.filters.reference = whole.cloned();
    }
    load_view(bam, &req.region, seq, label.to_owned(), &opts).map_err(|e| e.to_string())
}
