//! Comparing two alignment tracks (tumor/normal, before/after): positions where any allele's
//! fraction differs by more than a threshold, with enough depth in both.
//!
//! On screen this uses the per-base coverage both tracks already loaded, so it costs one pass
//! over the visible columns and no disk. Stepping beyond what is loaded runs a pileup over both
//! files on a background thread, a chunk at a time, until it finds the next difference.

use std::{
    path::PathBuf,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
        mpsc,
    },
};

use eframe::egui;
use ggv_layout::CoverageColumn;
use ggv_model::Region;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Settings {
    /// A difference must be more than this allele fraction, 0.01..1 (1: a full 100% difference).
    pub threshold: f64,
    /// Both tracks need at least this depth, 1..1000.
    pub min_depth: u32,
}

/// The settings the Compare window offers, and agents may set.
pub const THRESHOLDS: std::ops::RangeInclusive<f64> = 0.01..=1.0;
pub const MIN_DEPTHS: std::ops::RangeInclusive<u32> = 1..=1000;

impl Default for Settings {
    fn default() -> Self {
        Self { threshold: 0.2, min_depth: 10 }
    }
}

/// A, C, G, T, N, deletion: counts and depth.
pub type Alleles = ([u32; 6], u32);

pub const ALLELES: [&str; 6] = ["A", "C", "G", "T", "N", "DEL"];

pub fn from_column(c: &CoverageColumn) -> Alleles {
    let k = c.counts;
    ([k[0], k[1], k[2], k[3], k[4], c.deletions], c.depth)
}

fn from_stats(s: &ggv_query::PositionStats) -> Alleles {
    let b = &s.bases;
    ([b.a.total, b.c.total, b.g.total, b.t.total, b.n.total, s.deletions.total], s.depth)
}

/// The largest allele-fraction difference and its allele, if it reaches the threshold. A
/// gained allele and the one it displaced differ equally; the rarer one (the change, not the
/// background) is reported.
pub fn differs(a: Alleles, b: Alleles, s: Settings) -> Option<(f64, usize)> {
    if a.1 < s.min_depth || b.1 < s.min_depth {
        return None;
    }
    let (da, db) = (f64::from(a.1), f64::from(b.1));
    let mut best: Option<(f64, f64, usize)> = None;
    for k in 0..6 {
        let (fa, fb) = (f64::from(a.0[k]) / da, f64::from(b.0[k]) / db);
        let d = (fa - fb).abs();
        let better = match best {
            None => true,
            Some((bd, bsum, _)) => d > bd + 1e-12 || ((d - bd).abs() <= 1e-12 && fa + fb < bsum),
        };
        if better {
            best = Some((d, fa + fb, k));
        }
    }
    // Strictly more than the threshold, as MCP `compare` counts (the tolerance absorbs
    // floating-point noise); at 100%, where nothing is more, a full 100% difference.
    let at_least = (s.threshold + 1e-9).min(1.0 - 1e-9);
    best.filter(|b| b.0 > at_least).map(|(d, _, k)| (d, k))
}

/// A difference: 0-based position, the fraction difference and the allele.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Site {
    pub pos: u64,
    pub diff: f64,
    pub allele: usize,
}

/// Differences between two per-base coverages over `[from, to)`. Each coverage starts at its
/// own region start. Needs reference-aware columns (allele counts are only filled in where the
/// reference base is known).
pub fn sites(a: &[CoverageColumn], a_start: u64, b: &[CoverageColumn], b_start: u64, from: u64, to: u64, s: Settings) -> Vec<Site> {
    let lo = from.max(a_start).max(b_start);
    let hi = to.min(a_start + a.len() as u64).min(b_start + b.len() as u64);
    (lo..hi.max(lo))
        .filter_map(|p| {
            let (ca, cb) = (&a[(p - a_start) as usize], &b[(p - b_start) as usize]);
            differs(from_column(ca), from_column(cb), s).map(|(diff, allele)| Site { pos: p, diff, allele })
        })
        .collect()
}

/// Where Next or Previous goes from base `at`: a difference among `sites` (from the span both
/// tracks have checked, `[lo, hi)` in `checked`), else a scan of the files from the first
/// base not yet examined (forward: from it on; backward: everything before it).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Step {
    Site(u64),
    Scan(u64),
}

pub fn step(sites: &[Site], checked: Option<(u64, u64)>, at: u64, forward: bool) -> Step {
    match checked {
        // The bases after `at` up to `hi` were all checked.
        Some((lo, hi)) if forward && at + 1 >= lo && at < hi => match sites.iter().find(|s| s.pos > at && s.pos < hi) {
            Some(s) => Step::Site(s.pos),
            None => Step::Scan(hi),
        },
        Some((lo, hi)) if !forward && at > lo && at <= hi => match sites.iter().rev().find(|s| s.pos < at && s.pos >= lo) {
            Some(s) => Step::Site(s.pos),
            None => Step::Scan(lo),
        },
        // Nothing usable checked next to `at`: the files, from beside it.
        _ => Step::Scan(if forward { at + 1 } else { at }),
    }
}

/// What a scan found: a site on a contig (view naming), or none to the end of the genome.
pub enum Found {
    Site(String, Site),
    None,
    Error(String),
}

/// A background search for the next (or previous) difference.
pub struct Scan {
    rx: mpsc::Receiver<Found>,
    cancel: Arc<AtomicBool>,
    pub at: Arc<std::sync::Mutex<String>>,
}

const CHUNK: u64 = 250_000;

pub struct Files {
    pub a: PathBuf,
    pub b: PathBuf,
    pub reference: Option<PathBuf>,
    /// Each file's read filter, as its track draws it.
    pub opts: [ggv_query::PileupOptions; 2],
}

impl Scan {
    /// From `from` on the first of `contigs` (view naming, with lengths), onward in genome order
    /// (backward: the reverse). Forward examines `from` itself; backward the bases before it.
    pub fn start(files: Files, contigs: Vec<(String, u64)>, from: u64, forward: bool, settings: Settings, ctx: &egui::Context) -> Self {
        let (tx, rx) = mpsc::channel();
        let cancel = Arc::new(AtomicBool::new(false));
        let at = Arc::new(std::sync::Mutex::new(String::new()));
        let (c2, at2, ctx) = (cancel.clone(), at.clone(), ctx.clone());
        let _ = std::thread::Builder::new().name("ggv-diff".into()).spawn(move || {
            let found = run(&files, &contigs, from, forward, settings, &c2, &at2).unwrap_or_else(|e| Found::Error(e.to_string()));
            if !c2.load(Ordering::Relaxed) {
                let _ = tx.send(found);
            }
            ctx.request_repaint();
        });
        Self { rx, cancel, at }
    }

    pub fn poll(&self) -> Option<Found> {
        self.rx.try_recv().ok()
    }
}

impl Drop for Scan {
    fn drop(&mut self) {
        self.cancel.store(true, Ordering::Relaxed);
    }
}

fn run(files: &Files, contigs: &[(String, u64)], from: u64, forward: bool, s: Settings, cancel: &AtomicBool, at: &std::sync::Mutex<String>) -> anyhow::Result<Found> {
    let mut a = ggv_io::BamFile::open_with_reference(&files.a, files.reference.as_ref())?;
    let mut b = ggv_io::BamFile::open_with_reference(&files.b, files.reference.as_ref())?;
    for (ci, (name, len)) in contigs.iter().enumerate() {
        let (Ok(ra), Ok(rb)) = (a.resolve_str(name), b.resolve_str(name)) else { continue };
        let len = (*len).min(ra.end).min(rb.end);
        // The first contig starts at `from`; the rest are whole.
        let mut cursor = if ci == 0 { from.min(len) } else if forward { 0 } else { len };
        loop {
            if cancel.load(Ordering::Relaxed) {
                return Ok(Found::None);
            }
            let (start, end) = if forward { (cursor, (cursor + CHUNK).min(len)) } else { (cursor.saturating_sub(CHUNK), cursor) };
            if start >= end {
                break;
            }
            *at.lock().unwrap() = format!("{name}:{}", crate::util::commas(start + 1));
            let alleles = |f: &mut ggv_io::BamFile, contig: &str, opts: &ggv_query::PileupOptions| -> anyhow::Result<Vec<Alleles>> {
                let region = Region { contig: contig.to_owned(), start, end };
                let mut out = vec![([0u32; 6], 0u32); (end - start) as usize];
                ggv_query::pileup(f, &region, None, opts, |p| {
                    if let Some(slot) = out.get_mut((p.pos - 1 - start) as usize) {
                        *slot = from_stats(&p);
                    }
                    Ok(())
                })?;
                Ok(out)
            };
            let (xa, xb) = (alleles(&mut a, &ra.contig, &files.opts[0])?, alleles(&mut b, &rb.contig, &files.opts[1])?);
            let mut hits = (0..xa.len()).filter_map(|i| differs(xa[i], xb[i], s).map(|(diff, allele)| Site { pos: start + i as u64, diff, allele }));
            let hit = if forward { hits.next() } else { hits.last() };
            if let Some(site) = hit {
                return Ok(Found::Site(name.clone(), site));
            }
            cursor = if forward { end } else { start };
        }
    }
    Ok(Found::None)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn col(counts: [u32; 5], del: u32) -> CoverageColumn {
        CoverageColumn { depth: counts.iter().sum::<u32>() + del, counts, deletions: del, ref_base: Some(b'A'), alt_fraction: 0.0 }
    }

    #[test]
    fn finds_allele_fraction_differences() {
        let s = Settings::default();
        // 50% T in one, 0% in the other.
        assert_eq!(differs(from_column(&col([10, 0, 0, 10, 0], 0)), from_column(&col([20, 0, 0, 0, 0], 0)), s).map(|d| d.1), Some(3));
        // Exactly at the threshold (1/5 = 20% vs 0%): not counted, it must be more than.
        assert!(differs(from_column(&col([8, 0, 0, 2, 0], 0)), from_column(&col([10, 0, 0, 0, 0], 0)), s).is_none());
        assert!(differs(from_column(&col([79, 0, 0, 21, 0], 0)), from_column(&col([10, 0, 0, 0, 0], 0)), s).is_some());
        // 10% vs 0%: under the threshold.
        assert!(differs(from_column(&col([18, 0, 0, 2, 0], 0)), from_column(&col([20, 0, 0, 0, 0], 0)), s).is_none());
        // Deletions count as an allele.
        assert_eq!(differs(from_column(&col([10, 0, 0, 0, 0], 10)), from_column(&col([20, 0, 0, 0, 0], 0)), s).map(|d| d.1), Some(5));
        // Too shallow in one track.
        assert!(differs(from_column(&col([2, 0, 0, 2, 0], 0)), from_column(&col([20, 0, 0, 0, 0], 0)), s).is_none());
        // At 100% (the slider's top), a full difference counts and anything less doesn't.
        let full = Settings { threshold: 1.0, min_depth: 3 };
        assert_eq!(differs(from_column(&col([0, 0, 0, 10, 0], 0)), from_column(&col([10, 0, 0, 0, 0], 0)), full).map(|d| d.0), Some(1.0));
        assert!(differs(from_column(&col([1, 0, 0, 9, 0], 0)), from_column(&col([10, 0, 0, 0, 0], 0)), full).is_none());
    }

    #[test]
    fn sites_line_up_offset_coverages() {
        let same = col([20, 0, 0, 0, 0], 0);
        let het = col([10, 0, 0, 10, 0], 0);
        // a covers 100..110 with a het at 105; b covers 103..113, all reference.
        let mut a = vec![same; 10];
        a[5] = het;
        let b = vec![same; 10];
        let found = sites(&a, 100, &b, 103, 0, u64::MAX, Settings::default());
        assert_eq!(found.iter().map(|x| x.pos).collect::<Vec<_>>(), vec![105]);
        assert!(sites(&a, 100, &b, 103, 106, 200, Settings::default()).is_empty());
    }

    #[test]
    fn scan_matches_on_screen_differences() {
        // GGV_TEST_BAM against itself has no differences; against GGV_TEST_BAM2 (another sample)
        // the scan's first site must be one the in-memory comparison finds too.
        let (Ok(a), Ok(regions)) = (std::env::var("GGV_TEST_BAM"), std::env::var("GGV_TEST_REGIONS")) else { return };
        let region = regions.split(',').next().unwrap().to_owned();
        let bam = ggv_io::BamFile::open(&a).unwrap();
        let r = bam.resolve_str(&region).unwrap();
        let contigs = vec![(r.contig.clone(), r.end)];
        let cancel = AtomicBool::new(false);
        let at = std::sync::Mutex::new(String::new());
        let opts = ggv_query::PileupOptions::default();
        let files = Files { a: a.clone().into(), b: a.clone().into(), reference: None, opts: [opts.clone(), opts.clone()] };
        // Same file: no differences in the first chunk after the region start.
        let contigs1 = vec![(r.contig.clone(), (r.start + CHUNK).min(r.end))];
        let found = run(&files, &contigs1, r.start + 1, true, Settings::default(), &cancel, &at).unwrap();
        assert!(matches!(found, Found::None));
        let Ok(b) = std::env::var("GGV_TEST_BAM2") else { return };
        let files = Files { a: a.into(), b: b.into(), reference: None, opts: [opts.clone(), opts] };
        if let Found::Site(c, site) = run(&files, &contigs, r.start + 1, true, Settings::default(), &cancel, &at).unwrap() {
            eprintln!("first difference {c}:{} ({:.0}% {})", site.pos + 1, site.diff * 100.0, ALLELES[site.allele]);
            assert!(site.pos > r.start && site.diff >= 0.2);
        }
    }

    #[test]
    fn steps_start_at_the_first_unexamined_base() {
        let site = |pos| Site { pos, diff: 0.5, allele: 3 };
        let sites = [site(150), site(180)];
        // Inside the checked span: the next site there.
        assert_eq!(step(&sites, Some((100, 200)), 120, true), Step::Site(150));
        assert_eq!(step(&sites, Some((100, 200)), 170, false), Step::Site(150));
        // None left in it: the files from its edge, the edge itself included.
        assert_eq!(step(&sites, Some((100, 200)), 190, true), Step::Scan(200));
        assert_eq!(step(&sites, Some((100, 200)), 140, false), Step::Scan(100));
        // Unchecked bases between `at` and the span: the files from beside `at`.
        assert_eq!(step(&[site(50)], Some((100, 200)), 40, true), Step::Scan(41));
        assert_eq!(step(&[], Some((100, 200)), 99, true), Step::Scan(200));
        assert_eq!(step(&[site(250)], Some((100, 200)), 300, false), Step::Scan(300));
        // Nothing usable loaded (binned, no reference, another contig): from beside `at`.
        assert_eq!(step(&sites, None, 120, true), Step::Scan(121));
        assert_eq!(step(&sites, None, 120, false), Step::Scan(120));
    }

    #[test]
    fn scans_use_each_files_filter_and_start_where_asked() {
        let bam = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio/HG002.revio.slices.bam");
        let all = ggv_query::PileupOptions::default();
        let id = ggv_layout::expr::intern("hp == 1").unwrap();
        let hp1 = ggv_query::PileupOptions { read_filter: Some(ggv_layout::expr::predicate(ggv_layout::expr::get(id), "hp == 1")), ..Default::default() };
        let (cancel, at) = (AtomicBool::new(false), std::sync::Mutex::new(String::new()));
        let s = Settings { threshold: 0.2, min_depth: 3 };
        let contigs = vec![("chr20".to_owned(), 31_200_000)];
        let start = 31_150_000;
        let run_with = |opts: [&ggv_query::PileupOptions; 2], from, forward| {
            let files = Files { a: bam.clone(), b: bam.clone(), reference: None, opts: [opts[0].clone(), opts[1].clone()] };
            match run(&files, &contigs, from, forward, s, &cancel, &at).unwrap() {
                Found::Site(_, site) => Some(site.pos),
                Found::None => None,
                Found::Error(e) => panic!("{e}"),
            }
        };
        // The same file and filter: nothing.
        assert_eq!(run_with([&all, &all], start, true), None);
        // Only HP 1 reads in one: that haplotype's differences.
        let first = run_with([&all, &hp1], start, true).expect("a difference");
        // A scan starting at a site examines it; backward, from just after it.
        assert_eq!(run_with([&all, &hp1], first, true), Some(first));
        assert_eq!(run_with([&all, &hp1], first + 1, false), Some(first));
        assert!(run_with([&all, &hp1], first, false).is_none_or(|p| p < first));
        // Each filter applies to its own file: swapped, the same site.
        assert_eq!(run_with([&hp1, &all], start, true), Some(first));
    }
}
