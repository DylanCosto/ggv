//! Finding a read by name across alignment files.
//!
//! Names aren't indexed, so this scans: the current view first (instant if the read is on
//! screen), then the current contig, then every other contig, several contigs at once on as
//! many threads as there are cores, each with its own file handle. The first alignment found
//! carries its mate's position and its split parts (`SA`), so the scan stops there. Unmapped
//! reads without a position (at the end of a BAM) are not searched.

use std::{
    path::PathBuf,
    sync::{
        Arc, Mutex,
        atomic::{AtomicBool, AtomicUsize, Ordering},
        mpsc,
    },
};

use eframe::egui;
use ggv_model::Region;

/// Where an alignment of the read is, from its SAM line.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Found {
    pub file: PathBuf,
    pub contig: String,
    /// 0-based.
    pub pos: u64,
    pub flags: u16,
    pub mapq: u8,
    pub mate: Option<(String, u64)>,
    /// `SA:Z:` parts: contig, 0-based position.
    pub splits: Vec<(String, u64)>,
    pub sam: String,
}

/// Reads the fields we need from a SAM line.
pub fn from_sam(file: PathBuf, sam: String) -> Option<Found> {
    let f: Vec<&str> = sam.split('\t').collect();
    let flags: u16 = f.get(1)?.parse().ok()?;
    let contig = f.get(2)?.to_string();
    let pos: u64 = f.get(3)?.parse().ok()?;
    let mapq = f.get(4)?.parse().unwrap_or(255);
    let mate = match (f.get(6).copied(), f.get(7).and_then(|p| p.parse::<u64>().ok())) {
        (Some(m), Some(p)) if m != "*" && p > 0 && flags & 0x8 == 0 => Some((if m == "=" { contig.clone() } else { m.to_owned() }, p - 1)),
        _ => None,
    };
    let splits = f
        .iter()
        .skip(11)
        .find_map(|t| t.strip_prefix("SA:Z:"))
        .map(|sa| ggv_layout::parse_sa(sa).into_iter().map(|s| (s.contig, s.pos)).collect())
        .unwrap_or_default();
    (pos > 0).then(|| Found { file, contig, pos: pos - 1, flags, mapq, mate, splits, sam })
}

pub enum Msg {
    Found(Found),
    /// A file was fully searched without a hit.
    NotFound,
    Error(String),
}

pub struct Search {
    pub name: String,
    pub found: Option<Found>,
    pub done: bool,
    pub error: Option<String>,
    rx: mpsc::Receiver<Msg>,
    cancel: Arc<AtomicBool>,
    pending: usize,
    /// Contigs searched so far, out of `total`.
    progress: Arc<AtomicUsize>,
    pub total: usize,
    pub at: Arc<Mutex<String>>,
}

impl Search {
    /// Searches each file (with its contigs, view region first) for `name`.
    pub fn start(name: String, files: Vec<(PathBuf, Vec<Region>)>, reference: Option<PathBuf>, ctx: &egui::Context) -> Self {
        let (tx, rx) = mpsc::channel();
        let cancel = Arc::new(AtomicBool::new(false));
        let progress = Arc::new(AtomicUsize::new(0));
        let at = Arc::new(Mutex::new(String::new()));
        let total = files.iter().map(|f| f.1.len()).sum();
        let pending = files.len();
        let cores = std::thread::available_parallelism().map_or(4, |n| n.get()).min(8);
        for (path, regions) in files {
            // Workers per file share a queue of its regions, the view's first.
            let regions = Arc::new(regions);
            let next = Arc::new(AtomicUsize::new(0));
            let left = Arc::new(AtomicUsize::new(cores.min(regions.len()).max(1)));
            let hit = Arc::new(AtomicBool::new(false));
            for _ in 0..cores.min(regions.len()).max(1) {
                let (tx, cancel, progress, at, ctx, regions, next, left, hit, path, reference, name) = (
                    tx.clone(),
                    cancel.clone(),
                    progress.clone(),
                    at.clone(),
                    ctx.clone(),
                    regions.clone(),
                    next.clone(),
                    left.clone(),
                    hit.clone(),
                    path.clone(),
                    reference.clone(),
                    name.clone(),
                );
                let _ = std::thread::Builder::new().name("ggv-find-read".into()).spawn(move || {
                    let result = (|| -> Result<Option<String>, String> {
                        let mut bam = ggv_io::BamFile::open_with_reference(&path, reference.as_ref()).map_err(|e| e.to_string())?;
                        loop {
                            if cancel.load(Ordering::Relaxed) || hit.load(Ordering::Relaxed) {
                                return Ok(None);
                            }
                            let i = next.fetch_add(1, Ordering::Relaxed);
                            let Some(region) = regions.get(i) else { return Ok(None) };
                            *at.lock().unwrap() = format!("{} in {}", region.contig, crate::files::display_name(&path));
                            let found = bam.find_by_name(region, &name, &cancel).map_err(|e| e.to_string())?;
                            progress.fetch_add(1, Ordering::Relaxed);
                            ctx.request_repaint();
                            if found.is_some() {
                                hit.store(true, Ordering::Relaxed);
                                return Ok(found);
                            }
                        }
                    })();
                    let msg = match result {
                        Ok(Some(sam)) => from_sam(path.clone(), sam).map(Msg::Found),
                        Ok(None) => None,
                        Err(e) => Some(Msg::Error(e)),
                    };
                    let last = left.fetch_sub(1, Ordering::AcqRel) == 1;
                    if let Some(m) = msg {
                        let _ = tx.send(m);
                    }
                    if last {
                        // This file is done (any hit was sent first).
                        let _ = tx.send(Msg::NotFound);
                    }
                    ctx.request_repaint();
                });
            }
        }
        Self { name, found: None, done: pending == 0, error: None, rx, cancel, pending, progress, total, at }
    }

    /// Takes results; the first hit wins and stops every other search.
    pub fn poll(&mut self) -> bool {
        let mut changed = false;
        for m in self.rx.try_iter() {
            changed = true;
            match m {
                Msg::Found(f) => {
                    if self.found.is_none() {
                        self.found = Some(f);
                        self.cancel.store(true, Ordering::Relaxed);
                    }
                }
                Msg::Error(e) => self.error = Some(e),
                Msg::NotFound => {
                    self.pending = self.pending.saturating_sub(1);
                    if self.pending == 0 {
                        self.done = true;
                    }
                }
            }
        }
        if self.found.is_some() {
            self.done = true;
        }
        changed
    }

    pub fn progress(&self) -> (usize, usize) {
        (self.progress.load(Ordering::Relaxed), self.total)
    }
}

impl Drop for Search {
    fn drop(&mut self) {
        self.cancel.store(true, Ordering::Relaxed);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sam_fields() {
        let sam = "r1\t99\tchr1\t1000\t60\t100M\t=\t1200\t300\tACGT\tIIII\tSA:Z:chr5,500,+,50S50M,30,0;\tNM:i:0".to_owned();
        let f = from_sam("a.bam".into(), sam).unwrap();
        assert_eq!((f.contig.as_str(), f.pos, f.flags, f.mapq), ("chr1", 999, 99, 60));
        assert_eq!(f.mate, Some(("chr1".into(), 1199)));
        assert_eq!(f.splits, vec![("chr5".into(), 499)]);
        // Mate unmapped: no mate place.
        let g = from_sam("a.bam".into(), "r\t73\tchr2\t5\t0\t4M\t=\t5\t0\tACGT\t*".into()).unwrap();
        assert_eq!(g.mate, None);
        assert!(from_sam("a.bam".into(), "r\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*".into()).is_none());
    }

    #[test]
    fn finds_a_read_in_a_real_file() {
        // GGV_TEST_BAM + GGV_TEST_REGIONS: take a read from the first region and find it again,
        // searching the whole file.
        let (Ok(bam_path), Ok(regions)) = (std::env::var("GGV_TEST_BAM"), std::env::var("GGV_TEST_REGIONS")) else { return };
        let mut bam = ggv_io::BamFile::open(&bam_path).unwrap();
        let region = bam.resolve_str(regions.split(',').next().unwrap()).unwrap();
        let lines = bam.sam_lines(&region, 50, |_| true).unwrap();
        let name = lines.last().unwrap().split('\t').next().unwrap().to_owned();
        let all: Vec<Region> = bam.contigs().iter().map(|c| Region { contig: c.name.clone(), start: 0, end: c.length }).collect();
        let cancel = AtomicBool::new(false);
        let t = std::time::Instant::now();
        let mut hit = None;
        for r in &all {
            if let Some(s) = bam.find_by_name(r, &name, &cancel).unwrap() {
                hit = Some(s);
                break;
            }
        }
        let f = from_sam(bam_path.into(), hit.expect("found")).unwrap();
        eprintln!("found {name} at {}:{} in {:.2?}", f.contig, f.pos + 1, t.elapsed());
        assert_eq!(f.sam.split('\t').next(), Some(name.as_str()));
    }
}
