//! Sequence search: find a primer or oligo in the loaded reference, on both strands, allowing a
//! few mismatches, without going online.
//!
//! Matching is bit-parallel (shift-and with one state word per allowed mismatch), so each
//! reference base costs a handful of instructions however long the pattern is (up to 64 bases).
//! The reference is read in 4 Mb chunks spread over all cores; each thread has its own FASTA
//! handle, and memory stays at a chunk per thread.

use std::{
    path::PathBuf,
    sync::{
        Arc,
        atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering},
        mpsc,
    },
};

use eframe::egui;
use ggv_model::Region;

pub const MAX_LEN: usize = 64;
pub const MAX_MISMATCHES: u8 = 3;
/// Stop collecting after this many hits.
pub const MAX_HITS: usize = 5_000;
const CHUNK: u64 = 4 << 20;

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Hit {
    pub contig: String,
    /// 0-based start of the match on the forward strand.
    pub start: u64,
    pub mismatches: u8,
    /// Matches the reverse complement.
    pub reverse: bool,
}

/// IUPAC code → which of A, C, G, T it allows (bits 0..4).
fn allowed(b: u8) -> Option<u8> {
    Some(match b.to_ascii_uppercase() {
        b'A' => 0b0001,
        b'C' => 0b0010,
        b'G' => 0b0100,
        b'T' | b'U' => 0b1000,
        b'R' => 0b0101,
        b'Y' => 0b1010,
        b'S' => 0b0110,
        b'W' => 0b1001,
        b'K' => 0b1100,
        b'M' => 0b0011,
        b'B' => 0b1110,
        b'D' => 0b1101,
        b'H' => 0b1011,
        b'V' => 0b0111,
        b'N' => 0b1111,
        _ => return None,
    })
}

fn complement_set(s: u8) -> u8 {
    // A<->T, C<->G: reverse the four bits.
    ((s & 1) << 3) | ((s & 2) << 1) | ((s & 4) >> 1) | ((s & 8) >> 3)
}

/// A compiled pattern: per text code (A, C, G, T, other), the pattern positions it matches,
/// for the pattern and its reverse complement.
#[derive(Clone, Debug)]
pub struct Pattern {
    pub len: usize,
    fwd: [u64; 5],
    rev: [u64; 5],
    /// Both strands read the same (a palindrome): report forward hits only.
    palindrome: bool,
}

/// The bases of a typed sequence: spaces and hyphens dropped.
pub fn bases(text: &str) -> String {
    text.chars().filter(|c| !c.is_ascii_whitespace() && *c != '-').collect()
}

impl Pattern {
    pub fn new(text: &str) -> Result<Self, String> {
        let seq: Vec<u8> = bases(text).into_bytes();
        if seq.len() < 4 {
            return Err("enter at least 4 bases".into());
        }
        if seq.len() > MAX_LEN {
            return Err(format!("at most {MAX_LEN} bases"));
        }
        let sets: Vec<u8> = seq
            .iter()
            .map(|&b| allowed(b).ok_or_else(|| format!("'{}' is not a base (A C G T, or IUPAC codes)", b as char)))
            .collect::<Result<_, _>>()?;
        let rc: Vec<u8> = sets.iter().rev().map(|&s| complement_set(s)).collect();
        let masks = |sets: &[u8]| {
            let mut m = [0u64; 5];
            for (i, &s) in sets.iter().enumerate() {
                for (code, mask) in m.iter_mut().enumerate().take(4) {
                    if s & (1 << code) != 0 {
                        *mask |= 1 << i;
                    }
                }
                // Reference N (or other) matches only a pattern N.
                if s == 0b1111 {
                    m[4] |= 1 << i;
                }
            }
            m
        };
        Ok(Self { len: seq.len(), fwd: masks(&sets), rev: masks(&rc), palindrome: sets == rc })
    }
}

fn code(b: u8) -> usize {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

/// Matches of `p` with at most `k` mismatches in `seq`, as (offset of the match start in `seq`,
/// mismatches, reverse). Only matches starting before `limit` are reported.
pub fn scan(p: &Pattern, k: u8, seq: &[u8], limit: usize, mut hit: impl FnMut(usize, u8, bool)) {
    let k = usize::from(k.min(MAX_MISMATCHES));
    let high = 1u64 << (p.len - 1);
    let mut f = [0u64; MAX_MISMATCHES as usize + 1];
    let mut r = [0u64; MAX_MISMATCHES as usize + 1];
    for (i, &b) in seq.iter().enumerate() {
        let c = code(b);
        step(&mut f[..=k], p.fwd[c]);
        step(&mut r[..=k], p.rev[c]);
        if (f[k] | r[k]) & high != 0 && i + 1 >= p.len {
            let start = i + 1 - p.len;
            if start >= limit {
                continue;
            }
            if let Some(j) = (0..=k).find(|&j| f[j] & high != 0) {
                hit(start, j as u8, false);
            }
            if !p.palindrome
                && let Some(j) = (0..=k).find(|&j| r[j] & high != 0)
            {
                hit(start, j as u8, true);
            }
        }
    }
}

/// One base of shift-and with substitutions: state `j` holds prefixes matched with `j`
/// mismatches; a mismatch moves a prefix from `j - 1` to `j`.
#[inline(always)]
fn step(s: &mut [u64], mask: u64) {
    let mut prev_old = s[0];
    s[0] = ((s[0] << 1) | 1) & mask;
    for j in 1..s.len() {
        let old = s[j];
        s[j] = (((old << 1) | 1) & mask) | ((prev_old << 1) | 1);
        prev_old = old;
    }
}

enum Msg {
    Hits(Vec<Hit>),
    Error(String),
    Done,
}

/// A search running in the background.
pub struct Search {
    /// The bases searched for (spaces and hyphens dropped).
    pub pattern: String,
    /// Pattern length in bases: every hit covers `start..start + len`.
    pub len: u64,
    pub mismatches: u8,
    pub hits: Vec<Hit>,
    pub done: bool,
    pub error: Option<String>,
    pub truncated: bool,
    rx: mpsc::Receiver<Msg>,
    cancel: Arc<AtomicBool>,
    searched: Arc<AtomicU64>,
    pub total: u64,
    started: std::time::Instant,
    pub elapsed: Option<std::time::Duration>,
    workers: usize,
}

impl Search {
    /// Searches `regions` (any size; split into chunks) of the reference at `path`.
    pub fn start(path: PathBuf, pattern: &str, mismatches: u8, regions: Vec<Region>, ctx: &egui::Context) -> Result<Self, String> {
        let p = Arc::new(Pattern::new(pattern)?);
        let mismatches = mismatches.min(MAX_MISMATCHES).min((p.len / 3) as u8);
        let mut chunks = Vec::new();
        for r in regions {
            let mut s = r.start;
            while s < r.end {
                let e = (s + CHUNK).min(r.end);
                chunks.push((r.contig.clone(), s, e, r.end));
                s = e;
            }
        }
        let total: u64 = chunks.iter().map(|c| c.2 - c.1).sum();
        let chunks = Arc::new(chunks);
        let next = Arc::new(AtomicUsize::new(0));
        let cancel = Arc::new(AtomicBool::new(false));
        let searched = Arc::new(AtomicU64::new(0));
        let found = Arc::new(AtomicUsize::new(0));
        let (tx, rx) = mpsc::channel();
        let workers = std::thread::available_parallelism().map_or(4, |n| n.get()).min(chunks.len().max(1)).min(16);
        for _ in 0..workers {
            let (chunks, next, cancel, searched, found, tx, p, path, ctx) =
                (chunks.clone(), next.clone(), cancel.clone(), searched.clone(), found.clone(), tx.clone(), p.clone(), path.clone(), ctx.clone());
            std::thread::Builder::new()
                .name("ggv-seqsearch".into())
                .spawn(move || {
                    let mut fasta = match ggv_io::ReferenceFasta::open(&path) {
                        Ok(f) => f,
                        Err(e) => {
                            let _ = tx.send(Msg::Error(e.to_string()));
                            let _ = tx.send(Msg::Done);
                            return;
                        }
                    };
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        let Some((contig, s, e, contig_end)) = chunks.get(i) else { break };
                        if cancel.load(Ordering::Relaxed) {
                            break;
                        }
                        // Read past the chunk by the pattern length, so matches across the
                        // boundary are found by exactly one chunk (the one they start in).
                        let read_end = (*e + p.len as u64 - 1).min(*contig_end);
                        let region = Region { contig: contig.clone(), start: *s, end: read_end };
                        match fasta.fetch(&region) {
                            Ok(seq) => {
                                let mut hits = Vec::new();
                                scan(&p, mismatches, &seq, (*e - *s) as usize, |off, mm, reverse| {
                                    hits.push(Hit { contig: contig.clone(), start: *s + off as u64, mismatches: mm, reverse });
                                });
                                searched.fetch_add(*e - *s, Ordering::Relaxed);
                                if found.fetch_add(hits.len(), Ordering::Relaxed) + hits.len() > MAX_HITS {
                                    cancel.store(true, Ordering::Relaxed);
                                }
                                if !hits.is_empty() && tx.send(Msg::Hits(hits)).is_err() {
                                    break;
                                }
                            }
                            Err(e) => {
                                let _ = tx.send(Msg::Error(e.to_string()));
                            }
                        }
                        ctx.request_repaint();
                    }
                    let _ = tx.send(Msg::Done);
                    ctx.request_repaint();
                })
                .map_err(|e| e.to_string())?;
        }
        Ok(Self {
            pattern: bases(pattern),
            len: p.len as u64,
            mismatches,
            hits: Vec::new(),
            done: false,
            error: None,
            truncated: false,
            rx,
            cancel,
            searched,
            total,
            started: std::time::Instant::now(),
            elapsed: None,
            workers,
        })
    }

    /// Takes results from the workers; returns true if anything changed.
    pub fn poll(&mut self) -> bool {
        let mut changed = false;
        for m in self.rx.try_iter() {
            changed = true;
            match m {
                Msg::Hits(h) => self.hits.extend(h),
                Msg::Error(e) => self.error = Some(e),
                Msg::Done => {
                    self.workers -= 1;
                    if self.workers == 0 {
                        self.done = true;
                        self.elapsed = Some(self.started.elapsed());
                        self.truncated = self.cancel.load(Ordering::Relaxed) && self.hits.len() >= MAX_HITS;
                        self.hits.sort();
                        self.hits.truncate(MAX_HITS);
                    }
                }
            }
        }
        changed
    }

    pub fn progress(&self) -> f32 {
        if self.total == 0 { 1.0 } else { self.searched.load(Ordering::Relaxed) as f32 / self.total as f32 }
    }

    pub fn cancel(&self) {
        self.cancel.store(true, Ordering::Relaxed);
    }
}

impl Drop for Search {
    fn drop(&mut self) {
        self.cancel();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn find(pattern: &str, k: u8, text: &str) -> Vec<(usize, u8, bool)> {
        let p = Pattern::new(pattern).unwrap();
        let mut out = Vec::new();
        scan(&p, k, text.as_bytes(), usize::MAX, |o, m, r| out.push((o, m, r)));
        out
    }

    /// Hamming distance search, the slow way.
    fn naive(pattern: &str, k: u8, text: &str) -> Vec<(usize, u8, bool)> {
        let p: Vec<u8> = pattern.bytes().collect();
        let rc: Vec<u8> = p.iter().rev().map(|&b| match b { b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A', x => x }).collect();
        let t = text.as_bytes();
        let dist = |q: &[u8], w: &[u8]| q.iter().zip(w).filter(|(a, b)| !(**a == b'N' || a == b)).count() as u8;
        let mut out = Vec::new();
        for s in 0..=t.len().saturating_sub(p.len()) {
            let w = &t[s..s + p.len()];
            let (df, dr) = (dist(&p, w), dist(&rc, w));
            if df <= k {
                out.push((s, df, false));
            }
            if dr <= k && p != rc {
                out.push((s, dr, true));
            }
        }
        out
    }

    #[test]
    fn whole_genome_speed() {
        // GGV_BENCH_SEQSEARCH=<Mb>: a random genome of that size (uncompressed FASTA in the temp
        // folder), searched end to end like the app does, with 0 and 3 mismatches.
        let Ok(mb) = std::env::var("GGV_BENCH_SEQSEARCH") else { return };
        let mb: u64 = mb.parse().unwrap();
        let dir = std::env::temp_dir().join("ggv-seqsearch-bench");
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join(format!("random{mb}.fa"));
        let contigs = 4u64;
        let len = mb * 1_000_000 / contigs;
        if !fa.exists() {
            use std::io::Write;
            let mut f = std::io::BufWriter::new(std::fs::File::create(&fa).unwrap());
            let mut fai = String::new();
            let mut offset = 0u64;
            let mut x: u64 = 88172645463325252;
            for c in 0..contigs {
                let header = format!(">c{c}\n");
                f.write_all(header.as_bytes()).unwrap();
                offset += header.len() as u64;
                fai.push_str(&format!("c{c}\t{len}\t{offset}\t60\t61\n"));
                let mut line = Vec::with_capacity(61);
                for i in 0..len {
                    x ^= x << 13;
                    x ^= x >> 7;
                    x ^= x << 17;
                    line.push(b"ACGT"[(x >> 33) as usize % 4]);
                    if line.len() == 60 || i + 1 == len {
                        line.push(b'\n');
                        f.write_all(&line).unwrap();
                        offset += line.len() as u64;
                        line.clear();
                    }
                }
            }
            f.flush().unwrap();
            std::fs::write(fa.with_added_extension("fai"), fai).unwrap();
        }
        let regions: Vec<Region> = (0..contigs).map(|c| Region { contig: format!("c{c}"), start: 0, end: len }).collect();
        let ctx = egui::Context::default();
        for (pattern, k) in [("TCCCCAGAATCTGTCCCCATGG", 0), ("TCCCCAGAATCTGTCCCCATGG", 3)] {
            let t = std::time::Instant::now();
            let mut s = Search::start(fa.clone(), pattern, k, regions.clone(), &ctx).unwrap();
            while !s.done {
                s.poll();
                std::thread::sleep(std::time::Duration::from_millis(2));
            }
            let secs = t.elapsed().as_secs_f64();
            eprintln!("{mb} Mb, {k} mismatches: {} hits in {secs:.2} s ({:.0} Mb/s)", s.hits.len(), mb as f64 / secs);
            assert!(s.error.is_none(), "{:?}", s.error);
        }
    }

    #[test]
    fn exact_and_reverse_complement() {
        // GATTACA at 2; its reverse complement TGTAATC at 12.
        assert_eq!(find("GATTACA", 0, "CCGATTACAGGTGTAATCAA"), vec![(2, 0, false), (11, 0, true)]);
    }

    #[test]
    fn mismatches_agree_with_a_naive_search() {
        // A pseudo-random genome, with planted near-copies.
        let mut x: u64 = 12345;
        let mut text: Vec<u8> = (0..20_000)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x % 4) as usize]
            })
            .collect();
        let pat = "ACGTTGCAAGGCTTACCGAT";
        text[100..120].copy_from_slice(pat.as_bytes());
        text[5000..5020].copy_from_slice(b"ACGTTGCAAGCCTTACCGAT"); // 1 mismatch
        text[9000..9020].copy_from_slice(b"ATCGGTAAGCCTTGCAACGT"); // reverse complement, 0 mm
        let text = String::from_utf8(text).unwrap();
        for k in 0..=3 {
            let mut fast = find(pat, k, &text);
            let mut slow = naive(pat, k, &text);
            fast.sort();
            slow.sort();
            assert_eq!(fast, slow, "k = {k}");
        }
        assert!(find(pat, 1, &text).contains(&(5000, 1, false)));
        assert!(find(pat, 0, &text).contains(&(9000, 0, true)));
    }

    #[test]
    fn iupac_and_palindromes_and_limits() {
        // N in the pattern matches anything; an N in the reference matches only N.
        assert_eq!(find("GANTC", 0, "TTGACTCTT"), vec![(2, 0, false)]);
        assert!(find("GACTC", 0, "TTGANTCTT").is_empty());
        // EcoRI site reads the same on both strands: reported once.
        assert_eq!(find("GAATTC", 0, "AAGAATTCAA"), vec![(2, 0, false)]);
        assert!(Pattern::new("ACG").is_err());
        assert!(Pattern::new("ACGTX").is_err());
        assert!(Pattern::new(&"A".repeat(65)).is_err());
        assert_eq!(Pattern::new(&"A".repeat(64)).unwrap().len, 64);
        // Lower case and U are accepted.
        assert_eq!(find("gauuaca", 0, "CCGATTACAGG"), vec![(2, 0, false)]);
    }

    #[test]
    fn formatting_is_not_counted_in_hit_spans() {
        // A tiny indexed FASTA: GGAATGGAATGG at 2, its reverse complement at 20.
        let dir = std::env::temp_dir().join(format!("ggv-seqsearch-fmt-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("ref.fa");
        let seq = "CCGGAATGGAATGGTTTTTTCCATTCCATTCCAAAAAA";
        std::fs::write(&fa, format!(">c
{seq}
")).unwrap();
        std::fs::write(fa.with_added_extension("fai"), format!("c	{}	3	{}	{}
", seq.len(), seq.len(), seq.len() + 1)).unwrap();
        let ctx = egui::Context::default();
        let regions = vec![Region { contig: "c".into(), start: 0, end: seq.len() as u64 }];
        let run = |pattern: &str| {
            let mut s = Search::start(fa.clone(), pattern, 0, regions.clone(), &ctx).unwrap();
            while !s.done {
                s.poll();
                std::thread::sleep(std::time::Duration::from_millis(2));
            }
            s
        };
        let plain = run("GGAATGGAATGG");
        assert_eq!(plain.len, 12);
        assert_eq!(plain.hits.iter().map(|h| (h.start, h.reverse)).collect::<Vec<_>>(), [(2, false), (20, true)]);
        for formatted in ["GGAA TGGA ATGG", "GGAA-TGGA-ATGG", " GGAA	TGGA- ATGG
"] {
            let s = run(formatted);
            assert_eq!((s.len, s.pattern.as_str()), (12, "GGAATGGAATGG"), "{formatted:?}");
            assert_eq!(s.hits, plain.hits, "{formatted:?}");
        }
        // 64 bases with formatting is still 64 bases.
        let long = "ACGT ".repeat(16);
        assert_eq!(Pattern::new(&long).unwrap().len, 64);
        assert_eq!(bases(&long).len(), 64);
        let _ = std::fs::remove_dir_all(dir);
    }
}
