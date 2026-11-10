//! Signal tracks: numbers along the genome. bigWig (indexed, with zoom levels for wide views),
//! bedGraph (text, read whole) and SEG (copy-number segments, one row per sample).
//!
//! Every format answers the same question: the mean (and min / max) over each of `n` equal bins
//! of a region, weighting each value by how much of the bin it covers. bigWig answers wide
//! regions from the coarsest zoom level that is still finer than a bin, so a whole chromosome
//! costs a few small reads (as in IGV and UCSC). Format: <https://genome.ucsc.edu/goldenPath/help/bigWig.html>
//! and Kent et al. 2010 (BBI files).

use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, Read, Seek, SeekFrom},
    path::{Path, PathBuf},
};

use ggv_model::{Contig, contig_aliases};

/// One bin of a region.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Bin {
    /// Mean over the bases with data in the bin (weighted by overlap).
    pub mean: f64,
    pub min: f64,
    pub max: f64,
    /// Bases of the bin with data.
    pub covered: u64,
}

/// An interval with a value (0-based, half-open).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Interval {
    pub start: u64,
    pub end: u64,
    pub value: f64,
}

fn err(path: &Path, msg: impl Into<String>) -> crate::Error {
    crate::Error::Io { path: path.to_owned(), cause: std::io::Error::other(msg.into()) }
}

/// Accumulates intervals (or zoom summaries) into equal bins over `start..end`.
struct Binner {
    start: u64,
    end: u64,
    bins: Vec<(f64, f64, f64, f64)>, // (weighted sum, min, max, covered weight)
}

impl Binner {
    fn new(start: u64, end: u64, n: usize) -> Self {
        Self { start, end, bins: vec![(0.0, f64::INFINITY, f64::NEG_INFINITY, 0.0); n.max(1)] }
    }

    /// Adds a value over `[s, e)`; `mean` is the mean over the interval's covered bases,
    /// `covered_frac` the share of the interval with data (1 for raw data).
    fn add(&mut self, s: u64, e: u64, mean: f64, min: f64, max: f64, covered_frac: f64) {
        let (s, e) = (s.max(self.start), e.min(self.end));
        // NaN or infinity is no data: it poisoned the sum, and a NaN-only bin kept min = inf,
        // max = -inf.
        if s >= e || ![mean, min, max, covered_frac].iter().all(|v| v.is_finite()) {
            return;
        }
        let n = self.bins.len() as u64;
        let len = self.end - self.start;
        let first = ((s - self.start) * n / len) as usize;
        let last = (((e - self.start) * n).div_ceil(len) as usize).min(self.bins.len());
        for (k, b) in self.bins.iter_mut().enumerate().take(last).skip(first) {
            let bs = self.start + len * k as u64 / n;
            let be = self.start + len * (k as u64 + 1) / n;
            let overlap = e.min(be).saturating_sub(s.max(bs));
            if overlap == 0 {
                continue;
            }
            let w = overlap as f64 * covered_frac;
            b.0 += mean * w;
            b.1 = b.1.min(min);
            b.2 = b.2.max(max);
            // Summed unrounded: the mean divides by the same weights it was summed with, so
            // it stays within min and max.
            b.3 += w;
        }
    }

    fn finish(self) -> Vec<Option<Bin>> {
        self.bins
            .into_iter()
            .map(|(sum, min, max, w)| (w > 0.0).then(|| Bin { mean: (sum / w).max(min).min(max), min, max, covered: w.round() as u64 }))
            .collect()
    }
}

// ---------- bigWig ----------

struct Zoom {
    reduction: u32,
    index: u64,
}

pub struct BigWig {
    path: PathBuf,
    file: File,
    little: bool,
    /// Name → (id, length).
    chroms: Vec<(String, u32, u64)>,
    zooms: Vec<Zoom>,
    full_index: u64,
    compressed: bool,
}

impl BigWig {
    pub fn open(path: &Path) -> crate::Result<Self> {
        let io = |e: std::io::Error| crate::Error::Io { path: path.to_owned(), cause: e };
        let mut file = File::open(path).map_err(io)?;
        let mut h = [0u8; 64];
        file.read_exact(&mut h).map_err(io)?;
        let little = match h[..4] {
            [0x26, 0xFC, 0x8F, 0x88] => true,
            [0x88, 0x8F, 0xFC, 0x26] => false,
            _ => return Err(err(path, "not a bigWig file")),
        };
        let mut me = Self { path: path.to_owned(), file, little, chroms: Vec::new(), zooms: Vec::new(), full_index: 0, compressed: false };
        let zoom_levels = me.u16(&h[6..8]);
        let chrom_tree = me.u64(&h[8..16]);
        me.full_index = me.u64(&h[24..32]);
        me.compressed = me.u32(&h[52..56]) > 0;
        let mut zh = vec![0u8; 24 * zoom_levels as usize];
        me.file.read_exact(&mut zh).map_err(io)?;
        for z in zh.chunks_exact(24) {
            me.zooms.push(Zoom { reduction: me.u32(&z[0..4]), index: me.u64(&z[16..24]) });
        }
        me.zooms.sort_by_key(|z| z.reduction);
        me.read_chrom_tree(chrom_tree)?;
        Ok(me)
    }

    fn u16(&self, b: &[u8]) -> u16 {
        let a = [b[0], b[1]];
        if self.little { u16::from_le_bytes(a) } else { u16::from_be_bytes(a) }
    }
    fn u32(&self, b: &[u8]) -> u32 {
        let a = b[..4].try_into().unwrap();
        if self.little { u32::from_le_bytes(a) } else { u32::from_be_bytes(a) }
    }
    fn u64(&self, b: &[u8]) -> u64 {
        let a = b[..8].try_into().unwrap();
        if self.little { u64::from_le_bytes(a) } else { u64::from_be_bytes(a) }
    }
    fn f32(&self, b: &[u8]) -> f32 {
        f32::from_bits(self.u32(b))
    }

    fn read_at(&mut self, offset: u64, len: usize) -> crate::Result<Vec<u8>> {
        let mut buf = vec![0u8; len];
        self.file.seek(SeekFrom::Start(offset)).map_err(|e| crate::Error::Io { path: self.path.clone(), cause: e })?;
        self.file.read_exact(&mut buf).map_err(|e| crate::Error::Io { path: self.path.clone(), cause: e })?;
        Ok(buf)
    }

    fn read_chrom_tree(&mut self, offset: u64) -> crate::Result<()> {
        let h = self.read_at(offset, 32)?;
        if self.u32(&h[0..4]) != 0x78CA_8C91 {
            return Err(err(&self.path, "bigWig chromosome tree is damaged"));
        }
        let key_size = self.u32(&h[8..12]) as usize;
        let mut stack = vec![offset + 32];
        while let Some(node) = stack.pop() {
            let nh = self.read_at(node, 4)?;
            let (leaf, count) = (nh[0] == 1, self.u16(&nh[2..4]) as usize);
            let item = key_size + 8;
            let items = self.read_at(node + 4, item * count)?;
            for it in items.chunks_exact(item) {
                if leaf {
                    let name = String::from_utf8_lossy(&it[..key_size]).trim_end_matches('\0').to_owned();
                    self.chroms.push((name, self.u32(&it[key_size..key_size + 4]), u64::from(self.u32(&it[key_size + 4..key_size + 8]))));
                } else {
                    stack.push(self.u64(&it[key_size..key_size + 8]));
                }
            }
        }
        self.chroms.sort_by_key(|c| c.1);
        Ok(())
    }

    pub fn contigs(&self) -> Vec<Contig> {
        self.chroms.iter().map(|(n, _, l)| Contig { name: n.clone(), length: *l }).collect()
    }

    fn chrom(&self, name: &str) -> Option<(u32, u64)> {
        let find = |n: &str| self.chroms.iter().find(|c| c.0 == n).map(|c| (c.1, c.2));
        find(name).or_else(|| contig_aliases(name).iter().find_map(|a| find(a)))
    }

    /// Data blocks (offset, size) whose span overlaps `chrom:start-end`, from an R-tree.
    fn blocks(&mut self, index: u64, chrom: u32, start: u64, end: u64) -> crate::Result<Vec<(u64, u64)>> {
        let h = self.read_at(index, 48)?;
        if self.u32(&h[0..4]) != 0x2468_ACE0 {
            return Err(err(&self.path, "bigWig index is damaged"));
        }
        let mut out = Vec::new();
        let mut stack = vec![index + 48];
        let overlaps = |sc: u32, ss: u32, ec: u32, ee: u32| {
            let a = (u64::from(sc) << 32) | u64::from(ss);
            let b = (u64::from(ec) << 32) | u64::from(ee);
            let (qa, qb) = ((u64::from(chrom) << 32) | start, (u64::from(chrom) << 32) | end);
            a < qb && b > qa
        };
        while let Some(node) = stack.pop() {
            let nh = self.read_at(node, 4)?;
            let (leaf, count) = (nh[0] == 1, self.u16(&nh[2..4]) as usize);
            let item = if leaf { 32 } else { 24 };
            let items = self.read_at(node + 4, item * count)?;
            for it in items.chunks_exact(item) {
                let (sc, ss, ec, ee) = (self.u32(&it[0..4]), self.u32(&it[4..8]), self.u32(&it[8..12]), self.u32(&it[12..16]));
                if !overlaps(sc, ss, ec, ee) {
                    continue;
                }
                if leaf {
                    out.push((self.u64(&it[16..24]), self.u64(&it[24..32])));
                } else {
                    stack.push(self.u64(&it[16..24]));
                }
            }
        }
        out.sort_unstable();
        Ok(out)
    }

    fn block(&mut self, offset: u64, size: u64) -> crate::Result<Vec<u8>> {
        let raw = self.read_at(offset, size as usize)?;
        if !self.compressed {
            return Ok(raw);
        }
        let mut out = Vec::new();
        flate2::read::ZlibDecoder::new(&raw[..]).read_to_end(&mut out).map_err(|e| crate::Error::Io { path: self.path.clone(), cause: e })?;
        Ok(out)
    }

    /// The raw intervals overlapping `name:start-end`, clipped to it.
    pub fn intervals(&mut self, name: &str, start: u64, end: u64) -> crate::Result<Vec<Interval>> {
        let Some((id, len)) = self.chrom(name) else { return Ok(Vec::new()) };
        let end = end.min(len);
        let mut out = Vec::new();
        for (off, size) in self.blocks(self.full_index, id, start, end)? {
            let d = self.block(off, size)?;
            let mut p = 0;
            while p + 24 <= d.len() {
                let (cid, cstart, step, span, ty, count) = (self.u32(&d[p..]), self.u32(&d[p + 4..]), self.u32(&d[p + 12..]), self.u32(&d[p + 16..]), d[p + 20], self.u16(&d[p + 22..]) as usize);
                p += 24;
                let size = match ty {
                    1 => 12,
                    2 => 8,
                    3 => 4,
                    _ => return Err(err(&self.path, format!("bigWig section type {ty}"))),
                };
                for k in 0..count {
                    let it = &d[p + k * size..p + (k + 1) * size];
                    let (s, e, v) = match ty {
                        1 => (self.u32(&it[0..]), self.u32(&it[4..]), self.f32(&it[8..])),
                        2 => (self.u32(&it[0..]), self.u32(&it[0..]) + span, self.f32(&it[4..])),
                        _ => (cstart + k as u32 * step, cstart + k as u32 * step + span, self.f32(&it[0..])),
                    };
                    let (s, e) = (u64::from(s), u64::from(e));
                    if cid == id && e > start && s < end {
                        out.push(Interval { start: s.max(start), end: e.min(end), value: f64::from(v) });
                    }
                }
                p += count * size;
            }
        }
        out.sort_by_key(|i| i.start);
        Ok(out)
    }

    /// `n` bins over `name:start-end`, from the coarsest zoom level at most an eighth of a bin (summaries straddling bin edges are split by overlap) or the
    /// raw data).
    pub fn bins(&mut self, name: &str, start: u64, end: u64, n: usize) -> crate::Result<Vec<Option<Bin>>> {
        let Some((id, len)) = self.chrom(name) else { return Ok(vec![None; n.max(1)]) };
        // Bins always span the requested range (as for bedGraph and SEG), so callers can map
        // bin i to start + i * (end - start) / n; past the contig end they are just empty.
        // Clamping the range here instead stretched every bin and moved values.
        if start >= end || start >= len {
            return Ok(vec![None; n.max(1)]);
        }
        let per_bin = (end - start) / n.max(1) as u64;
        let zoom = self.zooms.iter().rev().find(|z| u64::from(z.reduction) * 8 <= per_bin).map(|z| z.index);
        let mut b = Binner::new(start, end, n);
        let end = end.min(len);
        match zoom {
            None => {
                for i in self.intervals(name, start, end)? {
                    b.add(i.start, i.end, i.value, i.value, i.value, 1.0);
                }
            }
            Some(index) => {
                for (off, size) in self.blocks(index, id, start, end)? {
                    let d = self.block(off, size)?;
                    for r in d.chunks_exact(32) {
                        if self.u32(&r[0..]) != id {
                            continue;
                        }
                        let (s, e, valid) = (u64::from(self.u32(&r[4..])), u64::from(self.u32(&r[8..])), u64::from(self.u32(&r[12..])));
                        if valid == 0 || e <= s {
                            continue;
                        }
                        let (min, max, sum) = (f64::from(self.f32(&r[16..])), f64::from(self.f32(&r[20..])), f64::from(self.f32(&r[24..])));
                        b.add(s, e, sum / valid as f64, min, max, (valid as f64 / (e - s) as f64).min(1.0));
                    }
                }
            }
        }
        Ok(b.finish())
    }
}

// ---------- bedGraph ----------

/// A bedGraph read whole: intervals per contig, in start order.
pub struct BedGraph {
    contigs: HashMap<String, Vec<Interval>>,
    /// Per contig, the furthest end of any row, `NaN` ones too (a contig of only those is in
    /// the file, with no data).
    ends: HashMap<String, u64>,
    /// Rows whose value is `NaN` or infinite, skipped as no data.
    pub non_finite: usize,
}

fn open_text(path: &Path) -> crate::Result<Box<dyn BufRead>> {
    let f = File::open(path).map_err(|e| crate::Error::Io { path: path.to_owned(), cause: e })?;
    let lower = path.to_string_lossy().to_ascii_lowercase();
    Ok(if lower.ends_with(".gz") || lower.ends_with(".bgz") {
        Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(f)))
    } else {
        Box::new(BufReader::new(f))
    })
}

impl BedGraph {
    pub fn open(path: &Path) -> crate::Result<Self> {
        let mut contigs: HashMap<String, Vec<Interval>> = HashMap::new();
        let mut ends: HashMap<String, u64> = HashMap::new();
        let mut non_finite = 0;
        for (n, line) in open_text(path)?.lines().enumerate() {
            let line = line.map_err(|e| crate::Error::Io { path: path.to_owned(), cause: e })?;
            let t = line.trim();
            if t.is_empty() || t.starts_with('#') || t.starts_with("track") || t.starts_with("browser") {
                continue;
            }
            let f: Vec<&str> = t.split_whitespace().collect();
            let parsed = (f.len() >= 4).then(|| Some((f[1].parse::<u64>().ok()?, f[2].parse::<u64>().ok()?, f[3].parse::<f64>().ok()?))).flatten();
            let Some((s, e, v)) = parsed else {
                return Err(err(path, format!("line {}: expected chrom, start, end, value", n + 1)));
            };
            let list = contigs.entry(f[0].to_owned()).or_default();
            let end = ends.entry(f[0].to_owned()).or_default();
            *end = (*end).max(e);
            // `f64::from_str` takes `NaN` and `inf` (ratio tracks over no coverage write them).
            if !v.is_finite() {
                non_finite += 1;
                continue;
            }
            list.push(Interval { start: s, end: e, value: v });
        }
        for v in contigs.values_mut() {
            v.sort_by_key(|i| i.start);
        }
        Ok(Self { contigs, ends, non_finite })
    }

    fn list(&self, name: &str) -> Option<&Vec<Interval>> {
        self.contigs.get(name).or_else(|| contig_aliases(name).iter().find_map(|a| self.contigs.get(a)))
    }

    pub fn contigs(&self) -> Vec<Contig> {
        let mut v: Vec<Contig> = self.ends.iter().map(|(n, &length)| Contig { name: n.clone(), length }).collect();
        v.sort_by(|a, b| a.name.cmp(&b.name));
        v
    }

    pub fn intervals(&self, name: &str, start: u64, end: u64) -> Vec<Interval> {
        let Some(l) = self.list(name) else { return Vec::new() };
        // Intervals don't nest in a bedGraph, so ends are in start order too.
        let first = l.partition_point(|i| i.end <= start);
        l[first..].iter().take_while(|i| i.start < end).map(|i| Interval { start: i.start.max(start), end: i.end.min(end), value: i.value }).collect()
    }

    pub fn bins(&self, name: &str, start: u64, end: u64, n: usize) -> Vec<Option<Bin>> {
        let mut b = Binner::new(start, end.max(start + 1), n);
        for i in self.intervals(name, start, end) {
            b.add(i.start, i.end, i.value, i.value, i.value, 1.0);
        }
        b.finish()
    }
}

// ---------- SEG ----------

/// Copy-number segments (`.seg`): ID, chrom, start, end, [markers,] value; one row per sample.
/// Starts are taken as 1-based (IGV's reading, matching the usual CBS output).
pub struct Seg {
    pub samples: Vec<String>,
    /// Per sample: contig → segments in start order.
    data: Vec<HashMap<String, Vec<Interval>>>,
    /// Per contig, the furthest end of any segment, `NaN` ones too.
    ends: HashMap<String, u64>,
    /// Segments whose value is `NaN` or infinite, skipped as no data.
    pub non_finite: usize,
}

impl Seg {
    pub fn open(path: &Path) -> crate::Result<Self> {
        let mut samples: Vec<String> = Vec::new();
        let mut data: Vec<HashMap<String, Vec<Interval>>> = Vec::new();
        let mut ends: HashMap<String, u64> = HashMap::new();
        // The header is the first row after any blank and `#` lines (IGV's `#track`, `#type=`).
        let mut first_row = true;
        let mut non_finite = 0;
        for (n, line) in open_text(path)?.lines().enumerate() {
            let line = line.map_err(|e| crate::Error::Io { path: path.to_owned(), cause: e })?;
            let f: Vec<&str> = line.split('\t').map(str::trim).collect();
            if line.trim().is_empty() || line.starts_with('#') || f.len() < 5 {
                continue;
            }
            let header = std::mem::replace(&mut first_row, false);
            let (Ok(s), Ok(e)) = (f[2].parse::<u64>(), f[3].parse::<u64>()) else {
                if header {
                    continue;
                }
                return Err(err(path, format!("line {}: expected ID, chrom, start, end, ..., value", n + 1)));
            };
            let Ok(v) = f[f.len() - 1].parse::<f64>() else { continue };
            // A sample (and contig) with only `NaN` segments keeps its row, with no data.
            let k = match samples.iter().position(|x| x == f[0]) {
                Some(k) => k,
                None => {
                    samples.push(f[0].to_owned());
                    data.push(HashMap::new());
                    samples.len() - 1
                }
            };
            let list = data[k].entry(f[1].to_owned()).or_default();
            let end = ends.entry(f[1].to_owned()).or_default();
            *end = (*end).max(e);
            if !v.is_finite() {
                non_finite += 1;
                continue;
            }
            list.push(Interval { start: s.saturating_sub(1), end: e, value: v });
        }
        for m in &mut data {
            for v in m.values_mut() {
                v.sort_by_key(|i| i.start);
            }
        }
        Ok(Self { samples, data, ends, non_finite })
    }

    /// Bins for one sample.
    pub fn bins(&self, sample: usize, name: &str, start: u64, end: u64, n: usize) -> Vec<Option<Bin>> {
        let mut b = Binner::new(start, end.max(start + 1), n);
        if let Some(m) = self.data.get(sample)
            && let Some(l) = m.get(name).or_else(|| contig_aliases(name).iter().find_map(|a| m.get(a)))
        {
            for i in l.iter().filter(|i| i.end > start && i.start < end) {
                b.add(i.start, i.end, i.value, i.value, i.value, 1.0);
            }
        }
        b.finish()
    }
}

// ---------- One interface ----------

pub enum Signal {
    BigWig(Box<BigWig>),
    BedGraph(BedGraph),
    Seg(Seg),
}

/// Whether a path names a signal file: `.bw`, `.bigwig`, `.bedgraph`, `.bg`, `.seg` (gzipped too).
pub fn is_signal(name: &str) -> bool {
    let n = name.to_ascii_lowercase();
    let b = n.trim_end_matches(".gz").trim_end_matches(".bgz");
    [".bw", ".bigwig", ".bedgraph", ".bg", ".seg"].iter().any(|e| b.ends_with(e))
}

/// The failure to read `path` as bedGraph, for a file whose name doesn't say it is a signal track.
fn not_a_signal(path: &Path, e: crate::Error) -> crate::Error {
    let why = match e {
        crate::Error::Io { cause, .. } => cause.to_string(),
        other => other.to_string(),
    };
    let cause = std::io::Error::other(format!("not a signal track by its name (bigWig .bw/.bigwig, bedGraph .bedgraph/.bg, SEG .seg), and it doesn't read as bedGraph text: {why}"));
    crate::Error::Io { path: path.to_owned(), cause }
}

impl Signal {
    pub fn open(path: &Path) -> crate::Result<Self> {
        let n = path.to_string_lossy().to_ascii_lowercase();
        let b = n.trim_end_matches(".gz").trim_end_matches(".bgz");
        if b.ends_with(".bw") || b.ends_with(".bigwig") {
            Ok(Self::BigWig(Box::new(BigWig::open(path)?)))
        } else if b.ends_with(".seg") {
            Ok(Self::Seg(Seg::open(path)?))
        } else {
            // Anything else is read as bedGraph text. When a file not named as a signal track
            // fails to parse, the name is the likelier mistake than the line it failed on.
            BedGraph::open(path).map(Self::BedGraph).map_err(|e| if is_signal(&n) { e } else { not_a_signal(path, e) })
        }
    }

    /// Contigs and lengths (a bigWig's own; for text formats, the furthest end seen).
    pub fn contigs(&self) -> Vec<Contig> {
        match self {
            Self::BigWig(b) => b.contigs(),
            Self::BedGraph(g) => g.contigs(),
            Self::Seg(s) => {
                let mut v: Vec<Contig> = s.ends.iter().map(|(name, &length)| Contig { name: name.clone(), length }).collect();
                v.sort_by(|a, b| a.name.cmp(&b.name));
                v
            }
        }
    }

    /// Resolves locus text (`chr7`, `chr7:1,000-2,000`) against this file's contigs.
    ///
    /// A text format knows only where its data ends, not how long a contig is: a range past the
    /// last value is kept as asked (empty bins there), not cut to it or refused. A bare contig
    /// is its data's extent.
    pub fn resolve(&self, text: &str) -> crate::Result<ggv_model::Region> {
        let contigs = self.contigs();
        if matches!(self, Self::BigWig(_)) {
            return crate::resolve_locus_str(&contigs, text);
        }
        let open: Vec<Contig> = contigs.iter().map(|c| Contig { name: c.name.clone(), length: u64::MAX }).collect();
        let mut r = crate::resolve_locus_str(&open, text)?;
        if r.end == u64::MAX {
            r.end = contigs.iter().find(|c| c.name == r.contig).map_or(0, |c| c.length);
        }
        Ok(r)
    }

    /// What was left out reading a text format: rows with a `NaN` or infinite value.
    pub fn warning(&self) -> Option<String> {
        let n = match self {
            Self::BigWig(_) => 0,
            Self::BedGraph(g) => g.non_finite,
            Self::Seg(s) => s.non_finite,
        };
        (n > 0).then(|| format!("{n} row{} with a NaN or infinite value skipped as no data", if n == 1 { "" } else { "s" }))
    }

    /// Rows: one for bigWig / bedGraph, one per sample for SEG.
    pub fn rows(&self) -> Vec<String> {
        match self {
            Self::Seg(s) => s.samples.clone(),
            _ => vec![String::new()],
        }
    }

    /// `n` bins of `row` over `name:start-end`.
    pub fn bins(&mut self, row: usize, name: &str, start: u64, end: u64, n: usize) -> crate::Result<Vec<Option<Bin>>> {
        Ok(match self {
            Self::BigWig(b) => b.bins(name, start, end, n)?,
            Self::BedGraph(g) => g.bins(name, start, end, n),
            Self::Seg(s) => s.bins(row, name, start, end, n),
        })
    }
}

/// Bins over `start..end` for `bins` asked: at most one per base, so none ends before it starts
/// (or is left empty between two that share a base).
pub fn bin_edges(start: u64, end: u64, bins: usize) -> Vec<(u64, u64)> {
    let len = end.saturating_sub(start).max(1);
    let n = (bins as u64).clamp(1, len);
    (0..n).map(|i| (start + len * i / n, start + len * (i + 1) / n)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bins_weight_by_overlap() {
        let mut b = Binner::new(0, 100, 2);
        b.add(0, 25, 4.0, 4.0, 4.0, 1.0); // a quarter of bin 0 at 4
        b.add(25, 50, 8.0, 8.0, 8.0, 1.0); // the other quarter at 8
        b.add(90, 120, 2.0, 2.0, 2.0, 1.0); // only 90..100 counts
        let r = b.finish();
        assert_eq!(r[0], Some(Bin { mean: 6.0, min: 4.0, max: 8.0, covered: 50 }));
        assert_eq!(r[1], Some(Bin { mean: 2.0, min: 2.0, max: 2.0, covered: 10 }));
    }

    #[test]
    fn sparse_zoom_mean_stays_within_min_max() {
        // Zoom summaries with a base or so of data each: their weights are fractions of a base.
        let mut b = Binner::new(0, 100, 1);
        b.add(0, 40, 1.0, 1.0, 1.0, 0.015);
        b.add(40, 80, 1.0, 1.0, 1.0, 0.015);
        let r = b.finish()[0].unwrap();
        assert_eq!((r.mean, r.min, r.max, r.covered), (1.0, 1.0, 1.0, 1));
    }

    #[test]
    fn bedgraph_and_seg() {
        let dir = std::env::temp_dir().join(format!("ggv-signal-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let bg = dir.join("t.bedGraph");
        std::fs::write(&bg, "track type=bedGraph\nchr1\t0\t10\t1.5\nchr1\t10\t20\t3\nchr2\t5\t6\t-1\n").unwrap();
        let mut s = Signal::open(&bg).unwrap();
        let r = s.bins(0, "1", 0, 20, 2).unwrap(); // `1` finds `chr1`
        assert_eq!((r[0].unwrap().mean, r[1].unwrap().mean), (1.5, 3.0));
        let seg = dir.join("t.seg");
        std::fs::write(&seg, "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\nA\t1\t1\t100\t10\t0.5\nB\t1\t51\t100\t5\t-0.8\n").unwrap();
        let mut s = Signal::open(&seg).unwrap();
        assert_eq!(s.rows(), ["A", "B"]);
        let b = s.bins(1, "chr1", 0, 100, 2).unwrap();
        assert_eq!((b[0], b[1].unwrap().mean), (None, -0.8));
        assert!(is_signal("x.bw") && is_signal("x.bedGraph.gz") && is_signal("x.seg") && !is_signal("x.bam"));
        // `#` and blank lines before the header (IGV's `#type=`) were refused at "line 2".
        let commented = dir.join("c.seg");
        std::fs::write(&commented, "#type=COPY_NUMBER\n\nID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\nA\t1\t1\t100\t10\t0.5\n").unwrap();
        let s = Signal::open(&commented).unwrap();
        assert_eq!(s.rows(), ["A"]);
        // A second non-numeric row is still an error.
        std::fs::write(&commented, "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\nA\t1\tx\t100\t10\t0.5\n").unwrap();
        assert!(Signal::open(&commented).is_err());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn non_finite_values_are_no_data() {
        // bughunt4: `NaN` gave `mean -inf  min inf  max -inf`; `inf` a bin of inf; a bin with
        // both averaged to the finite part's mean, counted as covered throughout.
        let dir = std::env::temp_dir().join(format!("ggv-signal-nan-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let bg = dir.join("n.bedGraph");
        std::fs::write(&bg, "20	0	10	1
20	10	20	NaN
20	20	30	3
20	30	40	inf
20	40	50	-Infinity
").unwrap();
        let mut s = Signal::open(&bg).unwrap();
        let b = s.bins(0, "20", 0, 40, 4).unwrap();
        assert_eq!(b[0], Some(Bin { mean: 1.0, min: 1.0, max: 1.0, covered: 10 }));
        assert_eq!((b[1], b[3]), (None, None));
        let one = s.bins(0, "20", 0, 40, 1).unwrap()[0].unwrap();
        assert_eq!((one.mean, one.min, one.max, one.covered), (2.0, 1.0, 3.0, 20));
        assert!(s.warning().unwrap().starts_with("3 rows with a NaN"), "{:?}", s.warning());
        let seg = dir.join("n.seg");
        std::fs::write(&seg, "ID	chrom	loc.start	loc.end	num.mark	seg.mean
A	20	1	10	5	1
A	20	11	20	5	NaN
").unwrap();
        let mut s = Signal::open(&seg).unwrap();
        let one = s.bins(0, "20", 0, 20, 1).unwrap()[0].unwrap();
        assert_eq!((one.mean, one.covered), (1.0, 10));
        assert!(s.warning().unwrap().starts_with("1 row with"));
        // And the binner itself, for bigWig values.
        let mut b = Binner::new(0, 10, 1);
        b.add(0, 5, f64::NAN, f64::NAN, f64::NAN, 1.0);
        assert_eq!(b.finish()[0], None);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn all_non_finite_contig_or_sample_is_there_with_no_data() {
        // bughunt5: skipping NaN rows before registering them made a contig of only NaN "not in
        // this file", and dropped a SEG sample of only NaN (the rows after it moved up).
        let dir = std::env::temp_dir().join(format!("ggv-signal-all-nan-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let bg = dir.join("t.bedGraph");
        std::fs::write(&bg, "chr1	0	100	1
chr1	100	200	NaN
chr2	0	500	NaN
").unwrap();
        let mut s = Signal::open(&bg).unwrap();
        let r = s.resolve("chr2").unwrap();
        assert_eq!((r.contig.as_str(), r.start, r.end), ("chr2", 0, 500));
        assert_eq!(s.contigs().iter().map(|c| (c.name.as_str(), c.length)).collect::<Vec<_>>(), [("chr1", 200), ("chr2", 500)]);
        assert_eq!(s.bins(0, "chr2", 0, 500, 2).unwrap(), [None, None]);
        let seg = dir.join("t.seg");
        std::fs::write(&seg, "ID	chrom	loc.start	loc.end	num.mark	seg.mean
A	chr1	1	100	5	0.5
B	chr1	1	100	5	NaN
C	chr1	1	100	5	-0.3
").unwrap();
        let mut s = Signal::open(&seg).unwrap();
        assert_eq!(s.rows(), ["A", "B", "C"]);
        assert_eq!(s.bins(1, "chr1", 0, 100, 2).unwrap(), [None, None]);
        assert_eq!(s.bins(2, "chr1", 0, 100, 1).unwrap()[0].unwrap().mean, -0.3);
        assert!(s.warning().unwrap().starts_with("1 row with"));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn text_formats_keep_the_region_past_their_last_value() {
        let dir = std::env::temp_dir().join(format!("ggv-signal-past-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let bg = dir.join("t.bedGraph");
        std::fs::write(&bg, "chr1\t0\t100\t5\nchr2\t0\t1000\t1\n").unwrap();
        let seg = dir.join("t.seg");
        std::fs::write(&seg, "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\nA\tchr1\t1\t100\t10\t5\n").unwrap();
        for path in [bg, seg] {
            let mut s = Signal::open(&path).unwrap();
            let span = |t: &str| s.resolve(t).map(|r| (r.start, r.end));
            // Was cut to 1-100, so two bins over chr1:1-200 were two bins over 1-100.
            assert_eq!(span("chr1:1-200").unwrap(), (0, 200));
            // Was "starts past the end of chr1 (length 100)".
            assert_eq!(span("chr1:150-300").unwrap(), (149, 300));
            assert_eq!(span("chr1").unwrap(), (0, 100));
            assert!(span("chr3").is_err());
            let b = s.bins(0, "chr1", 0, 200, 2).unwrap();
            assert_eq!((b[0].unwrap().mean, b[1]), (5.0, None));
        }
        std::fs::remove_dir_all(&dir).ok();
    }
}
