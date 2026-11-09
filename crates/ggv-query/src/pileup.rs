use std::io;

use ggv_io::{
    BamFile,
    pasta::pileup::{Entry, Overlaps, Pileup, Read},
};
use ggv_model::Region;
use serde::Serialize;

use crate::{Error, Result};

const FLAG_SECONDARY: u16 = 0x100;
const FLAG_QCFAIL: u16 = 0x200;
const FLAG_DUP: u16 = 0x400;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

const CIGAR_SOFT_CLIP: u32 = 4;
const CIGAR_HARD_CLIP: u32 = 5;

/// BAM stores quality 0xff for every base when qualities are missing (`*` in SAM).
const MISSING_QUALITY: u8 = 0xff;

/// A test each read must also pass, e.g. a filter expression: it gets the read and, when the
/// caller has it, the reference over the queried span (bases, 0-based position of the first).
pub type ReadTest = std::sync::Arc<dyn Fn(&Read, Option<(&[u8], u64)>) -> bool + Send + Sync>;

/// A [`ReadTest`] with the text it came from (for messages and provenance).
#[derive(Clone)]
pub struct ReadPredicate {
    pub text: String,
    pub test: ReadTest,
    /// The test looks at the reference over the whole read (a `mismatches` filter): it is
    /// handed the read's full extent, fetched through [`PileupOptions::reference`] when the
    /// query's own bases don't cover it.
    pub needs_reference: bool,
}

/// Fetches the reference bases of `[start, end)` (0-based) on a contig, or `None`.
pub type FetchReference = std::sync::Arc<dyn Fn(&str, u64, u64) -> Option<Vec<u8>> + Send + Sync>;

/// Reference for read tests that need a whole read: windows fetched as reads need them and
/// kept, so a query makes a few fetches rather than one per read. Chunks of a query on several
/// threads each move along their own window.
#[derive(Clone)]
pub struct ReadReference {
    fetch: FetchReference,
    windows: std::sync::Arc<std::sync::Mutex<Vec<Window>>>,
}

type Window = (std::sync::Arc<str>, u64, std::sync::Arc<Vec<u8>>);

/// Bases fetched past a read's end, so the reads after it reuse the window.
const REFERENCE_PAD: u64 = 64 * 1024;
/// Windows kept at once (about one per thread working on a query).
const REFERENCE_WINDOWS: usize = 32;

impl ReadReference {
    pub fn new(fetch: FetchReference) -> Self {
        Self { fetch, windows: Default::default() }
    }

    /// From an indexed FASTA or 2bit, opened on first use.
    pub fn from_path(path: impl Into<std::path::PathBuf>) -> Self {
        let path = path.into();
        let fasta: std::sync::Mutex<Option<ggv_io::ReferenceFasta>> = Default::default();
        Self::new(std::sync::Arc::new(move |contig, start, end| {
            let mut f = fasta.lock().unwrap_or_else(|e| e.into_inner());
            if f.is_none() {
                *f = Some(ggv_io::ReferenceFasta::open(&path).ok()?);
            }
            // `fetch` matches the FASTA's contig naming and clamps the end to the contig.
            f.as_mut()?.fetch(&Region { contig: contig.to_owned(), start, end }).ok().filter(|s| !s.is_empty())
        }))
    }

    /// Bases from a 0-based position that cover `[start, end)` on `contig` (as far as the
    /// reference goes), or `None` when the reference has none there.
    fn bases(&self, contig: &str, start: u64, end: u64) -> Option<(std::sync::Arc<Vec<u8>>, u64)> {
        let hit = |w: &Window| &*w.0 == contig && w.1 <= start && w.1 + w.2.len() as u64 >= end;
        {
            let ws = self.windows.lock().unwrap_or_else(|e| e.into_inner());
            if let Some(w) = ws.iter().find(|w| hit(w)) {
                return Some((w.2.clone(), w.1));
            }
        }
        let seq = (self.fetch)(contig, start, end.saturating_add(REFERENCE_PAD))?;
        let seq = std::sync::Arc::new(seq);
        let mut ws = self.windows.lock().unwrap_or_else(|e| e.into_inner());
        // Replace the window this thread has moved past: the nearest one starting before.
        let behind = ws.iter().enumerate().filter(|(_, w)| &*w.0 == contig && w.1 <= start).max_by_key(|(_, w)| w.1).map(|(i, _)| i);
        let w: Window = (contig.into(), start, seq.clone());
        match behind {
            Some(i) => ws[i] = w,
            None if ws.len() >= REFERENCE_WINDOWS => ws[0] = w,
            None => ws.push(w),
        }
        Some((seq, start))
    }
}

impl ReadReference {
    /// Bases covering all of `read` (from a 0-based position), when `have` doesn't already:
    /// for decoding a read that runs past the bases a caller fetched for its region.
    pub fn beyond(&self, read: &Read, contig: &str, have: Option<(&[u8], u64)>) -> Option<(std::sync::Arc<Vec<u8>>, u64)> {
        let (start, end) = read_extent(read);
        let covered = have.is_some_and(|(seq, from)| from <= start && from + seq.len() as u64 >= end);
        if covered { None } else { self.bases(contig, start, end) }
    }
}

impl std::fmt::Debug for ReadReference {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("ReadReference")
    }
}

/// The reference span a read's alignment covers, 0-based half-open.
fn read_extent(read: &Read) -> (u64, u64) {
    let start = read.position().max(0) as u64;
    let len: u64 = read
        .cigar()
        .iter()
        .filter(|c| matches!(*c & 0xf, 0 | 2 | 3 | 7 | 8))
        .map(|c| u64::from(c >> 4))
        .sum();
    (start, start + len)
}

impl std::fmt::Debug for ReadPredicate {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "ReadPredicate({:?})", self.text)
    }
}

/// Which reads and bases count. The defaults match what IGV shows and what
/// `samtools mpileup -B -Q 0 -x -d 0` counts: unmapped, secondary, QC-fail and duplicate
/// reads are dropped; everything else counts, and overlapping mates both count.
#[derive(Clone, Debug)]
pub struct PileupOptions {
    /// Skip reads with any of these flags set. Default `0x704`.
    pub exclude_flags: u16,
    /// Skip reads without all of these flags set. Default `0`.
    pub require_flags: u16,
    pub min_mapping_quality: u8,
    /// Bases below this quality are left out of `depth` and the base counts and tallied in
    /// `low_base_quality` instead. Deletions are never filtered by base quality.
    pub min_base_quality: u8,
    /// Also report positions in the region that no read covers (depth 0).
    pub include_zero_depth: bool,
    /// A further test every read must pass (a filter expression), applied after the flag and
    /// MAPQ filters.
    pub read_filter: Option<ReadPredicate>,
    /// Where a read test that needs the reference gets bases the query wasn't handed (reads
    /// reaching past the queried span, or queries given no bases at all).
    pub reference: Option<ReadReference>,
}

impl Default for PileupOptions {
    fn default() -> Self {
        Self {
            exclude_flags: 0x4 | FLAG_SECONDARY | FLAG_QCFAIL | FLAG_DUP,
            require_flags: 0,
            min_mapping_quality: 0,
            min_base_quality: 0,
            include_zero_depth: false,
            read_filter: None,
            reference: None,
        }
    }
}

/// A count split by strand.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub struct StrandCount {
    pub total: u32,
    pub fwd: u32,
    pub rev: u32,
}

impl StrandCount {
    fn add(&mut self, reverse: bool) {
        self.total += 1;
        if reverse {
            self.rev += 1;
        } else {
            self.fwd += 1;
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub struct BaseCounts {
    #[serde(rename = "A")]
    pub a: StrandCount,
    #[serde(rename = "C")]
    pub c: StrandCount,
    #[serde(rename = "G")]
    pub g: StrandCount,
    #[serde(rename = "T")]
    pub t: StrandCount,
    /// `N` and IUPAC ambiguity codes.
    #[serde(rename = "N")]
    pub n: StrandCount,
}

impl BaseCounts {
    fn slot(&mut self, base: u8) -> &mut StrandCount {
        match base {
            b'A' => &mut self.a,
            b'C' => &mut self.c,
            b'G' => &mut self.g,
            b'T' => &mut self.t,
            _ => &mut self.n,
        }
    }

    pub fn get(&self, base: u8) -> StrandCount {
        match base.to_ascii_uppercase() {
            b'A' => self.a,
            b'C' => self.c,
            b'G' => self.g,
            b'T' => self.t,
            _ => self.n,
        }
    }
}

/// The most common allele that differs from the reference (or, with no reference, from the
/// most common allele).
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct AltAllele {
    /// `A`, `C`, `G`, `T`, `N`, or `DEL`.
    pub allele: &'static str,
    pub count: StrandCount,
    /// `count.total / depth`, rounded to 3 places.
    pub fraction: f64,
}

/// Everything at one reference position.
#[derive(Clone, Debug, Default, PartialEq, Serialize)]
pub struct PositionStats {
    /// Shared by every position of a query, rather than copied into each.
    pub contig: std::sync::Arc<str>,
    /// 1-based.
    pub pos: u64,
    /// The reference base, when a reference FASTA was given.
    #[serde(rename = "ref")]
    pub ref_base: Option<char>,
    /// Reads with a base (passing `min_base_quality`) or a deletion here. Reference skips
    /// (`N` in the CIGAR, e.g. RNA-seq introns) are not counted.
    pub depth: u32,
    pub bases: BaseCounts,
    /// Reads with this base deleted.
    pub deletions: StrandCount,
    /// Reads with an insertion right after this base.
    pub insertions_after: StrandCount,
    pub ref_skips: u32,
    pub low_base_quality: u32,
    pub mean_mapping_quality: f64,
    pub mapping_quality_zero: u32,
    pub mean_base_quality: f64,
    /// Reads whose aligned part starts here, right after a soft clip.
    pub soft_clipped_left: u32,
    /// Reads whose aligned part ends here, right before a soft clip.
    pub soft_clipped_right: u32,
    /// Fraction of `depth` that is not the reference base (with no reference, or an `N`
    /// reference base: not the most common allele).
    pub alt_fraction: f64,
    pub top_alt: Option<AltAllele>,
    /// Reads in `depth` that aren't the reference allele (the count behind `alt_fraction`).
    #[serde(skip)]
    pub non_reference: u32,
}

impl PositionStats {
    fn empty(contig: &std::sync::Arc<str>, pos0: u64, ref_base: Option<char>) -> Self {
        Self {
            contig: contig.clone(),
            pos: pos0 + 1,
            ref_base,
            depth: 0,
            bases: BaseCounts::default(),
            deletions: StrandCount::default(),
            insertions_after: StrandCount::default(),
            ref_skips: 0,
            low_base_quality: 0,
            mean_mapping_quality: 0.0,
            mapping_quality_zero: 0,
            mean_base_quality: 0.0,
            soft_clipped_left: 0,
            soft_clipped_right: 0,
            alt_fraction: 0.0,
            top_alt: None,
            non_reference: 0,
        }
    }

    /// `alt_fraction` before rounding, for comparing with thresholds.
    pub fn exact_alt_fraction(&self) -> f64 {
        if self.depth == 0 { 0.0 } else { f64::from(self.non_reference) / f64::from(self.depth) }
    }

    fn finish(&mut self, mapq_sum: u64, mapq_n: u32, bq_sum: u64, bq_n: u32) {
        self.mean_mapping_quality = mean(mapq_sum, mapq_n);
        self.mean_base_quality = mean(bq_sum, bq_n);

        if self.depth == 0 {
            return;
        }

        let alleles = [
            ("A", self.bases.a),
            ("C", self.bases.c),
            ("G", self.bases.g),
            ("T", self.bases.t),
            ("N", self.bases.n),
            ("DEL", self.deletions),
        ];

        // An `N` (or IUPAC) reference base says nothing about what's expected, so it is treated
        // like having no reference: compare against the most common allele.
        let ref_allele: &str = match self.ref_base {
            Some('A') => "A",
            Some('C') => "C",
            Some('G') => "G",
            Some('T') => "T",
            _ => alleles.iter().max_by_key(|(_, c)| c.total).map_or("", |(a, _)| *a),
        };

        let ref_count = alleles
            .iter()
            .find(|(a, _)| *a == ref_allele)
            .map_or(0, |(_, c)| c.total);
        let depth = self.depth as f64;
        self.non_reference = self.depth - ref_count;
        self.alt_fraction = round3(f64::from(self.non_reference) / depth);

        self.top_alt = alleles
            .iter()
            .filter(|(a, c)| *a != ref_allele && c.total > 0)
            .max_by_key(|(_, c)| c.total)
            .map(|(a, c)| AltAllele {
                allele: a,
                count: *c,
                fraction: round3(c.total as f64 / depth),
            });
    }
}

/// How many reads the query saw and why some were left out.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub struct ReadCounts {
    pub seen: u64,
    pub used: u64,
    pub duplicate: u64,
    pub secondary: u64,
    pub supplementary: u64,
    pub qc_fail: u64,
    pub unmapped: u64,
    pub low_mapping_quality: u64,
    pub other_flags: u64,
    /// Rejected by the read filter expression.
    pub filtered_out: u64,
}

impl ReadCounts {
    /// Adds another tally (of different reads) to this one.
    pub fn merge(&mut self, o: &Self) {
        self.seen += o.seen;
        self.used += o.used;
        self.duplicate += o.duplicate;
        self.secondary += o.secondary;
        self.supplementary += o.supplementary;
        self.qc_fail += o.qc_fail;
        self.unmapped += o.unmapped;
        self.low_mapping_quality += o.low_mapping_quality;
        self.other_flags += o.other_flags;
        self.filtered_out += o.filtered_out;
    }

    /// Applies `opts`' read filters (flags, MAPQ) to `read` (on `contig`), counting it as seen
    /// and as used or rejected-and-why. Returns whether to use it. Every GGV layer filters through
    /// this, so queries, renders and the app always agree on which reads count.
    pub fn admit(&mut self, read: &Read, opts: &PileupOptions, contig: &str) -> bool {
        self.admit_in(read, opts, contig, None)
    }

    /// [`Self::admit`], giving the read test the reference the caller has (bases, 0-based
    /// position of the first). A test that needs the whole read gets it from
    /// `opts.reference` when these bases don't cover the read.
    pub fn admit_in(&mut self, read: &Read, opts: &PileupOptions, contig: &str, reference: Option<(&[u8], u64)>) -> bool {
        if !self.admit_fields(read.flags(), read.mapping_quality(), opts) {
            return false;
        }
        if let Some(p) = &opts.read_filter
            && !test_read(p, read, opts, contig, reference)
        {
            self.used -= 1;
            self.filtered_out += 1;
            return false;
        }
        true
    }

    /// [`Self::admit`] from the raw fields, for callers that skip decoding the record. It can't
    /// apply a read test, so callers must decode when `opts.read_filter` is set.
    pub fn admit_fields(&mut self, flags: u16, mapq: u8, opts: &PileupOptions) -> bool {
        self.seen += 1;
        if flags & opts.exclude_flags != 0 {
            self.reject(flags, opts.exclude_flags);
            return false;
        }
        if flags & opts.require_flags != opts.require_flags {
            self.other_flags += 1;
            return false;
        }
        if mapq < opts.min_mapping_quality {
            self.low_mapping_quality += 1;
            return false;
        }
        self.used += 1;
        true
    }

    fn reject(&mut self, flags: u16, excluded: u16) {
        let hit = flags & excluded;
        if hit & 0x4 != 0 {
            self.unmapped += 1;
        } else if hit & FLAG_DUP != 0 {
            self.duplicate += 1;
        } else if hit & FLAG_SECONDARY != 0 {
            self.secondary += 1;
        } else if hit & FLAG_SUPPLEMENTARY != 0 {
            self.supplementary += 1;
        } else if hit & FLAG_QCFAIL != 0 {
            self.qc_fail += 1;
        } else {
            self.other_flags += 1;
        }
    }
}

/// Runs `p` on `read`, with the reference over the whole read when the test needs it.
fn test_read(p: &ReadPredicate, read: &Read, opts: &PileupOptions, contig: &str, reference: Option<(&[u8], u64)>) -> bool {
    if !p.needs_reference {
        return (p.test)(read, reference);
    }
    let (start, end) = read_extent(read);
    let covered = reference.is_some_and(|(seq, from)| from <= start && from + seq.len() as u64 >= end);
    match &opts.reference {
        Some(source) if !covered => match source.bases(contig, start, end) {
            Some((seq, from)) => (p.test)(read, Some((&seq, from))),
            None => (p.test)(read, reference),
        },
        _ => (p.test)(read, reference),
    }
}

/// Runs a pileup over `region`, calling `emit` once per position in coordinate order.
///
/// `reference` is the reference sequence for exactly `region` (as from
/// [`ggv_io::ReferenceFasta::fetch`]); pass `None` to go without.
pub fn pileup(
    bam: &mut BamFile,
    region: &Region,
    reference: Option<&[u8]>,
    opts: &PileupOptions,
    emit: impl FnMut(PositionStats) -> io::Result<()>,
) -> Result<ReadCounts> {
    pileup_counting_from(bam, region, reference.map(|r| (r, region.start)), opts, 0, |_| {}, emit)
}

/// [`pileup`], also handing `on_read` every read it uses (after the filters), in the same pass:
/// for callers that want the reads themselves too, without decoding them twice.
pub fn pileup_with_reads(
    bam: &mut BamFile,
    region: &Region,
    reference: Option<&[u8]>,
    opts: &PileupOptions,
    on_read: impl FnMut(&Read),
    emit: impl FnMut(PositionStats) -> io::Result<()>,
) -> Result<ReadCounts> {
    pileup_counting_from(bam, region, reference.map(|r| (r, region.start)), opts, 0, on_read, emit)
}

/// [`pileup`], counting in the returned [`ReadCounts`] only reads that start at or after
/// `count_from` (all reads still count toward depth): for one chunk of a region split across
/// threads, so a read crossing a chunk boundary is counted once, in the chunk where it starts.
/// `reference` is the bases from a 0-based position on: it may cover more than `region` (the
/// whole of a split region, so read tests see every read whole) or less.
pub(crate) fn pileup_counting_from(
    bam: &mut BamFile,
    region: &Region,
    reference: Option<(&[u8], u64)>,
    opts: &PileupOptions,
    count_from: u64,
    mut on_read: impl FnMut(&Read),
    mut emit: impl FnMut(PositionStats) -> io::Result<()>,
) -> Result<ReadCounts> {
    let mut engine = Pileup::new(usize::MAX, Overlaps::Ignore);
    let mut counts = ReadCounts::default();
    let mut next_pos = region.start;
    let mut pileup_err: Option<Error> = None;
    let contig: std::sync::Arc<str> = region.contig.as_str().into();

    let ref_at = |pos0: u64| -> Option<char> {
        let (seq, from) = reference?;
        seq.get(usize::try_from(pos0.checked_sub(from)?).ok()?).map(|b| *b as char)
    };

    let mut drain = |engine: &mut Pileup, next_pos: &mut u64| -> Result<()> {
        while engine.next_position()? {
            let (_, pos0) = engine.position();
            let pos0 = pos0 as u64;
            if !region.contains(pos0) {
                continue;
            }
            if opts.include_zero_depth {
                while *next_pos < pos0 {
                    emit(PositionStats::empty(&contig, *next_pos, ref_at(*next_pos)))
                        .map_err(Error::Output)?;
                    *next_pos += 1;
                }
            }
            let stats = position_stats(engine, &contig, pos0, ref_at(pos0), opts);
            emit(stats).map_err(Error::Output)?;
            *next_pos = pos0 + 1;
        }
        Ok(())
    };

    bam.for_each_record(region, |record| {
        if pileup_err.is_some() {
            return Ok(());
        }
        let read = Read::new(record.as_ref().to_vec())?;
        let counted = (read.position().max(0) as u64) >= count_from;
        let mut scratch = ReadCounts::default();
        let tally = if counted { &mut counts } else { &mut scratch };
        if !tally.admit_in(&read, opts, &region.contig, reference) {
            return Ok(());
        }
        on_read(&read);

        let result = engine
            .push(Some(read))
            .map_err(Error::from)
            .and_then(|()| drain(&mut engine, &mut next_pos));
        if let Err(e) = result {
            pileup_err = Some(e);
        }
        Ok(())
    })?;

    if let Some(e) = pileup_err {
        return Err(e);
    }

    engine.push(None)?;
    drain(&mut engine, &mut next_pos)?;

    if opts.include_zero_depth {
        while next_pos < region.end {
            emit(PositionStats::empty(&contig, next_pos, ref_at(next_pos)))
                .map_err(Error::Output)?;
            next_pos += 1;
        }
    }

    Ok(counts)
}

fn position_stats(
    engine: &Pileup,
    contig: &std::sync::Arc<str>,
    pos0: u64,
    ref_base: Option<char>,
    opts: &PileupOptions,
) -> PositionStats {
    let mut s = PositionStats::empty(contig, pos0, ref_base);
    let (mut mapq_sum, mut mapq_n, mut bq_sum, mut bq_n) = (0u64, 0u32, 0u64, 0u32);

    for entry in engine.entries() {
        let read = engine.read(entry);

        if entry.is_refskip {
            s.ref_skips += 1;
            continue;
        }

        let reverse = read.is_reverse_complemented();

        if entry.is_del {
            s.deletions.add(reverse);
        } else {
            let qpos = entry.qpos as usize;
            let q = read.quality_scores().get(qpos).copied().unwrap_or(MISSING_QUALITY);
            let q_known = q != MISSING_QUALITY;
            if opts.min_base_quality > 0 && (!q_known || q < opts.min_base_quality) {
                s.low_base_quality += 1;
                continue;
            }
            s.bases.slot(base_letter(read, qpos)).add(reverse);
            if q_known {
                bq_sum += u64::from(q);
                bq_n += 1;
            }
            if entry.indel > 0 {
                s.insertions_after.add(reverse);
            }
        }

        s.depth += 1;
        let mapq = read.mapping_quality();
        if mapq != 0xff {
            mapq_sum += u64::from(mapq);
            mapq_n += 1;
        }
        if mapq == 0 {
            s.mapping_quality_zero += 1;
        }

        let (clip_left, clip_right) = soft_clips(read, entry);
        s.soft_clipped_left += u32::from(clip_left);
        s.soft_clipped_right += u32::from(clip_right);
    }

    s.finish(mapq_sum, mapq_n, bq_sum, bq_n);
    s
}

/// The letter for a read base: `A`, `C`, `G`, `T`, or `N` for anything else.
/// A read with no stored sequence (`SEQ *`, allowed for secondary alignments) or one shorter
/// than its CIGAR has no base here: `N`, as samtools counts it, rather than a panic.
fn base_letter(read: &Read, qpos: usize) -> u8 {
    if qpos >= read.sequence_len() {
        return b'N';
    }
    match read.base(qpos) {
        1 => b'A',
        2 => b'C',
        4 => b'G',
        8 => b'T',
        _ => b'N',
    }
}

/// Whether this position is a read's first/last aligned base with a soft clip beside it.
fn soft_clips(read: &Read, entry: &Entry) -> (bool, bool) {
    let op = |c: &u32| c & 0xf;
    let mut ops = read.cigar().iter().map(op).filter(|o| *o != CIGAR_HARD_CLIP);
    let left = entry.is_head && ops.next() == Some(CIGAR_SOFT_CLIP);
    let mut ops = read.cigar().iter().rev().map(op).filter(|o| *o != CIGAR_HARD_CLIP);
    let right = entry.is_tail && ops.next() == Some(CIGAR_SOFT_CLIP);
    (left, right)
}

fn mean(sum: u64, n: u32) -> f64 {
    if n == 0 { 0.0 } else { round3(sum as f64 / n as f64) }
}

fn round3(x: f64) -> f64 {
    (x * 1000.0).round() / 1000.0
}

#[cfg(test)]
mod tests {
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };

    use super::*;

    /// A raw BAM record: `len` bases aligned (`lenM`) at 0-based `pos`.
    fn read(pos: i32, len: u32) -> Read {
        let mut v = Vec::new();
        v.extend(0i32.to_le_bytes());
        v.extend(pos.to_le_bytes());
        v.extend([2u8, 60]);
        v.extend(0u16.to_le_bytes());
        v.extend(1u16.to_le_bytes());
        v.extend(0u16.to_le_bytes());
        v.extend(len.to_le_bytes());
        v.extend([(-1i32).to_le_bytes(), (-1i32).to_le_bytes(), 0i32.to_le_bytes()].concat());
        v.extend(b"r\0");
        v.extend((len << 4).to_le_bytes());
        v.extend(std::iter::repeat_n(0x11u8, len.div_ceil(2) as usize));
        v.extend(std::iter::repeat_n(30u8, len as usize));
        Read::new(v).unwrap()
    }

    /// Passes when it was handed reference over the whole read.
    fn whole_read_test() -> ReadPredicate {
        let test: ReadTest = Arc::new(|read: &Read, r: Option<(&[u8], u64)>| {
            let (s, e) = read_extent(read);
            r.is_some_and(|(seq, from)| from <= s && from + seq.len() as u64 >= e)
        });
        ReadPredicate { text: "mismatches > 0".into(), test, needs_reference: true }
    }

    #[test]
    fn read_tests_get_the_whole_read() {
        let fetches = Arc::new(AtomicUsize::new(0));
        let n = fetches.clone();
        let fetch: FetchReference = Arc::new(move |_: &str, s: u64, e: u64| {
            n.fetch_add(1, Ordering::Relaxed);
            Some(vec![b'A'; (e - s) as usize])
        });
        let opts = PileupOptions { read_filter: Some(whole_read_test()), reference: Some(ReadReference::new(fetch)), ..Default::default() };
        // The query's own bases cover only 109..110; reads run 100..250 and 120..270.
        let query = *b"G";
        let mut counts = ReadCounts::default();
        assert!(counts.admit_in(&read(100, 150), &opts, "chrQA", Some((&query, 109))));
        assert!(counts.admit_in(&read(120, 150), &opts, "chrQA", Some((&query, 109))));
        // Junction, repeat and SV passes give no bases at all.
        assert!(counts.admit(&read(130, 150), &opts, "chrQA"));
        assert_eq!(counts.used, 3);
        // One window serves the reads after the first.
        assert_eq!(fetches.load(Ordering::Relaxed), 1);
        // Covered by the query's own bases: no fetch.
        let wide = vec![b'A'; 400];
        assert!(counts.admit_in(&read(300, 50), &opts, "chrQA", Some((&wide, 200))));
        assert_eq!(fetches.load(Ordering::Relaxed), 1);
        // `beyond` gives bases only for a read the caller's bases don't cover.
        let source = opts.reference.as_ref().unwrap();
        assert!(source.beyond(&read(300, 50), "chrQA", Some((&wide, 200))).is_none());
        let (seq, from) = source.beyond(&read(100, 150), "chrQA", Some((&query, 109))).unwrap();
        assert!(from <= 100 && from + seq.len() as u64 >= 250);
    }

    #[test]
    fn thresholds_use_unrounded_fractions() {
        // 400 C of 2001 reads against A: 0.19990005, shown as 0.2.
        let mut s = PositionStats::empty(&"chrQ".into(), 100, Some('A'));
        s.depth = 2001;
        s.bases.a.total = 1601;
        s.bases.c.total = 400;
        s.finish(0, 0, 0, 0);
        assert_eq!(s.alt_fraction, 0.2);
        assert!(s.exact_alt_fraction() < 0.2);
    }
}
