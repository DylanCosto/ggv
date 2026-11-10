//! Tandem repeats: how long each read says a repeat is. A read that aligns across the repeat
//! and some flank on both sides gives the repeat's length on its haplotype directly: the read
//! bases between the repeat's two ends (insertions and deletions inside included). Long reads
//! span most repeats whole, so a repeat expansion is read off rather than inferred.

use ggv_io::{BamFile, pasta::pileup::Read};
use ggv_model::Region;
use serde::Serialize;

use crate::{PileupOptions, ReadCounts, Result};

/// Bases a read must align beyond each end of the repeat.
pub const FLANK: u64 = 10;

pub fn repeat_flank() -> u64 {
    FLANK
}

#[derive(Clone, Debug, Serialize)]
pub struct RepeatRead {
    /// Read bases between the repeat's ends.
    pub length: u64,
    pub haplotype: u8,
    pub reverse: bool,
}

#[derive(Clone, Debug, Serialize)]
pub struct RepeatAllele {
    /// Median length of its reads, and their shortest and longest.
    pub length: u64,
    pub min: u64,
    pub max: u64,
    pub reads: u32,
    /// Its reads by HP tag (none, 1, 2).
    pub by_haplotype: [u32; 3],
}

#[derive(Clone, Debug, Serialize)]
pub struct RepeatLengths {
    pub reference_length: u64,
    /// Reads that span the repeat with flank on both sides.
    pub spanning: u32,
    /// Reads that start or end inside it, or don't align across a flank (too short to
    /// measure it).
    pub partial: u32,
    /// Reads that skip the repeat or a flank (`N` in the CIGAR, e.g. an intron): no evidence.
    pub skipping: u32,
    /// One allele per cluster of read lengths (unphased reads), or one per haplotype with the
    /// untagged reads that match neither as their own (two or more); every spanning read is in
    /// one or in `unassigned`.
    pub alleles: Vec<RepeatAllele>,
    /// Untagged reads of phased reads that fit no haplotype allele and are too few (one) to be
    /// an allele of their own.
    pub unassigned: u32,
    /// Median length per HP tag, when reads carry one.
    pub haplotype_lengths: [Option<u64>; 3],
    pub reads: Vec<RepeatRead>,
    pub counts: ReadCounts,
}

const CIGAR_M: u32 = 0;
const CIGAR_I: u32 = 1;
const CIGAR_D: u32 = 2;
const CIGAR_N: u32 = 3;
const CIGAR_S: u32 = 4;
const CIGAR_EQ: u32 = 7;
const CIGAR_X: u32 = 8;

/// The read's query position at reference position `p` (insertions just before `p` count as
/// before it), and the read's aligned span.
fn query_at(read: &Read, p: u64) -> Option<u64> {
    let mut r = read.position().max(0) as u64;
    let mut q = 0u64;
    for &c in read.cigar() {
        let (op, len) = (c & 0xf, u64::from(c >> 4));
        match op {
            CIGAR_M | CIGAR_EQ | CIGAR_X => {
                if p < r + len {
                    return (p >= r).then(|| q + (p - r));
                }
                r += len;
                q += len;
            }
            CIGAR_D | CIGAR_N => {
                if p < r + len {
                    return (p >= r).then_some(q);
                }
                r += len;
            }
            CIGAR_I | CIGAR_S => q += len,
            _ => {}
        }
    }
    None
}

/// Reference bases in `[from, to)` the read aligns a base to, and whether a reference skip
/// (`N`) touches that span.
fn aligned_in(read: &Read, from: u64, to: u64) -> (u64, bool) {
    let mut r = read.position().max(0) as u64;
    let (mut aligned, mut skip) = (0, false);
    for &c in read.cigar() {
        let (op, len) = (c & 0xf, u64::from(c >> 4));
        let overlap = (r + len).min(to).saturating_sub(r.max(from));
        match op {
            CIGAR_M | CIGAR_EQ | CIGAR_X => aligned += overlap,
            CIGAR_N => skip |= overlap > 0,
            _ => {}
        }
        if matches!(op, CIGAR_M | CIGAR_D | CIGAR_N | CIGAR_EQ | CIGAR_X) {
            r += len;
        }
    }
    (aligned, skip)
}

#[derive(Debug, PartialEq, Eq)]
enum Span {
    /// Aligned across the repeat with flank on both sides: its length can be read off.
    Whole,
    /// Starts or ends too close, or doesn't align across a flank.
    Partial,
    /// Skips the repeat or a flank (`N`).
    Skip,
}

/// How a read covers the repeat `[start, end)`: spanning takes bases aligned on both flanks
/// (most of each: a deletion of the repeat may reach a base or two into one) and no reference
/// skip, so a read that skips or deletes the repeat with its flanks says nothing about it.
fn spans(read: &Read, start: u64, end: u64) -> Span {
    let (rs, re) = (read.position().max(0) as u64, aligned_end(read));
    if rs + FLANK > start || re < end + FLANK {
        return Span::Partial;
    }
    let (left, skip_left) = aligned_in(read, start - FLANK, start);
    let (right, skip_right) = aligned_in(read, end, end + FLANK);
    if skip_left || skip_right || aligned_in(read, start, end).1 {
        return Span::Skip;
    }
    if left < FLANK / 2 || right < FLANK / 2 {
        return Span::Partial;
    }
    Span::Whole
}

fn aligned_end(read: &Read) -> u64 {
    let span: u64 = read.cigar().iter().filter(|c| matches!(*c & 0xf, CIGAR_M | CIGAR_D | CIGAR_N | CIGAR_EQ | CIGAR_X)).map(|c| u64::from(c >> 4)).sum();
    read.position().max(0) as u64 + span
}

fn median(v: &mut [u64]) -> u64 {
    v.sort_unstable();
    v[v.len() / 2]
}

/// One allele from a set of reads.
fn allele(reads: &[&RepeatRead]) -> RepeatAllele {
    let mut lengths: Vec<u64> = reads.iter().map(|r| r.length).collect();
    let mut by_haplotype = [0u32; 3];
    for r in reads {
        by_haplotype[usize::from(r.haplotype.min(2))] += 1;
    }
    let length = median(&mut lengths);
    RepeatAllele { length, min: lengths[0], max: *lengths.last().expect("non-empty"), reads: reads.len() as u32, by_haplotype }
}

/// Splits read lengths into alleles, and counts the reads left out of them. Phased reads (most
/// with an HP tag, both haplotypes with two or more) give one allele per haplotype, merged when
/// they agree within 2% (at least 2 bases): read-length noise can hide the gap between two
/// close alleles that the haplotypes still tell apart. Untagged reads join the haplotype allele
/// whose lengths they fall among; those that fit neither (a mosaic expansion, an expanded
/// allele the tagger skipped) are split as unphased reads are, into allele(s) of their own,
/// when there are two or more: one such read is a noisy outlier more often than an allele, and
/// is counted as unassigned. Unphased reads are split by length (see [`best_split`]) into as many
/// alleles as their lengths cluster into, each of two reads or more.
fn alleles(reads: &[RepeatRead]) -> (Vec<RepeatAllele>, u32) {
    if reads.is_empty() {
        return (Vec::new(), 0);
    }
    let hp = |h: u8| reads.iter().filter(|r| r.haplotype == h).collect::<Vec<_>>();
    let (mut one, mut two) = (hp(1), hp(2));
    if (one.len() + two.len()) * 10 >= reads.len() * 7 && one.len() >= 2 && two.len() >= 2 {
        let (a, b) = (allele(&one), allele(&two));
        if a.length.abs_diff(b.length) <= 2.max(a.length.max(b.length) / 50) {
            return (vec![allele(&reads.iter().collect::<Vec<_>>())], 0);
        }
        let fits = |x: &RepeatAllele, l: u64| {
            let tolerance = 2.max(x.length / 50);
            l + tolerance >= x.min && l <= x.max + tolerance
        };
        let mut rest = Vec::new();
        for r in reads.iter().filter(|r| r.haplotype != 1 && r.haplotype != 2) {
            match (fits(&a, r.length), fits(&b, r.length)) {
                (true, true) if r.length.abs_diff(a.length) <= r.length.abs_diff(b.length) => one.push(r),
                (true, false) => one.push(r),
                (_, true) => two.push(r),
                (false, false) => rest.push(r.clone()),
            }
        }
        let mut out = vec![allele(&one), allele(&two)];
        // As unphased reads: an allele needs two reads.
        let unassigned = if rest.len() >= 2 { out.extend(alleles(&rest).0); 0 } else { rest.len() as u32 };
        out.sort_by_key(|x| x.length);
        return (out, unassigned);
    }
    let mut sorted: Vec<&RepeatRead> = reads.iter().collect();
    sorted.sort_by_key(|r| r.length);
    (clusters(&sorted), 0)
}

/// Two alleles are told apart when their medians differ by this many bases at least (or by a
/// tenth of the reads' median length, if that is more).
const MIN_ALLELE_SEPARATION: u64 = 3;

/// Unphased reads (sorted by length) as alleles: split where two groups of two reads or more
/// explain the lengths far better than one, and again inside each group.
fn clusters(sorted: &[&RepeatRead]) -> Vec<RepeatAllele> {
    match best_split(sorted) {
        Some(i) => {
            let mut out = clusters(&sorted[..i]);
            out.extend(clusters(&sorted[i..]));
            out
        }
        None => vec![allele(sorted)],
    }
}

/// Where to cut sorted read lengths into two alleles, if anywhere: the cut that leaves the least
/// spread around the two medians, among those with two reads or more each side whose medians are
/// [`MIN_ALLELE_SEPARATION`] apart and that at least halve the spread around one median.
///
/// The largest gap alone is not enough: a pair of stray long reads opens the widest gap and
/// would be called the second allele, with the real two (a few bases apart, joined by a noisy
/// read or two) merged into the first.
fn best_split(sorted: &[&RepeatRead]) -> Option<usize> {
    let n = sorted.len();
    if n < 4 {
        return None;
    }
    let mut prefix = Vec::with_capacity(n + 1);
    prefix.push(0u64);
    for r in sorted {
        prefix.push(prefix[prefix.len() - 1] + r.length);
    }
    let median_at = |from: usize, to: usize| from + (to - from) / 2;
    // Sum of distances from the median of `sorted[from..to]`.
    let spread = |from: usize, to: usize| {
        let mid = median_at(from, to);
        let m = sorted[mid].length;
        let below = m * (mid - from) as u64 - (prefix[mid] - prefix[from]);
        let above = (prefix[to] - prefix[mid]) - m * (to - mid) as u64;
        below + above
    };
    let whole = spread(0, n);
    let separation = MIN_ALLELE_SEPARATION.max(sorted[median_at(0, n)].length / 10);
    (2..=n - 2)
        .filter(|&i| sorted[i].length > sorted[i - 1].length)
        .filter(|&i| sorted[median_at(i, n)].length - sorted[median_at(0, i)].length >= separation)
        .map(|i| (spread(0, i) + spread(i, n), i))
        .filter(|&(cost, _)| cost * 2 <= whole)
        .min()
        .map(|(_, i)| i)
}

/// Measures the repeat `[start, end)` (0-based) on `contig` in every read that spans it.
pub fn repeat_lengths(bam: &mut BamFile, contig: &str, start: u64, end: u64, opts: &PileupOptions) -> Result<RepeatLengths> {
    let region = Region { contig: contig.to_owned(), start: start.saturating_sub(FLANK), end: end + FLANK };
    let mut counts = ReadCounts::default();
    let (mut reads, mut partial, mut skipping) = (Vec::new(), 0u32, 0u32);
    bam.for_each_record(&region, |record| {
        let read = Read::new(record.as_ref().to_vec())?;
        if !counts.admit(&read, opts, contig) {
            return Ok(());
        }
        match spans(&read, start, end) {
            Span::Whole => {}
            Span::Partial => {
                partial += 1;
                return Ok(());
            }
            Span::Skip => {
                skipping += 1;
                return Ok(());
            }
        }
        if let (Some(a), Some(b)) = (query_at(&read, start), query_at(&read, end)) {
            let haplotype = match crate::haplotype(&read) {
                Some(h @ 1..=2) => h as u8,
                _ => 0,
            };
            reads.push(RepeatRead { length: b.saturating_sub(a), haplotype, reverse: read.is_reverse_complemented() });
        }
        Ok(())
    })?;
    let mut haplotype_lengths = [None; 3];
    for (h, slot) in haplotype_lengths.iter_mut().enumerate() {
        let mut v: Vec<u64> = reads.iter().filter(|r| usize::from(r.haplotype) == h).map(|r| r.length).collect();
        if !v.is_empty() {
            *slot = Some(median(&mut v));
        }
    }
    let (alleles, unassigned) = alleles(&reads);
    Ok(RepeatLengths {
        reference_length: end - start,
        spanning: reads.len() as u32,
        partial,
        skipping,
        alleles,
        unassigned,
        haplotype_lengths,
        reads,
        counts,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A raw BAM record at 0-based `pos` with this CIGAR (`(len, op)`) and all-A bases.
    fn read(pos: i32, cigar: &[(u32, u32)]) -> Read {
        let qlen: u32 = cigar.iter().filter(|(_, op)| matches!(*op, CIGAR_M | CIGAR_I | CIGAR_S | CIGAR_EQ | CIGAR_X)).map(|(l, _)| l).sum();
        let mut v = Vec::new();
        v.extend(0i32.to_le_bytes());
        v.extend(pos.to_le_bytes());
        v.extend([2u8, 60]);
        v.extend(0u16.to_le_bytes());
        v.extend((cigar.len() as u16).to_le_bytes());
        v.extend(0u16.to_le_bytes());
        v.extend(qlen.to_le_bytes());
        v.extend([(-1i32).to_le_bytes(), (-1i32).to_le_bytes(), 0i32.to_le_bytes()].concat());
        v.extend(b"r ");
        for (len, op) in cigar {
            v.extend((len << 4 | op).to_le_bytes());
        }
        v.extend(std::iter::repeat_n(0x11u8, qlen.div_ceil(2) as usize));
        v.extend(std::iter::repeat_n(30u8, qlen as usize));
        Read::new(v).unwrap()
    }

    #[test]
    fn spanning_needs_aligned_flanks() {
        // Repeat at 0-based 149..155 (chrQA:150-155), reads from 100.
        let (s, e) = (149, 155);
        assert_eq!(spans(&read(100, &[(150, CIGAR_M)]), s, e), Span::Whole);
        // The repeat deleted, both flanks aligned: a real zero-length allele.
        assert_eq!(spans(&read(100, &[(49, CIGAR_M), (6, CIGAR_D), (95, CIGAR_M)]), s, e), Span::Whole);
        // An intron over the repeat and both flanks.
        assert_eq!(spans(&read(100, &[(25, CIGAR_M), (100, CIGAR_N), (25, CIGAR_M)]), s, e), Span::Skip);
        // A deletion swallowing the repeat and both flanks.
        assert_eq!(spans(&read(100, &[(30, CIGAR_M), (50, CIGAR_D), (30, CIGAR_M)]), s, e), Span::Partial);
        // Ends inside the right flank.
        assert_eq!(spans(&read(100, &[(60, CIGAR_M)]), s, e), Span::Partial);
    }

    fn r(length: u64, haplotype: u8) -> RepeatRead {
        RepeatRead { length, haplotype, reverse: false }
    }

    /// The alleles, when no read is left out.
    fn kept(reads: &[RepeatRead]) -> Vec<RepeatAllele> {
        let (a, unassigned) = alleles(reads);
        assert_eq!(unassigned, 0);
        a
    }

    #[test]
    fn one_or_two_alleles() {
        // Homozygous: tight around 60.
        let a = kept(&[r(60, 1), r(61, 1), r(59, 2), r(60, 2), r(60, 0)]);
        assert_eq!((a.len(), a[0].length, a[0].reads), (1, 60, 5));
        // Heterozygous expansion: 60 on HP 1, ~180 on HP 2.
        let a = kept(&[r(60, 1), r(61, 1), r(60, 1), r(178, 2), r(182, 2), r(180, 2)]);
        // Close alleles that noisy lengths overlap: the haplotypes still separate them.
        let close = kept(&[r(336, 2), r(345, 2), r(330, 2), r(369, 1), r(352, 1), r(372, 1)]);
        assert_eq!((close.len(), close[0].length, close[1].length), (2, 336, 369));
        assert_eq!(a.len(), 2);
        assert_eq!((a[0].length, a[0].by_haplotype, a[1].length, a[1].min, a[1].max), (60, [0, 3, 0], 180, 178, 182));
        // A single outlier read doesn't make an allele.
        assert_eq!(kept(&[r(60, 0), r(60, 0), r(61, 0), r(300, 0)]).len(), 1);
        // Phased and the same on both haplotypes: one allele.
        assert_eq!(kept(&[r(63, 1), r(63, 1), r(64, 2), r(63, 2)]).len(), 1);
    }

    #[test]
    fn stray_long_reads_do_not_hide_two_close_alleles() {
        // HTT in an unphased HG002 HiFi BAM: CAG tracts of 51 and 57 bases with three noisy reads
        // between them, and two reads at 70-71. The widest gap is before the two long reads;
        // the alleles are 51 and 57 (the read at 54, as near to either, goes with the later).
        let lengths = [vec![51; 9], vec![52, 53, 54], vec![57; 15], vec![70, 71]].concat();
        let reads: Vec<RepeatRead> = lengths.into_iter().map(|l| r(l, 0)).collect();
        let a = kept(&reads);
        assert_eq!(a.iter().map(|x| (x.length, x.reads, x.min, x.max)).collect::<Vec<_>>(), vec![(51, 11, 51, 53), (57, 16, 54, 57), (71, 2, 70, 71)]);
        // One stray read is never an allele, however far away, nor are two beside a noisy one.
        let one_stray: Vec<RepeatRead> = [vec![60; 20], vec![61, 300]].concat().into_iter().map(|l| r(l, 0)).collect();
        assert_eq!(kept(&one_stray).len(), 1);
        // Noise around one length stays one allele.
        let noisy: Vec<RepeatRead> = [58, 59, 60, 60, 60, 61, 62].into_iter().map(|l| r(l, 0)).collect();
        assert_eq!(kept(&noisy).len(), 1);
    }

    #[test]
    fn untagged_reads_stay_in_the_alleles() {
        // Mostly phased (25 on HP 1, 20 on HP 2), and four untagged reads at the reference's 40:
        // a third allele, not left out.
        let mut reads: Vec<RepeatRead> = (0..13).map(|_| r(25, 1)).chain((0..15).map(|_| r(20, 2))).collect();
        reads.extend((0..4).map(|_| r(40, 0)));
        let a = kept(&reads);
        assert_eq!(a.iter().map(|x| (x.length, x.reads)).collect::<Vec<_>>(), vec![(20, 15), (25, 13), (40, 4)]);
        assert_eq!(a.iter().map(|x| x.reads).sum::<u32>(), reads.len() as u32);
        // Untagged reads of a haplotype's length join it.
        reads.truncate(28);
        reads.extend([r(20, 0), r(26, 0)]);
        let a = kept(&reads);
        assert_eq!(a.iter().map(|x| (x.length, x.reads, x.by_haplotype)).collect::<Vec<_>>(), vec![(20, 16, [1, 0, 15]), (25, 14, [1, 13, 0])]);
    }

    #[test]
    fn one_untagged_outlier_is_not_an_allele() {
        // As chr20:31166646 in HG002: HP 2 at 51, HP 1 at 71 (to 81), and one untagged read at
        // 86 that fits neither: counted apart, not a third allele.
        let mut reads: Vec<RepeatRead> = (0..20).map(|_| r(51, 2)).chain((0..28).map(|i| r(66 + i % 16, 1))).collect();
        reads.extend([r(51, 0), r(71, 0), r(86, 0)]);
        let (a, unassigned) = alleles(&reads);
        assert_eq!(a.iter().map(|x| (x.length, x.reads)).collect::<Vec<_>>(), vec![(51, 21), (72, 29)]);
        assert_eq!(unassigned, 1);
        assert_eq!(a.iter().map(|x| x.reads).sum::<u32>() + unassigned, reads.len() as u32);
        // Two such reads are an allele, as two unphased reads past a gap are.
        reads.push(r(87, 0));
        let (a, unassigned) = alleles(&reads);
        assert_eq!((a.len(), a[2].length, a[2].reads, unassigned), (3, 87, 2, 0));
    }
}
