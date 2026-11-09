use ggv_model::Region;
use serde::Serialize;

use crate::{AlignedRead, Block, Indel};

/// One base of coverage, with allele counts for colouring variant columns.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize)]
pub struct CoverageColumn {
    /// Reads with a base or a deletion here (ref skips excluded), the same as
    /// `ggv query pileup`'s `depth` at default settings.
    pub depth: u32,
    /// A, C, G, T, N. Only filled in where the reference base is known (A/C/G/T).
    pub counts: [u32; 5],
    pub deletions: u32,
    pub ref_base: Option<u8>,
    /// Fraction of `depth` that is not the reference base; 0 where the reference is unknown.
    pub alt_fraction: f64,
}

impl CoverageColumn {
    /// Share of the reads here showing a base other than the reference (deletions not
    /// counted): what marks a column as a variant when drawing. `alt_fraction` counts deletions
    /// too, so under a heterozygous deletion every base would light up.
    pub fn mismatch_fraction(&self) -> f64 {
        let Some(r) = self.ref_base else { return 0.0 };
        let ref_k = [b'A', b'C', b'G', b'T'].iter().position(|&b| b == r.to_ascii_uppercase());
        let other: u32 = self.counts[..4].iter().enumerate().filter(|(k, _)| Some(*k) != ref_k).map(|(_, n)| *n).sum();
        if self.depth == 0 { 0.0 } else { f64::from(other) / f64::from(self.depth) }
    }
}

/// Coverage over a span of bases that shares one pixel column (or less).
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize)]
pub struct CoverageBin {
    pub mean: f64,
    pub max: u32,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub enum Coverage {
    /// One column per base of the region; computed from the decoded reads.
    PerBase(Vec<CoverageColumn>),
    /// Fixed-width bins across the region, for views too wide to load reads.
    Binned(Vec<CoverageBin>),
    /// Bins of depth estimated from the index (views too wide to count reads): read density
    /// from the index, scaled to depth by counting a few windows exactly. `max` equals `mean`.
    Estimated(Vec<CoverageBin>),
    /// The region is too wide even for binned coverage, and the index gave no estimate.
    Unavailable,
}

impl Coverage {
    pub fn max_depth(&self) -> u32 {
        match self {
            Self::PerBase(c) => c.iter().map(|c| c.depth).max().unwrap_or(0),
            Self::Binned(b) | Self::Estimated(b) => b.iter().map(|b| b.max).max().unwrap_or(0),
            Self::Unavailable => 0,
        }
    }
}

/// Adds +1 over every reference span of `items` that falls inside `region`.
fn add_spans(diff: &mut [i32], region: &Region, spans: impl Iterator<Item = (u64, u64)>) {
    for (s, e) in spans {
        let (s, e) = (s.max(region.start), e.min(region.end));
        if s < e {
            diff[(s - region.start) as usize] += 1;
            diff[(e - region.start) as usize] -= 1;
        }
    }
}

fn block_spans(blocks: &[Block]) -> impl Iterator<Item = (u64, u64)> + '_ {
    blocks.iter().map(|b| (b.ref_start, b.ref_start + u64::from(b.len)))
}

fn indel_spans(indels: &[Indel]) -> impl Iterator<Item = (u64, u64)> + '_ {
    indels.iter().map(|d| (d.pos, d.pos + u64::from(d.len)))
}

fn prefix_sum(diff: &[i32]) -> impl Iterator<Item = u32> + '_ {
    diff.iter().scan(0i32, |acc, d| {
        *acc += d;
        Some(*acc as u32)
    })
}

/// Per-base coverage from already-decoded reads: one pass over the BAM serves both the read
/// track and the coverage track, as in IGV. Bases that are not mismatches are counted as the
/// reference base, which is what a pileup reports when the reads were decoded against the same
/// reference.
pub fn coverage_from_reads<'a>(
    reads: impl IntoIterator<Item = &'a AlignedRead>,
    region: &Region,
    reference: Option<&[u8]>,
) -> Vec<CoverageColumn> {
    let len = region.len() as usize;
    let mut aligned = vec![0i32; len + 1];
    let mut deleted = vec![0i32; len + 1];
    let mut mismatch_counts = vec![[0u32; 5]; len];
    // Bases below the minimum base quality: aligned, but not observed (as in a pileup).
    let mut low = vec![0u32; len];

    for r in reads {
        add_spans(&mut aligned, region, block_spans(&r.blocks));
        add_spans(&mut deleted, region, indel_spans(&r.deletions));
        for m in &r.mismatches {
            if region.contains(m.pos) {
                mismatch_counts[(m.pos - region.start) as usize][base_index(m.base)] += 1;
            }
        }
        for &p in &r.low_quality {
            if region.contains(p) {
                low[(p - region.start) as usize] += 1;
            }
        }
    }

    prefix_sum(&aligned)
        .zip(prefix_sum(&deleted))
        .take(len)
        .enumerate()
        .map(|(i, (aligned, deletions))| {
            let aligned = aligned.saturating_sub(low[i]);
            let ref_base = reference.and_then(|r| r.get(i).copied());
            let depth = aligned + deletions;
            let mut col = CoverageColumn { depth, deletions, ref_base, ..Default::default() };
            if let Some(b @ (b'A' | b'C' | b'G' | b'T')) = ref_base {
                let mm = mismatch_counts[i];
                let mismatched: u32 = mm.iter().sum();
                col.counts = mm;
                col.counts[base_index(b)] += aligned.saturating_sub(mismatched);
                if depth > 0 {
                    let ref_count = col.counts[base_index(b)];
                    col.alt_fraction = round3(f64::from(depth - ref_count) / f64::from(depth));
                }
            }
            col
        })
        .collect()
}

/// Accumulates depth over a wide region without keeping reads, for binned coverage.
pub(crate) struct DepthAccumulator {
    region: Region,
    diff: Vec<i32>,
}

impl DepthAccumulator {
    pub(crate) fn new(region: &Region) -> Self {
        Self { region: region.clone(), diff: vec![0; region.len() as usize + 1] }
    }

    pub(crate) fn add(&mut self, read: &AlignedRead) {
        add_spans(&mut self.diff, &self.region, block_spans(&read.blocks));
        add_spans(&mut self.diff, &self.region, indel_spans(&read.deletions));
    }

    /// Adds a read's depth straight from its encoded CIGAR, with no decoding or allocation:
    /// aligned blocks (`M`, `=`, `X`) and deletions count, as in [`Self::add`].
    pub(crate) fn add_cigar(&mut self, pos: u64, cigar: impl Iterator<Item = u32>) {
        let mut ref_pos = pos;
        for c in cigar {
            let (op, len) = (c & 0xf, u64::from(c >> 4));
            match op {
                0 | 2 | 7 | 8 => {
                    add_spans(&mut self.diff, &self.region, std::iter::once((ref_pos, ref_pos + len)));
                    ref_pos += len;
                }
                3 => ref_pos += len,
                _ => {}
            }
        }
    }

    /// Adds `part`, an accumulator over a piece of this one's region.
    pub(crate) fn absorb(&mut self, part: &Self) {
        let offset = (part.region.start - self.region.start) as usize;
        for (i, d) in part.diff.iter().enumerate() {
            if let Some(slot) = self.diff.get_mut(offset + i) {
                *slot += d;
            }
        }
    }

    pub(crate) fn into_bins(self, bins: usize) -> Vec<CoverageBin> {
        let len = self.region.len() as usize;
        let bins = bins.clamp(1, len.max(1));
        let mut out = vec![CoverageBin::default(); bins];
        let mut sums = vec![0u64; bins];
        let mut sizes = vec![0u64; bins];
        for (i, d) in prefix_sum(&self.diff).take(len).enumerate() {
            let b = i * bins / len;
            sums[b] += u64::from(d);
            sizes[b] += 1;
            out[b].max = out[b].max.max(d);
        }
        for ((bin, sum), size) in out.iter_mut().zip(sums).zip(sizes) {
            bin.mean = if size == 0 { 0.0 } else { sum as f64 / size as f64 };
        }
        out
    }
}

fn base_index(b: u8) -> usize {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

fn round3(x: f64) -> f64 {
    (x * 1000.0).round() / 1000.0
}

#[cfg(test)]
mod mismatch_tests {
    use super::*;

    #[test]
    fn deletions_are_not_mismatches() {
        let del = CoverageColumn { depth: 10, counts: [5, 0, 0, 0, 0], deletions: 5, ref_base: Some(b'A'), alt_fraction: 0.5 };
        assert_eq!(del.mismatch_fraction(), 0.0);
        let snv = CoverageColumn { depth: 10, counts: [6, 0, 4, 0, 0], deletions: 0, ref_base: Some(b'a'), alt_fraction: 0.4 };
        assert!((snv.mismatch_fraction() - 0.4).abs() < 1e-9);
        assert_eq!(CoverageColumn { depth: 3, ref_base: None, ..Default::default() }.mismatch_fraction(), 0.0);
    }
}
