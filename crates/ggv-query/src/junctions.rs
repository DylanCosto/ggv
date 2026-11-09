//! Splice junctions from RNA-seq reads: every `N` in a read's CIGAR is an intron the read
//! skips. Counted straight from the records (nothing is decoded beyond the CIGAR and two tags),
//! with the same read filters as every other count, so a gene's junctions take a single pass.
//!
//! Strand and annotation are left to the caller, which has the reference (for the splice-site
//! motif) and the genes; [`ggv_model::splice_motif`] names a junction's motif and its strand.

use std::collections::BTreeMap;

use ggv_io::{BamFile, pasta::pileup::Read};
use ggv_model::Region;
use serde::Serialize;

use crate::{PileupOptions, ReadCounts, Result};

#[derive(Clone, Debug)]
pub struct JunctionOptions {
    pub reads: PileupOptions,
    /// A read counts for a junction only if it aligns at least this many bases on both sides
    /// of it (its anchor). Default 1.
    pub min_anchor: u32,
}

impl Default for JunctionOptions {
    fn default() -> Self {
        Self { reads: PileupOptions::default(), min_anchor: 1 }
    }
}

/// One intron and the reads that skip it. Positions are 0-based, half-open over the intron.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize)]
pub struct Junction {
    pub start: u64,
    pub end: u64,
    pub reads: u32,
    /// Reads with `NH` 1 (or no `NH`): mapped to one place.
    pub unique: u32,
    /// Reads on the forward and reverse strand (of the read, not the transcript).
    pub forward: u32,
    pub reverse: u32,
    /// Reads whose aligner gave the transcript's strand (`XS:A:+` / `XS:A:-`).
    pub xs_plus: u32,
    pub xs_minus: u32,
    /// The longest anchor among its reads (the shorter of the two flanking aligned blocks).
    pub max_anchor: u32,
}

#[derive(Clone, Debug, Serialize)]
pub struct Junctions {
    pub reads: ReadCounts,
    /// Reads (that passed the filters) with at least one junction.
    pub spliced_reads: u64,
    /// Junctions whose intron overlaps the region, in position order.
    pub junctions: Vec<Junction>,
}

const CIGAR_M: u32 = 0;
const CIGAR_I: u32 = 1;
const CIGAR_D: u32 = 2;
const CIGAR_N: u32 = 3;
const CIGAR_EQ: u32 = 7;
const CIGAR_X: u32 = 8;

/// The introns `read` skips with their anchors: `(start, end, anchor)`, 0-based half-open.
fn introns(read: &Read) -> Vec<(u64, u64, u32)> {
    let mut pos = read.position().max(0) as u64;
    // Aligned bases in each exon segment, and the introns between them.
    let mut segments = vec![0u32];
    let mut gaps = Vec::new();
    for &c in read.cigar() {
        let (op, len) = (c & 0xf, c >> 4);
        match op {
            CIGAR_M | CIGAR_EQ | CIGAR_X => {
                *segments.last_mut().expect("non-empty") += len;
                pos += u64::from(len);
            }
            CIGAR_D => pos += u64::from(len),
            CIGAR_N => {
                gaps.push((pos, pos + u64::from(len)));
                pos += u64::from(len);
                segments.push(0);
            }
            CIGAR_I => {}
            _ => {} // soft and hard clips, padding
        }
    }
    gaps.into_iter().enumerate().map(|(i, (s, e))| (s, e, segments[i].min(segments[i + 1]))).collect()
}

/// Counts the junctions of reads overlapping `region`: in one pass, or for a long region in
/// chunks on several threads (a read crossing a chunk boundary counts in the chunk it starts).
pub fn junctions(bam: &mut BamFile, region: &Region, opts: &JunctionOptions) -> Result<Junctions> {
    let parts = crate::chunks::run(bam, region, |b, chunk, count_from| junctions_pass(b, chunk, region, opts, count_from))?;
    let mut out = Junctions { reads: ReadCounts::default(), spliced_reads: 0, junctions: Vec::new() };
    let mut all: BTreeMap<(u64, u64), Junction> = BTreeMap::new();
    for p in parts {
        out.reads.merge(&p.reads);
        out.spliced_reads += p.spliced_reads;
        for j in p.junctions {
            let t = all.entry((j.start, j.end)).or_insert_with(|| Junction { start: j.start, end: j.end, ..Default::default() });
            t.reads += j.reads;
            t.unique += j.unique;
            t.forward += j.forward;
            t.reverse += j.reverse;
            t.xs_plus += j.xs_plus;
            t.xs_minus += j.xs_minus;
            t.max_anchor = t.max_anchor.max(j.max_anchor);
        }
    }
    out.junctions = all.into_values().collect();
    Ok(out)
}

/// One pass of [`junctions`] over `chunk`, keeping junctions that overlap `region` and counting
/// only reads that start at or after `count_from`.
fn junctions_pass(bam: &mut BamFile, chunk: &Region, region: &Region, opts: &JunctionOptions, count_from: u64) -> Result<Junctions> {
    let mut counts = ReadCounts::default();
    let mut spliced_reads = 0u64;
    let mut found: BTreeMap<(u64, u64), Junction> = BTreeMap::new();
    bam.for_each_record(chunk, |record| {
        let read = Read::new(record.as_ref().to_vec())?;
        if (read.position().max(0) as u64) < count_from {
            return Ok(());
        }
        if !counts.admit(&read, &opts.reads, &chunk.contig) {
            return Ok(());
        }
        let introns = introns(&read);
        if introns.is_empty() {
            return Ok(());
        }
        spliced_reads += 1;
        let unique = read.integer_tag(*b"NH").is_none_or(|n| n <= 1);
        let xs = read.tag(*b"XS").and_then(|(ty, v)| (ty == b'A').then(|| v.first().copied()).flatten());
        for (start, end, anchor) in introns {
            if anchor < opts.min_anchor || start >= region.end || end <= region.start {
                continue;
            }
            let j = found.entry((start, end)).or_insert_with(|| Junction { start, end, ..Default::default() });
            j.reads += 1;
            j.unique += u32::from(unique);
            if read.is_reverse_complemented() {
                j.reverse += 1;
            } else {
                j.forward += 1;
            }
            match xs {
                Some(b'+') => j.xs_plus += 1,
                Some(b'-') => j.xs_minus += 1,
                _ => {}
            }
            j.max_anchor = j.max_anchor.max(anchor);
        }
        Ok(())
    })?;
    Ok(Junctions { reads: counts, spliced_reads, junctions: found.into_values().collect() })
}

#[cfg(test)]
mod tests {
    use ggv_model::splice_motif as motif;

    #[test]
    fn motifs_and_strands() {
        assert_eq!(motif(*b"GT", *b"AG"), ("GT/AG", '+'));
        assert_eq!(motif(*b"ct", *b"ac"), ("GT/AG", '-'));
        assert_eq!(motif(*b"GC", *b"AG"), ("GC/AG", '+'));
        assert_eq!(motif(*b"AT", *b"AC"), ("AT/AC", '+'));
        assert_eq!(motif(*b"GT", *b"AT"), ("AT/AC", '-'));
        assert_eq!(motif(*b"AA", *b"GG"), ("non-canonical", '.'));
    }
}
