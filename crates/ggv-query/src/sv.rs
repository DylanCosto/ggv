//! Structural variants against the reads: how many reads carry a given deletion, insertion or
//! rearrangement, how many show the reference across it, and (for finding what a call set
//! missed) clusters of large events in the reads themselves.
//!
//! Long reads carry an SV in their CIGAR (a long `D` or `I`, often split into a few nearby
//! pieces by the aligner) or as a split alignment (`SA`) that jumps across it; short reads show
//! it mostly as soft clips and split reads. A read supports a deletion or insertion when its
//! events near the breakpoint add up to the SV's size (within 30%), or a split alignment bridges
//! it; it supports the reference when it spans the site without one.

use ggv_io::{BamFile, pasta::pileup::Read};
use ggv_model::Region;
use serde::Serialize;

use crate::{PileupOptions, ReadCounts, Result};

const CIGAR_M: u32 = 0;
const CIGAR_I: u32 = 1;
const CIGAR_D: u32 = 2;
const CIGAR_N: u32 = 3;
const CIGAR_S: u32 = 4;
const CIGAR_EQ: u32 = 7;
const CIGAR_X: u32 = 8;
/// CIGAR events shorter than this are ignored (alignment noise, small indels) when looking for
/// SVs; testing a site, so are those under half its size (see [`min_piece`]).
const MIN_PIECE: u32 = 20;
/// Soft clips this long count as evidence of a breakpoint.
const MIN_CLIP: u32 = 50;
/// A read must align this far past each side of a site to show it.
const FLANK: u64 = 50;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum SvKind {
    Deletion,
    Insertion,
    /// A tandem duplication of `start..end`: long reads carry the extra copy as an insertion of
    /// its size in or beside it; split reads jump from its end back to its start.
    Duplication,
    /// Inversions, breakends: two breakpoints joined by split reads.
    Rearrangement,
}

/// An SV to test: 0-based breakpoints (`end == start` for an insertion) and its length.
#[derive(Clone, Copy, Debug)]
pub struct SvSite {
    pub kind: SvKind,
    pub start: u64,
    pub end: u64,
    pub len: u64,
}

impl SvSite {
    /// How far from a breakpoint an event may sit and still count: aligners place an SV in a
    /// repeat anywhere along it.
    pub fn slack(&self) -> u64 {
        (100 + self.len / 2).min(2000)
    }
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct SvSupport {
    pub supporting: u32,
    pub reference: u32,
    /// Reads with a different allele here: an event of the same kind, at least half the
    /// SV's size, that doesn't match it (the other haplotype of a multi-allelic site).
    pub other_allele: u32,
    /// Reads soft-clipped at a breakpoint that neither support nor span it.
    pub clipped: u32,
    /// Of the supporting reads: on the forward and reverse strand, and by HP tag (none, 1, 2).
    pub forward: u32,
    pub reverse: u32,
    pub support_by_haplotype: [u32; 3],
    pub reference_by_haplotype: [u32; 3],
    /// Supporting reads that show it as one CIGAR event (or a few nearby) vs a split alignment.
    pub in_cigar: u32,
    pub split: u32,
    /// The median size of the event in the supporting reads that carry it in their CIGAR.
    pub median_event_len: Option<u64>,
    pub reads: ReadCounts,
}

impl SvSupport {
    pub fn fraction(&self) -> Option<f64> {
        let n = self.supporting + self.reference + self.other_allele;
        (n > 0).then(|| (f64::from(self.supporting) / f64::from(n) * 1000.0).round() / 1000.0)
    }

    /// 0/0, 0/1 or 1/1 from the fraction, or ./. with fewer than 5 informative reads.
    pub fn genotype(&self) -> &'static str {
        match self.fraction() {
            _ if self.supporting + self.reference + self.other_allele < 5 => "./.",
            Some(f) if f < 0.2 => "0/0",
            Some(f) if f > 0.8 => "1/1",
            _ => "0/1",
        }
    }
}

/// One read's alignment, as far as SVs go.
struct Shape {
    start: u64,
    end: u64,
    /// `(ref position, length)` of deletions and insertions of at least the minimum piece.
    deletions: Vec<(u64, u32)>,
    insertions: Vec<(u64, u32)>,
    clip_left: u32,
    clip_right: u32,
    /// Other alignments of the read on the same contig (from `SA`): start and end.
    others: Vec<(u64, u64)>,
    reverse: bool,
    haplotype: usize,
}

fn ref_span(cigar: impl Iterator<Item = (u32, u32)>) -> u64 {
    cigar.filter(|(op, _)| matches!(*op, CIGAR_M | CIGAR_D | CIGAR_N | CIGAR_EQ | CIGAR_X)).map(|(_, n)| u64::from(n)).sum()
}

/// The smallest CIGAR event that counts at `site`: `MIN_PIECE`, or half the SV's size when that
/// is smaller, so a 15 bp deletion in a read at a 20 bp SV is that read's allele (supporting it
/// or another allele), not the reference. A rearrangement has no size: `MIN_PIECE`.
fn min_piece(site: &SvSite) -> u32 {
    if site.kind == SvKind::Rearrangement {
        return MIN_PIECE;
    }
    MIN_PIECE.min((site.len / 2).max(1) as u32)
}

fn shape(read: &Read, contig: &str, min_piece: u32) -> Shape {
    let start = read.position().max(0) as u64;
    let mut pos = start;
    let (mut deletions, mut insertions) = (Vec::new(), Vec::new());
    let ops: Vec<(u32, u32)> = read.cigar().iter().map(|c| (c & 0xf, c >> 4)).collect();
    for &(op, len) in &ops {
        match op {
            CIGAR_M | CIGAR_EQ | CIGAR_X | CIGAR_N => pos += u64::from(len),
            CIGAR_D => {
                if len >= min_piece {
                    deletions.push((pos, len));
                }
                pos += u64::from(len);
            }
            CIGAR_I if len >= min_piece => insertions.push((pos, len)),
            _ => {}
        }
    }
    let clip = |o: Option<&(u32, u32)>| o.filter(|(op, _)| *op == CIGAR_S).map_or(0, |(_, n)| *n);
    // SA: contig,pos,strand,CIGAR,mapq,NM; ...
    let others = read
        .string_tag(*b"SA")
        .map(|sa| {
            String::from_utf8_lossy(sa)
                .split(';')
                .filter_map(|part| {
                    let f: Vec<&str> = part.split(',').collect();
                    if f.len() < 4 || f[0] != contig {
                        return None;
                    }
                    let s = f[1].parse::<u64>().ok()?.checked_sub(1)?;
                    let cig = parse_cigar(f[3])?;
                    Some((s, s + ref_span(cig.into_iter())))
                })
                .collect()
        })
        .unwrap_or_default();
    Shape {
        start,
        end: pos,
        deletions,
        insertions,
        clip_left: clip(ops.first()),
        clip_right: clip(ops.last()),
        others,
        reverse: read.is_reverse_complemented(),
        haplotype: match crate::haplotype(read) {
            Some(1) => 1,
            Some(2) => 2,
            _ => 0,
        },
    }
}

fn parse_cigar(text: &str) -> Option<Vec<(u32, u32)>> {
    let mut out = Vec::new();
    let mut n = 0u32;
    for c in text.bytes() {
        if c.is_ascii_digit() {
            n = n.checked_mul(10)?.checked_add(u32::from(c - b'0'))?;
        } else {
            let op = match c {
                b'M' => CIGAR_M,
                b'I' => CIGAR_I,
                b'D' => CIGAR_D,
                b'N' => CIGAR_N,
                b'S' => CIGAR_S,
                b'=' => CIGAR_EQ,
                b'X' => CIGAR_X,
                _ => 5, // H, P: no reference
            };
            out.push((op, n));
            n = 0;
        }
    }
    Some(out)
}

/// How far apart (on the reference) pieces of one event split by the aligner may be.
const PIECE_GAP: u64 = 100;

/// The size of `site`'s event in this read, if its CIGAR carries one near the site: a deletion
/// (or insertion; for a tandem duplication, an insertion in or beside it) within the slack, or
/// a run of pieces each within 100 bp of the last, whose sum is 70–130% of the SV's length (the
/// closest such). Only a run is summed: a small indel elsewhere near the site isn't part of it.
fn cigar_event(s: &Shape, site: &SvSite) -> Option<u64> {
    let w = site.slack();
    let near = |p: u64, a: u64, b: u64| p + w >= a && p <= b + w;
    // (start, end on the reference, length), in reference order.
    let pieces: Vec<(u64, u64, u64)> = match site.kind {
        SvKind::Deletion => s
            .deletions
            .iter()
            .map(|&(p, n)| (p, p + u64::from(n), u64::from(n)))
            .filter(|&(p, e, _)| near(p, site.start, site.end) && near(e, site.start, site.end))
            .collect(),
        SvKind::Insertion | SvKind::Duplication => {
            s.insertions.iter().filter(|(p, _)| near(*p, site.start, site.end)).map(|&(p, n)| (p, p, u64::from(n))).collect()
        }
        SvKind::Rearrangement => return None,
    };
    let len = site.len.max(1);
    let fits = |t: u64| t * 10 >= len * 7 && t * 10 <= len * 13;
    let mut best: Option<u64> = None;
    for i in 0..pieces.len() {
        let mut total = 0;
        for j in i..pieces.len() {
            if j > i && pieces[j].0 > pieces[j - 1].1 + PIECE_GAP {
                break;
            }
            total += pieces[j].2;
            if fits(total) && best.is_none_or(|b| total.abs_diff(len) < b.abs_diff(len)) {
                best = Some(total);
            }
        }
    }
    best
}

/// Whether one alignment of the read ends at one breakpoint and another starts at the other:
/// across a deletion, or (a duplication) from its end back to its start.
fn bridges(s: &Shape, site: &SvSite) -> bool {
    let w = site.slack();
    let close = |a: u64, b: u64| a.abs_diff(b) <= w;
    let parts: Vec<(u64, u64)> = std::iter::once((s.start, s.end)).chain(s.others.iter().copied()).collect();
    let skips = || parts.iter().any(|a| parts.iter().any(|b| a != b && close(a.1, site.start) && close(b.0, site.end)));
    let repeats = || parts.iter().any(|a| parts.iter().any(|b| a != b && close(a.0, site.start) && close(b.1, site.end)));
    match site.kind {
        SvKind::Insertion => false,
        SvKind::Deletion => skips(),
        SvKind::Duplication => repeats(),
        SvKind::Rearrangement => skips() || repeats(),
    }
}

/// Read support for `site`, from reads overlapping its breakpoints.
pub fn sv_support(bam: &mut BamFile, contig: &str, site: &SvSite, opts: &PileupOptions) -> Result<SvSupport> {
    let w = site.slack();
    let region = Region { contig: contig.to_owned(), start: site.start.saturating_sub(w), end: site.end + w + 1 };
    let mut out = SvSupport::default();
    let mut lengths = Vec::new();
    let mut counts = ReadCounts::default();
    // A rearrangement's breakpoints may be far apart: read around each, not between.
    let regions: Vec<Region> = if matches!(site.kind, SvKind::Rearrangement | SvKind::Duplication) && site.end - site.start > 4 * w {
        vec![
            Region { contig: contig.to_owned(), start: site.start.saturating_sub(w), end: site.start + w },
            Region { contig: contig.to_owned(), start: site.end.saturating_sub(w), end: site.end + w },
        ]
    } else {
        vec![region]
    };
    // Each read once: its supplementary alignments (found at the other breakpoint, or in the
    // same window) are already in its `SA`.
    let mut seen = std::collections::HashSet::new();
    for r in &regions {
        bam.for_each_record(r, |record| {
            let read = Read::new(record.as_ref().to_vec())?;
            if !counts.admit(&read, opts, contig) || !seen.insert((read.name().to_vec(), read.flags() & 0xc0)) {
                return Ok(());
            }
            let s = shape(&read, contig, min_piece(site));
            let event = cigar_event(&s, site);
            if event.is_some() || bridges(&s, site) {
                out.supporting += 1;
                out.support_by_haplotype[s.haplotype] += 1;
                if s.reverse {
                    out.reverse += 1;
                } else {
                    out.forward += 1;
                }
                match event {
                    Some(n) => {
                        out.in_cigar += 1;
                        lengths.push(n);
                    }
                    None => out.split += 1,
                }
                return Ok(());
            }
            let spans = |a: u64, b: u64| s.start + FLANK <= a && s.end >= b + FLANK;
            // A large CIGAR event over or at breakpoint `b`: the read doesn't show the
            // reference there (a deletion called as a breakend pair, say).
            let jumps = |b: u64| {
                s.deletions.iter().any(|&(p, n)| p <= b + FLANK && p + u64::from(n) + FLANK >= b) || s.insertions.iter().any(|(p, _)| p.abs_diff(b) <= w)
            };
            // An event of the same kind at least half the SV's size: another allele (smaller
            // ones are alignment noise, common in long reads).
            let other = match site.kind {
                SvKind::Deletion => s.deletions.iter().any(|(p, n)| *p + w >= site.start && *p <= site.end + w && u64::from(*n) * 2 >= site.len),
                SvKind::Insertion => s.insertions.iter().any(|(p, n)| p.abs_diff(site.start) <= w && u64::from(*n) * 2 >= site.len),
                SvKind::Duplication => s.insertions.iter().any(|(p, n)| *p + w >= site.start && *p <= site.end + w && u64::from(*n) * 2 >= site.len),
                SvKind::Rearrangement => jumps(site.start) || jumps(site.end),
            };
            if other {
                out.other_allele += 1;
                return Ok(());
            }
            let clean = match site.kind {
                SvKind::Deletion => spans(site.start, site.end),
                SvKind::Insertion => spans(site.start, site.start),
                // Each copy aligns cleanly over one junction: only a read across the whole
                // span without the extra copy shows the reference.
                SvKind::Duplication => spans(site.start, site.end),
                SvKind::Rearrangement => spans(site.start, site.start) || spans(site.end, site.end),
            };
            let clipped_at = |p: u64| (s.clip_left >= MIN_CLIP && s.start.abs_diff(p) <= w) || (s.clip_right >= MIN_CLIP && s.end.abs_diff(p) <= w);
            if clean {
                out.reference += 1;
                out.reference_by_haplotype[s.haplotype] += 1;
            } else if clipped_at(site.start) || clipped_at(site.end) {
                out.clipped += 1;
            }
            Ok(())
        })?;
    }
    lengths.sort_unstable();
    out.median_event_len = lengths.get(lengths.len() / 2).copied();
    out.reads = counts;
    Ok(out)
}

/// A cluster of large events in the reads: a candidate SV.
#[derive(Clone, Debug, Serialize)]
pub struct SvEvent {
    pub kind: SvKind,
    /// 0-based position of the breakpoint (the deletion's start, the insertion point).
    pub pos: u64,
    pub len: u64,
    pub reads: u32,
    pub by_haplotype: [u32; 3],
}

/// Clusters the deletions and insertions of at least `min_len` that reads in `region` carry
/// (summing pieces within 100 bp in a read), keeping clusters of `min_reads` reads or more.
pub fn sv_events(bam: &mut BamFile, region: &Region, min_len: u32, min_reads: u32, opts: &PileupOptions) -> Result<(Vec<SvEvent>, ReadCounts)> {
    let parts = crate::chunks::run(bam, region, |b, chunk, count_from| events_pass(b, chunk, region, min_len, opts, count_from))?;
    let mut counts = ReadCounts::default();
    let mut raw = Vec::new();
    for (events, c) in parts {
        raw.extend(events);
        counts.merge(&c);
    }
    Ok((cluster(raw, min_reads), counts))
}

/// One read's large event: kind, position, length, haplotype.
type ReadEvent = (SvKind, u64, u64, usize);

/// Every read's large deletions and insertions in `chunk` that lie in `region`, from reads
/// starting at or after `count_from`.
fn events_pass(bam: &mut BamFile, chunk: &Region, region: &Region, min_len: u32, opts: &PileupOptions, count_from: u64) -> Result<(Vec<ReadEvent>, ReadCounts)> {
    let mut counts = ReadCounts::default();
    let mut raw: Vec<ReadEvent> = Vec::new();
    bam.for_each_record(chunk, |record| {
        let read = Read::new(record.as_ref().to_vec())?;
        if (read.position().max(0) as u64) < count_from || !counts.admit(&read, opts, &chunk.contig) {
            return Ok(());
        }
        let s = shape(&read, &region.contig, MIN_PIECE);
        for (kind, events) in [(SvKind::Deletion, &s.deletions), (SvKind::Insertion, &s.insertions)] {
            // Merge pieces of one event split by the aligner: within 100 bp of where the last
            // piece ends on the reference (an insertion's end is where it is).
            let mut merged: Vec<(u64, u64, u64)> = Vec::new();
            for &(p, n) in events {
                let end = if kind == SvKind::Deletion { p + u64::from(n) } else { p };
                match merged.last_mut() {
                    Some((_, mn, me)) if p <= *me + 100 => {
                        *mn += u64::from(n);
                        *me = end;
                    }
                    _ => merged.push((p, u64::from(n), end)),
                }
            }
            for (p, n, _) in merged {
                if n >= u64::from(min_len) && p >= region.start && p < region.end {
                    raw.push((kind, p, n, s.haplotype));
                }
            }
        }
        Ok(())
    })?;
    Ok((raw, counts))
}

/// Clusters read events of one kind near each other with similar sizes; keeps clusters of
/// `min_reads` reads or more.
fn cluster(mut raw: Vec<ReadEvent>, min_reads: u32) -> Vec<SvEvent> {
    // Every field, so the clusters don't depend on the order the reads came in.
    raw.sort_by_key(|e| (e.0 as u8, e.1, e.2, e.3));
    // Per cluster: its kind and its reads' (position, length, haplotype).
    type Cluster = (SvKind, Vec<(u64, u64, usize)>);
    let mut clusters: Vec<Cluster> = Vec::new();
    for (kind, p, n, hp) in raw {
        let fits = |c: &Cluster| {
            let (lp, ln, _) = *c.1.last().expect("non-empty");
            c.0 == kind && p.abs_diff(lp) <= (100 + ln / 2).min(2000) && n.abs_diff(ln) * 10 <= ln.max(n) * 3
        };
        match clusters.iter_mut().rev().take(8).find(|c| fits(c)) {
            Some(c) => c.1.push((p, n, hp)),
            None => clusters.push((kind, vec![(p, n, hp)])),
        }
    }
    let median = |mut v: Vec<u64>| {
        v.sort_unstable();
        v[v.len() / 2]
    };
    let mut out: Vec<SvEvent> = clusters
        .into_iter()
        .filter(|c| c.1.len() as u32 >= min_reads)
        .map(|(kind, ev)| {
            let mut by_haplotype = [0u32; 3];
            for e in &ev {
                by_haplotype[e.2] += 1;
            }
            SvEvent { kind, pos: median(ev.iter().map(|e| e.0).collect()), len: median(ev.iter().map(|e| e.1).collect()), reads: ev.len() as u32, by_haplotype }
        })
        .collect();
    out.sort_by_key(|e| e.pos);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cigars_and_slack() {
        assert_eq!(parse_cigar("10S100M2000D50M5H"), Some(vec![(CIGAR_S, 10), (CIGAR_M, 100), (CIGAR_D, 2000), (CIGAR_M, 50), (5, 5)]));
        assert_eq!(ref_span(parse_cigar("10S100M2000D50M").unwrap().into_iter()), 2150);
        let site = SvSite { kind: SvKind::Deletion, start: 1000, end: 1300, len: 300 };
        assert_eq!(site.slack(), 250);
        let shape = |dels: Vec<(u64, u32)>| Shape { start: 0, end: 5000, deletions: dels, insertions: vec![], clip_left: 0, clip_right: 0, others: vec![], reverse: false, haplotype: 0 };
        // One deletion of the right size, or two pieces near the site adding up to it.
        assert_eq!(cigar_event(&shape(vec![(1010, 295)]), &site), Some(295));
        assert_eq!(cigar_event(&shape(vec![(1000, 150), (1160, 160)]), &site), Some(310));
        // Too small, or too far away.
        assert_eq!(cigar_event(&shape(vec![(1000, 100)]), &site), None);
        assert_eq!(cigar_event(&shape(vec![(3000, 300)]), &site), None);
        // A small SV keeps the pieces half its size: a read's 15 bp deletion at a 20 bp one.
        let small = SvSite { kind: SvKind::Deletion, start: 1000, end: 1020, len: 20 };
        assert_eq!((min_piece(&small), min_piece(&site)), (10, MIN_PIECE));
        assert_eq!(cigar_event(&shape(vec![(1000, 15)]), &small), Some(15));
        // The SV plus an unrelated small deletion near it: the SV's piece alone, not the sum.
        assert_eq!(cigar_event(&shape(vec![(1000, 20), (1060, 12)]), &small), Some(20));
        assert_eq!(cigar_event(&shape(vec![(950, 12), (1010, 295)]), &site), Some(295));
        // Pieces far apart on the reference aren't one event.
        assert_eq!(cigar_event(&shape(vec![(1000, 150), (1400, 150)]), &site), None);
    }

    #[test]
    fn duplications_are_insertions_of_their_size() {
        // A 586 bp tandem duplication, carried as an insertion at its start or end.
        let dup = SvSite { kind: SvKind::Duplication, start: 10_000, end: 10_586, len: 586 };
        let ins = |at: Vec<(u64, u32)>| Shape { start: 9000, end: 12_000, deletions: vec![], insertions: at, clip_left: 0, clip_right: 0, others: vec![], reverse: false, haplotype: 0 };
        assert_eq!(cigar_event(&ins(vec![(10_001, 586)]), &dup), Some(586));
        assert_eq!(cigar_event(&ins(vec![(10_590, 580)]), &dup), Some(580));
        assert_eq!(cigar_event(&ins(vec![(10_001, 100)]), &dup), None);
        // A split read from its end back to its start; a deletion's split doesn't.
        let split = |others| Shape { start: 9000, end: 10_586, others, ..ins(vec![]) };
        assert!(bridges(&split(vec![(10_000, 11_500)]), &dup));
        assert!(!bridges(&split(vec![(10_000, 11_500)]), &SvSite { kind: SvKind::Deletion, ..dup }));
        assert_eq!(min_piece(&SvSite { kind: SvKind::Rearrangement, start: 0, end: 10, len: 0 }), MIN_PIECE);
    }
}
