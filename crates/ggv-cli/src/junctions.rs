//! Splice junctions across one or more RNA-seq files, side by side: each junction's motif and
//! strand (from the reference), whether the annotation knows it, and per file its reads and its
//! share of the reads leaving its donor and entering its acceptor. A change in those shares
//! between samples is what aberrant splicing looks like. Used by `ggv query junctions` and the
//! MCP `junctions` tool.

use std::collections::{BTreeMap, HashMap};

use ggv_annot::{Annotation, Strand};
use ggv_query::{Junction, Junctions};
use serde::Serialize;

#[derive(Clone, Debug, Serialize)]
pub struct Row {
    /// 1-based, inclusive: the intron's first and last base.
    pub start: u64,
    pub end: u64,
    pub length: u64,
    pub strand: char,
    pub motif: &'static str,
    /// known, novel (known sites), novel donor, novel acceptor, novel; or unannotated when no
    /// annotation was given. Only annotation on the junction's strand counts; a junction of
    /// unknown strand (`.`) is matched against both.
    pub annotation: &'static str,
    pub genes: Vec<String>,
    /// One per file, in the order given.
    pub samples: Vec<Count>,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct Count {
    pub reads: u32,
    pub unique: u32,
    /// This junction's share of the reads of every junction from its donor (its 5' splice site
    /// on the transcript), and into its acceptor; `None` with no strand to tell them apart.
    pub donor_share: Option<f64>,
    pub acceptor_share: Option<f64>,
    pub max_anchor: u32,
}

#[derive(Clone, Debug, Serialize)]
pub struct Report {
    pub junctions: Vec<Row>,
    /// Junctions before `min_reads` and `max` were applied.
    pub total: usize,
}

fn round3(x: f64) -> f64 {
    (x * 1000.0).round() / 1000.0
}

/// Builds the table. `base(p)` is the reference base at 0-based `p` on `contig` (for motifs);
/// a junction is kept when any file has at least `min_reads` reads for it; at most `max` are
/// listed, the most-read first when there are more, then in position order.
pub fn report(
    contig: &str,
    files: &[Junctions],
    base: impl Fn(u64) -> Option<u8>,
    ann: Option<&Annotation>,
    min_reads: u32,
    max: usize,
) -> Report {
    let mut merged: BTreeMap<(u64, u64), Vec<Option<&Junction>>> = BTreeMap::new();
    for (i, f) in files.iter().enumerate() {
        for j in &f.junctions {
            merged.entry((j.start, j.end)).or_insert_with(|| vec![None; files.len()])[i] = Some(j);
        }
    }
    let total = merged.len();
    let (lo, hi) = merged.keys().fold((u64::MAX, 0), |(lo, hi), &(s, e)| (lo.min(s), hi.max(e)));
    let introns = ann.map(|a| a.introns(contig, lo.saturating_sub(1), hi + 1)).unwrap_or_default();

    let strand_of = |(s, e): (u64, u64), per: &[Option<&Junction>]| -> (&'static str, char) {
        let two = |p: u64| Some([base(p)?, base(p + 1)?]);
        let (motif, strand) = match (two(s), two(e.saturating_sub(2))) {
            (Some(a), Some(b)) => ggv_query::motif(a, b),
            _ => ("unknown", '.'),
        };
        if strand != '.' {
            return (motif, strand);
        }
        // Non-canonical: the aligner's XS votes, else the annotation's strand (when its
        // transcripts there agree on one).
        let (plus, minus) = per.iter().flatten().fold((0, 0), |(p, m), j| (p + j.xs_plus, m + j.xs_minus));
        let mut annotated = introns.iter().filter(|i| i.start == s && i.end == e && i.strand != Strand::Unknown).map(|i| if i.strand == Strand::Reverse { '-' } else { '+' });
        let first = annotated.next();
        let by_annotation = first.filter(|f| annotated.all(|o| o == *f));
        (motif, if plus > minus { '+' } else if minus > plus { '-' } else { by_annotation.unwrap_or('.') })
    };
    // Annotation on the other strand doesn't describe a junction; with no strand, either does.
    let same_strand = |junction: char, t: Strand| match junction {
        '+' => t != Strand::Reverse,
        '-' => t != Strand::Forward,
        _ => true,
    };

    let mut rows: Vec<(Row, (u64, u64))> = merged
        .iter()
        .map(|(&(s, e), per)| {
            let (motif, strand) = strand_of((s, e), per);
            let here: Vec<&ggv_annot::Intron> = introns.iter().filter(|i| same_strand(strand, i.strand)).collect();
            let starts = here.iter().any(|i| i.start == s);
            let ends = here.iter().any(|i| i.end == e);
            let annotation = match (ann.is_some(), here.iter().any(|i| i.start == s && i.end == e), starts, ends) {
                (false, ..) => "unannotated",
                (true, true, ..) => "known",
                (true, false, true, true) => "novel (known sites)",
                // The donor is the intron's start on +, its end on −.
                (true, false, true, false) => if strand == '-' { "novel donor" } else { "novel acceptor" },
                (true, false, false, true) => if strand == '-' { "novel acceptor" } else { "novel donor" },
                _ => "novel",
            };
            let mut genes: Vec<String> = here.iter().filter(|i| i.start == s || i.end == e).map(|i| i.gene.to_string()).collect();
            if genes.is_empty()
                && let Some(a) = ann
            {
                // Genes at either end (not every gene a long intron spans), on its strand.
                genes = a
                    .query(contig, s.saturating_sub(1), s + 1)
                    .chain(a.query(contig, e.saturating_sub(1), e + 1))
                    .filter(|t| same_strand(strand, t.strand))
                    .map(|t| t.gene_name.to_string())
                    .collect();
            }
            genes.sort();
            genes.dedup();
            let samples = per
                .iter()
                .map(|j| j.map_or_else(Count::default, |j| Count { reads: j.reads, unique: j.unique, max_anchor: j.max_anchor, ..Default::default() }))
                .collect();
            (Row { start: s + 1, end: e, length: e - s, strand, motif, annotation, genes, samples }, (s, e))
        })
        .collect();

    // Shares: per file, the reads of every junction leaving each donor and entering each
    // acceptor (sites keyed by position and strand).
    for f in 0..files.len() {
        let mut donors: HashMap<(u64, char), u32> = HashMap::new();
        let mut acceptors: HashMap<(u64, char), u32> = HashMap::new();
        let sites = |r: &Row, (s, e): (u64, u64)| match r.strand {
            '+' => Some(((s, '+'), (e, '+'))),
            '-' => Some(((e, '-'), (s, '-'))),
            _ => None,
        };
        for (r, k) in &rows {
            if let Some((d, a)) = sites(r, *k) {
                *donors.entry(d).or_default() += r.samples[f].reads;
                *acceptors.entry(a).or_default() += r.samples[f].reads;
            }
        }
        for (r, k) in &mut rows {
            if let Some((d, a)) = sites(r, *k) {
                let n = f64::from(r.samples[f].reads);
                let share = |total: u32| (total > 0).then(|| round3(n / f64::from(total)));
                r.samples[f].donor_share = share(donors[&d]);
                r.samples[f].acceptor_share = share(acceptors[&a]);
            }
        }
    }

    let mut rows: Vec<Row> = rows.into_iter().map(|(r, _)| r).filter(|r| r.samples.iter().any(|c| c.reads >= min_reads)).collect();
    if rows.len() > max {
        rows.sort_by_key(|r| std::cmp::Reverse(r.samples.iter().map(|c| c.reads).max().unwrap_or(0)));
        rows.truncate(max);
        rows.sort_by_key(|r| (r.start, r.end));
    }
    Report { junctions: rows, total }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn j(start: u64, end: u64, reads: u32) -> Junction {
        Junction { start, end, reads, unique: reads, forward: reads, max_anchor: 20, ..Default::default() }
    }

    #[test]
    fn shares_strands_and_annotation() {
        // Forward genome: GT at 100, AG ending at 200 and at 300 (two acceptors, one donor).
        let base = |p: u64| match p {
            100 => Some(b'G'),
            101 => Some(b'T'),
            198 | 298 => Some(b'A'),
            199 | 299 => Some(b'G'),
            _ => Some(b'C'),
        };
        let a = Junctions { reads: Default::default(), spliced_reads: 10, junctions: vec![j(100, 200, 9), j(100, 300, 1)] };
        let b = Junctions { reads: Default::default(), spliced_reads: 10, junctions: vec![j(100, 200, 5), j(100, 300, 5)] };
        let r = report("chr1", &[a, b], base, None, 1, 100);
        assert_eq!(r.total, 2);
        let (first, second) = (&r.junctions[0], &r.junctions[1]);
        assert_eq!((first.start, first.end, first.length, first.strand, first.motif), (101, 200, 100, '+', "GT/AG"));
        assert_eq!(first.annotation, "unannotated");
        // One donor used by both: 9 of 10 reads in the first file, 5 of 10 in the second.
        assert_eq!((first.samples[0].donor_share, first.samples[1].donor_share), (Some(0.9), Some(0.5)));
        assert_eq!((second.samples[0].donor_share, second.samples[0].acceptor_share), (Some(0.1), Some(1.0)));
        // min_reads applies to the best file; max keeps the most-read.
        assert_eq!(report("chr1", &[Junctions { reads: Default::default(), spliced_reads: 0, junctions: vec![j(100, 200, 9), j(100, 300, 1)] }], base, None, 2, 100).junctions.len(), 1);
    }

    fn gtf(lines: &[(&str, char, &str)]) -> Annotation {
        // (gene, strand, exon start-end, 1-based)
        let text: String = lines
            .iter()
            .map(|(g, strand, span)| {
                let (a, b) = span.split_once('-').unwrap();
                format!("chr1	X	exon	{a}	{b}	.	{strand}	.	gene_id \"{g}\"; transcript_id \"T{g}\"; gene_name \"{g}\";
")
            })
            .collect();
        Annotation::parse(text.as_bytes(), ggv_annot::Format::Gtf).unwrap()
    }

    // GT at 0-based 100, AG ending at 200: a + strand intron 101-200 (1-based).
    fn gt_ag(p: u64) -> Option<u8> {
        Some(match p {
            100 => b'G',
            101 => b'T',
            198 => b'A',
            199 => b'G',
            _ => b'C',
        })
    }

    #[test]
    fn shared_introns_keep_every_gene() {
        let one = Junctions { reads: Default::default(), spliced_reads: 5, junctions: vec![j(100, 200, 5)] };
        for order in [["GENE_A", "GENE_B"], ["GENE_B", "GENE_A"]] {
            let ann = gtf(&[(order[0], '+', "51-100"), (order[0], '+', "201-250"), (order[1], '+', "51-100"), (order[1], '+', "201-250")]);
            let r = report("chr1", std::slice::from_ref(&one), gt_ag, Some(&ann), 1, 100);
            assert_eq!(r.junctions[0].annotation, "known");
            assert_eq!(r.junctions[0].genes, ["GENE_A", "GENE_B"], "{order:?}");
        }
    }

    #[test]
    fn annotation_on_the_other_strand_does_not_count() {
        let one = Junctions { reads: Default::default(), spliced_reads: 5, junctions: vec![j(100, 200, 5)] };
        let plus = gtf(&[("G", '+', "51-100"), ("G", '+', "201-250")]);
        let minus = gtf(&[("G", '-', "51-100"), ("G", '-', "201-250")]);
        let r = report("chr1", std::slice::from_ref(&one), gt_ag, Some(&plus), 1, 100);
        assert_eq!((r.junctions[0].strand, r.junctions[0].annotation, r.junctions[0].genes.clone()), ('+', "known", vec!["G".to_owned()]));
        let r = report("chr1", std::slice::from_ref(&one), gt_ag, Some(&minus), 1, 100);
        assert_eq!((r.junctions[0].strand, r.junctions[0].annotation), ('+', "novel"));
        assert!(r.junctions[0].genes.is_empty(), "{:?}", r.junctions[0].genes);
        // No motif and no XS: the strand is unknown, and either strand's annotation matches.
        let r = report("chr1", std::slice::from_ref(&one), |_| Some(b'C'), Some(&minus), 1, 100);
        assert_eq!((r.junctions[0].strand, r.junctions[0].annotation), ('-', "known"), "takes the annotation's strand");
    }
}
