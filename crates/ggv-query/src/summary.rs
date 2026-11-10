use std::collections::BTreeMap;

use ggv_io::BamFile;
use ggv_model::Region;
use serde::Serialize;

use crate::{PileupOptions, PositionStats, ReadCounts, Result};

#[derive(Clone, Debug)]
pub struct SummaryOptions {
    pub pileup: PileupOptions,
    /// A position is "variable" when its most common non-reference allele is at least this
    /// fraction of the reads there. Default 0.2.
    pub min_alt_fraction: f64,
    /// …and its depth is at least this. Default 8.
    pub min_variable_depth: u32,
    /// The most variable sites to list (highest alt fraction first). Default 25.
    pub max_variable_sites: usize,
}

impl Default for SummaryOptions {
    fn default() -> Self {
        Self {
            pileup: PileupOptions::default(),
            min_alt_fraction: 0.2,
            min_variable_depth: 8,
            max_variable_sites: 25,
        }
    }
}

/// A one-screen answer to "what does this region look like?".
#[derive(Clone, Debug, Serialize)]
pub struct RegionSummary {
    pub region: String,
    pub length: u64,
    pub reads: ReadCounts,
    pub mean_depth: f64,
    pub median_depth: u32,
    pub min_depth: u32,
    pub max_depth: u32,
    /// Fraction of positions with depth ≥ 1, 10, 20, 30.
    pub breadth: BTreeMap<String, f64>,
    pub variable_site_count: u64,
    pub variable_sites: Vec<VariableSite>,
}

/// A variable position, or a run of adjacent positions deleted in enough reads, listed once.
#[derive(Clone, Debug, Serialize)]
pub struct VariableSite {
    /// The position (1-based), or a deletion run's first deleted base.
    pub pos: u64,
    /// The reference base at `pos`.
    #[serde(rename = "ref")]
    pub ref_base: Option<char>,
    /// Depth, alt reads and fraction at the position; for a deletion run, at the base of the
    /// run where the most reads carry the deletion.
    pub depth: u32,
    pub alt: String,
    pub alt_count: u32,
    pub alt_fwd: u32,
    pub alt_rev: u32,
    pub alt_fraction: f64,
    /// Bases the site covers: 1, or the length of a deletion run.
    #[serde(skip_serializing_if = "is_one")]
    pub length: u64,
}

/// For serde's `skip_serializing_if`, which passes a reference.
fn is_one(n: &u64) -> bool {
    *n == 1
}

/// The pileup's name for the deleted allele.
const DELETION: &str = "DEL";

impl VariableSite {
    /// The change as people read it: `C>T`, `G>DEL` for one deleted base, `DEL 15 bp` for a run.
    #[must_use]
    pub fn change(&self) -> String {
        if self.length > 1 { format!("{} {} bp", self.alt, self.length) } else { format!("{}>{}", self.ref_base.unwrap_or('.'), self.alt) }
    }

    fn is_deletion(&self) -> bool {
        self.alt == DELETION
    }

    /// The position just past the site (1-based).
    fn end(&self) -> u64 {
        self.pos + self.length
    }

    /// Adds `next`, a deletion that starts where this one ends, keeping the counts of whichever
    /// base more reads delete.
    fn extend(&mut self, next: &VariableSite) {
        debug_assert_eq!(self.end(), next.pos, "runs are extended by the adjacent site");
        let length = self.length + next.length;
        if next.alt_fraction > self.alt_fraction {
            *self = VariableSite { pos: self.pos, ref_base: self.ref_base, ..next.clone() };
        }
        self.length = length;
    }
}

/// The listed allele of a variable site: deep enough, and its share of the reads at least the
/// minimum, compared unrounded (400 of 2001 reads is below 0.2, though shown as 0.2).
fn variable_alt<'a>(s: &'a PositionStats, opts: &SummaryOptions) -> Option<&'a crate::AltAllele> {
    let alt = s.top_alt.as_ref()?;
    (s.depth >= opts.min_variable_depth && s.depth > 0 && f64::from(alt.count.total) / f64::from(s.depth) >= opts.min_alt_fraction).then_some(alt)
}

/// The variable sites at one position: the deletion there when enough reads carry it, and the
/// most common other allele when it passes too.
///
/// The deletion is kept apart from the top allele so that a long deletion stays one run: where
/// the reads without it carry a common mismatch, that base outnumbers the deletion and was the
/// only site listed, cutting a 6.8 kb heterozygous deletion into 28 rows. Both are listed there
/// now, the deletion as part of its run.
fn sites_at(s: &PositionStats, opts: &SummaryOptions) -> (Option<VariableSite>, Option<VariableSite>) {
    let site = |allele: &str, count: &crate::StrandCount| VariableSite {
        pos: s.pos,
        ref_base: s.ref_base,
        depth: s.depth,
        alt: allele.to_owned(),
        alt_count: count.total,
        alt_fwd: count.fwd,
        alt_rev: count.rev,
        alt_fraction: (f64::from(count.total) / f64::from(s.depth.max(1)) * 1000.0).round() / 1000.0,
        length: 1,
    };
    let passes = |count: u32| s.depth >= opts.min_variable_depth && s.depth > 0 && f64::from(count) / f64::from(s.depth) >= opts.min_alt_fraction;
    let deletion = passes(s.deletions.total).then(|| site(DELETION, &s.deletions));
    // The listed allele must pass on its own: a site whose non-reference reads are split
    // among several alleles, none of them common, is noise rather than a variant.
    let other = variable_alt(s, opts).filter(|alt| alt.allele != DELETION).map(|alt| site(alt.allele, &alt.count));
    (deletion, other)
}

/// What one pass over (part of) a region adds up.
#[derive(Default)]
struct Parts {
    histogram: BTreeMap<u32, u64>,
    depth_sum: u64,
    /// The most variable sites, already cut to the number asked for.
    variable: Vec<VariableSite>,
    variable_count: u64,
    /// A deletion run that starts at the pass's first position, and one that ends at its last:
    /// neither is in `variable` or counted yet, since the neighbouring pass may hold the rest.
    head: Option<VariableSite>,
    tail: Option<VariableSite>,
    /// The head run reaches the pass's last position too (a deletion over the whole chunk).
    head_reaches_end: bool,
    reads: ReadCounts,
}

impl Parts {
    fn add(&mut self, site: VariableSite) {
        self.variable_count += 1;
        self.variable.push(site);
    }

    /// Files a finished deletion run of the pass over `region`: at an edge when it may go on in
    /// the next pass, else in the list.
    fn file(&mut self, site: VariableSite, region: &Region) {
        debug_assert!(site.is_deletion(), "other sites are added as they are met");
        if site.pos == region.start + 1 {
            self.head_reaches_end = site.end() == region.end + 1;
            self.head = Some(site);
        } else if site.end() == region.end + 1 {
            self.tail = Some(site);
        } else {
            self.add(site);
        }
    }
}

fn top_sites(variable: &mut Vec<VariableSite>, max: usize) {
    variable.sort_by(|a, b| b.alt_fraction.total_cmp(&a.alt_fraction).then(b.depth.cmp(&a.depth)).then(a.pos.cmp(&b.pos)));
    variable.truncate(max);
}

fn accumulate(bam: &mut BamFile, region: &Region, reference: Option<(&[u8], u64)>, opts: &SummaryOptions, count_from: u64) -> Result<Parts> {
    let mut pileup_opts = opts.pileup.clone();
    pileup_opts.include_zero_depth = true;
    let mut p = Parts::default();
    // The deletion run still open: a deletion is one site however many bases it spans, so a
    // deleted base right after the last one extends it instead of starting another.
    let mut open: Option<VariableSite> = None;
    p.reads = crate::pileup::pileup_counting_from(bam, region, reference, &pileup_opts, count_from, |_| {}, |s: PositionStats| {
        *p.histogram.entry(s.depth).or_default() += 1;
        p.depth_sum += u64::from(s.depth);
        let (deletion, other) = sites_at(&s, opts);
        match (&mut open, deletion) {
            (Some(run), Some(site)) if run.end() == site.pos => run.extend(&site),
            (_, deletion) => {
                if let Some(done) = std::mem::replace(&mut open, deletion) {
                    p.file(done, region);
                }
            }
        }
        if let Some(site) = other {
            p.add(site);
        }
        Ok(())
    })?;
    if let Some(done) = open {
        p.file(done, region);
    }
    top_sites(&mut p.variable, opts.max_variable_sites);
    Ok(p)
}

/// The region in one pass, or, when it is long, in chunks on several threads (each with its own
/// handle on the file), added up to the same answer.
fn parts(bam: &mut BamFile, region: &Region, reference: Option<&[u8]>, opts: &SummaryOptions) -> Result<Parts> {
    // Every chunk gets the whole region's reference: read tests then see reads that cross a
    // chunk boundary whole, and a reference shorter than the region (a FASTA contig shorter
    // than the BAM's) just runs out, as in one pass.
    let reference = reference.map(|r| (r, region.start));
    let results = crate::chunks::run(bam, region, |b, chunk, count_from| accumulate(b, chunk, reference, opts, count_from))?;
    let mut total = Parts::default();
    // A deletion run that reached the end of the passes so far: the next pass may continue it.
    let mut carried: Option<VariableSite> = None;
    for p in results {
        for (d, c) in p.histogram {
            *total.histogram.entry(d).or_default() += c;
        }
        total.depth_sum += p.depth_sum;
        total.variable.extend(p.variable);
        total.variable_count += p.variable_count;
        total.reads.merge(&p.reads);
        let head = match (carried.take(), p.head) {
            (Some(mut run), Some(head)) if run.end() == head.pos => {
                run.extend(&head);
                Some(run)
            }
            (before, head) => {
                before.into_iter().for_each(|site| total.add(site));
                head
            }
        };
        if p.head_reaches_end {
            carried = head;
        } else {
            head.into_iter().for_each(|site| total.add(site));
            carried = p.tail;
        }
    }
    carried.into_iter().for_each(|site| total.add(site));
    top_sites(&mut total.variable, opts.max_variable_sites);
    Ok(total)
}

pub fn summarize(
    bam: &mut BamFile,
    region: &Region,
    reference: Option<&[u8]>,
    opts: &SummaryOptions,
) -> Result<RegionSummary> {
    let Parts { histogram, depth_sum, mut variable, variable_count, reads, .. } = parts(bam, region, reference, opts)?;
    variable.sort_by_key(|v| v.pos);

    let n = region.len().max(1);
    let at_least = |min: u32| -> f64 {
        let covered: u64 = histogram.range(min..).map(|(_, c)| c).sum();
        (covered as f64 / n as f64 * 1000.0).round() / 1000.0
    };

    let median_depth = {
        let half = n.div_ceil(2);
        let mut seen = 0;
        histogram
            .iter()
            .find(|(_, c)| {
                seen += **c;
                seen >= half
            })
            .map_or(0, |(d, _)| *d)
    };

    Ok(RegionSummary {
        region: region.to_string(),
        length: region.len(),
        reads,
        mean_depth: (depth_sum as f64 / n as f64 * 100.0).round() / 100.0,
        median_depth,
        min_depth: histogram.keys().next().copied().unwrap_or(0),
        max_depth: histogram.keys().next_back().copied().unwrap_or(0),
        breadth: [1, 10, 20, 30]
            .into_iter()
            .map(|m| (format!(">={m}"), at_least(m)))
            .collect(),
        variable_site_count: variable_count,
        variable_sites: variable,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{AltAllele, StrandCount};

    #[test]
    fn variable_sites_compare_unrounded() {
        let site = |alt: u32, depth: u32| PositionStats {
            depth,
            top_alt: Some(AltAllele { allele: "C", count: StrandCount { total: alt, fwd: alt, rev: 0 }, fraction: 0.2 }),
            ..Default::default()
        };
        let opts = SummaryOptions::default();
        assert!(variable_alt(&site(400, 2001), &opts).is_none());
        assert!(variable_alt(&site(400, 2000), &opts).is_some());
        // An exact boundary still counts.
        assert!(variable_alt(&site(3, 10), &SummaryOptions { min_alt_fraction: 0.3, ..Default::default() }).is_some());
    }

    fn deletion(pos: u64, length: u64, alt_fraction: f64) -> VariableSite {
        VariableSite { pos, ref_base: Some('A'), depth: 30, alt: DELETION.into(), alt_count: 20, alt_fwd: 10, alt_rev: 10, alt_fraction, length }
    }

    #[test]
    fn a_deletion_run_is_one_site_with_its_deepest_base() {
        let mut run = deletion(100, 1, 0.5);
        run.extend(&deletion(101, 1, 0.8));
        run.extend(&deletion(102, 1, 0.6));
        assert_eq!((run.pos, run.length, run.end(), run.alt_fraction), (100, 3, 103, 0.8));
        assert_eq!(run.change(), "DEL 3 bp");
    }

    /// A base the intact reads mismatch at, inside a heterozygous deletion: the mismatch is the
    /// top allele, and the deletion is still seen, so its run isn't cut there.
    #[test]
    fn a_deletion_under_a_commoner_mismatch_is_still_a_deletion_site() {
        let at = |top: (&'static str, u32), deleted: u32| PositionStats {
            pos: 500,
            depth: 33,
            deletions: StrandCount { total: deleted, fwd: deleted, rev: 0 },
            top_alt: Some(AltAllele { allele: top.0, count: StrandCount { total: top.1, fwd: top.1, rev: 0 }, fraction: 0.0 }),
            ..Default::default()
        };
        let opts = SummaryOptions::default();
        let (deletion, other) = sites_at(&at(("A", 15), 11), &opts);
        assert_eq!(deletion.map(|d| (d.alt, d.alt_count, d.alt_fraction)), Some((DELETION.to_owned(), 11, 0.333)));
        assert_eq!(other.map(|o| (o.alt, o.alt_count)), Some(("A".to_owned(), 15)));
        // The deletion as the top allele is one site, not two.
        let (deletion, other) = sites_at(&at((DELETION, 11), 11), &opts);
        assert!(deletion.is_some() && other.is_none());
        // Too few reads delete the base: only the mismatch.
        let (deletion, other) = sites_at(&at(("A", 15), 3), &opts);
        assert!(deletion.is_none() && other.is_some());
    }

    #[test]
    fn runs_at_a_pass_edge_wait_for_the_next_pass() {
        // A pass over 0-based [100, 200): positions 101..=200.
        let region = Region { contig: "c".into(), start: 100, end: 200 };
        let mut p = Parts::default();
        p.file(deletion(101, 5, 0.5), &region);
        p.file(deletion(150, 5, 0.5), &region);
        p.file(deletion(196, 5, 0.5), &region);
        assert_eq!((p.head.as_ref().map(|s| s.pos), p.tail.as_ref().map(|s| s.pos), p.variable_count, p.head_reaches_end), (Some(101), Some(196), 1, false));
        let mut whole = Parts::default();
        whole.file(deletion(101, 100, 0.5), &region);
        assert!(whole.head_reaches_end && whole.tail.is_none() && whole.variable_count == 0);
    }
}
