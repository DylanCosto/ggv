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

#[derive(Clone, Debug, Serialize)]
pub struct VariableSite {
    pub pos: u64,
    #[serde(rename = "ref")]
    pub ref_base: Option<char>,
    pub depth: u32,
    pub alt: String,
    pub alt_count: u32,
    pub alt_fwd: u32,
    pub alt_rev: u32,
    pub alt_fraction: f64,
}

/// The listed allele of a variable site: deep enough, and its share of the reads at least the
/// minimum, compared unrounded (400 of 2001 reads is below 0.2, though shown as 0.2).
fn variable_alt<'a>(s: &'a PositionStats, opts: &SummaryOptions) -> Option<&'a crate::AltAllele> {
    let alt = s.top_alt.as_ref()?;
    (s.depth >= opts.min_variable_depth && s.depth > 0 && f64::from(alt.count.total) / f64::from(s.depth) >= opts.min_alt_fraction).then_some(alt)
}

/// What one pass over (part of) a region adds up.
#[derive(Default)]
struct Parts {
    histogram: BTreeMap<u32, u64>,
    depth_sum: u64,
    /// The most variable sites, already cut to the number asked for.
    variable: Vec<VariableSite>,
    variable_count: u64,
    reads: ReadCounts,
}

fn top_sites(variable: &mut Vec<VariableSite>, max: usize) {
    variable.sort_by(|a, b| b.alt_fraction.total_cmp(&a.alt_fraction).then(b.depth.cmp(&a.depth)).then(a.pos.cmp(&b.pos)));
    variable.truncate(max);
}

fn accumulate(bam: &mut BamFile, region: &Region, reference: Option<(&[u8], u64)>, opts: &SummaryOptions, count_from: u64) -> Result<Parts> {
    let mut pileup_opts = opts.pileup.clone();
    pileup_opts.include_zero_depth = true;
    let mut p = Parts::default();
    p.reads = crate::pileup::pileup_counting_from(bam, region, reference, &pileup_opts, count_from, |_| {}, |s: PositionStats| {
        *p.histogram.entry(s.depth).or_default() += 1;
        p.depth_sum += u64::from(s.depth);
        // The listed allele must pass on its own: a site whose non-reference reads are split
        // among several alleles, none of them common, is noise rather than a variant.
        if let Some(alt) = variable_alt(&s, opts) {
            p.variable_count += 1;
            p.variable.push(VariableSite {
                pos: s.pos,
                ref_base: s.ref_base,
                depth: s.depth,
                alt: alt.allele.to_owned(),
                alt_count: alt.count.total,
                alt_fwd: alt.count.fwd,
                alt_rev: alt.count.rev,
                alt_fraction: alt.fraction,
            });
        }
        Ok(())
    })?;
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
    for p in results {
        for (d, c) in p.histogram {
            *total.histogram.entry(d).or_default() += c;
        }
        total.depth_sum += p.depth_sum;
        total.variable.extend(p.variable);
        total.variable_count += p.variable_count;
        total.reads.merge(&p.reads);
    }
    top_sites(&mut total.variable, opts.max_variable_sites);
    Ok(total)
}

pub fn summarize(
    bam: &mut BamFile,
    region: &Region,
    reference: Option<&[u8]>,
    opts: &SummaryOptions,
) -> Result<RegionSummary> {
    let Parts { histogram, depth_sum, mut variable, variable_count, reads } = parts(bam, region, reference, opts)?;
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
}
