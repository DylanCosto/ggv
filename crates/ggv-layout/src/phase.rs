//! Phasing reads by the alleles they carry, for data without `HP` tags.
//!
//! Heterozygous sites are found from the reads themselves: positions where one alternate base
//! makes up 25–75% of the reads covering it (at least 6 reads, 3 of them alternate). Sites are
//! then oriented one at a time, left to right: each is lined up with the haplotype that the
//! reads covering it were already assigned to by earlier sites (a greedy single-pass phaser,
//! the core of what WhatsHap-style tools refine). Reads end up in haplotype A or B by the
//! balance of their votes; reads with no informative site, or an even split, in neither.
//!
//! Sites come from mismatches, so they need reads decoded against a reference.

use std::collections::HashMap;

use crate::AlignedRead;

/// Minimum base quality for a base to count towards a site.
const MIN_BASE_QUALITY: u8 = 10;
const MIN_DEPTH: u32 = 6;
const MIN_ALT: u32 = 3;
const MIN_FRACTION: f64 = 0.25;
const MAX_FRACTION: f64 = 0.75;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Phasing {
    /// Per read index: `Some(true)` haplotype A, `Some(false)` B, `None` unassigned. Indices
    /// not among the phased reads are `None`.
    pub hap: Vec<Option<bool>>,
    /// The heterozygous sites used (0-based positions), with their alternate base.
    pub sites: Vec<(u64, u8)>,
    /// Stretches of sites linked by shared reads; haplotypes A and B are only consistent
    /// within one.
    pub blocks: usize,
    /// Per site (as in `sites`): its block, numbered from 0 left to right.
    pub site_blocks: Vec<usize>,
    /// Per read index (as in `hap`): the block of the sites it covers, if any.
    pub read_blocks: Vec<Option<usize>>,
}

/// Whether the read has an observed base at `pos` (aligned, and not below the minimum base
/// quality the reads were decoded with).
fn covers(r: &AlignedRead, pos: u64) -> bool {
    // Blocks are in reference order: a binary search, not a scan (a nanopore read has hundreds
    // of blocks and spans thousands of candidate sites).
    let i = r.blocks.partition_point(|b| b.ref_start + u64::from(b.len) <= pos);
    r.blocks.get(i).is_some_and(|b| b.ref_start <= pos) && !r.is_low_quality(pos)
}

fn usable(q: u8) -> bool {
    q == 0xff || q >= MIN_BASE_QUALITY
}

/// Phases `members` (indices into `reads`).
pub fn phase(reads: &[AlignedRead], members: &[usize]) -> Phasing {
    // Candidate positions: enough mismatches of one base.
    let mut alt_counts: HashMap<u64, [u32; 5]> = HashMap::new();
    for &i in members {
        for m in reads[i].mismatches.iter().filter(|m| usable(m.quality)) {
            let k = match m.base {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4,
            };
            alt_counts.entry(m.pos).or_default()[k] += 1;
        }
    }
    let mut candidates: Vec<(u64, u8, u32)> = alt_counts
        .into_iter()
        .filter_map(|(pos, c)| {
            let (k, &n) = c[..4].iter().enumerate().max_by_key(|(_, n)| **n)?;
            (n >= MIN_ALT).then_some((pos, b"ACGT"[k], n))
        })
        .collect();
    candidates.sort_unstable_by_key(|c| c.0);

    // Depth at each candidate, then keep the heterozygous-looking ones.
    let mut depth = vec![0u32; candidates.len()];
    for &i in members {
        let r = &reads[i];
        let from = candidates.partition_point(|c| c.0 < r.start);
        for (j, c) in candidates[from..].iter().enumerate().take_while(|(_, c)| c.0 < r.end) {
            if covers(r, c.0) {
                depth[from + j] += 1;
            }
        }
    }
    let sites: Vec<(u64, u8)> = candidates
        .iter()
        .zip(&depth)
        .filter(|((_, _, alt), d)| {
            let f = f64::from(*alt) / f64::from((**d).max(1));
            **d >= MIN_DEPTH && (MIN_FRACTION..=MAX_FRACTION).contains(&f)
        })
        .map(|((p, b, _), _)| (*p, *b))
        .collect();

    // Each read's alleles at the sites it covers: +1 alternate, −1 reference.
    let mut obs: Vec<Vec<(usize, i32)>> = vec![Vec::new(); members.len()];
    for (mi, &i) in members.iter().enumerate() {
        let r = &reads[i];
        let from = sites.partition_point(|s| s.0 < r.start);
        for (j, &(pos, alt)) in sites[from..].iter().enumerate().take_while(|(_, s)| s.0 < r.end) {
            if !covers(r, pos) {
                continue;
            }
            // Mismatches are in position order.
            let at = r.mismatches.binary_search_by_key(&pos, |m| m.pos).ok().map(|k| &r.mismatches[k]);
            match at {
                Some(m) if m.base == alt && usable(m.quality) => obs[mi].push((from + j, 1)),
                Some(_) => {} // another base, or too low quality: no evidence
                None => obs[mi].push((from + j, -1)),
            }
        }
    }

    // Site by site: orient it to agree with the votes reads already carry, then add its votes.
    let mut by_site: Vec<Vec<(usize, i32)>> = vec![Vec::new(); sites.len()];
    for (mi, o) in obs.iter().enumerate() {
        for &(s, a) in o {
            by_site[s].push((mi, a));
        }
    }
    let mut votes = vec![0i32; members.len()];
    // Blocks: a site joins the block of the phased reads covering it (merging blocks those
    // reads join up), or starts one. Union-find over block numbers.
    let mut parent: Vec<usize> = Vec::new();
    fn root(parent: &mut [usize], mut b: usize) -> usize {
        while parent[b] != b {
            parent[b] = parent[parent[b]];
            b = parent[b];
        }
        b
    }
    let mut read_block: Vec<Option<usize>> = vec![None; members.len()];
    let mut site_block = Vec::with_capacity(sites.len());
    for readings in &by_site {
        let score: i32 = readings.iter().filter(|(mi, _)| votes[*mi] != 0).map(|(mi, a)| votes[*mi].signum() * a).sum();
        let mut block = None;
        for &(mi, _) in readings.iter().filter(|(mi, _)| votes[*mi] != 0) {
            let Some(b) = read_block[mi] else { continue };
            let r = root(&mut parent, b);
            match block {
                None => block = Some(r),
                Some(k) if k != r => parent[r.max(k)] = r.min(k),
                Some(_) => {}
            }
            if let Some(k) = block {
                block = Some(root(&mut parent, k));
            }
        }
        let block = block.unwrap_or_else(|| {
            parent.push(parent.len());
            parent.len() - 1
        });
        site_block.push(block);
        let orient = if score < 0 { -1 } else { 1 };
        for &(mi, a) in readings {
            votes[mi] += orient * a;
            read_block[mi] = Some(block);
        }
    }
    // Number the blocks 0, 1, ... from the left.
    let mut number: HashMap<usize, usize> = HashMap::new();
    let site_blocks: Vec<usize> = site_block
        .iter()
        .map(|&b| {
            let r = root(&mut parent, b);
            let next = number.len();
            *number.entry(r).or_insert(next)
        })
        .collect();
    let blocks = number.len();

    let mut hap = vec![None; reads.len()];
    let mut read_blocks = vec![None; reads.len()];
    for (mi, &i) in members.iter().enumerate() {
        hap[i] = match votes[mi].signum() {
            1 => Some(true),
            -1 => Some(false),
            _ => None,
        };
        read_blocks[i] = read_block[mi].map(|b| number[&root(&mut parent, b)]);
    }
    Phasing { hap, sites, blocks, site_blocks, read_blocks }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Block, Mismatch};

    fn read(start: u64, end: u64, alts: &[(u64, u8)]) -> AlignedRead {
        AlignedRead {
            name: String::new(),
            start,
            end,
            reverse: false,
            mapq: 60,
            flags: 0,
            blocks: vec![Block { ref_start: start, query_start: 0, len: (end - start) as u32 }],
            deletions: vec![],
            ref_skips: vec![],
            insertions: vec![],
            mismatches: alts.iter().map(|&(pos, base)| Mismatch { pos, base, quality: 30 }).collect(),
            soft_clip_left: 0,
            soft_clip_right: 0,
            clip_mismatches: vec![],
            haplotype: None,
            read_group: None,
            template_length: 0,
            mate: None,
            supplementary: vec![],
            mods: vec![],
            low_quality: vec![],
            tag: None,
            bases: None,
        }
    }

    #[test]
    fn two_haplotypes_chained_through_overlapping_reads() {
        // Sites at 100 (T), 300 (G), 500 (C). Haplotype 1 carries T, ref, C; haplotype 2
        // carries ref, G, ref. Reads span one or two sites, so phase passes along the chain.
        let mut reads = Vec::new();
        for _ in 0..5 {
            reads.push(read(50, 350, &[(100, b'T')])); // h1: 100 alt, 300 ref
            reads.push(read(50, 350, &[(300, b'G')])); // h2: 100 ref, 300 alt
            reads.push(read(250, 550, &[(500, b'C')])); // h1: 300 ref, 500 alt
            reads.push(read(250, 550, &[(300, b'G')])); // h2: 300 alt, 500 ref
        }
        // A read with no site, and one with a sequencing error only.
        reads.push(read(600, 700, &[]));
        reads.push(read(50, 150, &[(120, b'A')]));
        let members: Vec<usize> = (0..reads.len()).collect();
        let p = phase(&reads, &members);
        assert_eq!(p.sites, vec![(100, b'T'), (300, b'G'), (500, b'C')]);
        assert_eq!(p.blocks, 1);
        let h1 = p.hap[0];
        assert!(h1.is_some());
        for k in 0..5 {
            assert_eq!(p.hap[4 * k], h1, "h1 read at 100");
            assert_eq!(p.hap[4 * k + 2], h1, "h1 read at 500: linked through site 300");
            assert_eq!(p.hap[4 * k + 1], h1.map(|h| !h));
            assert_eq!(p.hap[4 * k + 3], h1.map(|h| !h));
        }
        assert_eq!(p.hap[20], None, "no site in the read");
        assert_eq!(p.hap[21], h1.map(|h| !h), "the error at 120 is no site; the read is ref at 100: h2");
    }

    #[test]
    fn homozygous_and_rare_changes_are_not_sites() {
        let mut reads = Vec::new();
        for k in 0..10 {
            // 100: every read alt (homozygous). 200: one read alt (noise).
            let mut alts = vec![(100, b'T')];
            if k == 0 {
                alts.push((200, b'C'));
            }
            reads.push(read(0, 300, &alts));
        }
        let p = phase(&reads, &(0..10).collect::<Vec<_>>());
        assert!(p.sites.is_empty());
        assert!(p.hap.iter().all(Option::is_none));
        assert_eq!(p.blocks, 0);
    }

    #[test]
    fn unlinked_sites_make_separate_blocks() {
        let mut reads = Vec::new();
        for _ in 0..4 {
            reads.push(read(0, 200, &[(100, b'T')]));
            reads.push(read(0, 200, &[]));
            reads.push(read(1000, 1200, &[(1100, b'G')]));
            reads.push(read(1000, 1200, &[]));
        }
        let p = phase(&reads, &(0..16).collect::<Vec<_>>());
        assert_eq!(p.sites.len(), 2);
        assert_eq!(p.blocks, 2);
        assert_eq!(p.site_blocks, vec![0, 1]);
        for k in 0..4 {
            assert_eq!(p.read_blocks[4 * k], Some(0));
            assert_eq!(p.read_blocks[4 * k + 2], Some(1));
        }
    }

    #[test]
    fn linked_sites_share_a_block() {
        // Sites at 110 and 130 joined by the same reads; 300 on its own.
        let mut reads = Vec::new();
        for _ in 0..5 {
            reads.push(read(100, 140, &[(110, b'T'), (130, b'G')]));
            reads.push(read(100, 140, &[]));
            reads.push(read(290, 310, &[(300, b'C')]));
            reads.push(read(290, 310, &[]));
        }
        reads.push(read(500, 600, &[]));
        let p = phase(&reads, &(0..reads.len()).collect::<Vec<_>>());
        assert_eq!(p.sites.iter().map(|s| s.0).collect::<Vec<_>>(), vec![110, 130, 300]);
        assert_eq!(p.blocks, 2);
        assert_eq!(p.site_blocks, vec![0, 0, 1]);
        assert_eq!(p.read_blocks[0], Some(0));
        assert_eq!(p.read_blocks[3], Some(1));
        assert_eq!(p.read_blocks[20], None, "covers no site");
    }
}
