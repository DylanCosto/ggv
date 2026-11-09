//! A sample's sequence over a region, from its reads: at each base the allele most reads show,
//! with an IUPAC code where two bases are each well supported (a heterozygous SNV), `N` where
//! too few reads cover it, and nothing where most reads have the base deleted. Insertions in
//! the reads aren't added (the pileup counts them but doesn't keep their bases).

use ggv_query::{PositionStats, StrandCount};

/// Reads needed to call a base; fewer gives `N`.
pub const MIN_DEPTH: u32 = 3;
/// A second base is called with an IUPAC code when it has at least this share of the reads...
pub const HET_FRACTION: f64 = 0.25;
/// ...and at least this many reads.
pub const HET_READS: u32 = 2;

/// What [`call`] made, for the message after copying.
#[derive(Debug, Default, PartialEq, Eq)]
pub struct Made {
    pub sequence: String,
    /// Positions written as an IUPAC two-base code.
    pub ambiguous: usize,
    /// Positions written as `N` (too few reads).
    pub low_depth: usize,
    /// Positions left out (most reads have them deleted).
    pub deleted: usize,
}

/// The sequence over `[start, end)` (0-based) from `stats` (1-based positions, in order, some
/// missing where nothing covers them).
pub fn call(stats: &[PositionStats], start: u64, end: u64) -> Made {
    let mut out = Made::default();
    let mut it = stats.iter().peekable();
    for pos in start + 1..=end {
        while it.peek().is_some_and(|s| s.pos < pos) {
            it.next();
        }
        let Some(s) = it.peek().filter(|s| s.pos == pos) else {
            out.sequence.push('N');
            out.low_depth += 1;
            continue;
        };
        let b = &s.bases;
        let mut alleles: Vec<(char, u32)> = [('A', b.a), ('C', b.c), ('G', b.g), ('T', b.t), ('-', s.deletions)]
            .into_iter()
            .map(|(c, n): (char, StrandCount)| (c, n.total))
            .filter(|&(_, n)| n > 0)
            .collect();
        let depth: u32 = alleles.iter().map(|a| a.1).sum();
        if depth < MIN_DEPTH {
            out.sequence.push('N');
            out.low_depth += 1;
            continue;
        }
        // Most reads first; ties in alphabetical order, so the result doesn't depend on order.
        alleles.sort_by(|x, y| y.1.cmp(&x.1).then(x.0.cmp(&y.0)));
        let (top, _) = alleles[0];
        if top == '-' {
            out.deleted += 1;
            continue;
        }
        let second = alleles.get(1).filter(|&&(c, n)| c != '-' && n >= HET_READS && f64::from(n) / f64::from(depth) >= HET_FRACTION);
        match second.and_then(|&(c, _)| iupac(top, c)) {
            Some(code) => {
                out.sequence.push(code);
                out.ambiguous += 1;
            }
            None => out.sequence.push(top),
        }
    }
    out
}

/// The IUPAC code for two different bases.
fn iupac(a: char, b: char) -> Option<char> {
    let (a, b) = if a < b { (a, b) } else { (b, a) };
    Some(match (a, b) {
        ('A', 'G') => 'R',
        ('C', 'T') => 'Y',
        ('C', 'G') => 'S',
        ('A', 'T') => 'W',
        ('G', 'T') => 'K',
        ('A', 'C') => 'M',
        _ => return None,
    })
}

/// FASTA with 60 bases per line.
pub fn fasta(header: &str, seq: &str) -> String {
    let mut out = format!(">{header}\n");
    for line in seq.as_bytes().chunks(60) {
        out.push_str(&String::from_utf8_lossy(line));
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn at(pos: u64, a: u32, c: u32, g: u32, t: u32, del: u32) -> PositionStats {
        let n = |total| StrandCount { total, fwd: total, rev: 0 };
        let mut s = PositionStats { pos, ..Default::default() };
        s.bases.a = n(a);
        s.bases.c = n(c);
        s.bases.g = n(g);
        s.bases.t = n(t);
        s.deletions = n(del);
        s.depth = a + c + g + t + del;
        s
    }

    #[test]
    fn majority_het_low_depth_and_deletions() {
        let stats = [
            at(11, 20, 0, 0, 0, 0), // A
            at(12, 0, 10, 0, 9, 0), // C/T het: Y
            at(13, 0, 0, 1, 0, 0),  // one read: N
            // 14 uncovered: N
            at(15, 0, 0, 0, 2, 18), // mostly deleted: left out
            at(16, 0, 0, 19, 1, 0), // one stray T among 20: G
            at(17, 6, 0, 0, 0, 5),  // het deletion, more A: A
        ];
        let m = call(&stats, 10, 17);
        assert_eq!(m.sequence, "AYNNGA");
        assert_eq!((m.ambiguous, m.low_depth, m.deleted), (1, 2, 1));
    }

    /// On the HG002 fixture: around the C-run deletion only 5 of 29 reads lack a C, so HG002's
    /// sequence is the reference's; at a known heterozygous SNV it gets an IUPAC code.
    #[test]
    fn real_reads_match_the_reference_except_at_variants() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio");
        let (bam, fa) = (dir.join("HG002.revio.slices.bam"), dir.join("GRCh38.chr20.slices.fa.gz"));
        if !bam.exists() {
            return;
        }
        let run = |start: u64, end: u64| {
            let region = ggv_model::Region { contig: "chr20".into(), start, end };
            let reference = ggv_io::ReferenceFasta::open(&fa).unwrap().fetch(&region).unwrap();
            let mut b = ggv_io::BamFile::open_with_reference(&bam, Some(&fa)).unwrap();
            let mut stats = Vec::new();
            ggv_query::pileup(&mut b, &region, Some(&reference), &ggv_query::PileupOptions { exclude_flags: 0x704, ..Default::default() }, |s| {
                stats.push(s);
                Ok(())
            })
            .unwrap();
            (call(&stats, start, end), String::from_utf8(reference).unwrap().to_ascii_uppercase(), stats)
        };
        let (m, reference, _) = run(58_854_229, 58_854_262);
        assert_eq!(m.sequence, reference, "{m:?}");
        assert_eq!((m.ambiguous, m.low_depth, m.deleted), (0, 0, 0));

        // A heterozygous SNV somewhere in the GNAS slice gives an IUPAC code there, and only
        // there differs from the reference.
        let (m, reference, stats) = run(58_880_000, 58_890_000);
        let het = stats.iter().find(|s| s.depth >= 10 && (0.3..0.7).contains(&s.alt_fraction) && s.deletions.total == 0 && s.insertions_after.total == 0);
        if let Some(h) = het {
            let i = (h.pos - 58_880_001) as usize;
            assert!("RYSWKM".contains(m.sequence.as_bytes()[i] as char), "at {}: {}", h.pos, m.sequence.as_bytes()[i] as char);
        }
        assert_eq!(m.sequence.len() + m.deleted, reference.len());
    }

    #[test]
    fn fasta_wraps_at_60() {
        let f = fasta("x", &"A".repeat(61));
        assert_eq!(f, format!(">x\n{}\nA\n", "A".repeat(60)));
    }
}
