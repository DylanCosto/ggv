//! GGV_TEST_PHASED_BAM + GGV_TEST_REF + GGV_TEST_PHASED_VCF + GGV_TEST_PHASED_REGIONS: reads phased
//! by their own alleles must agree (up to swapping A and B) with the haplotype each read gets
//! from an independently phased VCF, e.g. phase from transmission in a trio, which uses no reads
//! at all.

use ggv_layout::{AlignedRead, ViewOptions, load_view, phase::phase};

/// The read's haplotype from phased SNVs (`a|b`, first sample): +1 votes for the first
/// haplotype, −1 for the second.
fn truth(r: &AlignedRead, phased: &[(u64, u8, u8, u8)]) -> Option<bool> {
    let mut vote = 0i32;
    for &(pos, _, alt, first) in phased.iter().filter(|p| r.start <= p.0 && p.0 < r.end) {
        let covered = r.blocks.iter().any(|b| b.ref_start <= pos && pos < b.ref_start + u64::from(b.len));
        if !covered {
            continue;
        }
        let allele = match r.mismatches.iter().find(|m| m.pos == pos) {
            Some(m) if m.base == alt => 1,
            Some(_) => continue,
            None => 0,
        };
        vote += if allele == first { 1 } else { -1 };
    }
    (vote != 0).then_some(vote > 0)
}

#[test]
fn allele_phasing_agrees_with_trio_phase() {
    let vars = ["GGV_TEST_PHASED_BAM", "GGV_TEST_REF", "GGV_TEST_PHASED_VCF", "GGV_TEST_PHASED_REGIONS"].map(std::env::var);
    let [Ok(bam), Ok(reference), Ok(vcf), Ok(regions)] = vars else { return };
    let mut b = ggv_io::BamFile::open(&bam).unwrap();
    let mut fasta = ggv_io::ReferenceFasta::open(&reference).unwrap();
    let mut calls = ggv_io::VariantFile::open(&vcf).unwrap();
    for locus in regions.split(',') {
        let region = b.resolve_str(locus).unwrap();
        let seq = fasta.fetch(&region).unwrap();
        let view = load_view(&mut b, &region, Some(seq), String::new(), &ViewOptions::default()).unwrap();
        let (reads, _) = view.reads.as_ref().unwrap();
        // Phased heterozygous SNVs: (pos, ref, alt, allele on the first haplotype).
        let phased: Vec<(u64, u8, u8, u8)> = calls
            .query(&region, false)
            .unwrap()
            .into_iter()
            .filter(|v| v.reference.len() == 1 && v.alts.len() == 1 && v.alts[0].len() == 1)
            .filter_map(|v| {
                let gt = v.samples.first()?.gt.as_bytes().to_vec();
                (gt.len() == 3 && gt[1] == b'|' && gt[0] != gt[2]).then(|| (v.start, v.reference.as_bytes()[0], v.alts[0].as_bytes()[0], gt[0] - b'0'))
            })
            .collect();
        let members: Vec<usize> = (0..reads.len()).collect();
        let t = std::time::Instant::now();
        let p = phase(reads, &members);
        let took = t.elapsed();
        let (mut same, mut swapped) = (0, 0);
        for (i, r) in reads.iter().enumerate() {
            if let (Some(want), Some(got)) = (truth(r, &phased), p.hap[i]) {
                if want == got { same += 1 } else { swapped += 1 }
            }
        }
        let both = same + swapped;
        let agree = same.max(swapped) as f64 / both.max(1) as f64;
        eprintln!(
            "{locus}: {} reads, {} phased by alleles ({} sites, {} blocks) in {took:.2?}; {} trio-phased SNVs; agreement {:.1}% over {both} reads placed by both",
            reads.len(),
            p.hap.iter().filter(|h| h.is_some()).count(),
            p.sites.len(),
            p.blocks,
            phased.len(),
            agree * 100.0
        );
        if p.blocks == 1 && both >= 20 {
            assert!(agree >= 0.95, "{locus}: {agree}");
        }
    }
}
