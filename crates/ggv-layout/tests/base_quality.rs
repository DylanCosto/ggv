//! GGV_TEST_BAM + GGV_TEST_REF + GGV_TEST_REGIONS: with a minimum base quality, coverage from
//! decoded reads (what the view draws, and what MCP `reads`, `phase` and `render` use) must give
//! the same depth and allele counts as the pileup engine (`evidence`, `summarize`) at every
//! position. Regression test: MCP coverage, reads, phase and render once ignored the minimum base
//! quality.

use ggv_layout::{Coverage, ViewOptions, load_view, phase::phase};
use ggv_query::PileupOptions;

#[test]
fn read_coverage_matches_pileup_under_a_base_quality_minimum() {
    let vars = ["GGV_TEST_BAM", "GGV_TEST_REF", "GGV_TEST_REGIONS"].map(std::env::var);
    let [Ok(bam), Ok(reference), Ok(regions)] = vars else { return };
    let mut b = ggv_io::BamFile::open(&bam).unwrap();
    let mut fasta = ggv_io::ReferenceFasta::open(&reference).unwrap();
    for locus in regions.split(',') {
        let region = b.resolve_str(locus).unwrap();
        let seq = fasta.fetch(&region).unwrap();
        for min_bq in [0u8, 20, 30, 40] {
            let filters = PileupOptions { min_base_quality: min_bq, ..Default::default() };
            let view = load_view(&mut b, &region, Some(seq.clone()), String::new(), &ViewOptions { filters: filters.clone(), ..Default::default() }).unwrap();
            let Coverage::PerBase(cols) = &view.coverage else { panic!("{locus}: not per base") };
            let mut checked = 0;
            let mut low = 0u64;
            ggv_query::pileup(&mut b, &region, Some(&seq), &filters, |p| {
                let c = &cols[(p.pos - 1 - region.start) as usize];
                assert_eq!(c.depth, p.depth, "{locus} q{min_bq} depth at {}", p.pos);
                if matches!(c.ref_base, Some(b'A' | b'C' | b'G' | b'T')) {
                    let want = [p.bases.a.total, p.bases.c.total, p.bases.g.total, p.bases.t.total, p.bases.n.total];
                    assert_eq!(c.counts, want, "{locus} q{min_bq} alleles at {}", p.pos);
                }
                low += u64::from(p.low_base_quality);
                checked += 1;
                Ok(())
            })
            .unwrap();
            // Phasing sees only observed bases: no read is counted at a site where its base is
            // below the minimum.
            let (reads, _) = view.reads.as_ref().unwrap();
            let p = phase(reads, &(0..reads.len()).collect::<Vec<_>>());
            for &(pos, _) in &p.sites {
                let observed = reads.iter().filter(|r| r.start <= pos && pos < r.end && !r.is_low_quality(pos)).count();
                assert!(observed > 0, "{locus}: site {pos} without an observed base");
            }
            eprintln!("{locus} q{min_bq}: {checked} positions agree; {low} bases below the minimum; {} het sites", p.sites.len());
        }
    }
}
