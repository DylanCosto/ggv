//! Checks that coverage computed from decoded reads (what the renderer and app draw) matches
//! the pileup engine (what `ggv query pileup` reports, itself checked against samtools).
//!
//! Needs real data, so it runs only when these are set (otherwise it passes trivially):
//!   GGV_TEST_BAM     indexed BAM
//!   GGV_TEST_REF     indexed FASTA (optional)
//!   GGV_TEST_REGIONS comma-separated loci, e.g. `chr20:144000-147000,chr20:1890000-1940000`

use ggv_io::{BamFile, ReferenceFasta};
use ggv_layout::{Coverage, ViewOptions, load_view};
use ggv_query::PileupOptions;

#[test]
fn coverage_from_reads_matches_pileup() {
    let (Ok(bam_path), Ok(regions)) = (std::env::var("GGV_TEST_BAM"), std::env::var("GGV_TEST_REGIONS"))
    else {
        eprintln!("GGV_TEST_BAM / GGV_TEST_REGIONS not set; skipping");
        return;
    };
    let ref_path = std::env::var("GGV_TEST_REF").ok();

    for locus in regions.split(',') {
        let mut bam = BamFile::open_with_reference(&bam_path, std::env::var("GGV_TEST_REF").ok()).unwrap();
        let region = bam.resolve_str(locus).unwrap();
        let reference = ref_path
            .as_ref()
            .map(|p| ReferenceFasta::open(p).unwrap().fetch(&region).unwrap());

        let view = load_view(&mut bam, &region, reference.clone(), String::new(), &ViewOptions::default()).unwrap();
        let Coverage::PerBase(columns) = &view.coverage else { panic!("{locus}: expected per-base coverage") };

        let mut checked = 0;
        ggv_query::pileup(&mut bam, &region, reference.as_deref(), &PileupOptions::default(), |s| {
            let c = &columns[(s.pos - 1 - region.start) as usize];
            assert_eq!(c.depth, s.depth, "{locus} pos {}: depth", s.pos);
            assert_eq!(c.deletions, s.deletions.total, "{locus} pos {}: deletions", s.pos);
            if matches!(s.ref_base, Some('A' | 'C' | 'G' | 'T')) {
                let b = &s.bases;
                assert_eq!(
                    c.counts,
                    [b.a.total, b.c.total, b.g.total, b.t.total, b.n.total],
                    "{locus} pos {}: allele counts",
                    s.pos
                );
                assert!((c.alt_fraction - s.alt_fraction).abs() < 1e-9, "{locus} pos {}: alt fraction", s.pos);
            }
            checked += 1;
            Ok(())
        })
        .unwrap();
        eprintln!("{locus}: {checked} positions match");
        assert!(checked > 0, "{locus}: no covered positions");
    }
}

#[test]
fn binned_fast_path_matches_per_base() {
    let (Ok(bam_path), Ok(regions)) = (std::env::var("GGV_TEST_BAM"), std::env::var("GGV_TEST_REGIONS")) else { return };
    for locus in regions.split(',') {
        let mut bam = BamFile::open_with_reference(&bam_path, std::env::var("GGV_TEST_REF").ok()).unwrap();
        let region = bam.resolve_str(locus).unwrap();
        let per_base = load_view(&mut bam, &region, None, String::new(), &ViewOptions::default()).unwrap();
        let Coverage::PerBase(cols) = &per_base.coverage else { panic!("expected per-base") };
        let bins = 500;
        let opts = ViewOptions { bins, max_read_window: 0, ..Default::default() };
        let binned = load_view(&mut bam, &region, None, String::new(), &opts).unwrap();
        let Coverage::Binned(got) = &binned.coverage else { panic!("expected binned") };

        let len = cols.len();
        for (b, bin) in got.iter().enumerate() {
            // Base i is in bin floor(i * bins / len), so bin b is [ceil(b*len/bins), ceil((b+1)*len/bins)).
            let slice = &cols[(b * len).div_ceil(bins)..((b + 1) * len).div_ceil(bins)];
            let max = slice.iter().map(|c| c.depth).max().unwrap_or(0);
            let mean = slice.iter().map(|c| f64::from(c.depth)).sum::<f64>() / slice.len().max(1) as f64;
            assert_eq!(bin.max, max, "{locus} bin {b}: max");
            assert!((bin.mean - mean).abs() < 1e-9, "{locus} bin {b}: mean {} vs {mean}", bin.mean);
        }
        assert_eq!(per_base.read_counts, binned.read_counts, "{locus}: read filter counts");
        eprintln!("{locus}: {bins} bins match");
    }
}

#[test]
fn dump_mod_calls() {
    // GGV_TEST_MODS_BAM + GGV_TEST_MODS_REGION + GGV_TEST_MODS_OUT: writes `pos code count prob_sum`
    // per (1-based position, code) for comparison with `samtools mpileup --output-mods`.
    let (Ok(bam), Ok(locus), Ok(out)) =
        (std::env::var("GGV_TEST_MODS_BAM"), std::env::var("GGV_TEST_MODS_REGION"), std::env::var("GGV_TEST_MODS_OUT"))
    else {
        return;
    };
    let mut bam = BamFile::open(&bam).unwrap();
    let region = bam.resolve_str(&locus).unwrap();
    let view = load_view(&mut bam, &region, None, String::new(), &ViewOptions::default()).unwrap();
    let mut tally: std::collections::BTreeMap<(u64, u8), (u32, u64)> = Default::default();
    for r in &view.reads.as_ref().unwrap().0 {
        for m in r.mods.iter().filter(|m| region.contains(m.pos(r.start))) {
            let e = tally.entry((m.pos(r.start) + 1, m.code)).or_default();
            e.0 += 1;
            e.1 += u64::from(m.prob);
        }
    }
    let text: String = tally.iter().map(|((p, c), (n, s))| format!("{p}\t{}\t{n}\t{s}\n", *c as char)).collect();
    std::fs::write(&out, text).unwrap();
    eprintln!("{} (position, code) entries", tally.len());
}

/// The index estimate against counted depth: over each region (a few Mb of a whole-genome BAM),
/// 50 bins of estimated depth should track the counted binned depth. Prints both and checks
/// the overall level is within 25% and the bins correlate.
#[test]
fn index_estimate_tracks_counted_depth() {
    let (Ok(bam_path), Ok(regions)) = (std::env::var("GGV_TEST_WGS_BAM"), std::env::var("GGV_TEST_WGS_REGIONS")) else { return };
    for locus in regions.split(',') {
        let mut bam = BamFile::open(&bam_path).unwrap();
        let region = bam.resolve_str(locus).unwrap();
        let opts = ViewOptions { bins: 50, max_read_window: 0, max_coverage_window: u64::MAX, ..Default::default() };
        let t = std::time::Instant::now();
        let counted = load_view(&mut bam, &region, None, String::new(), &opts).unwrap();
        let t_count = t.elapsed();
        let t = std::time::Instant::now();
        let est = load_view(&mut bam, &region, None, String::new(), &ViewOptions { max_coverage_window: 0, ..opts }).unwrap();
        let t_est = t.elapsed();
        let (Coverage::Binned(c), Coverage::Estimated(e)) = (&counted.coverage, &est.coverage) else { panic!("{locus}: wrong kinds") };
        let c: Vec<f64> = c.iter().map(|b| b.mean).collect();
        let e: Vec<f64> = e.iter().map(|b| b.mean).collect();
        let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
        let (mc, me) = (mean(&c), mean(&e));
        let cov: f64 = c.iter().zip(&e).map(|(a, b)| (a - mc) * (b - me)).sum();
        let sd = |v: &[f64], m: f64| v.iter().map(|x| (x - m).powi(2)).sum::<f64>().sqrt();
        let r = cov / (sd(&c, mc) * sd(&e, me));
        eprintln!("{locus}: counted mean {mc:.1} ({t_count:.0?}), estimated mean {me:.1} ({t_est:.0?}), r = {r:.3}");
        for (a, b) in c.iter().zip(&e).step_by(5) {
            eprintln!("    {a:7.1} {b:7.1}");
        }
        assert!((me / mc - 1.0).abs() < 0.25, "{locus}: level off");
        assert!(r > 0.8, "{locus}: shape off (r = {r})");
    }
}

/// Per-bin error of the index estimate against counted depth, in 1 Mb bins, with and without a
/// MAPQ filter (the estimate once overshot counted depth by a wide margin). Bins with a counted
/// mean under 2 are skipped (a relative error there says little).
///
/// Unfiltered, nine bins in ten must be within 8%. With a MAPQ filter the typical bin must be
/// within 3% and nine in ten within 15%, but single bins in repeats can still be far off: the index can't see MAPQ, and
/// the share of low-MAPQ reads can change inside a bin. Those are printed, not asserted.
#[test]
fn index_estimate_per_bin_error() {
    let (Ok(bam_path), Ok(regions)) = (std::env::var("GGV_TEST_WGS_BAM"), std::env::var("GGV_TEST_WGS_REGIONS")) else { return };
    for min_mapq in [0u8, 21] {
        let mut errors = Vec::new();
        for locus in regions.split(',') {
            let mut bam = BamFile::open(&bam_path).unwrap();
            let region = bam.resolve_str(locus).unwrap();
            let bins = (region.len() / 1_000_000).max(1) as usize;
            let filters = ggv_query::PileupOptions { min_mapping_quality: min_mapq, ..Default::default() };
            let opts = ViewOptions { filters, bins, max_read_window: 0, max_coverage_window: u64::MAX, ..Default::default() };
            let counted = load_view(&mut bam, &region, None, String::new(), &opts).unwrap();
            let t = std::time::Instant::now();
            let est = load_view(&mut bam, &region, None, String::new(), &ViewOptions { max_coverage_window: 0, ..opts }).unwrap();
            let t_est = t.elapsed();
            let (Coverage::Binned(c), Coverage::Estimated(e)) = (&counted.coverage, &est.coverage) else { panic!("{locus}: wrong kinds") };
            for (k, (a, b)) in c.iter().zip(e).enumerate() {
                if a.mean >= 2.0 {
                    let err = b.mean / a.mean - 1.0;
                    errors.push(err.abs());
                    if err.abs() > 0.3 {
                        eprintln!("  q{min_mapq} {locus} bin {k}: counted {:.1} estimated {:.1} ({:+.0}%)", a.mean, b.mean, err * 100.0);
                    }
                }
            }
            eprintln!("  q{min_mapq} {locus}: estimate in {t_est:.0?}");
        }
        errors.sort_by(f64::total_cmp);
        let at = |q: f64| errors[((errors.len() - 1) as f64 * q) as usize] * 100.0;
        eprintln!("MAPQ >= {min_mapq}: {} bins, |error| median {:.1}%, 90th {:.1}%, worst {:.1}%", errors.len(), at(0.5), at(0.9), at(1.0));
        if min_mapq == 0 {
            assert!(at(0.9) < 8.0, "unfiltered: 90th-percentile bin error {:.1}%", at(0.9));
        } else {
            assert!(at(0.5) < 3.0, "MAPQ >= {min_mapq}: median bin error {:.1}%", at(0.5));
            assert!(at(0.9) < 15.0, "MAPQ >= {min_mapq}: 90th-percentile bin error {:.1}%", at(0.9));
        }
    }
}
