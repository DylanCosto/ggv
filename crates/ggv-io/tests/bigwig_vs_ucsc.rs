//! GGV_TEST_BIGWIG + GGV_TEST_BIGWIG_JSON: a bigWig and UCSC's API answer for the same file over
//! whole chromosomes (getData/track ... maxItemsOutput=-1, saved as JSON; comma-separated for
//! several). Every raw interval must match, and bins over a whole chromosome (served from zoom
//! levels) must agree with bins computed from the raw values.

#[test]
fn bigwig_matches_ucsc() {
    let (Ok(bw), Ok(files)) = (std::env::var("GGV_TEST_BIGWIG"), std::env::var("GGV_TEST_BIGWIG_JSON")) else { return };
    let mut b = ggv_io::signal::BigWig::open(std::path::Path::new(&bw)).unwrap();
    for f in files.split(',') {
        let j: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(f).unwrap()).unwrap();
        let (chrom, end) = (j["chrom"].as_str().unwrap(), j["end"].as_u64().unwrap());
        let track = j["track"].as_str().unwrap();
        let want: Vec<(u64, u64, f64)> = j[track].as_array().unwrap().iter().map(|i| (i["start"].as_u64().unwrap(), i["end"].as_u64().unwrap(), i["value"].as_f64().unwrap())).collect();
        let got: Vec<(u64, u64, f64)> = b.intervals(chrom, 0, end).unwrap().iter().map(|i| (i.start, i.end, i.value)).collect();
        assert_eq!(got.len(), want.len(), "{chrom}: interval count");
        for (g, w) in got.iter().zip(&want) {
            assert!(g.0 == w.0 && g.1 == w.1 && (g.2 - w.2).abs() < 1e-4, "{chrom}: {g:?} vs {w:?}");
        }
        // Whole chromosome in 100 bins: zoom levels against the raw data, bin by bin.
        let zoomed = b.bins(chrom, 0, end, 100).unwrap();
        let mut raw = vec![(0.0f64, 0u64); 100];
        for &(s, e, v) in &want {
            for k in 0..100u64 {
                let (bs, be) = (end * k / 100, end * (k + 1) / 100);
                let o = e.min(be).saturating_sub(s.max(bs));
                raw[k as usize].0 += v * o as f64;
                raw[k as usize].1 += o;
            }
        }
        let mut worst: f64 = 0.0;
        for (z, (sum, n)) in zoomed.iter().zip(&raw) {
            if *n > 0 {
                let want_mean = sum / *n as f64;
                worst = worst.max((z.unwrap().mean - want_mean).abs());
            }
        }
        // Zoom summaries straddling bin edges are split by overlap, as if even: allow 3% of the
        // track's range.
        let range = want.iter().map(|w| w.2).fold(f64::NEG_INFINITY, f64::max) - want.iter().map(|w| w.2).fold(f64::INFINITY, f64::min);
        eprintln!("{chrom}: {} intervals identical; whole-chromosome bins from zoom levels within {worst:.3} of the raw means (range {range})", got.len());
        assert!(worst < 0.03 * range, "{chrom}: zoom bins off by {worst}");
    }
}
