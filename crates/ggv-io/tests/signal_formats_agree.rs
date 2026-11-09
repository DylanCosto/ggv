//! The same values as bigWig and bedGraph bin identically, including windows that run
//! past the contig end (chrQA is 1,000 bases; the desktop pads its requests by a view width).

use std::path::PathBuf;

use ggv_io::signal::Signal;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/fixtures").join(name)
}

#[test]
fn bigwig_matches_bedgraph_at_the_contig_end() {
    let (bw, bg) = (fixture("signal.bw"), fixture("signal.bedGraph"));
    assert!(bw.exists() && bg.exists(), "{} or {} is missing", bw.display(), bg.display());
    let mut bw = Signal::open(&bw).unwrap();
    let mut bg = Signal::open(&bg).unwrap();
    for (start, end, n) in [(0, 1000, 100), (0, 3000, 300), (200, 1400, 120), (600, 1000, 40), (200, 1400, 1200), (0, 1000, 7)] {
        let a = bw.bins(0, "chrQA", start, end, n).unwrap();
        let b = bg.bins(0, "chrQA", start, end, n).unwrap();
        assert_eq!(a.len(), n, "{start}-{end}");
        for (i, (x, y)) in a.iter().zip(&b).enumerate() {
            match (x, y) {
                (Some(x), Some(y)) => {
                    assert!((x.mean - y.mean).abs() < 1e-6, "{start}-{end}/{n} bin {i}: bigWig {} bedGraph {}", x.mean, y.mean);
                    assert!((x.max - y.max).abs() < 1e-6 && (x.min - y.min).abs() < 1e-6, "{start}-{end}/{n} bin {i}");
                }
                (None, None) => {}
                _ => panic!("{start}-{end}/{n} bin {i}: bigWig {x:?} bedGraph {y:?}"),
            }
        }
    }
    // chrQA:601-1000 is +4 throughout; past the end there is nothing.
    let bins = bw.bins(0, "chrQA", 600, 1400, 80).unwrap();
    assert!(bins[..40].iter().all(|b| b.is_some_and(|b| b.mean == 4.0)));
    assert!(bins[40..].iter().all(Option::is_none));
}
