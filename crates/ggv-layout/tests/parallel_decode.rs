//! Reads decoded on several threads (wide spans) are the same reads, in the same order, with
//! the same counts and coverage, as one pass.
//!
//! Needs real data, so it runs only when these are set (otherwise it passes trivially):
//!   GGV_TEST_BAM     indexed BAM
//!   GGV_TEST_REF     indexed FASTA (optional)
//!   GGV_TEST_WIDE    a locus of a megabase or more, e.g. `chr20:30000000-32000000`

use ggv_io::{BamFile, ReferenceFasta};
use ggv_layout::{Coverage, ViewOptions, load_view};

#[test]
fn chunked_decoding_matches_one_pass() {
    let (Ok(bam_path), Ok(locus)) = (std::env::var("GGV_TEST_BAM"), std::env::var("GGV_TEST_WIDE")) else {
        eprintln!("GGV_TEST_BAM / GGV_TEST_WIDE not set; skipping");
        return;
    };
    let mut bam = BamFile::open_with_reference(&bam_path, std::env::var("GGV_TEST_REF").ok()).unwrap();
    let region = bam.resolve_str(&locus).unwrap();
    let reference = std::env::var("GGV_TEST_REF").ok().map(|p| ReferenceFasta::open(p).unwrap().fetch(&region).unwrap());
    let opts = ViewOptions { max_read_window: u64::MAX, max_rows: 0, ..Default::default() };
    let chunked = load_view(&mut bam, &region, reference.clone(), String::new(), &ViewOptions { read_chunk: 100_000, ..opts.clone() }).unwrap();
    let single = load_view(&mut bam, &region, reference, String::new(), &ViewOptions { read_chunk: u64::MAX, ..opts }).unwrap();
    let (a, b) = (&chunked.reads.as_ref().unwrap().0, &single.reads.as_ref().unwrap().0);
    assert!(!b.is_empty(), "{locus}: no reads to compare");
    assert_eq!(a.len(), b.len(), "{locus}: read count");
    for (x, y) in a.iter().zip(b) {
        assert_eq!((x.start, &x.name, x.flags, x.mismatches.len()), (y.start, &y.name, y.flags, y.mismatches.len()), "{locus}: read order and content");
    }
    assert_eq!(chunked.read_counts, single.read_counts, "{locus}: counts");
    let (Coverage::PerBase(ca), Coverage::PerBase(cb)) = (&chunked.coverage, &single.coverage) else { panic!("per-base coverage") };
    assert!(ca == cb, "{locus}: coverage");
    eprintln!("{locus}: {} reads identical", a.len());
}
