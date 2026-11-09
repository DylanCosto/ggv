//! GGV_TEST_2BIT + GGV_TEST_2BIT_FASTA: a UCSC .2bit and the same genome as indexed FASTA (e.g.
//! sacCer3.2bit and sacCer3.fa from hgdownload) must give the same bases everywhere: every
//! contig whole, and a thousand random spans.

use ggv_model::Region;

#[test]
fn twobit_matches_fasta() {
    let (Ok(two), Ok(fa)) = (std::env::var("GGV_TEST_2BIT"), std::env::var("GGV_TEST_2BIT_FASTA")) else { return };
    let mut t = ggv_io::ReferenceFasta::open(&two).unwrap();
    let mut f = ggv_io::ReferenceFasta::open(&fa).unwrap();
    let contigs = t.contigs().to_vec();
    assert_eq!(contigs, f.contigs().to_vec(), "contig names and lengths");
    let mut bases = 0u64;
    for c in &contigs {
        let r = Region { contig: c.name.clone(), start: 0, end: c.length };
        assert_eq!(t.fetch(&r).unwrap(), f.fetch(&r).unwrap(), "{}", c.name);
        bases += c.length;
    }
    let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
    for _ in 0..1000 {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        let c = &contigs[(x % contigs.len() as u64) as usize];
        let start = (x >> 8) % c.length;
        let end = (start + 1 + (x >> 40) % 5000).min(c.length);
        let r = Region { contig: c.name.clone(), start, end };
        assert_eq!(t.fetch(&r).unwrap(), f.fetch(&r).unwrap(), "{}:{start}-{end}", c.name);
    }
    eprintln!("{} contigs, {bases} bases identical", contigs.len());
}
