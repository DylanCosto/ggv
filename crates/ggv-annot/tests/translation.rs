//! GGV_TEST_REF + GGV_TEST_ANNOTATIONS: every coding transcript whose coding bases are all known
//! in the reference must translate, through `codon_positions`, to a protein that starts with Met
//! (mostly), ends with a stop, and has no stop in between. That exercises the coding-index
//! mapping on both strands, across every exon boundary, against real annotations.

use ggv_annot::{Annotation, Strand, coding};
use ggv_model::Region;

#[test]
fn coding_sequences_translate_cleanly() {
    let (Ok(reference), Ok(annotations)) = (std::env::var("GGV_TEST_REF"), std::env::var("GGV_TEST_ANNOTATIONS")) else { return };
    let mut fasta = ggv_io::ReferenceFasta::open(&reference).unwrap();
    for path in annotations.split(',') {
        let a = Annotation::read(path).unwrap();
        let (mut checked, mut met, mut by_strand) = (0, 0, [0, 0]);
        for contig in a.contig_names().to_vec() {
            let Ok(r) = ggv_io::resolve_locus_str(fasta.contigs(), &contig) else { continue };
            for t in a.query(&contig, 0, u64::MAX).filter(|t| t.cds.is_some() && t.biotype.contains("protein_coding")) {
                let Ok(seq) = fasta.fetch(&Region { contig: r.contig.clone(), start: t.start, end: t.end.min(r.end) }) else { continue };
                let base = |p: u64| seq.get((p - t.start) as usize).copied();
                let len = t.cds_len();
                if len % 3 != 0 || len < 6 {
                    continue; // incomplete annotations (cds_start_NF / end_NF)
                }
                let mut protein = Vec::new();
                for k in 0..len / 3 {
                    let pos = coding::codon_positions(t, k).expect("codon inside the CDS");
                    let mut c = pos.map(|p| base(p).unwrap_or(b'N'));
                    if t.strand == Strand::Reverse {
                        c = c.map(coding::complement);
                    }
                    protein.push(coding::amino_acid(t, k, c));
                }
                if protein.contains(&b'X') {
                    continue; // masked reference
                }
                let text = String::from_utf8_lossy(&protein).into_owned();
                assert_eq!(protein.last(), Some(&b'*'), "{} ({:?}) does not end in a stop: {text}", t.id, t.strand);
                assert!(!protein[..protein.len() - 1].contains(&b'*'), "{} ({:?}) has an internal stop: {text}", t.id, t.strand);
                met += usize::from(protein[0] == b'M');
                by_strand[usize::from(t.strand == Strand::Reverse)] += 1;
                checked += 1;
            }
        }
        eprintln!("{path}: {checked} complete coding sequences translate cleanly ({} +, {} -), {met} start with Met", by_strand[0], by_strand[1]);
        assert!(checked > 0, "{path}: nothing to check (all masked?)");
        assert!(met * 10 >= checked * 9, "most should start with Met: {met} of {checked}");
    }
}
