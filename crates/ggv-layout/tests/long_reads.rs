//! Long-read features on real PacBio Revio HiFi reads of HG002 (tests/data/hg002-revio, see its
//! README): HP haplotype tags, MM/ML methylation calls, split reads (SA) and a heterozygous
//! 6.8 kb deletion. The fixture is in the repository, so these run everywhere, CI included.

use std::path::PathBuf;

use ggv_layout::{AlignedRead, ViewOptions, load_view, phase::phase};
use ggv_model::Region;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio").join(name)
}

/// Reads over `locus`, decoded against the reference, with the region and its bases.
fn reads(locus: &str) -> (Vec<AlignedRead>, Region, Vec<u8>) {
    let mut bam = ggv_io::BamFile::open(data("HG002.revio.slices.bam")).unwrap();
    let mut fasta = ggv_io::ReferenceFasta::open(data("GRCh38.chr20.slices.fa.gz")).unwrap();
    let region = bam.resolve_str(locus).unwrap();
    let seq = fasta.fetch(&region).unwrap();
    let opts = ViewOptions { max_read_window: 200_000, ..Default::default() };
    let view = load_view(&mut bam, &region, Some(seq.clone()), String::new(), &opts).unwrap();
    (view.reads.unwrap().0.to_vec(), region, seq)
}

/// pbsv.DEL.567, 0/1: chr20:31,170,827-31,177,629 (6,802 bp).
const DEL: (u64, u64) = (31_170_827, 31_177_629);

fn carries_deletion(r: &AlignedRead) -> bool {
    r.deletions.iter().any(|d| d.len >= 6000 && d.pos < DEL.1 && d.pos + u64::from(d.len) > DEL.0)
}

#[test]
fn deletion_is_on_one_haplotype() {
    let (rs, _, _) = reads("chr20:31160000-31190000");
    let carriers: Vec<&AlignedRead> = rs.iter().filter(|r| carries_deletion(r)).collect();
    let intact: Vec<&AlignedRead> = rs.iter().filter(|r| r.start < DEL.0 - 500 && r.end > DEL.1 + 500 && !carries_deletion(r)).collect();
    assert!(carriers.len() >= 10, "{} reads carry the deletion", carriers.len());
    assert!(carriers.iter().all(|r| r.haplotype == Some(1)), "deletion reads: {:?}", carriers.iter().map(|r| r.haplotype).collect::<Vec<_>>());
    let tagged: Vec<_> = intact.iter().filter_map(|r| r.haplotype).collect();
    assert!(tagged.len() >= 8 && tagged.iter().all(|&h| h == 2), "intact reads' HP: {tagged:?}");
}

#[test]
fn split_reads_list_their_other_parts() {
    let (rs, _, _) = reads("chr20:31160000-31190000");
    let supplementary: Vec<&AlignedRead> = rs.iter().filter(|r| r.flags & 0x800 != 0).collect();
    assert!(supplementary.len() >= 50, "{} supplementary records", supplementary.len());
    for r in &supplementary {
        assert!(!r.supplementary.is_empty(), "{}: supplementary record without SA parts", r.name);
        assert!(r.supplementary.iter().all(|s| s.mapq <= 60 && !s.cigar.is_empty()), "{}: {:?}", r.name, r.supplementary);
    }
}

#[test]
fn allele_phasing_matches_hp_tags() {
    // GGV's phasing from the reads' own alleles against PacBio's haplotagging (WhatsHap), up to
    // which haplotype is called first.
    let (rs, _, _) = reads("chr20:58830000-58930000");
    let members: Vec<usize> = (0..rs.len()).collect();
    let p = phase(&rs, &members);
    let (mut same, mut swapped) = (0, 0);
    for (r, a) in rs.iter().zip(&p.hap) {
        if let (Some(h), Some(a)) = (r.haplotype, a) {
            if (h == 1) == *a { same += 1 } else { swapped += 1 }
        }
    }
    let both = same + swapped;
    assert!(both >= 40, "only {both} reads both tagged and phased");
    assert!(same.max(swapped) as f64 >= 0.98 * both as f64, "agree {same} / disagree {swapped}");
}

#[test]
fn methylation_calls_sit_on_cpgs() {
    // 5mC calls (MM `C+m`, ML) are made at CpGs in the read; mapped to the reference they must
    // land on a C followed by G (forward reads) or on the G of a CpG (reverse reads).
    let (rs, region, seq) = reads("chr20:58800000-58950000");
    let at = |p: u64| (p >= region.start && p < region.end).then(|| seq[(p - region.start) as usize]);
    let (mut cpg, mut all) = (0, 0);
    for r in &rs {
        for m in r.mods.iter().filter(|m| m.pos(r.start) > region.start && m.pos(r.start) + 1 < region.end) {
            assert_eq!((m.code, m.canonical), (b'm', b'C'));
            all += 1;
            let pos = m.pos(r.start);
            let on_cpg = if r.reverse { at(pos) == Some(b'G') && at(pos - 1) == Some(b'C') } else { at(pos) == Some(b'C') && at(pos + 1) == Some(b'G') };
            cpg += usize::from(on_cpg);
        }
    }
    assert!(all > 50_000, "{all} calls");
    assert!(cpg as f64 >= 0.97 * all as f64, "{cpg} of {all} calls on CpGs");
}

#[test]
fn gnas_ab_region_is_methylated_on_one_haplotype() {
    // Imprinting: at GNAS A/B (chr20:58,888,000-58,890,000) one parental copy is methylated and
    // the other not, so methylation splits by haplotype.
    let (rs, _, _) = reads("chr20:58880000-58900000");
    let mut by = [(0u32, 0u32); 2];
    for r in &rs {
        let Some(h @ 1..=2) = r.haplotype else { continue };
        for m in r.mods.iter().filter(|m| (58_888_000..58_890_000).contains(&m.pos(r.start))) {
            let e = &mut by[usize::from(h - 1)];
            e.1 += 1;
            e.0 += u32::from(m.prob >= 128);
        }
    }
    let frac = by.map(|(m, n)| { assert!(n >= 500, "{n} calls"); f64::from(m) / f64::from(n) });
    let (lo, hi) = (frac[0].min(frac[1]), frac[0].max(frac[1]));
    assert!(lo < 0.3 && hi > 0.75, "methylated fraction by haplotype: {frac:?}");
}

#[test]
fn any_tag_groups_like_its_dedicated_grouping() {
    // Loaded with the HP tag, grouping by `tag HP` gives the same groups (same reads) as
    // grouping by haplotype; PS values come through as the phase-set numbers.
    let mut bam = ggv_io::BamFile::open(data("HG002.revio.slices.bam")).unwrap();
    let region = bam.resolve_str("chr20:31165000-31184000").unwrap();
    let opts = ViewOptions { tag: Some(*b"HP"), ..Default::default() };
    let view = load_view(&mut bam, &region, None, String::new(), &opts).unwrap();
    let (reads, _) = view.reads.as_ref().unwrap();
    let pack = |group| ggv_layout::pack_with(reads, &ggv_layout::PackOptions { max_rows: 1000, min_gap: 2, group, sort: None, link: false });
    let (by_tag, by_hp) = (pack(ggv_layout::GroupBy::Tag(*b"HP")), pack(ggv_layout::GroupBy::Haplotype));
    let counts = |p: &ggv_layout::Packed| p.groups.iter().map(|g| g.reads).collect::<Vec<_>>();
    assert_eq!(counts(&by_tag), counts(&by_hp));
    assert_eq!(by_tag.groups.iter().map(|g| g.label.as_str()).collect::<Vec<_>>(), ["HP: 1", "HP: 2", "no HP"]);
    assert!(reads.iter().all(|r| r.tag.as_deref() == r.haplotype.map(|h| h.to_string()).as_deref()));

    let opts = ViewOptions { tag: Some(*b"PS"), ..Default::default() };
    let gnas = bam.resolve_str("chr20:58880000-58900000").unwrap();
    let view = load_view(&mut bam, &gnas, None, String::new(), &opts).unwrap();
    let ps: std::collections::BTreeSet<_> = view.reads.as_ref().unwrap().0.iter().filter_map(|r| r.tag.as_deref().map(str::to_owned)).collect();
    assert!(ps.contains("58891317"), "{ps:?}");
}

/// Two loads as the app makes them at an 89 kb view (100 kb spans, the read window), before and
/// after panning 12 kb: how many reads shown in both views changed rows.
fn rows_changed_after_pan(sticky: bool) -> (usize, usize) {
    let mut bam = ggv_io::BamFile::open(data("HG002.revio.slices.bam")).unwrap();
    let opts = ViewOptions::default();
    let load = |bam: &mut ggv_io::BamFile, s: u64| {
        let region = bam.resolve_str(&format!("chr20:{}-{}", s + 1, s + 100_000)).unwrap();
        load_view(bam, &region, None, String::new(), &opts).unwrap()
    };
    let (a, b) = (load(&mut bam, 58_811_000), load(&mut bam, 58_823_000));
    let (ra, rb) = (&a.reads.as_ref().unwrap().0, &b.reads.as_ref().unwrap().0);
    let po = ggv_layout::PackOptions { max_rows: 10_000, min_gap: 2, group: ggv_layout::GroupBy::None, sort: None, link: false };
    let pa = ggv_layout::pack_with(ra, &po);
    let pb = if sticky {
        ggv_layout::pack_sticky(rb, (0..rb.len()).collect(), &po, &ggv_layout::row_hints(ra, &pa))
    } else {
        ggv_layout::pack_with(rb, &po)
    };
    let row_of = |p: &ggv_layout::Packed| -> std::collections::HashMap<usize, usize> { p.rows.iter().enumerate().flat_map(|(r, row)| row.iter().map(move |&i| (i, r))).collect() };
    let (rows_a, rows_b) = (row_of(&pa), row_of(&pb));
    // Reads on screen in both views (the overlap of two 89 kb views 12 kb apart).
    let (vs, ve) = (58_828_500u64, 58_905_500u64);
    let mut both = 0;
    let mut moved = 0;
    for (i, r) in ra.iter().enumerate() {
        if r.end <= vs || r.start >= ve {
            continue;
        }
        if let Some(j) = rb.iter().position(|x| x.name == r.name && x.start == r.start && x.flags == r.flags) {
            both += 1;
            moved += usize::from(rows_a[&i] != rows_b[&j]);
        }
    }
    (moved, both)
}

#[test]
fn panning_keeps_reads_in_their_rows() {
    let (moved, both) = rows_changed_after_pan(false);
    eprintln!("plain packing: {moved} of {both} reads on screen changed rows after a 12 kb pan");
    let (moved, both) = rows_changed_after_pan(true);
    eprintln!("sticky packing: {moved} of {both} changed rows");
    assert_eq!(moved, 0, "reads on screen must keep their rows");
}
