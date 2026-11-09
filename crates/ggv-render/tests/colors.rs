//! Read colours in headless pictures on real long reads (tests/data/hg002-revio): methylation
//! (MM/ML) and haplotype (HP), the same colours as the app.

use std::path::PathBuf;

use ggv_layout::{ViewOptions, load_view};
use ggv_render::{ColorBy, RenderOptions, build_scene};

fn svg(locus: &str, color_by: ColorBy) -> String {
    let data = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio");
    let mut bam = ggv_io::BamFile::open(data.join("HG002.revio.slices.bam")).unwrap();
    let mut fasta = ggv_io::ReferenceFasta::open(data.join("GRCh38.chr20.slices.fa.gz")).unwrap();
    let region = bam.resolve_str(locus).unwrap();
    let seq = fasta.fetch(&region).unwrap();
    let view = load_view(&mut bam, &region, Some(seq), "HG002".into(), &ViewOptions::default()).unwrap();
    build_scene(&view, &RenderOptions { color_by, ..Default::default() }).to_svg()
}

const METHYLATED: &str = "rgb(220,20,20)";
const UNMETHYLATED: &str = "rgb(20,20,220)";
const HP1: &str = "rgb(230,150,150)";
const HP2: &str = "rgb(150,150,230)";

#[test]
fn methylation_and_haplotype_colours() {
    // GNAS A/B: methylated on one haplotype, not on the other, so both colours appear.
    let locus = "chr20:58888000-58890000";
    let mods = svg(locus, ColorBy::BaseMods);
    assert!(mods.contains(METHYLATED) && mods.contains(UNMETHYLATED), "methylation colours missing");
    let hap = svg(locus, ColorBy::Haplotype);
    assert!(hap.contains(HP1) && hap.contains(HP2), "haplotype colours missing");
    let plain = svg(locus, ColorBy::None);
    for c in [METHYLATED, UNMETHYLATED, HP1, HP2] {
        assert!(!plain.contains(c), "{c} drawn without being asked for");
    }
}

#[test]
fn grouped_pictures_label_each_group() {
    let data = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio");
    let mut bam = ggv_io::BamFile::open(data.join("HG002.revio.slices.bam")).unwrap();
    let region = bam.resolve_str("chr20:31165000-31184000").unwrap();
    let mut view = load_view(&mut bam, &region, None, "HG002".into(), &ViewOptions::default()).unwrap();
    ggv_layout::repack(&mut view, 100, ggv_layout::GroupBy::Haplotype, None, false);
    let svg = build_scene(&view, &RenderOptions::default()).to_svg();
    // 38 reads on each haplotype here (samtools agrees), the rest untagged.
    for label in ["HP 1 (38 reads)", "HP 2 (38 reads)", "no HP ("] {
        assert!(svg.contains(label), "{label} missing");
    }
}
