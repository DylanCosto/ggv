use super::*;

const GTF: &str = "\
##description: test
chr1\tHAVANA\tgene\t1001\t5000\t.\t+\t.\tgene_id \"ENSG01.3\"; gene_name \"ALPHA\";
chr1\tHAVANA\ttranscript\t1001\t5000\t.\t+\t.\tgene_id \"ENSG01.3\"; transcript_id \"ENST01.1\"; gene_name \"ALPHA\"; transcript_type \"protein_coding\"; tag \"basic\";
chr1\tHAVANA\texon\t1001\t1200\t.\t+\t.\tgene_id \"ENSG01.3\"; transcript_id \"ENST01.1\";
chr1\tHAVANA\texon\t4001\t5000\t.\t+\t.\tgene_id \"ENSG01.3\"; transcript_id \"ENST01.1\";
chr1\tHAVANA\tCDS\t1101\t1200\t.\t+\t0\tgene_id \"ENSG01.3\"; transcript_id \"ENST01.1\";
chr1\tHAVANA\tCDS\t4001\t4497\t.\t+\t0\tgene_id \"ENSG01.3\"; transcript_id \"ENST01.1\";
chr1\tHAVANA\tstop_codon\t4498\t4500\t.\t+\t0\tgene_id \"ENSG01.3\"; transcript_id \"ENST01.1\";
chr1\tHAVANA\tSelenocysteine\t4101\t4103\t.\t+\t.\tgene_id \"ENSG01.3\"; transcript_id \"ENST01.1\";
chr1\tHAVANA\ttranscript\t1001\t3000\t.\t+\t.\tgene_id \"ENSG01.3\"; transcript_id \"ENST02.1\"; gene_name \"ALPHA\"; tag \"MANE_Select\";
chr1\tHAVANA\texon\t1001\t1100\t.\t+\t.\tgene_id \"ENSG01.3\"; transcript_id \"ENST02.1\";
chr1\tHAVANA\texon\t2901\t3000\t.\t+\t.\tgene_id \"ENSG01.3\"; transcript_id \"ENST02.1\";
chr1\tHAVANA\texon\t20001\t20500\t.\t-\t.\tgene_id \"ENSG02.1\"; transcript_id \"ENST03.1\"; gene_name \"BETA\"; transcript_type \"lncRNA\";
chr2\tHAVANA\texon\t11\t20\t.\t-\t.\tgene_id \"ENSG03.1\"; transcript_id \"ENST04.1\"; gene_name \"GAMMA\";
";

const GFF3: &str = "\
##gff-version 3
chr20\tensembl\tgene\t142590\t145751\t.\t+\t.\tID=gene:ENSG00000125788;Name=DEFB126;biotype=protein_coding;gene_id=ENSG00000125788
chr20\tensembl\tmRNA\t142590\t145751\t.\t+\t.\tID=transcript:ENST00000382398;Parent=gene:ENSG00000125788;Name=DEFB126-201;biotype=protein_coding;tag=basic,Ensembl_canonical,MANE_Select;transcript_id=ENST00000382398
chr20\tensembl\texon\t142590\t142686\t.\t+\t.\tParent=transcript:ENST00000382398;Name=ENSE1
chr20\tensembl\texon\t145415\t145751\t.\t+\t.\tParent=transcript:ENST00000382398;Name=ENSE2
chr20\tensembl\tCDS\t142629\t142686\t.\t+\t0\tID=CDS:ENSP1;Parent=transcript:ENST00000382398
chr20\tensembl\tCDS\t145415\t145692\t.\t+\t2\tID=CDS:ENSP1;Parent=transcript:ENST00000382398
chr20\tensembl\tgene\t200000\t201000\t.\t-\t.\tID=gene:X1;Name=A%3BB;biotype=misc
chr20\tensembl\texon\t200000\t201000\t.\t-\t.\tParent=gene:X1
##FASTA
>chr20
ACGT
";

const BED: &str = "\
track name=test
chr3\t1000\t5000\tTX1\t0\t+\t1500\t4500\t0\t3\t100,200,300,\t0,1000,3700,
chr3\t7000\t8000\tPLAIN
";

fn gtf() -> Annotation {
    Annotation::parse(GTF.as_bytes(), Format::Gtf).unwrap()
}

#[test]
fn gtf_builds_transcripts_with_zero_based_coordinates() {
    let a = gtf();
    assert_eq!(a.len(), 4);
    assert_eq!(a.gene_count(), 3);
    let t = a.query("chr1", 0, 10_000).find(|t| &*t.id == "ENST01.1").unwrap();
    assert_eq!((t.start, t.end), (1000, 5000));
    assert_eq!(t.exons, vec![(1000, 1200), (4000, 5000)]);
    // The stop codon counts as coding (Ensembl-style GTFs leave it out of CDS).
    assert_eq!(t.cds, Some((1100, 4500)));
    assert_eq!(t.translation, Translation { phase: 0, table: 1, selenocysteine: vec![(4100, 4103)] });
    assert!(!t.exons_from_span);
    assert_eq!(t.strand, Strand::Forward);
    assert_eq!(&*t.gene_name, "ALPHA");
    assert_eq!(&*t.biotype, "protein_coding");
    assert_eq!(t.tags, TAG_BASIC);
}

#[test]
fn mane_select_beats_longer_cds() {
    let a = gtf();
    let reps: Vec<&str> = a.query("chr1", 0, 10_000).filter(|t| t.representative).map(|t| &*t.id).collect();
    assert_eq!(reps, vec!["ENST02.1"]);
}

#[test]
fn query_finds_overlaps_only_and_handles_long_transcripts() {
    let a = gtf();
    // 3500 is inside ENST01 (1000-5000) but past ENST02 (1000-3000).
    let ids: Vec<&str> = a.query("chr1", 3500, 3600).map(|t| &*t.id).collect();
    assert_eq!(ids, vec!["ENST01.1"]);
    assert_eq!(a.query("chr1", 5000, 20_000).count(), 0, "half-open ends");
    assert_eq!(a.query("chr1", 20_499, 20_500).count(), 1);
    // Aliases: `1` finds `chr1`.
    assert_eq!(a.query("1", 0, 10_000).count(), 2);
    assert_eq!(a.query("chrX", 0, 10_000).count(), 0);
}

#[test]
fn find_by_name_or_id_case_insensitive_without_version() {
    let a = gtf();
    let g = a.find("alpha").unwrap();
    assert_eq!((g.contig.as_str(), g.start, g.end), ("chr1", 1000, 5000));
    assert_eq!(a.find("ENSG01").unwrap().name, "ALPHA");
    assert_eq!(a.find("enst03.1").unwrap().name, "BETA");
    assert!(a.find("DELTA").is_none());
}

#[test]
fn gff3_resolves_parents_names_and_tags() {
    let a = Annotation::parse(GFF3.as_bytes(), Format::Gff3).unwrap();
    assert_eq!(a.len(), 2);
    let t = a.query("chr20", 142_000, 146_000).next().unwrap();
    assert_eq!(&*t.id, "ENST00000382398");
    assert_eq!(&*t.gene_name, "DEFB126");
    assert_eq!(&*t.gene_id, "ENSG00000125788");
    assert_eq!(t.exons, vec![(142_589, 142_686), (145_414, 145_751)]);
    assert_eq!(t.cds, Some((142_628, 145_692)));
    assert_eq!(t.tags, TAG_BASIC | TAG_CANONICAL | TAG_MANE_SELECT);
    // Exons straight under a gene; percent-encoded name.
    let g = a.query("chr20", 200_000, 200_001).next().unwrap();
    assert_eq!(&*g.gene_name, "A;B");
    assert_eq!(g.strand, Strand::Reverse);
}

#[test]
fn bed12_blocks_and_thick_range() {
    let a = Annotation::parse(BED.as_bytes(), Format::Bed).unwrap();
    assert_eq!(a.len(), 2);
    let t = a.query("chr3", 0, 6000).next().unwrap();
    assert_eq!(t.exons, vec![(1000, 1100), (2000, 2200), (4700, 5000)]);
    assert_eq!(t.cds, Some((1500, 4500)));
    let plain = a.query("chr3", 7500, 7501).next().unwrap();
    assert_eq!(plain.exons, vec![(7000, 8000)]);
    assert_eq!(plain.cds, None);
    assert_eq!(a.find("tx1").unwrap().start, 1000);
}

#[test]
fn contigs_for_browsing_without_a_bam() {
    let a = gtf();
    let c = a.contigs();
    assert_eq!(c.iter().map(|c| (c.name.as_str(), c.length)).collect::<Vec<_>>(), vec![("chr1", 20_500), ("chr2", 20)]);
}

#[test]
fn reads_real_files_when_given() {
    // GGV_TEST_ANNOTATIONS=a.gtf,b.gff3.gz checks that real files parse into something sane.
    let Ok(paths) = std::env::var("GGV_TEST_ANNOTATIONS") else { return };
    for p in paths.split(',') {
        let t = std::time::Instant::now();
        let a = Annotation::read(p).unwrap();
        eprintln!("{p}: {} transcripts, {} genes in {:.2?}", a.len(), a.gene_count(), t.elapsed());
        assert!(a.len() > 0);
    }
}

#[test]
fn cache_round_trip_matches_the_parse() {
    let dir = std::env::temp_dir().join(format!("ggv-annot-cache-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let source = dir.join("genes.gtf");
    std::fs::write(&source, GTF).unwrap();
    let cache_file = dir.join("genes.ggvannot");
    let key = cache::SourceKey::of(&source).unwrap();
    let parsed = Annotation::read(&source).unwrap();
    cache::write(&parsed, &key, &cache_file).unwrap();
    let cached = cache::open(&source, &cache_file).expect("a fresh cache opens");
    assert!(cached.is_cached() && !parsed.is_cached());
    assert!(cached.parts.iter().all(|p| p.loaded.get().is_none()), "nothing decoded until asked");

    assert_eq!((cached.len(), cached.gene_count()), (parsed.len(), parsed.gene_count()));
    assert_eq!(cached.contigs(), parsed.contigs());
    assert_eq!(cached.genes(), parsed.genes());
    for (contig, s, e) in [("chr1", 0, 10_000), ("chr1", 2_950, 2_960), ("1", 19_000, 21_000), ("chr2", 0, 100), ("chr9", 0, 100)] {
        let a: Vec<_> = parsed.query(contig, s, e).collect();
        let b: Vec<_> = cached.query(contig, s, e).collect();
        assert_eq!(a, b, "{contig}:{s}-{e}");
    }
    assert!(cached.parts[0].loaded.get().is_some() && cached.contig_index("chr2").is_some());
    for name in ["alpha", "ENSG01", "ENSG01.3", "ENST02", "enst03.1", "GAMMA", "nope"] {
        assert_eq!(cached.find(name), parsed.find(name), "{name}");
    }

    // A changed source no longer matches its cache.
    std::fs::write(&source, format!("{GTF}chr3\tX\texon\t1\t10\t.\t+\t.\tgene_id \"G9\"; transcript_id \"T9\";\n")).unwrap();
    assert!(cache::open(&source, &cache_file).is_none());
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn a_shared_intron_keeps_every_gene() {
    // Two genes with the same exons: the intron 125-225 (0-based) belongs to both.
    let gtf = "chrQA\tX\texon\t101\t125\t.\t+\t.\tgene_id \"GA\"; transcript_id \"TA\"; gene_name \"GENE_A\";\n\
               chrQA\tX\texon\t226\t250\t.\t+\t.\tgene_id \"GA\"; transcript_id \"TA\"; gene_name \"GENE_A\";\n\
               chrQA\tX\texon\t101\t125\t.\t+\t.\tgene_id \"GB\"; transcript_id \"TB\"; gene_name \"GENE_B\";\n\
               chrQA\tX\texon\t226\t250\t.\t+\t.\tgene_id \"GB\"; transcript_id \"TB\"; gene_name \"GENE_B\";\n";
    let a = Annotation::parse(gtf.as_bytes(), Format::Gtf).unwrap();
    let mut genes: Vec<String> = a.introns("chrQA", 0, 1000).iter().map(|i| (i.start, i.end, i.gene.to_string())).inspect(|i| assert_eq!((i.0, i.1), (125, 225))).map(|i| i.2).collect();
    genes.sort();
    assert_eq!(genes, ["GENE_A", "GENE_B"]);
}

#[test]
fn coding_frame_genetic_code_and_selenocysteine() {
    const GTF: &str = "\
chrM\tqa\texon\t1\t9\t.\t+\t.\tgene_id \"MITO\"; transcript_id \"MITO_tx\";
chrM\tqa\tCDS\t1\t9\t.\t+\t0\tgene_id \"MITO\"; transcript_id \"MITO_tx\";
chrP\tqa\texon\t1\t10\t.\t+\t.\tgene_id \"P\"; transcript_id \"P_tx\";
chrP\tqa\tCDS\t1\t4\t.\t+\t1\tgene_id \"P\"; transcript_id \"P_tx\";
chrP\tqa\tCDS\t7\t10\t.\t+\t0\tgene_id \"P\"; transcript_id \"P_tx\";
chrP\tqa\texon\t21\t30\t.\t-\t.\tgene_id \"R\"; transcript_id \"R_tx\";
chrP\tqa\tCDS\t21\t24\t.\t-\t0\tgene_id \"R\"; transcript_id \"R_tx\";
chrP\tqa\tCDS\t27\t30\t.\t-\t2\tgene_id \"R\"; transcript_id \"R_tx\";
";
    let a = Annotation::parse(GTF.as_bytes(), Format::Gtf).unwrap();
    let get = |id: &str| a.query("chrM", 0, 100).chain(a.query("chrP", 0, 100)).find(|t| &*t.id == id).unwrap().clone();
    assert_eq!(get("MITO_tx").translation.table, 2, "mitochondrial code on chrM");
    assert_eq!(get("P_tx").translation, Translation { phase: 1, table: 1, selenocysteine: vec![] });
    // On the minus strand the 5'-most CDS row is the highest one.
    assert_eq!(get("R_tx").translation.phase, 2);
    assert!(is_mitochondrial("MT") && is_mitochondrial("chrM") && !is_mitochondrial("chrMX"));

    // GFF3: phase, and selenocysteine as a transl_except on the CDS.
    const GFF: &str = "\
##gff-version 3
chr1\tqa\tgene\t1\t12\t.\t+\t.\tID=g1;Name=SEL
chr1\tqa\tmRNA\t1\t12\t.\t+\t.\tID=t1;Parent=g1
chr1\tqa\texon\t1\t12\t.\t+\t.\tParent=t1
chr1\tqa\tCDS\t1\t12\t.\t+\t1\tParent=t1;transl_except=(pos:4..6%2Caa:Sec)
";
    let a = Annotation::parse(GFF.as_bytes(), Format::Gff3).unwrap();
    let t = a.query("chr1", 0, 12).next().unwrap();
    assert_eq!(t.translation, Translation { phase: 1, table: 1, selenocysteine: vec![(3, 6)] });
}

#[test]
fn gff3_transcripts_without_exons_are_kept_as_spans() {
    const GFF: &str = "\
##gff-version 3
chrQA\tqa\tgene\t101\t150\t.\t+\t.\tID=GENE_A;Name=GENE_A
chrQA\tqa\tmRNA\t101\t150\t.\t+\t.\tID=TX_A;Parent=GENE_A
chrQA\tqa\tfive_prime_UTR\t101\t110\t.\t+\t.\tID=U_A;Parent=TX_A
chrQA\tqa\tgene\t201\t250\t.\t+\t.\tID=GENE_B;Name=GENE_B
chrQA\tqa\tmRNA\t201\t250\t.\t+\t.\tID=TX_B;Parent=GENE_B
chrQA\tqa\texon\t201\t250\t.\t+\t.\tID=EX_B;Parent=TX_B
chrQA\tqa\tregion\t1\t1000\t.\t.\t.\tID=chrQA_region
";
    let a = Annotation::parse(GFF.as_bytes(), Format::Gff3).unwrap();
    let ts: Vec<_> = a.query("chrQA", 100, 250).collect();
    assert_eq!(ts.iter().map(|t| &*t.id).collect::<Vec<_>>(), vec!["TX_A", "TX_B"]);
    assert_eq!((ts[0].exons.clone(), &*ts[0].gene_name, ts[0].exons_from_span), (vec![(100, 150)], "GENE_A", true));
    assert!(!ts[1].exons_from_span);
    assert!(a.find("TX_A").is_some());
    // Its consequences say the exon is only the span.
    let c = coding::consequence(ts[0], 120, "A", "G", |_| None).unwrap();
    assert!(c.warnings.iter().any(|w| w.contains("no exons")));
}

/// GTF rows like UCSC's ncbiRefSeq (no MANE tags): per transcript of `gene` on `contig`, one
/// exon over 1001-4000 and a CDS from 1101 to `cds_end` (none if 0).
fn refseq_gtf(contig: &str, gene: &str, txs: &[(&str, u32)]) -> String {
    let row = |kind: &str, s: u32, e: u32, tx: &str| {
        format!("{contig}	ncbiRefSeq	{kind}	{s}	{e}	.	-	.	gene_id \"{gene}\"; transcript_id \"{tx}\"; gene_name \"{gene}\";
")
    };
    let mut gtf = String::new();
    for &(tx, cds_end) in txs {
        gtf += &row("transcript", 1001, 4000, tx);
        gtf += &row("exon", 1001, 4000, tx);
        if cds_end > 0 {
            gtf += &row("CDS", 1101, cds_end, tx);
        }
    }
    gtf
}

fn mains(a: &Annotation, contig: &str) -> Vec<String> {
    a.query(contig, 0, 10_000).filter(|t| t.representative).map(|t| t.id.to_string()).collect()
}

#[test]
fn representative_is_mane_from_the_shipped_list() {
    // BRAF: NM_004333 is MANE Select; NM_001374258 (Plus Clinical) and an XM_ model are longer.
    let braf = [("XM_017012559.2", 3000), ("NM_001374258.1", 2900), ("NM_004333.6", 2500), ("NR_027676.2", 0)];
    let a = Annotation::parse(refseq_gtf("chr7", "BRAF", &braf).as_bytes(), Format::Gtf).unwrap();
    assert_eq!(mains(&a, "chr7"), vec!["NM_004333.6"]);
    let tags = |id: &str| a.query("chr7", 0, 10_000).find(|t| &*t.id == id).unwrap().tags;
    assert_eq!(tags("NM_004333.6"), TAG_MANE_SELECT, "listed transcripts are tagged");
    assert_eq!(tags("NM_001374258.1"), TAG_MANE_PLUS_CLINICAL);
    assert_eq!(tags("XM_017012559.2"), 0);
    // Another version of the Select accession still beats Plus Clinical, but isn't tagged.
    let older = [("NM_001374258.1", 2900), ("NM_004333.5", 2500)];
    let a = Annotation::parse(refseq_gtf("chr7", "BRAF", &older).as_bytes(), Format::Gtf).unwrap();
    assert_eq!(mains(&a, "chr7"), vec!["NM_004333.5"]);
    assert_eq!(a.query("chr7", 0, 10_000).find(|t| &*t.id == "NM_004333.5").unwrap().tags, 0);
    // Plus Clinical over a longer curated one.
    let plus = [("NM_999999001.1", 3000), ("NM_001374258.1", 2900)];
    let a = Annotation::parse(refseq_gtf("chr7", "BRAF", &plus).as_bytes(), Format::Gtf).unwrap();
    assert_eq!(mains(&a, "chr7"), vec!["NM_001374258.1"]);
    // UCSC copies (`_2`) and Ensembl ids are matched too.
    let copies = [("NM_999999001.1", 3000), ("NM_004333.6_2", 2500)];
    let a = Annotation::parse(refseq_gtf("chr7", "BRAF", &copies).as_bytes(), Format::Gtf).unwrap();
    assert_eq!(mains(&a, "chr7"), vec!["NM_004333.6_2"]);
    let ens = [("ENST00000999999.1", 3000), ("ENST00000646891.2", 2500)];
    let a = Annotation::parse(refseq_gtf("chr7", "BRAF", &ens).as_bytes(), Format::Gtf).unwrap();
    assert_eq!(mains(&a, "chr7"), vec!["ENST00000646891.2"]);
}

#[test]
fn representative_without_mane_is_the_longest_curated_cds() {
    // Not in MANE: the longest coding curated transcript, not the oldest accession; a longer
    // predicted model still loses.
    let txs = [("XM_000000003.1", 3500), ("NM_999999001.1", 2000), ("NM_999999002.1", 3000), ("NR_999999000.1", 0)];
    let a = Annotation::parse(refseq_gtf("chr7", "NOMANE", &txs).as_bytes(), Format::Gtf).unwrap();
    assert_eq!(mains(&a, "chr7"), vec!["NM_999999002.1"]);
    // Only models: the longest coding one.
    let models = [("XM_000001.1", 1500), ("XM_000002.1", 1800)];
    let a = Annotation::parse(refseq_gtf("chr8", "MODEL", &models).as_bytes(), Format::Gtf).unwrap();
    assert_eq!(mains(&a, "chr8"), vec!["XM_000002.1"]);
    assert_eq!(curated_number("NM_000451.4_2"), Some(451));
    assert_eq!(curated_number("ENST0001.1"), None);
}

#[test]
fn transcript_ids_find_the_primary_copy_without_the_ucsc_suffix() {
    // UCSC keeps the plain id on the fix-patch copy and names the chr10 one `_2`.
    let gtf = refseq_gtf("chr10_KQ090021v1_fix", "PTEN", &[("NM_000314.8", 2000)]) + &refseq_gtf("chr10", "PTEN", &[("NM_000314.8_2", 2000)]);
    let a = Annotation::parse(gtf.as_bytes(), Format::Gtf).unwrap();
    for id in ["NM_000314.8", "nm_000314", "NM_000314.8_2"] {
        assert_eq!(a.find(id).unwrap().contig, "chr10", "{id}");
        let t = a.find_transcript(id).unwrap();
        assert_eq!((&*t.id, t.contig), ("NM_000314.8_2", a.contig_index("chr10").unwrap() as u32), "{id}");
    }
    assert_eq!(a.find("PTEN").unwrap().contig, "chr10");
    // Consequences name the NCBI accession.
    let t = a.find_transcript("NM_000314.8").unwrap();
    assert_eq!(coding::consequence(t, 1500, "A", "G", |_| None).unwrap().transcript, "NM_000314.8");
}

#[test]
fn pseudo_autosomal_genes_find_their_chrx_copy() {
    // UCSC lists chrY first, with the plain id there and `_2` on chrX; a GENCODE-style file
    // has chrX first. Either way a name finds chrX (chrY's PARs are masked in analysis sets),
    // and `find_all` still lists both.
    let y = refseq_gtf("chrY", "SHOX", &[("NM_000451.4", 2000)]);
    let x = refseq_gtf("chrX", "SHOX", &[("NM_000451.4_2", 2000)]);
    for gtf in [y.clone() + &x, x + &y] {
        let a = Annotation::parse(gtf.as_bytes(), Format::Gtf).unwrap();
        for name in ["SHOX", "shox", "NM_000451.4", "NM_000451", "NM_000451.4_2"] {
            assert_eq!(a.find(name).unwrap().contig, "chrX", "{name}");
        }
        assert_eq!(a.find_transcript("NM_000451.4").unwrap().contig, a.contig_index("chrX").unwrap() as u32);
        assert_eq!(a.find_all("SHOX").iter().map(|g| g.contig.as_str()).collect::<Vec<_>>(), ["chrX", "chrY"]);
    }
    assert_eq!((copy_rank("chrX"), copy_rank("Y"), copy_rank("chrY_KI270740v1_random")), (0, 1, 2));
}

#[test]
fn representative_is_mane_in_real_refseq_when_given() {
    // GGV_TEST_REFSEQ_GTF=hg38.ncbiRefSeq.gtf.gz: the main transcript of well-known genes is
    // their MANE Select one, though the file has no MANE tags (on the primary contig, for genes
    // UCSC also places on fix patches).
    let Ok(path) = std::env::var("GGV_TEST_REFSEQ_GTF") else { return };
    if !Path::new(&path).exists() {
        return;
    }
    let a = Annotation::read(&path).unwrap();
    let genes = [
        ("BRAF", "NM_004333.6"),
        ("GNAS", "NM_000516.7"),
        ("KRAS", "NM_004985.5"),
        ("TP53", "NM_000546.6"),
        ("BRCA1", "NM_007294.4"),
        ("EGFR", "NM_005228.5"),
        ("TTN", "NM_001267550.2"),
        ("DMD", "NM_004006.3"),
        ("SCN1A", "NM_001165963.4"),
        ("NF1", "NM_001042492.3"),
        ("LMNA", "NM_170707.4"),
        ("RET", "NM_020975.6"),
        ("PTEN", "NM_000314.8"),
    ];
    for (gene, mane) in genes {
        let g = a.find(gene).unwrap();
        let main: Vec<_> = a.query(&g.contig, g.start, g.end).filter(|t| t.representative && &*t.gene_name == gene).collect();
        assert!(main.len() == 1 && without_copy_suffix(&main[0].id) == mane, "{gene} on {}: {:?}", g.contig, main.iter().map(|t| &t.id).collect::<Vec<_>>());
        assert!(main[0].tags & TAG_MANE_SELECT != 0, "{gene}");
        eprintln!("{gene}: {} on {}", main[0].id, g.contig);
    }
    // Non-coding MANE Select genes (NR_): UCSC's file has older versions of some, still the
    // gene's main transcript.
    let noncoding = [
        ("H19", "NR_002196"),
        ("GAS5", "NR_002578"),
        ("MIAT", "NR_033320"),
        ("PCA3", "NR_132312"),
        ("IFNG-AS1", "NR_104125"),
        ("LINC00632", "NR_173140"),
        ("PCAT7", "NR_121567"),
    ];
    for (gene, mane) in noncoding {
        let g = a.find(gene).unwrap();
        let main: Vec<_> = a.query(&g.contig, g.start, g.end).filter(|t| t.representative && &*t.gene_name == gene).collect();
        assert!(main.len() == 1 && strip_version(without_copy_suffix(&main[0].id)) == mane, "{gene} on {}: {:?}", g.contig, main.iter().map(|t| &t.id).collect::<Vec<_>>());
        eprintln!("{gene}: {} on {}", main[0].id, g.contig);
    }
}

/// GNAS in the HG002 fixture: its MANE Select transcript starts at the last promoter, 37 kb
/// after the XLαs first exon, which is coding in the curated NM_080425.4.
#[test]
fn one_per_gene_is_the_main_transcript_or_the_best_one_where_it_does_not_reach() {
    let gtf = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio/ncbiRefSeq.slices.gtf");
    let a = Annotation::open(&gtf).unwrap();
    let at = |pos: u64| -> Vec<String> {
        let here = a.query("chr20", pos, pos + 1).filter(|t| t.start <= pos && pos < t.end);
        one_per_gene(here).into_iter().filter(|t| &*t.gene_name == "GNAS").map(|t| without_copy_suffix(&t.id).to_owned()).collect()
    };
    // Inside the MANE transcript: that one, though a dozen others are there too.
    assert_eq!(at(58_900_000), ["NM_000516.7"]);
    // In the upstream exon it lacks (1-based 58,854,246): the curated coding transcript with
    // the longest coding sequence, not a predicted XM_ model and not nothing.
    assert_eq!(at(58_854_245), ["NM_080425.4"]);
    assert!(!a.query("chr20", 58_854_245, 58_854_246).any(|t| t.representative && &*t.gene_name == "GNAS"));
}
