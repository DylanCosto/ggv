//! `ggv batch` and `ggv query` on small files written here: BAMs from SAM text, indexed FASTAs.

use std::path::{Path, PathBuf};

use clap::Parser;

use crate::batch::{self, Batch};

/// A temporary folder, removed when dropped.
struct Dir(PathBuf);

impl Dir {
    fn new(name: &str) -> Self {
        static N: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        let n = N.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let d = std::env::temp_dir().join(format!("ggv-cli-{name}-{}-{n}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        Self(d)
    }

    fn file(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// A made-up sequence of `len` bases (a prefix of any longer one).
fn bases(len: usize) -> String {
    let mut x: u32 = 12345;
    (0..len)
        .map(|_| {
            x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
            b"ACGT"[(x >> 16) as usize % 4] as char
        })
        .collect()
}

fn fasta(d: &Dir, name: &str, contigs: &[(&str, &str)]) -> PathBuf {
    let p = d.file(name);
    let text: String = contigs.iter().map(|(n, s)| format!(">{n}\n{s}\n")).collect();
    std::fs::write(&p, text).unwrap();
    ggv_io::build_index(&p).unwrap();
    p
}

/// A SAM line: forward strand, MAPQ 60, base qualities 40.
fn read(name: &str, contig: &str, pos: u64, cigar: &str, seq: &str) -> String {
    format!("{name}\t0\t{contig}\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{}", "I".repeat(seq.len()))
}

/// Reads of 40 bases tiling `seq` (1-based from 1) every 10 bases.
fn tiling(contig: &str, seq: &str) -> Vec<String> {
    (0..seq.len().saturating_sub(40)).step_by(10).map(|i| read(&format!("{contig}-{i}"), contig, i as u64 + 1, "40M", &seq[i..i + 40])).collect()
}

/// An indexed BAM of `reads` (SAM lines, sorted here).
fn bam(d: &Dir, name: &str, contigs: &[(&str, u64)], reads: &[String]) -> PathBuf {
    use ggv_io::pasta::{bam, sam};
    use sam::alignment::io::Write as _;
    let mut lines = reads.to_vec();
    let key = |l: &String| {
        let f: Vec<&str> = l.split('\t').collect();
        (contigs.iter().position(|c| c.0 == f[2]).unwrap(), f[3].parse::<u64>().unwrap())
    };
    lines.sort_by_key(key);
    let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for (n, len) in contigs {
        text += &format!("@SQ\tSN:{n}\tLN:{len}\n");
    }
    for l in &lines {
        text += l;
        text.push('\n');
    }
    let mut reader = sam::io::Reader::new(text.as_bytes());
    let header = reader.read_header().unwrap();
    let p = d.file(name);
    let mut w = bam::io::Writer::new(std::fs::File::create(&p).unwrap());
    w.write_alignment_header(&header).unwrap();
    for r in reader.record_bufs(&header) {
        w.write_alignment_record(&header, &r.unwrap()).unwrap();
    }
    w.try_finish().unwrap();
    drop(w);
    ggv_io::build_index(&p).unwrap();
    p
}

fn batch(d: &Dir, list: &str, list_name: &str, bams: &[&Path], reference: Option<&Path>, jobs: usize) -> (batch::Outcome, String) {
    batch_genes(d, list, list_name, bams, reference, None, jobs)
}

/// [`batch`] with a gene annotation.
fn batch_genes(d: &Dir, list: &str, list_name: &str, bams: &[&Path], reference: Option<&Path>, genes: Option<&Path>, jobs: usize) -> (batch::Outcome, String) {
    batch_vcfs(d, list, list_name, bams, reference, genes, &[], jobs)
}

/// [`batch_genes`] with `-v` files.
#[allow(clippy::too_many_arguments)]
fn batch_vcfs(d: &Dir, list: &str, list_name: &str, bams: &[&Path], reference: Option<&Path>, genes: Option<&Path>, vcfs: &[&Path], jobs: usize) -> (batch::Outcome, String) {
    let path = d.file(list_name);
    std::fs::write(&path, list).unwrap();
    let out = d.file(&format!("out-{list_name}"));
    let b = Batch {
        list: path,
        bams: bams.iter().map(|p| p.to_path_buf()).collect(),
        reference: reference.map(Path::to_path_buf),
        vcfs: vcfs.iter().map(|p| p.to_path_buf()).collect(),
        genes: genes.map(Path::to_path_buf),
        out: out.clone(),
        width: 300,
        scale: 1.0,
        flank: 20,
        max_rows: 10,
        jobs,
        pileup: ggv_query::PileupOptions::default(),
        color_by: ggv_render::ColorBy::None,
        group: ggv_layout::GroupBy::None,
        tag: None,
        link: false,
    };
    let outcome = batch::run(&b).unwrap();
    (outcome, std::fs::read_to_string(out.join("summary.tsv")).unwrap())
}

/// summary.tsv's rows (not `#` lines), split into columns.
fn rows(tsv: &str) -> Vec<Vec<String>> {
    tsv.lines().filter(|l| !l.starts_with('#')).map(|l| l.split('\t').map(String::from).collect()).collect()
}

#[test]
fn batch_survives_a_reference_shorter_than_the_contig() {
    // A GRCh37 BAM with an hg38 reference: the place runs past the reference's chr1. Slicing the
    // reference there panicked, and the whole run was lost.
    let d = Dir::new("short-ref");
    let seq = bases(300);
    let bam = bam(&d, "a.bam", &[("chrQ", 300)], &tiling("chrQ", &seq));
    let fa = fasta(&d, "ref.fa", &[("chrQ", &seq[..150])]);
    let (o, tsv) = batch(&d, "chrQ\t100\t250\tend-of-reference\nchrQ\t10\t30\tok\n", "ends.bed", &[&bam], Some(&fa), 1);
    assert_eq!((o.written, o.failed), (2, 0), "{tsv}");
}

#[test]
fn batch_reads_hla_contigs_and_a_byte_order_mark() {
    // `HLA-A*01:01:01` exists too: the bare name used to land there, at position 1.
    let d = Dir::new("hla");
    let seq = bases(100);
    let bam = bam(&d, "hla.bam", &[("HLA-A*01:01:01:01", 100), ("HLA-A*01:01:01", 100)], &tiling("HLA-A*01:01:01:01", &seq));
    let list = "\u{feff}HLA-A*01:01:01:01\t9\t29\tbedrow\r\nHLA-A*01:01:01:01\r\nHLA-A*01:01:01:01:15 locusline\r\n";
    let (o, tsv) = batch(&d, list, "hla.txt", &[&bam], None, 1);
    assert_eq!((o.written, o.failed), (3, 0), "{tsv}");
    let r = rows(&tsv);
    let got: Vec<(&str, &str)> = r.iter().map(|c| (c[1].as_str(), c[2].as_str())).collect();
    assert_eq!(got, [("HLA-A*01:01:01:01:10-29", "bedrow"), ("HLA-A*01:01:01:01:1-100", ""), ("HLA-A*01:01:01:01:15", "locusline")]);
}

#[test]
fn batch_finds_a_contig_in_any_alignment_file() {
    // The first file lacks chrQ; the place used to fail outright, depending on `-b` order.
    let d = Dir::new("any-bam");
    let seq = bases(200);
    let r_only = bam(&d, "r.bam", &[("chrR", 200)], &tiling("chrR", &seq));
    let q = bam(&d, "q.bam", &[("chrQ", 200)], &tiling("chrQ", &seq));
    let (o, tsv) = batch(&d, "chrQ:105\nchrZ:5\n", "mix.txt", &[&r_only, &q], None, 1);
    assert_eq!((o.written, o.failed), (1, 1), "{tsv}");
    let r = rows(&tsv);
    assert_eq!((r[0][1].as_str(), r[0][4].as_str()), ("chrQ:105", ""));
    assert!(r[0][5].starts_with("depth "), "{tsv}");
    assert!(r[1][8].contains("in none of the alignment files"), "{tsv}");
}

#[test]
fn batch_keeps_insertion_sites_and_reports_unreadable_lines() {
    let d = Dir::new("bad-lines");
    let seq = bases(200);
    let bam = bam(&d, "q.bam", &[("chrQ", 200)], &tiling("chrQ", &seq));
    let (o, tsv) = batch(&d, "chrQ\t105\t105\tinsertion-site\nchrQ 10x\nchrQ:150\n", "l.txt", &[&bam], None, 2);
    // The unreadable line counts as failed (exit 3), in list order.
    assert_eq!((o.written, o.failed), (2, 1), "{tsv}");
    let r = rows(&tsv);
    assert_eq!((r[0][0].as_str(), r[0][1].as_str(), r[0][2].as_str()), ("1", "chrQ:105", "insertion-site"));
    assert_eq!((r[1][0].as_str(), r[1][1].as_str()), ("2", "chrQ 10x"));
    assert!(r[1][7].starts_with("line 2 of the list is not a place"), "{tsv}");
    assert_eq!(r[2][1], "chrQ:150");
}

#[test]
fn batch_records_a_truncated_file_once() {
    let d = Dir::new("trunc");
    let seq = bases(200);
    let whole = bam(&d, "q.bam", &[("chrQ", 200)], &tiling("chrQ", &seq));
    // The same file without its 28-byte end-of-file marker, and the same index.
    let cut = d.file("cut");
    std::fs::create_dir_all(&cut).unwrap();
    let bytes = std::fs::read(&whole).unwrap();
    let t = cut.join("q.bam");
    std::fs::write(&t, &bytes[..bytes.len() - 28]).unwrap();
    for e in std::fs::read_dir(&d.0).unwrap().flatten() {
        let n = e.file_name().to_string_lossy().into_owned();
        if n.starts_with("q.bam.") {
            std::fs::copy(e.path(), cut.join(&n)).unwrap();
        }
    }
    let (o, tsv) = batch(&d, "chrQ:50\nchrQ:100\nchrQ:150\n", "t.txt", &[&t], None, 3);
    assert_eq!((o.written, o.failed), (3, 0), "{tsv}");
    let warnings: Vec<&str> = tsv.lines().filter(|l| l.starts_with("## warning: ")).collect();
    assert_eq!(warnings.len(), 1, "{tsv}");
    assert!(warnings[0].contains("no BGZF end-of-file marker"));
    let index = std::fs::read_to_string(d.file("out-t.txt").join("index.html")).unwrap();
    assert!(index.contains("no BGZF end-of-file marker"));
    // And stderr says it once per file, however often it is opened.
    assert!(crate::first_time("ggv-cli test warning"));
    assert!(!crate::first_time("ggv-cli test warning"));
}

#[test]
fn junctions_refuse_a_reference_without_the_contig() {
    // A `mismatches` filter used to count none, and every motif was "unknown", with exit 0.
    let d = Dir::new("junctions");
    let seq = bases(200);
    let spliced = read("s", "chrQ", 30, "10M50N10M", &format!("{}{}", &seq[29..39], &seq[89..99]));
    let bam = bam(&d, "q.bam", &[("chrQ", 200), ("chrR", 50)], &[spliced]);
    let good = fasta(&d, "ref.fa", &[("chrQ", &seq), ("chrR", &seq[..50])]);
    let partial = fasta(&d, "ref2.fa", &[("chrR", &seq[..50])]);
    let run = |reference: &Path| {
        let cli = crate::Cli::try_parse_from(["ggv", "query", "junctions", &bam.display().to_string(), "chrQ", "-f", &reference.display().to_string(), "--filter", "mismatches == 0", "--format", "jsonl"]).unwrap();
        crate::run(cli).map_err(|e| format!("{e:#}"))
    };
    assert!(run(&good).is_ok());
    let e = run(&partial).unwrap_err();
    assert!(e.contains("chrQ") && e.contains("not in this file"), "{e}");
}

#[test]
fn summary_reads_line_says_what_the_filter_removed() {
    let r = ggv_query::ReadCounts { seen: 10, used: 6, duplicate: 1, filtered_out: 3, ..Default::default() };
    assert!(crate::reads_line(&r).ends_with(" other=0 filtered=3"), "{}", crate::reads_line(&r));
}

#[test]
fn signal_bins_never_end_before_they_start() {
    // Five bins over three bases gave `chr1 1 0` rows.
    let e = crate::bin_edges(0, 3, 5);
    assert_eq!(e, [(0, 1), (1, 2), (2, 3)]);
    assert_eq!(crate::bin_edges(0, 200, 2), [(0, 100), (100, 200)]);
    assert_eq!(crate::bin_edges(10, 10, 4), [(10, 11)]);
}

/// `base` changed to another base (`k` picks which).
fn other(base: u8, k: usize) -> char {
    let b = b"ACGT";
    let i = b.iter().position(|&x| x == base).unwrap();
    b[(i + 1 + k % 3) % 4] as char
}

#[test]
fn batch_names_only_credible_changes_every_gene_and_insertions() {
    // bughunt3: one error read became `p.` missense, reads without sequence an `N` change, four
    // isoforms of one gene hid the antisense gene, and insertion places reported a noise base.
    let d = Dir::new("consequence");
    let seq = bases(300);
    let s = seq.as_bytes();
    let mut reads = Vec::new();
    for (k, start) in (0..=256usize).step_by(4).enumerate() {
        let mut r: Vec<u8> = s[start..start + 40].to_vec();
        // chrQ:100 — one read of ten with another base.
        if start == 80 {
            r[99 - start] = other(s[99], 0) as u8;
        }
        // chrQ:150 — every other read.
        if start <= 149 && start + 40 > 149 && k % 2 == 0 {
            r[149 - start] = other(s[149], 0) as u8;
        }
        let r = String::from_utf8(r).unwrap();
        // chrQ:200 — an insertion after the base in every read spanning it.
        if start < 195 && start + 40 > 205 {
            let at = 200 - start;
            reads.push(read(&format!("r{start}"), "chrQ", start as u64 + 1, &format!("{at}M2I{}M", 40 - at), &format!("{}AC{}", &r[..at], &r[at..])));
        } else {
            reads.push(read(&format!("r{start}"), "chrQ", start as u64 + 1, "40M", &r));
        }
    }
    // chrQ:60 — five reads without sequence.
    for start in [25u64, 27, 29, 31, 33] {
        reads.push(format!("n{start}\t0\tchrQ\t{}\t60\t40M\t*\t0\t0\t*\t*", start + 1));
    }
    let bam = bam(&d, "q.bam", &[("chrQ", 300)], &reads);
    let fa = fasta(&d, "ref.fa", &[("chrQ", &seq)]);
    // G1: a curated transcript and four predicted ones; G2: a non-coding antisense gene.
    let mut gtf = String::new();
    let row = |kind: &str, a: u64, b: u64, strand: char, gene: &str, tx: &str| format!("chrQ\tt\t{kind}\t{a}\t{b}\t.\t{strand}\t{}\tgene_id \"{gene}\"; transcript_id \"{tx}\"; gene_name \"{gene}\";\n", if kind == "CDS" { "0" } else { "." });
    for tx in ["XM_1.1", "XM_2.1", "XM_3.1", "XM_4.1", "NM_9.1"] {
        gtf += &row("transcript", 1, 300, '+', "G1", tx);
        gtf += &row("exon", 1, 300, '+', "G1", tx);
        gtf += &row("CDS", 1, 297, '+', "G1", tx);
    }
    gtf += &row("transcript", 120, 180, '-', "G2", "NR_5.1");
    gtf += &row("exon", 120, 180, '-', "G2", "NR_5.1");
    let genes = d.file("g.gtf");
    std::fs::write(&genes, gtf).unwrap();
    let vcf_alt = other(s[149], 1);
    let list = format!("chrQ:100\tnoise\nchrQ:60\tno-sequence\nchrQ:150\treal\nchrQ\t200\t200\tins-site\nchrQ\t150\tv1\t{}\t{vcf_alt}\t.\tPASS\t.\n", s[149] as char);
    let (o, tsv) = batch_genes(&d, &list, "c.txt", &[&bam], Some(&fa), Some(&genes), 1);
    assert_eq!((o.written, o.failed), (5, 0), "{tsv}");
    let r = rows(&tsv);
    // No `c.` change from one read in ten, nor from reads without sequence.
    assert!(!r[0][3].contains("c.100"), "{tsv}");
    assert!(!r[1][3].contains('>'), "{tsv}");
    // Half the reads: a change, named on the curated transcript of the coding gene.
    let read_alt = format!("c.150{}>{}", s[149] as char, other(s[149], 0));
    assert!(r[2][3].starts_with("G1 ") && r[2][3].contains(&read_alt), "{tsv}");
    let report = std::fs::read_to_string(d.file("out-c.txt").join(&r[2][6])).unwrap();
    assert!(report.contains("NM_9.1") && report.contains("G2 NR_5.1"), "{report}");
    // At the insertion site, the insertion, not a base.
    assert!(r[3][4].contains(" INS ") && r[3][4].ends_with('%') && !r[3][4].ends_with(" 0%"), "{tsv}");
    // A VCF row's alleles, not the reads'.
    assert!(r[4][3].contains(&format!("c.150{}>{vcf_alt}", s[149] as char)) && !r[4][3].contains(&read_alt), "{tsv}");
}

#[test]
fn batch_links_work_for_hash_contigs() {
    // PanSN names: `#` in a file name started the link's URL fragment.
    let d = Dir::new("pansn");
    let seq = bases(200);
    let bam = bam(&d, "p.bam", &[("HG002#1#chrQ", 200)], &tiling("HG002#1#chrQ", &seq));
    let (o, tsv) = batch(&d, "HG002#1#chrQ:15\nchrZ:5\n", "p.txt", &[&bam], None, 1);
    assert_eq!((o.written, o.failed), (1, 1), "{tsv}");
    let r = rows(&tsv);
    assert_eq!(r[0][5], "001-HG002_1_chrQ-15.png");
    let index = std::fs::read_to_string(d.file("out-p.txt").join("index.html")).unwrap();
    assert!(index.contains("href=\"001-HG002_1_chrQ-15.html\"") && !index.contains("href=\"\""), "{index}");
}

#[test]
fn demo_suggests_commands_ggv_accepts() {
    // `ggv mcp <bam>` was printed, which clap rejects: alignments take `-b`.
    // A folder that breaks unquoted commands: a space, parentheses, an apostrophe and `&`.
    let dir = PathBuf::from("C:/Program Files (x86)/O'Brien & co/v1");
    let d = ggv_io::demo::Demo {
        bam: dir.join("a.bam"),
        vcfs: vec![dir.join("a.vcf.gz"), dir.join("b.vcf.gz")],
        genes: dir.join("g.gtf"),
        reference: dir.join("r.fa.gz"),
        dir,
    };
    // Shell words: double quotes group.
    let words = |c: &str| {
        let (mut out, mut w, mut quoted) = (Vec::new(), String::new(), false);
        for ch in c.chars() {
            match ch {
                '"' => quoted = !quoted,
                ' ' if !quoted => out.extend((!w.is_empty()).then(|| std::mem::take(&mut w))),
                ch => w.push(ch),
            }
        }
        out.extend((!w.is_empty()).then_some(w));
        out
    };
    let commands = crate::demo_commands(&d);
    let ggv: Vec<&String> = commands.iter().filter(|c| c.starts_with("ggv ")).collect();
    assert_eq!(ggv.len(), 2, "{commands:?}");
    for c in ggv {
        let cli = crate::Cli::try_parse_from(words(c)).unwrap_or_else(|e| panic!("{c}: {e}"));
        if let crate::Command::Mcp { bams, vcfs, genes, reference, .. } = cli.command {
            assert_eq!((bams, vcfs.len(), genes.len(), reference), (vec![d.bam.clone()], 2, 1, Some(d.reference.clone())));
        }
    }
}

#[test]
fn info_reads_a_bam_without_an_index() {
    let d = Dir::new("noindex");
    let seq = bases(200);
    let indexed = bam(&d, "q.bam", &[("chrQ", 200)], &tiling("chrQ", &seq));
    let bare = d.file("bare");
    std::fs::create_dir_all(&bare).unwrap();
    let p = bare.join("noidx.bam");
    std::fs::copy(&indexed, &p).unwrap();
    let info = ggv_io::BamFile::header_info(&p).unwrap();
    assert_eq!((info.index.as_str(), info.contigs.len()), ("none", 1));
    let run = |args: &[&str]| crate::run(crate::Cli::try_parse_from(args).unwrap()).map_err(|e| format!("{e:#}"));
    let p = p.display().to_string();
    run(&["ggv", "info", &p, "--format", "json"]).unwrap();
    // Region queries still need the index.
    assert!(run(&["ggv", "query", "pileup", &p, "chrQ:5"]).unwrap_err().contains("no index found"));
}

#[test]
fn batch_genome_fetches_only_what_f_and_g_leave_out() {
    // bughunt4: `--genome hg38 -f ref -g genes` downloaded hg38 although nothing of it was used.
    let args = |extra: &[&str]| {
        let mut v = vec!["ggv", "batch", "l.txt", "-b", "a.bam", "-o", "out", "--genome", "hg38"];
        v.extend_from_slice(extra);
        match crate::Cli::try_parse_from(v).unwrap().command {
            crate::Command::Batch(a) => a,
            _ => unreachable!(),
        }
    };
    let local = || ggv_io::genomes::Local { reference: "cache/hg38.2bit".into(), genes: "cache/hg38.gtf.gz".into() };
    let mut both = args(&["-f", "my.fa", "-g", "my.gtf"]);
    crate::batch_genome(&mut both, |_| panic!("fetched a genome -f and -g replace")).unwrap();
    assert_eq!((both.reference.unwrap(), both.genes.unwrap()), ("my.fa".into(), "my.gtf".into()));
    // One missing: fetched, and the given one still wins.
    let mut genes_only = args(&["-g", "my.gtf"]);
    crate::batch_genome(&mut genes_only, |_| Ok(local())).unwrap();
    assert_eq!((genes_only.reference.unwrap(), genes_only.genes.unwrap()), ("cache/hg38.2bit".into(), "my.gtf".into()));
    let mut neither = args(&[]);
    crate::batch_genome(&mut neither, |n| if n == "hg38" { Ok(local()) } else { panic!("{n}") }).unwrap();
    assert_eq!(neither.genes.unwrap(), PathBuf::from("cache/hg38.gtf.gz"));
}

/// A bgzipped, tabix-indexed VCF of `body` with one sample, `S1`.
fn vcf(d: &Dir, name: &str, body: &str) -> PathBuf {
    use std::io::Write;
    let p = d.file(name);
    let text = format!(
        "##fileformat=VCFv4.2
##contig=<ID=chrQ,length=300>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"t\">
##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
{body}"
    );
    let mut w = ggv_io::pasta::bgzf::io::Writer::new(std::fs::File::create(&p).unwrap());
    w.write_all(text.as_bytes()).unwrap();
    w.finish().unwrap();
    ggv_io::build_index(&p).unwrap();
    p
}

#[test]
fn batch_summary_is_the_most_severe_change_at_the_place() {
    // bughunt4: the summary named whichever overlapping gene came first (an intron of GA over
    // a coding change in GB); any VCF record overlapping the place (a 150 bp <DUP>, a RefCall,
    // a deletion starting upstream) hid the change every read shows; an MNV starting a base
    // upstream was the summary of a list line at the SNV's base; and `chrQ:150<TAB>30<TAB>99`
    // was a BED row on a contig `chrQ:150`.
    let d = Dir::new("severity");
    let seq = bases(300);
    let s = seq.as_bytes();
    let read_alt = other(s[149], 0);
    let mut alt_seq = s.to_vec();
    alt_seq[149] = read_alt as u8;
    let alt_seq = String::from_utf8(alt_seq).unwrap();
    let reads: Vec<String> = (0..=256usize).step_by(4).map(|i| read(&format!("r{i}"), "chrQ", i as u64 + 1, "40M", &alt_seq[i..i + 40])).collect();
    let bam = bam(&d, "q.bam", &[("chrQ", 300)], &reads);
    let fa = fasta(&d, "ref.fa", &[("chrQ", &seq)]);
    // GA (first in the file) has chrQ:150 in an intron; GB codes there (c.51).
    let row = |kind: &str, a: u64, b: u64, gene: &str| format!("chrQ	t	{kind}	{a}	{b}	.	+	{}	gene_id \"{gene}\"; transcript_id \"NM_{gene}.1\"; gene_name \"{gene}\";
", if kind == "CDS" { "0" } else { "." });
    let gtf = [row("transcript", 1, 300, "GA"), row("exon", 1, 60, "GA"), row("exon", 240, 300, "GA"), row("CDS", 1, 60, "GA"), row("CDS", 240, 297, "GA"), row("transcript", 100, 200, "GB"), row("exon", 100, 200, "GB"), row("CDS", 100, 198, "GB")].concat();
    let genes = d.file("g.gtf");
    std::fs::write(&genes, gtf).unwrap();
    let r150 = s[149] as char;
    let overlapping = vcf(
        &d,
        "over.vcf.gz",
        &format!(
            "chrQ	100	dup	N	<DUP>	.	PASS	SVTYPE=DUP;END=250	GT	0/1
chrQ	145	del	{}	{}	.	PASS	.	GT	0/1
chrQ	150	rc	{r150}	{}	.	RefCall	.	GT	./.
chrQ	150	homref	{r150}	{}	.	PASS	.	GT	0/0
",
            &seq[144..154],
            &seq[144..145],
            other(s[149], 1),
            other(s[149], 2)
        ),
    );
    let (o, tsv) = batch_vcfs(&d, "chrQ:150	reads
chrQ:150	30	99
", "l.txt", &[&bam], Some(&fa), Some(&genes), &[&overlapping], 1);
    assert_eq!((o.written, o.failed), (2, 0), "{tsv}");
    let r = rows(&tsv);
    assert!(r[0][3].starts_with(&format!("GB c.51{r150}>{read_alt}")), "{tsv}");
    assert_eq!((r[1][1].as_str(), r[1][2].as_str()), ("chrQ:150", "30 99"), "{tsv}");
    // A called SNV at the base, after an MNV starting a base upstream: the SNV is the summary.
    let snv = other(s[149], 1);
    let mnv = vcf(&d, "mnv.vcf.gz", &format!("chrQ	149	mnv	{}	{}{}	.	PASS	.	GT	0/1
chrQ	150	snv	{r150}	{snv}	.	PASS	.	GT	0/1
", &seq[148..150], other(s[148], 0), other(s[149], 2)));
    let (_, tsv) = batch_vcfs(&d, "chrQ:150	hot
", "m.txt", &[&bam], Some(&fa), Some(&genes), &[&mnv], 1);
    assert!(rows(&tsv)[0][3].starts_with(&format!("GB c.51{r150}>{snv}")), "{tsv}");
}

#[test]
fn batch_counts_insertions_at_a_symbolic_ins_call() {
    // bughunt5: a PASS `<INS>` carried 1/1 at the place was no call (no sequence ALT), so the
    // per-file column showed the top base, not the insertion the reads carry.
    let d = Dir::new("symbolic-ins");
    let seq = bases(300);
    let s = seq.as_bytes();
    let reads: Vec<String> = (0..=256usize)
        .step_by(4)
        .map(|start| {
            let r = &seq[start..start + 40];
            if start < 195 && start + 40 > 205 {
                let at = 200 - start;
                read(&format!("r{start}"), "chrQ", start as u64 + 1, &format!("{at}M2I{}M", 40 - at), &format!("{}AC{}", &r[..at], &r[at..]))
            } else {
                read(&format!("r{start}"), "chrQ", start as u64 + 1, "40M", r)
            }
        })
        .collect();
    let bam = bam(&d, "q.bam", &[("chrQ", 300)], &reads);
    let fa = fasta(&d, "ref.fa", &[("chrQ", &seq)]);
    let ins = vcf(&d, "ins.vcf.gz", &format!("chrQ	200	ins	{}	<INS>	.	PASS	SVTYPE=INS;END=200	GT	1/1\n", s[199] as char));
    let (o, tsv) = batch_vcfs(&d, "chrQ:200\tsite\n", "i.txt", &[&bam], Some(&fa), None, &[&ins], 1);
    assert_eq!((o.written, o.failed), (1, 0), "{tsv}");
    let r = rows(&tsv);
    assert!(r[0][4].contains(" INS ") && !r[0][4].ends_with(" 0%"), "{tsv}");
}

#[test]
fn batch_summary_without_a_change_names_the_first_of_equal_genes() {
    // bughunt5: at an intron of a coding and a non-coding gene, both their main transcripts,
    // the summary with no change named the last (non-coding) gene; with a change, the first.
    let d = Dir::new("tie");
    let seq = bases(300);
    let bam = bam(&d, "q.bam", &[("chrQ", 300)], &tiling("chrQ", &seq));
    let fa = fasta(&d, "ref.fa", &[("chrQ", &seq)]);
    let row = |kind: &str, a: u64, b: u64, gene: &str, tx: &str| format!("chrQ	t	{kind}	{a}	{b}	.	+	{}	gene_id \"{gene}\"; transcript_id \"{tx}\"; gene_name \"{gene}\";\n", if kind == "CDS" { "0" } else { "." });
    // The non-coding gene first in the file: the coding one still comes first.
    let gtf = [
        row("transcript", 80, 220, "BBB", "NR_2.1"),
        row("exon", 80, 100, "BBB", "NR_2.1"),
        row("exon", 200, 220, "BBB", "NR_2.1"),
        row("transcript", 1, 300, "AAA", "NM_1.1"),
        row("exon", 1, 60, "AAA", "NM_1.1"),
        row("exon", 240, 300, "AAA", "NM_1.1"),
        row("CDS", 1, 60, "AAA", "NM_1.1"),
        row("CDS", 240, 299, "AAA", "NM_1.1"),
    ]
    .concat();
    let genes = d.file("g.gtf");
    std::fs::write(&genes, gtf).unwrap();
    let (o, tsv) = batch_genes(&d, "chrQ:150\tnochange\n", "t.txt", &[&bam], Some(&fa), Some(&genes), 1);
    assert_eq!((o.written, o.failed), (1, 0), "{tsv}");
    assert_eq!(rows(&tsv)[0][3], "AAA intron", "{tsv}");
}
