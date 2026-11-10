//! The MCP tools on small files written here: BAMs from SAM text, an indexed FASTA, a GTF.

use std::path::{Path, PathBuf};

use serde_json::{Value, json};

use super::{
    Server,
    tools::{Defaults, State},
};

/// A temporary folder, removed when dropped.
struct Dir(PathBuf);

impl Dir {
    fn new(name: &str) -> Self {
        static N: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        let n = N.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let d = std::env::temp_dir().join(format!("ggv-mcp-{name}-{}-{n}", std::process::id()));
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

/// A made-up sequence of `len` bases.
fn bases(len: usize) -> String {
    let mut x: u32 = 12345;
    (0..len)
        .map(|_| {
            x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
            b"ACGT"[(x >> 16) as usize % 4] as char
        })
        .collect()
}

fn other(b: char) -> char {
    if b == 'A' { 'C' } else { 'A' }
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

fn server(defaults: Defaults) -> Server {
    Server { state: State::new(defaults) }
}

/// A tool call through the protocol: the answer, or the error text.
fn call(s: &mut Server, tool: &str, args: Value) -> Result<Value, String> {
    let msg = json!({"jsonrpc": "2.0", "id": 1, "method": "tools/call", "params": {"name": tool, "arguments": args}});
    let reply: Value = serde_json::from_str(&s.handle(&msg.to_string()).unwrap()).unwrap();
    let text = reply["result"]["content"][0]["text"].as_str().unwrap_or_else(|| panic!("{reply}")).to_owned();
    if reply["result"]["isError"] == true { Err(text) } else { Ok(serde_json::from_str(&text).unwrap()) }
}

fn s(p: &Path) -> String {
    p.display().to_string()
}

/// `n` reads of `len` bases from 1-based `pos`; with `alt: (k, at)`, the first `k` have another
/// base at 1-based `at`.
fn reads(prefix: &str, reference: &str, contig: &str, pos: u64, len: usize, n: usize, alt: Option<(usize, u64)>) -> Vec<String> {
    (0..n)
        .map(|k| {
            let mut seq: Vec<char> = reference[pos as usize - 1..pos as usize - 1 + len].chars().collect();
            if let Some((alt, at)) = alt
                && k < alt
            {
                let i = (at - pos) as usize;
                seq[i] = other(seq[i]);
            }
            read(&format!("{prefix}{k}"), contig, pos, &format!("{len}M"), &seq.iter().collect::<String>())
        })
        .collect()
}

#[test]
fn compare_reads_the_region_against_its_own_files() {
    let d = Dir::new("compare");
    let r = bases(600);
    let fa = fasta(&d, "ref.fa", &[("chrA", &r)]);
    let a = bam(&d, "a.bam", &[("chrA", 600)], &reads("a", &r, "chrA", 451, 50, 10, None));
    let b = bam(&d, "b.bam", &[("chrA", 600)], &reads("b", &r, "chrA", 451, 50, 10, Some((5, 470))));
    let unrelated = bam(&d, "q.bam", &[("chrQ", 600)], &[read("q", "chrQ", 1, "10M", "ACGTACGTAC")]);
    let short = bam(&d, "short.bam", &[("chrA", 300)], &[read("s", "chrA", 1, "10M", &r[..10])]);
    let args = json!({"region": "chrA:1-600", "a": s(&a), "b": s(&b), "reference": s(&fa), "threshold": 0.2, "min_depth": 10});
    for default in [None, Some(unrelated), Some(short)] {
        let mut srv = server(Defaults { bams: default.clone().into_iter().collect(), ..Default::default() });
        let out = call(&mut srv, "compare", args.clone()).unwrap_or_else(|e| panic!("{default:?}: {e}"));
        assert_eq!(out["differences"], 1, "{default:?}: {out}");
        assert_eq!(out["shown"][0]["pos"], 470);
        assert_eq!(out["query"]["compared"], "chrA:1-600", "{default:?}: not shortened");
    }
}

#[test]
fn compare_reports_insertions() {
    let d = Dir::new("compare-ins");
    let r = bases(400);
    let fa = fasta(&d, "ref.fa", &[("chrQA", &r)]);
    let plain = reads("r", &r, "chrQA", 100, 50, 10, None);
    let mut with_ins = reads("p", &r, "chrQA", 100, 50, 5, None);
    for k in 0..5 {
        let seq = format!("{}TT{}", &r[99..109], &r[109..149]);
        with_ins.push(read(&format!("i{k}"), "chrQA", 100, "10M2I40M", &seq));
    }
    let a = bam(&d, "reference.bam", &[("chrQA", 400)], &plain);
    let b = bam(&d, "insertion.bam", &[("chrQA", 400)], &with_ins);
    let mut srv = server(Defaults { reference: Some(fa), ..Default::default() });
    let out = call(&mut srv, "compare", json!({"region": "chrQA:100-149", "a": s(&a), "b": s(&b)})).unwrap();
    assert_eq!(out["differences"], 1, "{out}");
    let site = &out["shown"][0];
    assert_eq!((site["pos"].as_u64(), site["allele"].as_str()), (Some(109), Some("INS")), "{out}");
    assert_eq!((site["a_fraction"].as_f64(), site["b_fraction"].as_f64()), (Some(0.0), Some(0.5)));
}

#[test]
fn supplied_arguments_of_the_wrong_kind_are_refused() {
    let d = Dir::new("args");
    let r = bases(300);
    let fa = fasta(&d, "ref.fa", &[("chrQA", &r)]);
    let one = bam(&d, "one.bam", &[("chrQA", 300)], &reads("r", &r, "chrQA", 101, 30, 9, None));
    let mut srv = server(Defaults { bams: vec![one.clone()], reference: Some(fa), ..Default::default() });
    let ev = |srv: &mut Server, extra: Value| {
        let mut args = json!({"locus": "chrQA:110"});
        args.as_object_mut().unwrap().extend(extra.as_object().unwrap().clone());
        call(srv, "evidence", args)
    };
    let depth = |v: Value| v["results"][0]["tracks"][0]["depth"].as_u64().unwrap();
    assert_eq!(depth(ev(&mut srv, json!({})).unwrap()), 9);
    assert_eq!(depth(ev(&mut srv, json!({"min_mapq": 255})).unwrap()), 0);
    assert_eq!(depth(ev(&mut srv, json!({"filter": "mapq >= 255"})).unwrap()), 0);
    for (bad, name) in [
        (json!({"min_mapq": "255"}), "min_mapq"),
        (json!({"min_mapq": 256}), "min_mapq"),
        (json!({"filter": {"expression": "mapq >= 255"}}), "filter"),
        (json!({"bams": [true]}), "bams[0]"),
        (json!({"min_base_quality": -1}), "min_base_quality"),
        (json!({"include_duplicates": "yes"}), "include_duplicates"),
    ] {
        let e = ev(&mut srv, bad.clone()).expect_err(&bad.to_string());
        assert!(e.contains(&format!("`{name}`")), "{bad}: {e}");
        // The same through the tools' own argument handling, without the schema check.
        let mut args = json!({"locus": "chrQA:110"});
        args.as_object_mut().unwrap().extend(bad.as_object().unwrap().clone());
        let e = srv.state.call("evidence", &args).err().unwrap_or_else(|| panic!("{bad} accepted"));
        assert!(e.to_string().contains(name.trim_end_matches("[0]")), "{bad}: {e}");
    }
    let cmp = |extra: Value| {
        let mut args = json!({"region": "chrQA:101-130", "a": s(&one), "b": s(&one)});
        args.as_object_mut().unwrap().extend(extra.as_object().unwrap().clone());
        args
    };
    assert_eq!(call(&mut srv, "compare", cmp(json!({"min_depth": 1}))).unwrap()["differences"], 0);
    for (bad, name) in [(json!({"threshold": -0.01}), "threshold"), (json!({"threshold": "0.2"}), "threshold"), (json!({"min_depth": 4_294_967_296u64}), "min_depth")] {
        let e = call(&mut srv, "compare", cmp(bad.clone())).expect_err(&bad.to_string());
        assert!(e.contains(&format!("`{name}`")), "{bad}: {e}");
        let e = srv.state.call("compare", &cmp(bad.clone())).err().unwrap_or_else(|| panic!("{bad} accepted"));
        assert!(e.to_string().contains(name), "{bad}: {e}");
    }
}

#[test]
fn coverage_bins_agree_on_their_bounds() {
    let d = Dir::new("bins");
    let r = bases(300);
    let one = bam(&d, "reads.bam", &[("chrQA", 300)], &reads("r", &r, "chrQA", 110, 10, 9, None));
    let mut srv = server(Defaults { bams: vec![one], ..Default::default() });
    for q in [0, 1] {
        for bins in [3, 6, 7] {
            let out = call(&mut srv, "coverage", json!({"region": "chrQA:110-119", "bins": bins, "min_base_quality": q})).unwrap();
            for b in out["results"][0]["depth"].as_array().unwrap() {
                assert_eq!(b[0].as_f64(), Some(9.0), "q{q} bins {bins}: {out}");
            }
        }
    }
}

#[test]
fn coverage_says_when_a_file_ends_before_the_region() {
    let d = Dir::new("short");
    let r = bases(1000);
    let long = bam(&d, "reads.bam", &[("chrQA", 1000)], &reads("l", &r, "chrQA", 101, 20, 10, None));
    let short = bam(&d, "short.bam", &[("chrQA", 120)], &reads("s", &r, "chrQA", 101, 20, 10, None));
    let mut srv = server(Defaults { bams: vec![long, short], ..Default::default() });
    for q in [0, 1] {
        let out = call(&mut srv, "coverage", json!({"region": "chrQA:100-200", "bins": 5, "min_base_quality": q})).unwrap();
        let (a, b) = (&out["results"][0], &out["results"][1]);
        assert_eq!(a["region"], "chrQA:100-200");
        assert!(a.get("clipped").is_none(), "{a}");
        assert_eq!(b["region"], "chrQA:100-120", "q{q}: {b}");
        assert!((b["bin_size"].as_f64().unwrap() - 4.2).abs() < 1e-9, "q{q}: {b}");
        assert!(b["clipped"].as_str().unwrap().contains("ends at 120"), "{b}");
    }
}

#[test]
fn a_named_reference_that_cannot_be_read_is_an_error() {
    let d = Dir::new("ref");
    let r = bases(300);
    let fa = fasta(&d, "ref.fa", &[("chrQ", &r)]);
    let one = bam(&d, "reads.bam", &[("chrQ", 300)], &reads("r", &r, "chrQ", 91, 20, 10, Some((8, 101))));
    let absent = d.file("absent.fa");
    let gtf = d.file("g.gtf");
    std::fs::write(&gtf, "chrQ\tX\texon\t51\t150\t.\t+\t.\tgene_id \"G\"; transcript_id \"T\"; gene_name \"G\";\nchrQ\tX\tCDS\t51\t150\t.\t+\t0\tgene_id \"G\"; transcript_id \"T\"; gene_name \"G\";\n").unwrap();
    let mut srv = server(Defaults { bams: vec![one.clone()], ..Default::default() });
    let out = call(&mut srv, "summarize", json!({"region": "chrQ:101", "reference": s(&fa)})).unwrap();
    assert_eq!(out["results"][0]["variable_sites"][0]["ref"], r[100..101]);
    let e = call(&mut srv, "summarize", json!({"region": "chrQ:101", "reference": s(&absent)})).unwrap_err();
    assert!(e.contains("absent.fa"), "{e}");
    // Left out: reference-free, as before.
    assert!(call(&mut srv, "summarize", json!({"region": "chrQ:101"})).unwrap()["results"][0]["variable_sites"][0]["ref"].is_null());
    let change = json!({"locus": "chrQ:101", "ref": &r[100..101], "alt": other(r.as_bytes()[100] as char).to_string(), "genes": [s(&gtf)]});
    assert!(call(&mut srv, "consequence", change.clone()).is_ok());
    let mut bad = change.clone();
    bad["reference"] = json!(s(&absent));
    assert!(call(&mut srv, "consequence", bad).unwrap_err().contains("absent.fa"));
    // A default reference that can't be read is reported too.
    let mut srv = server(Defaults { bams: vec![one], reference: Some(absent), ..Default::default() });
    assert!(call(&mut srv, "summarize", json!({"region": "chrQ:101"})).unwrap_err().contains("absent.fa"));
}

#[test]
fn wide_summaries_keep_the_reference() {
    // Wider than the old 20 Mb cut-off, past which the reference was dropped.
    let d = Dir::new("wide");
    let len = 20_000_001;
    let r = "C".repeat(len);
    let fa = fasta(&d, "ref.fa", &[("chrQ", &r)]);
    let one = bam(&d, "reads.bam", &[("chrQ", len as u64)], &reads("r", &r, "chrQ", 91, 20, 10, Some((8, 101))));
    let mut srv = server(Defaults { bams: vec![one], reference: Some(fa), ..Default::default() });
    for region in ["chrQ:1-20000000", "chrQ:1-20000001"] {
        let out = call(&mut srv, "summarize", json!({"region": region})).unwrap();
        let site = &out["results"][0]["variable_sites"][0];
        assert_eq!((site["ref"].as_str(), site["alt"].as_str(), site["alt_fraction"].as_f64()), (Some("C"), Some("A"), Some(0.8)), "{region}: {out}");
    }
}

#[test]
fn phase_says_which_block_each_site_is_in() {
    let d = Dir::new("phase");
    let r = bases(400);
    let fa = fasta(&d, "ref.fa", &[("chrQA", &r)]);
    let mut all = Vec::new();
    for k in 0..10 {
        let mut seq: Vec<char> = r[100..140].chars().collect();
        if k < 5 {
            seq[9] = other(seq[9]); // 110
            seq[29] = other(seq[29]); // 130
        }
        all.push(read(&format!("x{k}"), "chrQA", 101, "40M", &seq.iter().collect::<String>()));
    }
    all.extend(reads("y", &r, "chrQA", 291, 20, 10, Some((5, 300))));
    let one = bam(&d, "blocks.bam", &[("chrQA", 400)], &all);
    let mut srv = server(Defaults { bams: vec![one], reference: Some(fa), ..Default::default() });
    let out = call(&mut srv, "phase", json!({"region": "chrQA:100-329"})).unwrap();
    assert_eq!(out["blocks"], 2, "{out}");
    let sites: Vec<(u64, u64)> = out["sites"].as_array().unwrap().iter().map(|s| (s["pos"].as_u64().unwrap(), s["block"].as_u64().unwrap())).collect();
    assert_eq!(sites, [(110, 1), (130, 1), (300, 2)]);
    let blocks = out["phase_blocks"].as_array().unwrap();
    assert_eq!((blocks[0]["start"].as_u64(), blocks[0]["end"].as_u64(), blocks[0]["sites"].as_u64()), (Some(110), Some(130), Some(2)));
    for b in blocks {
        assert_eq!(b["phase_A"].as_u64().unwrap() + b["phase_B"].as_u64().unwrap(), 10, "{b}");
        assert_eq!((b["phase_A"].as_u64(), b["phase_B"].as_u64()), (Some(5), Some(5)), "{b}");
    }
}

#[test]
fn genes_say_how_many_there_are() {
    let d = Dir::new("genes");
    let gtf = |n: usize| {
        let p = d.file(&format!("{n}.gtf"));
        let mut text = String::new();
        for i in 0..n {
            let attrs = format!("gene_id \"g{i:03}\"; transcript_id \"t{i:03}\"; gene_name \"g{i:03}\";");
            for f in ["transcript", "exon"] {
                text += &format!("chrQ\tqa\t{f}\t1\t100\t.\t+\t.\t{attrs}\n");
            }
        }
        std::fs::write(&p, text).unwrap();
        p
    };
    let (hundred, more) = (gtf(100), gtf(101));
    let mut srv = server(Defaults::default());
    for all in [true, false] {
        let out = call(&mut srv, "genes", json!({"region": "chrQ:1-100", "all": all, "genes": [s(&hundred)]})).unwrap();
        let f = &out["results"][0];
        assert_eq!((f["total"].as_u64(), f["truncated"].as_bool(), f["transcripts"].as_array().unwrap().len()), (Some(100), Some(false), 100), "{all}");
        let out = call(&mut srv, "genes", json!({"region": "chrQ:1-100", "all": all, "genes": [s(&more)]})).unwrap();
        let f = &out["results"][0];
        assert_eq!((f["total"].as_u64(), f["truncated"].as_bool(), f["transcripts"].as_array().unwrap().len()), (Some(101), Some(true), 100), "{all}");
        let out = call(&mut srv, "genes", json!({"region": "chrQ:1-100", "all": all, "max": 500, "genes": [s(&more)]})).unwrap();
        let f = &out["results"][0];
        assert_eq!((f["truncated"].as_bool(), f["transcripts"].as_array().unwrap().len()), (Some(false), 101), "{all}");
    }
    let out = call(&mut srv, "genes", json!({"name": "g100", "all": true, "genes": [s(&more)]})).unwrap();
    assert_eq!((out["results"][0]["total"].as_u64(), out["results"][0]["truncated"].as_bool()), (Some(1), Some(false)));
}

/// A bgzipped, tabix-indexed VCF of `body` on `contig`.
fn vcf(d: &Dir, name: &str, contig: (&str, u64), body: &str) -> PathBuf {
    use std::io::Write;
    let p = d.file(name);
    let text = format!("##fileformat=VCFv4.2\n##contig=<ID={},length={}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{body}", contig.0, contig.1);
    let mut w = ggv_io::pasta::bgzf::io::Writer::new(std::fs::File::create(&p).unwrap());
    w.write_all(text.as_bytes()).unwrap();
    w.finish().unwrap();
    ggv_io::build_index(&p).unwrap();
    p
}

/// GTF rows for one transcript: its span as one exon, coding from `cds`.
fn gtf_tx(contig: &str, gene: &str, tx: &str, (s, e): (u64, u64), cds: Option<(u64, u64)>) -> String {
    let attrs = format!("gene_id \"{gene}\"; transcript_id \"{tx}\"; gene_name \"{gene}\";");
    let mut t = format!("{contig}\tqa\ttranscript\t{s}\t{e}\t.\t+\t.\t{attrs}\n{contig}\tqa\texon\t{s}\t{e}\t.\t+\t.\t{attrs}\n");
    if let Some((cs, ce)) = cds {
        t += &format!("{contig}\tqa\tCDS\t{cs}\t{ce}\t.\t+\t0\t{attrs}\n");
    }
    t
}

#[test]
fn variants_and_genes_resolve_gene_names_and_refuse_unknown_contigs() {
    let d = Dir::new("resolve");
    let calls = vcf(&d, "calls.vcf.gz", ("chrQ", 1000), "chrQ\t150\t.\tA\tC\t.\tPASS\t.\nchrQ\t700\t.\tA\tG\t.\tPASS\t.\n");
    let gtf = d.file("g.gtf");
    std::fs::write(&gtf, gtf_tx("chrQ", "GQ", "TQ", (101, 200), None)).unwrap();
    // No BAM or reference: the VCF's and the annotation's own contigs.
    let mut srv = server(Defaults { vcfs: vec![calls], genes: vec![gtf], ..Default::default() });
    let out = call(&mut srv, "variants", json!({"region": "GQ"})).unwrap();
    assert_eq!((out["results"][0]["records"].as_u64(), out["query"]["searched"].as_str()), (Some(1), Some("chrQ:101-200")), "{out}");
    assert_eq!(call(&mut srv, "variants", json!({"region": "chrQ"})).unwrap()["results"][0]["records"], 2);
    for region in ["chrO:1-100", "NOPE"] {
        let e = call(&mut srv, "variants", json!({"region": region})).unwrap_err();
        assert!(e.contains(region.split(':').next().unwrap()), "{region}: {e}");
        assert!(call(&mut srv, "genes", json!({"region": region})).is_err(), "{region}");
    }
    assert_eq!(call(&mut srv, "genes", json!({"region": "GQ"})).unwrap()["results"][0]["total"], 1);
    // Past the last gene is empty, not out of bounds.
    assert_eq!(call(&mut srv, "genes", json!({"region": "chrQ:500-900"})).unwrap()["results"][0]["total"], 0);
}

#[test]
fn genes_by_transcript_id_and_on_several_contigs() {
    let d = Dir::new("genes-by-id");
    let gtf = d.file("g.gtf");
    let mut text = gtf_tx("chr7", "BRAF", "NM_001374258.1", (101, 400), Some((101, 399)));
    text += &gtf_tx("chr7", "BRAF", "NM_004333.6", (101, 400), Some((101, 300)));
    text += &gtf_tx("chrX", "SHOX", "NM_000451.4_2", (51, 150), Some((51, 149)));
    text += &gtf_tx("chrY", "SHOX", "NM_000451.4", (51, 150), Some((51, 149)));
    std::fs::write(&gtf, text).unwrap();
    let mut srv = server(Defaults { genes: vec![gtf], ..Default::default() });
    let ids = |out: &Value| out["results"][0]["transcripts"].as_array().unwrap().iter().map(|t| format!("{}:{}", t["contig"].as_str().unwrap(), t["transcript"].as_str().unwrap())).collect::<Vec<_>>();
    // The main transcript is the MANE Select one, and a transcript id gets that transcript.
    assert_eq!(ids(&call(&mut srv, "genes", json!({"name": "BRAF"})).unwrap()), ["chr7:NM_004333.6"]);
    for name in ["NM_001374258.1", "nm_001374258", "NM_004333.6"] {
        let out = call(&mut srv, "genes", json!({"name": name})).unwrap();
        assert_eq!(ids(&out), [format!("chr7:{}", if name.contains("4333") { "NM_004333.6" } else { "NM_001374258.1" })], "{name}: {out}");
    }
    let out = call(&mut srv, "genes", json!({"name": "SHOX"})).unwrap();
    // UCSC's copy suffix (`_2`) isn't an accession: the contig tells the copies apart.
    assert_eq!(ids(&out), ["chrX:NM_000451.4", "chrY:NM_000451.4"], "{out}");
    assert_eq!(out["results"][0]["transcripts"][0]["mane_select"], true, "{out}");
    assert!(out["results"][0]["note"].as_str().unwrap().contains("2 contigs"), "{out}");
    // A version the file doesn't have: the one it has, said so.
    let out = call(&mut srv, "genes", json!({"name": "NM_004333.99"})).unwrap();
    assert_eq!(ids(&out), ["chr7:NM_004333.6"], "{out}");
    assert!(out["results"][0]["note"].as_str().unwrap().contains("NM_004333.99 is not in g.gtf: it has NM_004333.6"), "{out}");
}

#[test]
fn pseudo_autosomal_gene_regions_are_on_chrx_whatever_the_row_order() {
    let d = Dir::new("par");
    let gtf = d.file("g.gtf");
    // UCSC's order: chrY first, the plain id there.
    let text = gtf_tx("chrY", "SHOX", "NM_000451.4", (51, 150), Some((51, 149))) + &gtf_tx("chrX", "SHOX", "NM_000451.4_2", (61, 160), Some((61, 159)));
    std::fs::write(&gtf, text).unwrap();
    let reference = fasta(&d, "r.fa", &[("chrY", &bases(300)), ("chrX", &bases(300))]);
    let mut srv = server(Defaults { genes: vec![gtf], reference: Some(reference), ..Default::default() });
    for name in ["SHOX", "NM_000451.4"] {
        let out = call(&mut srv, "genes", json!({"region": name})).unwrap();
        let t = &out["results"][0]["transcripts"][0];
        assert_eq!((out["results"][0]["contig"].as_str(), t["start"].as_u64()), (Some("chrX"), Some(61)), "{name}: {out}");
    }
}

#[test]
fn open_is_idempotent_and_compare_says_when_a_and_b_are_one_file() {
    let d = Dir::new("open-twice");
    let r = bases(300);
    let one = bam(&d, "one.bam", &[("chrQA", 300)], &reads("r", &r, "chrQA", 101, 30, 9, None));
    let vcf = vcf(&d, "calls.vcf.gz", ("chrQA", 300), "chrQA\t110\t.\tA\tC\t.\tPASS\t.\n");
    let mut srv = server(Defaults { bams: vec![one.clone()], ..Default::default() });
    let other_spelling = d.0.join(".").join("one.bam");
    for files in [json!([s(&one), s(&vcf)]), json!([s(&one), s(&vcf)]), json!([s(&other_spelling)])] {
        let out = call(&mut srv, "open", json!({"files": files})).unwrap();
        assert_eq!((out["defaults"]["bams"].as_array().unwrap().len(), out["defaults"]["vcfs"].as_array().unwrap().len()), (1, 1), "{out}");
    }
    let out = call(&mut srv, "compare", json!({"region": "chrQA:101-130", "a": s(&one), "b": s(&other_spelling), "min_depth": 1})).unwrap();
    assert!(out["warning"].as_str().unwrap().contains("same file"), "{out}");
}

#[test]
fn compare_warns_about_a_truncated_file() {
    let d = Dir::new("compare-trunc");
    let r = bases(300);
    let good = bam(&d, "good.bam", &[("chrQA", 300)], &reads("r", &r, "chrQA", 101, 30, 9, None));
    let cut = d.file("cut.bam");
    let bytes = std::fs::read(&good).unwrap();
    std::fs::write(&cut, &bytes[..bytes.len() - 28]).unwrap();
    std::fs::copy(good.with_added_extension("bai"), cut.with_added_extension("bai")).unwrap();
    // Also with a default BAM the call doesn't read: only the pair is warned about.
    let elsewhere = bam(&d, "elsewhere.bam", &[("chrQA", 300)], &reads("e", &r, "chrQA", 101, 30, 9, None));
    let mut srv = server(Defaults { bams: vec![elsewhere], ..Default::default() });
    let out = call(&mut srv, "compare", json!({"region": "chrQA:101-130", "a": s(&good), "b": s(&cut), "min_depth": 1})).unwrap();
    let w = out["file_warnings"].as_array().unwrap_or_else(|| panic!("{out}"));
    assert!(w.len() == 1 && w[0].as_str().unwrap().contains("cut.bam: no BGZF end-of-file marker"), "{out}");
}

#[test]
fn files_of_different_assemblies_are_warned_about() {
    let d = Dir::new("assembly");
    let r = bases(600);
    let fa = fasta(&d, "ref.fa", &[("chr1", &r)]);
    // chr1's length says hg19; the reference (600 bp) and the hg38 genes disagree.
    let hg19 = bam(&d, "sample.bam", &[("chr1", 249_250_621)], &reads("r", &r, "chr1", 101, 30, 9, None));
    let same = bam(&d, "same.bam", &[("chr1", 600)], &reads("r", &r, "chr1", 101, 30, 9, None));
    let genes = d.file("hg38.ncbiRefSeq.gtf");
    std::fs::write(&genes, gtf_tx("chr1", "G", "NM_1.1", (101, 200), None)).unwrap();
    let warnings = |srv: &mut Server, args: Value| -> Vec<String> {
        let out = call(srv, "coverage", args).unwrap();
        out["file_warnings"].as_array().map(|w| w.iter().map(|x| x.as_str().unwrap().to_owned()).collect()).unwrap_or_default()
    };
    let mut srv = server(Defaults { bams: vec![hg19.clone()], reference: Some(fa.clone()), ..Default::default() });
    let w = warnings(&mut srv, json!({"region": "chr1:101-130"}));
    assert!(w.len() == 1 && w[0].contains("assembly mismatch") && w[0].contains("ref.fa and sample.bam are different assemblies"), "{w:?}");
    let mut srv = server(Defaults { bams: vec![hg19], genes: vec![genes.clone()], ..Default::default() });
    let w = warnings(&mut srv, json!({"region": "chr1:101-130"}));
    assert!(w.len() == 1 && w[0].contains("hg38.ncbiRefSeq.gtf is hg38"), "{w:?}");
    // Files that agree: nothing.
    let mut srv = server(Defaults { bams: vec![same], reference: Some(fa), ..Default::default() });
    assert!(warnings(&mut srv, json!({"region": "chr1:101-130"})).is_empty());
}

#[test]
fn consequence_refuses_alt_equal_to_ref() {
    let d = Dir::new("ref-is-alt");
    let gtf = d.file("g.gtf");
    std::fs::write(&gtf, gtf_tx("chrQ", "G", "T", (51, 150), Some((51, 149)))).unwrap();
    let mut srv = server(Defaults { genes: vec![gtf], ..Default::default() });
    for (r, a) in [("A", "A"), ("ac", "AC")] {
        let e = call(&mut srv, "consequence", json!({"locus": "chrQ:101", "ref": r, "alt": a})).unwrap_err();
        assert!(e.contains("not a change"), "{e}");
    }
    assert!(call(&mut srv, "consequence", json!({"locus": "chrQ:101", "ref": "A", "alt": "C"})).is_ok());
    // Several ALTs: said so, not called symbolic.
    let e = call(&mut srv, "consequence", json!({"locus": "chrQ:101", "ref": "A", "alt": "C,T"})).unwrap_err();
    assert!(e.contains("multi-allelic") && !e.contains("symbolic"), "{e}");
    // A range must be REF's span.
    let e = call(&mut srv, "consequence", json!({"locus": "chrQ:101-110", "ref": "A", "alt": "C"})).unwrap_err();
    assert!(e.contains("10 bases") && e.contains("has 1"), "{e}");
    assert!(call(&mut srv, "consequence", json!({"locus": "chrQ:101-102", "ref": "AC", "alt": "A"})).is_ok());
}

#[test]
fn transcript_ids_resolve_to_the_primary_copy_and_are_named_without_the_ucsc_suffix() {
    let d = Dir::new("ucsc-copies");
    let r = bases(400);
    let fa = fasta(&d, "ref.fa", &[("chr10", &r), ("chr10_KQ090021v1_fix", &r)]);
    let gtf = d.file("g.gtf");
    let mut text = gtf_tx("chr10_KQ090021v1_fix", "PTEN", "NM_000314.8", (101, 300), Some((101, 299)));
    text += &gtf_tx("chr10", "PTEN", "NM_000314.8_2", (101, 300), Some((101, 299)));
    std::fs::write(&gtf, text).unwrap();
    let one = bam(&d, "one.bam", &[("chr10", 400)], &reads("r", &r, "chr10", 101, 30, 5, None));
    let mut srv = server(Defaults { bams: vec![one], reference: Some(fa), genes: vec![gtf], ..Default::default() });
    let out = call(&mut srv, "genes", json!({"name": "NM_000314.8"})).unwrap();
    let t = &out["results"][0]["transcripts"][0];
    assert_eq!((t["contig"].as_str(), t["transcript"].as_str()), (Some("chr10"), Some("NM_000314.8")), "{out}");
    // As a region.
    let out = call(&mut srv, "coverage", json!({"region": "NM_000314.8", "bins": 1})).unwrap();
    assert!(out.to_string().contains("chr10:101-300"), "{out}");
    // HGVS names the accession.
    let base = r.as_bytes()[149] as char;
    let out = call(&mut srv, "consequence", json!({"locus": "chr10:150", "ref": base.to_string(), "alt": other(base).to_string()})).unwrap();
    assert_eq!(out["consequences"][0]["transcript"], "NM_000314.8", "{out}");
}

#[test]
fn info_and_open_show_signal_tracks_and_signal_bins_cover_every_base() {
    let d = Dir::new("signal");
    let bg = d.file("t.bedGraph");
    std::fs::write(&bg, "chr1\t0\t100\t5\nchr1\t100\t103\t7\nchr2\t0\t1000\t1\n").unwrap();
    let mut srv = server(Defaults { signals: vec![bg.clone()], ..Default::default() });
    let out = call(&mut srv, "info", json!({})).unwrap();
    assert_eq!(out["defaults"]["signals"][0], s(&bg), "{out}");
    assert_eq!((out["files"][0]["kind"].as_str(), out["files"][0]["contigs"].as_u64()), (Some("signal"), Some(2)), "{out}");
    let mut srv = server(Defaults::default());
    let out = call(&mut srv, "open", json!({"files": [s(&bg)]})).unwrap();
    assert_eq!(out["defaults"]["signals"][0], s(&bg), "{out}");
    // Six bins over three bases: three of one base each, none empty.
    let out = call(&mut srv, "signal", json!({"region": "chr1:101-103", "bins": 6})).unwrap();
    let f = &out["results"][0];
    assert_eq!((f["bins"].as_u64(), f["bin_size"].as_f64()), (Some(3), Some(1.0)), "{out}");
    assert_eq!(f["rows"][0]["mean"], json!([7.0, 7.0, 7.0]), "{out}");
    // Past the contig's end (a text file doesn't know it; the reference does): cut there, said so.
    let reference = fasta(&d, "r.fa", &[("chr1", &bases(150)), ("chr2", &bases(1000))]);
    let out = call(&mut srv, "signal", json!({"region": "chr1:101-300", "bins": 10, "reference": s(&reference)})).unwrap();
    let f = &out["results"][0];
    assert_eq!((f["region"].as_str(), f["bin_size"].as_f64()), (Some("chr1:101-150"), Some(5.0)), "{out}");
    assert!(f["note"].as_str().unwrap().contains("chr1 ends at 150"), "{out}");
    assert!(call(&mut srv, "signal", json!({"region": "chr1:201-300", "reference": s(&reference)})).unwrap()["results"][0]["error"].is_string());
}

#[test]
fn render_stops_at_the_contig_end() {
    let d = Dir::new("render-end");
    let r = bases(300);
    let one = bam(&d, "one.bam", &[("chrQA", 300)], &reads("r", &r, "chrQA", 251, 50, 5, None));
    let mut srv = server(Defaults { bams: vec![one], ..Default::default() });
    let msg = json!({"jsonrpc": "2.0", "id": 1, "method": "tools/call", "params": {"name": "render", "arguments": {"locus": "chrQA:290"}}});
    let reply = srv.handle(&msg.to_string()).unwrap();
    assert!(reply.contains("chrQA:240-300 (61 bp)"), "{reply}");
}

#[test]
fn svs_tests_duplications_and_breakend_pairs_as_the_reads_carry_them() {
    // The fixtures' own events, written as other callers write them (see the pbsv records):
    // a 6.8 kb deletion as <DEL> and as a BND pair; a 586 bp tandem dup (pbsv.INS.1051) as
    // <DUP>; a 94 bp insertion (pbsv.INS.574) as <INS> without SVLEN.
    let data = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio");
    if !data.join("HG002.revio.slices.bam").exists() {
        return;
    }
    let d = Dir::new("svs-dup-bnd");
    let calls = d.file("syn.vcf.gz");
    {
        use std::io::Write;
        let text = "##fileformat=VCFv4.2\n##contig=<ID=chr20,length=64444167>\n\
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"t\">\n##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">\n\
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"l\">\n##INFO=<ID=MATEID,Number=.,Type=String,Description=\"m\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr20\t31170827\tsymDEL\tT\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=31177629;SVLEN=-6802\n\
chr20\t31170827\tbnd1\tT\tT[chr20:31177630[\t.\tPASS\tSVTYPE=BND;MATEID=bnd2\n\
chr20\t31177630\tbnd2\tG\t]chr20:31170827]G\t.\tPASS\tSVTYPE=BND;MATEID=bnd1\n\
chr20\t31191350\tinsNoLen\tG\t<INS>\t.\tPASS\tSVTYPE=INS;END=31191350\n\
chr20\t58830532\tdupT\tA\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=58831118;SVLEN=586\n";
        let mut w = ggv_io::pasta::bgzf::io::Writer::new(std::fs::File::create(&calls).unwrap());
        w.write_all(text.as_bytes()).unwrap();
        w.finish().unwrap();
        ggv_io::build_index(&calls).unwrap();
    }
    let mut srv = server(Defaults { bams: vec![data.join("HG002.revio.slices.bam")], vcfs: vec![calls], ..Default::default() });
    let rec = |out: &Value, id: &str| out["records"].as_array().unwrap().iter().find(|r| r["id"] == id).cloned().unwrap_or_else(|| panic!("{id}: {out}"));
    let reads = |r: &Value| (r["reads"][0]["supporting"].as_u64().unwrap(), r["reads"][0]["reference"].as_u64().unwrap(), r["reads"][0]["genotype"].as_str().unwrap().to_owned());

    let out = call(&mut srv, "svs", json!({"region": "chr20:31160000-31195000", "min_len": 20})).unwrap();
    let (del, bnd) = (rec(&out, "symDEL"), rec(&out, "bnd1"));
    // The pair is tested once, as the deletion it is, with the same reads as <DEL>.
    assert_eq!(reads(&bnd), reads(&del), "{bnd}");
    assert_eq!((bnd["tested_as"].as_str(), bnd["mate"]["id"].as_str(), bnd["pos"].as_u64()), (Some("deletion"), Some("bnd2"), Some(31170827)), "{bnd}");
    assert!(out["records"].as_array().unwrap().iter().all(|r| r["id"] != "bnd2"), "{out}");
    let (sup, reference, gt) = reads(&del);
    assert!(sup >= 10 && reference <= sup + 2 && gt == "0/1", "{del}");
    // The unsized insertion is named, not dropped.
    assert!(out["warnings"].as_array().unwrap().iter().any(|w| w.as_str().unwrap().contains("31191350") && w.as_str().unwrap().contains("no SVLEN")), "{out}");
    assert_eq!(out["records_in_region"], 4, "{out}");

    let out = call(&mut srv, "svs", json!({"region": "chr20:58820000-58860000", "min_len": 20})).unwrap();
    let dup = rec(&out, "dupT");
    let (sup, reference, gt) = reads(&dup);
    assert!(sup >= 15 && reference <= 2 && gt == "1/1", "{dup}");
    // Its insertion is explained by it.
    assert!(out["unexplained"].as_array().unwrap().iter().all(|e| e["pos"].as_u64().unwrap().abs_diff(58830533) > 700), "{out}");
}

#[test]
fn open_names_a_missing_fasta_index() {
    let d = Dir::new("noindex");
    let fa = d.file("small.fa");
    std::fs::write(&fa, ">c\nACGT\n").unwrap();
    let mut srv = server(Defaults::default());
    let e = call(&mut srv, "open", json!({"files": [s(&fa)]})).unwrap_err();
    assert!(e.contains("no FASTA index (looked for small.fa.fai)") && !e.contains("cannot find"), "{e}");
    assert_eq!(e.matches("small.fa:").count(), 1, "the path once: {e}");
}

#[test]
fn svs_counts_a_smaller_deletion_as_the_reads_allele() {
    // pbsv.DEL.568: a 20 bp deletion in a GGAAT repeat, 1/1 (HP2 reads delete 20 bp, HP1 reads
    // 15 bp). The 15 bp reads were counted as reference, making it 0/1.
    let data = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio");
    if !data.join("HG002.revio.slices.bam").exists() {
        return;
    }
    let mut srv = server(Defaults {
        bams: vec![data.join("HG002.revio.slices.bam")],
        reference: Some(data.join("GRCh38.chr20.slices.fa.gz")),
        vcfs: vec![data.join("HG002.revio.pbsv.vcf.gz")],
        ..Default::default()
    });
    let out = call(&mut srv, "svs", json!({"region": "chr20:31,178,000-31,179,000", "min_len": 20})).unwrap();
    let text = out.to_string();
    let rec = out["records"].as_array().unwrap_or_else(|| panic!("{text}")).iter().find(|r| r["id"] == "pbsv.DEL.568").unwrap_or_else(|| panic!("{text}"));
    let reads = &rec["reads"][0];
    // 15 bp is within the 70-130% a supporting event may be, so those reads support it.
    assert_eq!((reads["genotype"].as_str(), reads["supporting"].as_u64(), reads["reference"].as_u64()), (Some("1/1"), Some(28), Some(4)), "{rec}");
    assert_eq!(reads["by_haplotype"]["reference"]["1"], 0, "{rec}");
}

#[test]
fn answers_without_a_reference_say_what_that_changes() {
    let d = Dir::new("no-ref");
    let r = bases(300);
    let fa = fasta(&d, "ref.fa", &[("chrQA", &r)]);
    // Every read carries another base at 120: a homozygous change.
    let one = bam(&d, "hom.bam", &[("chrQA", 300)], &reads("r", &r, "chrQA", 101, 50, 10, Some((10, 120))));
    let mut srv = server(Defaults { bams: vec![one], ..Default::default() });
    let bare = call(&mut srv, "evidence", json!({"locus": "chrQA:120"})).unwrap();
    assert!(bare["results"][0]["ref"].is_null() && bare["results"][0]["no_reference"].as_str().unwrap().contains("most common one may be it"), "{bare}");
    let summary = call(&mut srv, "summarize", json!({"region": "chrQA:101-150"})).unwrap();
    assert_eq!(summary["results"][0]["variable_site_count"], 0, "{summary}");
    assert!(summary["no_reference"].as_str().unwrap().contains("homozygous variants are not among them"), "{summary}");
    // With the reference, the change is found and nothing is noted.
    let with = call(&mut srv, "evidence", json!({"locus": "chrQA:120", "reference": s(&fa)})).unwrap();
    assert!(with["results"][0].get("no_reference").is_none() && with["results"][0]["ref"].is_string(), "{with}");
    let summary = call(&mut srv, "summarize", json!({"region": "chrQA:101-150", "reference": s(&fa)})).unwrap();
    assert!(summary.get("no_reference").is_none() && summary["results"][0]["variable_site_count"] == 1, "{summary}");
}

#[test]
fn a_misspelt_argument_is_refused_not_dropped() {
    let d = Dir::new("typo");
    let r = bases(300);
    let one = bam(&d, "a.bam", &[("chrQA", 300)], &reads("r", &r, "chrQA", 101, 50, 10, None));
    let mut srv = server(Defaults { bams: vec![one.clone()], ..Default::default() });
    // The misspelling used to be dropped: the call ran on the default and said nothing.
    let e = call(&mut srv, "summarize", json!({"region": "chrQA:101-150", "fliter": "mapq >= 60"})).unwrap_err();
    assert!(e.contains("`fliter` is not an argument") && e.contains("did you mean `filter`?"), "{e}");
    // The files every tool finds its contigs in are arguments of every tool.
    for tool in ["summarize", "coverage"] {
        assert!(call(&mut srv, tool, json!({"region": "chrQA:101-150", "bam": s(&one)})).is_ok(), "{tool} takes `bam`");
    }
}

#[test]
fn phase_lists_as_many_sites_as_asked_and_counts_them_all() {
    let d = Dir::new("phase-max");
    let r = bases(400);
    let fa = fasta(&d, "ref.fa", &[("chrQA", &r)]);
    let mut all = Vec::new();
    for k in 0..10 {
        let mut seq: Vec<char> = r[100..140].chars().collect();
        if k < 5 {
            for at in [9, 19, 29] {
                seq[at] = other(seq[at]); // 110, 120, 130
            }
        }
        all.push(read(&format!("x{k}"), "chrQA", 101, "40M", &seq.iter().collect::<String>()));
    }
    let one = bam(&d, "sites.bam", &[("chrQA", 400)], &all);
    let mut srv = server(Defaults { bams: vec![one], reference: Some(fa), ..Default::default() });
    let out = call(&mut srv, "phase", json!({"region": "chrQA:100-150", "max_sites": 2})).unwrap();
    assert_eq!((out["het_sites"].as_u64(), out["sites_shown"].as_u64(), out["sites"].as_array().map(Vec::len)), (Some(3), Some(2), Some(2)), "{out}");
    assert!(out["sites_truncated"].as_str().unwrap().contains("2 of 3"), "{out}");
    assert_eq!(out["phase_blocks"][0]["sites"], 3, "blocks count every site: {out}");
    let all_sites = call(&mut srv, "phase", json!({"region": "chrQA:100-150"})).unwrap();
    assert!(all_sites["sites_truncated"].is_null() && all_sites["sites_shown"] == 3, "{all_sites}");
}

#[test]
fn long_alleles_are_cut_and_qual_is_not_widened() {
    let d = Dir::new("long-allele");
    let deleted = bases(200);
    let body = format!("chrQA\t100\tsv1\t{deleted}\t{}\t45.1\tPASS\t.\nchrQA\t400\t.\tA\tC\t33.3\tPASS\t.\n", &deleted[..1]);
    let calls = vcf(&d, "calls.vcf.gz", ("chrQA", 1000), &body);
    let mut srv = server(Defaults { vcfs: vec![calls], ..Default::default() });
    let out = call(&mut srv, "variants", json!({"region": "chrQA:1-1000"})).unwrap();
    let (sv, snv) = (&out["results"][0]["shown"][0], &out["results"][0]["shown"][1]);
    assert_eq!(sv["ref"], format!("{}… (200 bp)", &deleted[..20]), "{sv}");
    assert_eq!((sv["ref_length"].as_u64(), sv["alt_lengths"][0].as_u64(), sv["alt"][0].as_str()), (Some(200), Some(1), Some(&deleted[..1])), "{sv}");
    assert_eq!(sv["qual"], 45.1, "the file's float, not 45.099998…: {sv}");
    assert!(snv.get("ref_length").is_none() && snv["ref"] == "A" && snv["qual"] == 33.3, "{snv}");
}

#[test]
fn compare_counts_each_file_under_its_own_filters() {
    let d = Dir::new("compare-own");
    let r = bases(400);
    let fa = fasta(&d, "ref.fa", &[("chrQA", &r)]);
    // Ten reads from 100, half with another base at 120; ten more from 101, all reference.
    let mut all = reads("alt", &r, "chrQA", 100, 50, 10, Some((5, 120)));
    all.extend(reads("ref", &r, "chrQA", 101, 50, 10, None));
    let one = bam(&d, "mixed.bam", &[("chrQA", 400)], &all);
    let mut srv = server(Defaults { reference: Some(fa), ..Default::default() });
    let same = call(&mut srv, "compare", json!({"region": "chrQA:100-149", "a": s(&one), "b": s(&one)})).unwrap();
    assert_eq!(same["differences"], 0, "{same}");
    assert!(same["warning"].as_str().unwrap().contains("same file"), "{same}");
    // The same file split by where its reads start: the alt is in one half only.
    let split = call(&mut srv, "compare", json!({"region": "chrQA:100-149", "a": s(&one), "b": s(&one), "a_filters": {"filter": "start == 100"}, "b_filters": {"filter": "start == 101"}})).unwrap();
    assert_eq!((split["differences"].as_u64(), split["shown"][0]["pos"].as_u64()), (Some(1), Some(120)), "{split}");
    assert_eq!((split["shown"][0]["a_fraction"].as_f64(), split["shown"][0]["b_fraction"].as_f64()), (Some(0.5), Some(0.0)), "{split}");
    assert_eq!(split["query"]["filters"]["a"]["expression"], "start == 100", "{split}");
    assert!(split["warning"].is_null(), "one file under two filters is a real comparison: {split}");
}

#[test]
fn repeats_warn_when_the_region_is_narrower_than_the_tract() {
    let d = Dir::new("repeat-narrow");
    let r = bases(400);
    // Eight reads delete 20 bases at 150-169, which the aligner might as well have put at
    // 170-189: asked for 171-185 alone, they have no aligned left flank.
    let deleting: Vec<String> = (0..8).map(|k| read(&format!("d{k}"), "chrQA", 101, "49M20D51M", &format!("{}{}", &r[100..149], &r[169..220]))).collect();
    let whole = reads("w", &r, "chrQA", 101, 120, 2, None);
    let one = bam(&d, "tract.bam", &[("chrQA", 400)], &[deleting, whole].concat());
    let mut srv = server(Defaults { bams: vec![one], ..Default::default() });
    let narrow = call(&mut srv, "repeats", json!({"region": "chrQA:171-185"})).unwrap();
    let f = &narrow["files"][0];
    assert_eq!((f["spanning"].as_u64(), f["partial"].as_u64()), (Some(2), Some(8)), "{narrow}");
    assert!(f["warning"].as_str().unwrap().contains("give the whole tract"), "{narrow}");
    let tract = call(&mut srv, "repeats", json!({"region": "chrQA:140-195"})).unwrap();
    assert!(tract["files"][0]["warning"].is_null() && tract["files"][0]["spanning"] == 10, "{tract}");
    assert_eq!(tract["files"][0]["alleles"].as_array().map(Vec::len), Some(2), "{tract}");
}
