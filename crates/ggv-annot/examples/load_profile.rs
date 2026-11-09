//! Times loading an annotation and reports memory: `cargo run --release -p ggv-annot --example load_profile -- genes.gtf.gz`.
use std::time::Instant;

fn rss_mb() -> f64 {
    let out = std::process::Command::new("ps").args(["-o", "rss=", "-p", &std::process::id().to_string()]).output().unwrap();
    String::from_utf8_lossy(&out.stdout).trim().parse::<f64>().unwrap_or(0.0) / 1024.0
}

fn main() {
    let path = std::env::args().nth(1).expect("path");
    let before = rss_mb();
    let t = Instant::now();
    let a = ggv_annot::Annotation::open(&path).unwrap();
    let load = t.elapsed();
    println!("open: {:.2?}  cached {}  transcripts {}  genes {}  rss {:.0} MB (was {:.0})", load, a.is_cached(), a.len(), a.gene_count(), rss_mb(), before);
    let t = Instant::now();
    let n17 = a.query("chr17", 0, u64::MAX).count();
    println!("first chr17 query (decodes the contig): {:.2?} ({n17} transcripts), rss {:.0} MB", t.elapsed(), rss_mb());
    let t = Instant::now();
    let exons: usize = a.query("chr1", 0, u64::MAX).count();
    println!("first chr1 query: {:.2?}, rss {:.0} MB", t.elapsed(), rss_mb());
    let t = Instant::now();
    for _ in 0..1000 {
        std::hint::black_box(a.query("chr17", 43_000_000, 43_200_000).count());
    }
    println!("query x1000: {:.2?} ({} chr1 transcripts)", t.elapsed(), exons);
    let t = Instant::now();
    for g in ["BRCA1", "TP53", "ENSG00000139618", "ENST00000003084.11"] {
        std::hint::black_box(a.find(g));
    }
    println!("find x4: {:.2?}", t.elapsed());
}
