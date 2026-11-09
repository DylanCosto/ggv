//! How big one decoded read span is: `span_bytes <bam> <fasta> <locus>` (`MODS=0` leaves out
//! the base-modification calls, as the app does unless colouring by them).
use std::time::Instant;
use ggv_layout::{AlignedRead, ViewOptions, load_view};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let mut bam = ggv_io::BamFile::open_with_reference(&a[1], Some(&a[2])).unwrap();
    let mut fasta = ggv_io::ReferenceFasta::open(&a[2]).unwrap();
    let region = bam.resolve_str(&a[3]).unwrap();
    let seq = fasta.fetch(&region).unwrap();
    let t = Instant::now();
    let opts = ViewOptions { max_read_window: 20_000_000, mods: std::env::var("MODS").as_deref() != Ok("0"), ..Default::default() };
    let view = load_view(&mut bam, &region, Some(seq), String::new(), &opts).unwrap();
    let rs: &[AlignedRead] = &view.reads.as_ref().unwrap().0;
    let n = |f: &dyn Fn(&AlignedRead) -> usize| rs.iter().map(f).sum::<usize>();
    println!("{} reads in {:.1}s", rs.len(), t.elapsed().as_secs_f64());
    println!("blocks {} dels {} ins {} mism {} clipmm {} mods {} lowq {} sa {}",
        n(&|r| r.blocks.len()), n(&|r| r.deletions.len()), n(&|r| r.insertions.len()), n(&|r| r.mismatches.len()),
        n(&|r| r.clip_mismatches.len()), n(&|r| r.mods.len()), n(&|r| r.low_quality.len()), n(&|r| r.supplementary.len()));
    println!("spare capacity: mism {} mods {} blocks {}", n(&|r| r.mismatches.capacity() - r.mismatches.len()),
        n(&|r| r.mods.capacity() - r.mods.len()), n(&|r| r.blocks.capacity() - r.blocks.len()));
    breakdown(rs);
    println!("ModCall {} Mismatch {} Block {} AlignedRead {}", size_of::<ggv_layout::ModCall>(), size_of::<ggv_layout::Mismatch>(), size_of::<ggv_layout::Block>(), size_of::<AlignedRead>());
}

fn breakdown(rs: &[AlignedRead]) {
    let mut by = std::collections::BTreeMap::<(char, char), [u64; 4]>::new();
    for r in rs { for m in &r.mods { by.entry((m.code as char, m.canonical as char)).or_default()[usize::from(m.prob / 64)] += 1; } }
    for (k, v) in by { println!("{k:?} by ML quarter {v:?}"); }
    let bases: u64 = rs.iter().map(|r| r.end - r.start).sum();
    println!("aligned bases {bases}");
}
