//! Splitting a long region across threads: each chunk gets its own handle on the file, and the
//! caller adds the parts up. A read is fetched by every chunk it overlaps, so passes count it
//! only in the chunk where it starts (`count_from`; the first chunk counts every read it meets,
//! as a single pass over the whole region would).

use ggv_io::BamFile;
use ggv_model::Region;

use crate::Result;

/// Regions at least twice this long are split.
pub const CHUNK: u64 = 500_000;

/// Runs `pass(handle, chunk, count_from)` over `region`: once on `bam` when the region is short
/// (or the file can't be reopened), else on up to 8 chunks in parallel. Results are in order.
pub fn run<T: Send>(bam: &mut BamFile, region: &Region, pass: impl Fn(&mut BamFile, &Region, u64) -> Result<T> + Sync) -> Result<Vec<T>> {
    run_sized(bam, region, CHUNK, pass)
}

/// [`run`] with chunks of at least `chunk` bases (regions at least twice that are split).
pub fn run_sized<T: Send>(bam: &mut BamFile, region: &Region, chunk: u64, pass: impl Fn(&mut BamFile, &Region, u64) -> Result<T> + Sync) -> Result<Vec<T>> {
    let threads = std::thread::available_parallelism().map_or(1, |n| n.get()).min(8) as u64;
    let n = (region.len() / chunk.max(1)).clamp(1, threads);
    let handles: Option<Vec<BamFile>> = (n > 1).then(|| (0..n).map(|_| bam.reopen().ok()).collect()).flatten();
    let Some(mut handles) = handles else { return Ok(vec![pass(bam, region, 0)?]) };
    let step = region.len().div_ceil(n);
    let chunks: Vec<Region> = (0..n)
        .map(|i| Region { contig: region.contig.clone(), start: region.start + i * step, end: (region.start + (i + 1) * step).min(region.end) })
        .filter(|c| c.start < c.end)
        .collect();
    let pass = &pass;
    std::thread::scope(|scope| {
        let jobs: Vec<_> = chunks
            .iter()
            .zip(handles.iter_mut())
            .enumerate()
            .map(|(i, (chunk, b))| scope.spawn(move || pass(b, chunk, if i == 0 { 0 } else { chunk.start })))
            .collect();
        jobs.into_iter().map(|j| j.join().expect("chunk worker")).collect()
    })
}
