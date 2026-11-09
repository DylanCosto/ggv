//! The demo: two slices of public PacBio Revio HiFi data for GIAB HG002 on GRCh38 chr20 (reads
//! with haplotype tags and methylation, DeepVariant and pbsv calls, RefSeq genes and the
//! reference), about 4.3 MB. It isn't built into the binaries: [`fetch`] downloads it once, from
//! the `demo-v1` GitHub release (`GGV_DEMO_URL` overrides), and keeps it in the cache.
//!
//! The files are `tests/data/hg002-revio` (its README says where they come from); the release
//! workflow publishes them.

use std::path::{Path, PathBuf};

/// Where the files are published, one URL per file name.
const BASE_URL: &str = "https://github.com/DylanCosto/ggv/releases/download/demo-v1";

/// Every file, by name, and the whole download's size (for progress before it starts).
pub const FILES: &[&str] = &[
    "HG002.revio.slices.bam",
    "HG002.revio.slices.bam.bai",
    "HG002.revio.deepvariant.phased.vcf.gz",
    "HG002.revio.deepvariant.phased.vcf.gz.tbi",
    "HG002.revio.pbsv.vcf.gz",
    "HG002.revio.pbsv.vcf.gz.tbi",
    "GRCh38.chr20.slices.fa.gz",
    "GRCh38.chr20.slices.fa.gz.fai",
    "GRCh38.chr20.slices.fa.gz.gzi",
    "ncbiRefSeq.slices.gtf",
];
pub const SIZE: u64 = 4_306_830;

/// The demo's files on disk.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Demo {
    pub dir: PathBuf,
    pub bam: PathBuf,
    pub vcfs: Vec<PathBuf>,
    pub genes: PathBuf,
    pub reference: PathBuf,
}

/// The places worth seeing, in tour order: (locus, what it shows).
pub const STOPS: &[(&str, &str)] = &[
    (
        "chr20:31,165,000-31,184,000",
        "A heterozygous 6.8 kb deletion: grouped by haplotype, every read carrying it is HP 1 and every intact read spanning it is HP 2. Split reads (supplementary alignments) jump the gap.",
    ),
    ("chr20:58,854,246", "A one-base deletion in a run of Cs: Evidence counts the reads with and without it, by strand, with quality and warnings."),
    (
        "chr20:58,888,000-58,890,000",
        "GNAS imprinting: coloured by base modifications, one haplotype's CpGs are about 17% methylated and the other's 88%.",
    ),
];

impl Demo {
    fn at(dir: PathBuf) -> Self {
        Self {
            bam: dir.join("HG002.revio.slices.bam"),
            vcfs: vec![dir.join("HG002.revio.deepvariant.phased.vcf.gz"), dir.join("HG002.revio.pbsv.vcf.gz")],
            genes: dir.join("ncbiRefSeq.slices.gtf"),
            reference: dir.join("GRCh38.chr20.slices.fa.gz"),
            dir,
        }
    }

    /// Every file to open, alignments first.
    pub fn files(&self) -> Vec<PathBuf> {
        let mut v = vec![self.bam.clone()];
        v.extend(self.vcfs.iter().cloned());
        v.push(self.genes.clone());
        v.push(self.reference.clone());
        v
    }
}

/// Where the demo is kept: `<cache>/ggv/demo/v1` (GGV_DEMO_DIR overrides).
pub fn dir() -> Option<PathBuf> {
    if let Some(d) = crate::folders::env_var_os("DEMO_DIR") {
        return Some(PathBuf::from(d));
    }
    let genomes = crate::genomes::cache_dir()?;
    Some(genomes.parent()?.join("demo").join("v1"))
}

/// The demo, if it has been downloaded.
pub fn cached() -> Option<Demo> {
    let dir = dir()?;
    FILES.iter().all(|f| dir.join(f).is_file()).then(|| Demo::at(dir))
}

/// Downloads whatever of the demo is missing (resuming partial files) and returns it.
pub fn fetch(progress: crate::genomes::Progress) -> std::io::Result<Demo> {
    let dir = dir().ok_or_else(|| std::io::Error::other("no cache folder (set GGV_DEMO_DIR)"))?;
    let base = crate::folders::env_var_os("DEMO_URL").and_then(|u| u.into_string().ok()).unwrap_or_else(|| BASE_URL.to_owned());
    fetch_from(base.trim_end_matches('/'), &dir, progress)
}

fn fetch_from(base: &str, dir: &Path, progress: crate::genomes::Progress) -> std::io::Result<Demo> {
    std::fs::create_dir_all(dir)?;
    for name in FILES {
        let dest = dir.join(name);
        if !dest.is_file() {
            crate::genomes::download(&format!("{base}/{name}"), &dest, progress)?;
        }
    }
    let d = Demo::at(dir.to_owned());
    // Opening them proves the download is whole.
    crate::BamFile::open(&d.bam).map_err(|e| std::io::Error::other(format!("{}: {e}", d.bam.display())))?;
    crate::ReferenceFasta::open(&d.reference).map_err(|e| std::io::Error::other(format!("{}: {e}", d.reference.display())))?;
    Ok(d)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn size_matches_the_files() {
        let src = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio");
        let total: u64 = FILES.iter().map(|f| std::fs::metadata(src.join(f)).unwrap().len()).sum();
        // Text files may differ by line endings on a Windows checkout.
        assert!(total.abs_diff(SIZE) < 10_000, "{total} vs {SIZE}");
    }

    #[test]
    fn downloads_and_opens() {
        // From the source tree through curl's file:// support: the same path as a real download.
        let src = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio").canonicalize().unwrap();
        let mut url = src.to_string_lossy().replace('\\', "/");
        if let Some(rest) = url.strip_prefix("//?/") {
            url = rest.to_owned();
        }
        let base = format!("file:///{}", url.trim_start_matches('/'));
        let dir = std::env::temp_dir().join(format!("ggv-demo-test-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let d = match fetch_from(&base, &dir, &mut |_, _, _| {}) {
            Ok(d) => d,
            // No curl on this machine: nothing to test.
            Err(e) if e.to_string().contains("couldn't run curl") => return,
            Err(e) => panic!("{e}"),
        };
        for f in d.files() {
            assert!(f.exists(), "{}", f.display());
        }
        let _ = std::fs::remove_dir_all(dir);
    }
}
