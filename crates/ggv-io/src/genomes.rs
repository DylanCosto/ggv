//! Built-in genomes: a name (`hg38`) instead of finding a reference and gene annotation.
//!
//! Each genome is UCSC's `.2bit` (read at any position without unpacking; see
//! [`crate::twobit`]) and its NCBI RefSeq gene annotation, downloaded once into a cache folder
//! (`GGV_GENOME_DIR`, else `%LOCALAPPDATA%\ggv\genomes`, `~/Library/Caches/ggv/genomes` or
//! `$XDG_CACHE_HOME/ggv/genomes`). Downloads use the system's `curl` (Windows 10+, macOS and
//! Linux ship it), into a `.part` file that resumes if interrupted.

use std::{
    path::{Path, PathBuf},
    process::{Command, Stdio},
    time::Duration,
};

pub struct Genome {
    /// What to type: `hg38`.
    pub id: &'static str,
    pub name: &'static str,
    /// Other names it answers to (case-insensitive).
    pub aliases: &'static [&'static str],
    /// Approximate download size of the reference and genes, in bytes (for the list).
    pub size: u64,
}

const UCSC: &str = "https://hgdownload.soe.ucsc.edu/goldenPath";

pub const GENOMES: &[Genome] = &[
    Genome { id: "hg38", name: "Human GRCh38 (hg38)", aliases: &["grch38", "hg38", "human"], size: 877_000_000 },
    Genome { id: "hg19", name: "Human GRCh37 (hg19)", aliases: &["grch37", "hg19", "b37"], size: 836_000_000 },
    Genome { id: "hs1", name: "Human T2T-CHM13 v2.0 (hs1)", aliases: &["t2t", "chm13", "t2t-chm13", "hs1"], size: 852_000_000 },
    Genome { id: "mm39", name: "Mouse GRCm39 (mm39)", aliases: &["grcm39", "mm39", "mouse"], size: 745_000_000 },
    Genome { id: "mm10", name: "Mouse GRCm38 (mm10)", aliases: &["grcm38", "mm10"], size: 743_000_000 },
    Genome { id: "sacCer3", name: "Yeast S. cerevisiae R64 (sacCer3)", aliases: &["saccer3", "yeast", "r64"], size: 3_400_000 },
];

impl Genome {
    pub fn reference_url(&self) -> String {
        format!("{UCSC}/{0}/bigZips/{0}.2bit", self.id)
    }

    pub fn genes_url(&self) -> String {
        format!("{UCSC}/{0}/bigZips/genes/{0}.ncbiRefSeq.gtf.gz", self.id)
    }
}

/// A genome by id or alias (`GRCh38`, `T2T`, `mouse`...).
pub fn find(name: &str) -> Option<&'static Genome> {
    let n = name.trim().to_ascii_lowercase();
    GENOMES.iter().find(|g| g.id.eq_ignore_ascii_case(&n) || g.aliases.contains(&n.as_str()))
}

/// Where downloaded genomes are kept.
pub fn cache_dir() -> Option<PathBuf> {
    if let Some(d) = crate::folders::env_var_os("GENOME_DIR") {
        return Some(PathBuf::from(d));
    }
    Some(crate::folders::cache_dir()?.join("genomes"))
}

/// A genome's files on disk (they exist once [`fetch`] has run).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Local {
    pub reference: PathBuf,
    pub genes: PathBuf,
}

pub fn local(g: &Genome) -> Option<Local> {
    let dir = cache_dir()?.join(g.id);
    Some(Local { reference: dir.join(format!("{}.2bit", g.id)), genes: dir.join(format!("{}.ncbiRefSeq.gtf.gz", g.id)) })
}

/// The genome's files if both are already downloaded.
pub fn cached(g: &Genome) -> Option<Local> {
    local(g).filter(|l| l.reference.is_file() && l.genes.is_file())
}

/// Download progress: which file, bytes so far, bytes in all (0 if unknown).
pub type Progress<'a> = &'a mut dyn FnMut(&str, u64, u64);

/// Downloads whatever of the genome is missing (resuming partial files) and returns its files.
pub fn fetch(g: &Genome, progress: Progress) -> std::io::Result<Local> {
    let l = local(g).ok_or_else(|| std::io::Error::other("no cache folder (set GGV_GENOME_DIR)"))?;
    if let Some(dir) = l.reference.parent() {
        std::fs::create_dir_all(dir)?;
    }
    for (url, dest) in [(g.reference_url(), &l.reference), (g.genes_url(), &l.genes)] {
        if !dest.is_file() {
            download(&url, dest, progress)?;
        }
    }
    // Opening it once proves the reference is whole (a cut-off .2bit fails here, not later).
    crate::ReferenceFasta::open(&l.reference).map_err(|e| std::io::Error::other(format!("{}: {e}", l.reference.display())))?;
    Ok(l)
}

fn content_length(url: &str) -> u64 {
    let out = Command::new("curl").args(["-sIL", url]).output();
    out.ok()
        .map(|o| String::from_utf8_lossy(&o.stdout).into_owned())
        .and_then(|h| h.lines().rev().find_map(|l| l.to_ascii_lowercase().strip_prefix("content-length:").and_then(|v| v.trim().parse().ok())))
        .unwrap_or(0)
}

/// `url` to `dest` with curl, through `dest.part` (resumed when present).
pub(crate) fn download(url: &str, dest: &Path, progress: Progress) -> std::io::Result<()> {
    let part = dest.with_extension(format!("{}.part", dest.extension().and_then(|e| e.to_str()).unwrap_or("")));
    let name = dest.file_name().map_or_else(String::new, |n| n.to_string_lossy().into_owned());
    let total = content_length(url);
    let mut child = Command::new("curl")
        .args(["-sSfL", "--retry", "3", "-C", "-", "-o"])
        .arg(&part)
        .arg(url)
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| std::io::Error::other(format!("couldn't run curl to download {url}: {e}")))?;
    loop {
        let done = std::fs::metadata(&part).map_or(0, |m| m.len());
        progress(&name, done, total);
        if let Some(status) = child.try_wait()? {
            if !status.success() {
                let mut msg = String::new();
                if let Some(mut e) = child.stderr.take() {
                    use std::io::Read;
                    let _ = e.read_to_string(&mut msg);
                }
                return Err(std::io::Error::other(format!("downloading {url} failed: {}", msg.trim())));
            }
            break;
        }
        std::thread::sleep(Duration::from_millis(250));
    }
    let done = std::fs::metadata(&part)?.len();
    if total > 0 && done != total {
        return Err(std::io::Error::other(format!("{url}: got {done} of {total} bytes; run again to resume")));
    }
    progress(&name, done, total.max(done));
    std::fs::rename(&part, dest)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_and_aliases() {
        assert_eq!(find("GRCh38").map(|g| g.id), Some("hg38"));
        assert_eq!(find("hg19").map(|g| g.id), Some("hg19"));
        assert_eq!(find("T2T").map(|g| g.id), Some("hs1"));
        assert_eq!(find("mouse").map(|g| g.id), Some("mm39"));
        assert_eq!(find("SACCER3").map(|g| g.id), Some("sacCer3"));
        assert!(find("hg99").is_none());
        assert_eq!(find("hg38").unwrap().reference_url(), "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit");
        assert!(find("hs1").unwrap().genes_url().ends_with("hs1/bigZips/genes/hs1.ncbiRefSeq.gtf.gz"));
    }
}
