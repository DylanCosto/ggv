//! What a file is, from its name: shared by the command line, drag-and-drop and recent files.

use std::path::{Path, PathBuf};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Signal,
    Alignments,
    Variants,
    Genes,
    Reference,
    Session,
}

pub fn kind_of(name: &str) -> Option<Kind> {
    let lower = name.to_ascii_lowercase();
    let base = lower.trim_end_matches(".gz").trim_end_matches(".bgz");
    if ggv_io::signal::is_signal(&lower) {
        Some(Kind::Signal)
    } else if base.ends_with(".bam") || base.ends_with(".cram") {
        Some(Kind::Alignments)
    } else if base.ends_with(".vcf") || base.ends_with(".bcf") {
        Some(Kind::Variants)
    } else if [".gtf", ".gff", ".gff3", ".bed", ".bed12"].iter().any(|e| base.ends_with(e)) {
        Some(Kind::Genes)
    } else if [".fa", ".fasta", ".fna", ".2bit"].iter().any(|e| base.ends_with(e)) {
        Some(Kind::Reference)
    } else if lower.ends_with(".json") {
        Some(Kind::Session)
    } else {
        None
    }
}

/// Files sorted by kind; anything unrecognised is listed in `unknown`.
#[derive(Debug, Default, PartialEq, Eq)]
pub struct Sorted {
    pub bams: Vec<PathBuf>,
    pub vcfs: Vec<PathBuf>,
    pub genes: Vec<PathBuf>,
    pub signals: Vec<PathBuf>,
    pub reference: Option<PathBuf>,
    pub session: Option<PathBuf>,
    pub unknown: Vec<PathBuf>,
}

pub fn sort(paths: impl IntoIterator<Item = PathBuf>) -> Sorted {
    let mut out = Sorted::default();
    for p in paths {
        match kind_of(&p.to_string_lossy()) {
            Some(Kind::Alignments) => out.bams.push(p),
            Some(Kind::Variants) => out.vcfs.push(p),
            Some(Kind::Genes) => out.genes.push(p),
            Some(Kind::Signal) => out.signals.push(p),
            Some(Kind::Reference) => out.reference = Some(p),
            Some(Kind::Session) => out.session = Some(p),
            None => out.unknown.push(p),
        }
    }
    out
}

/// Recently opened files, newest first, kept next to the last session.
pub fn recent_path() -> Option<PathBuf> {
    crate::session::last_session_path().map(|p| p.with_file_name("recent.json"))
}

pub fn load_recent() -> Vec<PathBuf> {
    recent_path()
        .and_then(|p| std::fs::read_to_string(p).ok())
        .and_then(|t| serde_json::from_str(&t).ok())
        .unwrap_or_default()
}

/// Adds `opened` to the front of the recent list (deduplicated, at most 15) and saves it.
pub fn remember(recent: &mut Vec<PathBuf>, opened: &[PathBuf]) {
    for p in opened.iter().rev() {
        recent.retain(|r| r != p);
        recent.insert(0, p.clone());
    }
    recent.truncate(15);
    if let Some(path) = recent_path() {
        if let Some(dir) = path.parent() {
            let _ = std::fs::create_dir_all(dir);
        }
        let _ = std::fs::write(path, serde_json::to_string_pretty(recent).unwrap_or_default());
    }
}

pub fn display_name(p: &Path) -> String {
    p.file_name().map_or_else(|| p.display().to_string(), |n| n.to_string_lossy().into_owned())
}

/// `p` rebuilt from its components: one separator throughout (`C:/x\y/z` on Windows is
/// `C:\x\y\z`), without `.` parts. Reads nothing from disk.
pub fn normalize(p: &Path) -> PathBuf {
    p.components().collect()
}

/// Whether `a` and `b` name the same file: for existing files any spelling of it (separators,
/// relative, `..`, case on Windows, links), else the same path part by part.
pub fn same_file(a: &Path, b: &Path) -> bool {
    match (std::fs::canonicalize(a), std::fs::canonicalize(b)) {
        (Ok(x), Ok(y)) => x == y,
        _ => normalize(a) == normalize(b),
    }
}

/// Whether `want` (as an agent names an open file: a path in any spelling, or its file name)
/// names `path`.
pub fn names(path: &Path, want: &str) -> bool {
    path.file_name().is_some_and(|n| n.to_string_lossy() == want) || same_file(path, Path::new(want))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kinds_by_extension() {
        assert_eq!(kind_of("a.BAM"), Some(Kind::Alignments));
        assert_eq!(kind_of("a.cram"), Some(Kind::Alignments));
        assert_eq!(kind_of("calls.vcf.gz"), Some(Kind::Variants));
        assert_eq!(kind_of("x.bcf"), Some(Kind::Variants));
        assert_eq!(kind_of("genes.gff3.gz"), Some(Kind::Genes));
        assert_eq!(kind_of("hg38.fa.gz"), Some(Kind::Reference));
        assert_eq!(kind_of("review.ggv.json"), Some(Kind::Session));
        assert_eq!(kind_of("chr1:100"), None);
    }

    #[test]
    fn sorting_dropped_files() {
        let s = sort(["a.bam", "b.vcf.gz", "g.gtf", "r.fasta", "notes.txt"].map(PathBuf::from));
        assert_eq!(s.bams, vec![PathBuf::from("a.bam")]);
        assert_eq!(s.vcfs, vec![PathBuf::from("b.vcf.gz")]);
        assert_eq!(s.genes, vec![PathBuf::from("g.gtf")]);
        assert_eq!(s.reference, Some(PathBuf::from("r.fasta")));
        assert_eq!(s.unknown, vec![PathBuf::from("notes.txt")]);
    }

    #[test]
    fn paths_name_files_in_any_spelling() {
        let dir = std::env::temp_dir().join(format!("ggv-same-file-{}", std::process::id()));
        std::fs::create_dir_all(dir.join("sub")).unwrap();
        let f = dir.join("calls.vcf.gz");
        std::fs::write(&f, "x").unwrap();
        let slashes = PathBuf::from(f.to_string_lossy().replace('\\', "/"));
        assert!(same_file(&f, &slashes));
        assert!(same_file(&f, &dir.join("sub").join("..").join("calls.vcf.gz")));
        assert!(names(&f, &slashes.to_string_lossy()) && names(&f, "calls.vcf.gz"));
        assert!(!names(&f, "other.vcf.gz"));
        // Missing files compare part by part.
        assert!(same_file(Path::new("/no/such/./x.bam"), Path::new("/no/such/x.bam")));
        assert!(!same_file(Path::new("/no/such/x.bam"), Path::new("/no/such/y.bam")));
        #[cfg(windows)]
        assert_eq!(normalize(Path::new("C:/data\\a/b.bam")), PathBuf::from(r"C:\data\a\b.bam"));
        let _ = std::fs::remove_dir_all(dir);
    }
}
