//! Whether files are the same assembly: each file's build and how GGV knows it, and where the
//! files disagree, so an hg19 BAM against an hg38 reference or annotation is explained before
//! anyone trusts its mismatches or consequence labels. Shared by the desktop and `ggv mcp`.
//!
//! A build is *declared* when a file names it (a VCF's `##reference`, a file name such as
//! `hg38.ncbiRefSeq.gtf.gz`), *inferred* when only chromosome 1's length matches a known build,
//! and *unknown* otherwise. An inferred build is a hint, not a verification: two assemblies can
//! share chr1's length, so every contig the files share is compared too. Custom assemblies are
//! fine: unknown is not an error.

use std::{
    collections::BTreeMap,
    path::{Path, PathBuf},
};

/// How GGV knows a file's build.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum How {
    Declared,
    Inferred,
    Unknown,
}

impl How {
    pub fn text(self) -> &'static str {
        match self {
            How::Declared => "declared",
            How::Inferred => "inferred from chr1's length (not verified)",
            How::Unknown => "unknown",
        }
    }
}

/// One file's assembly.
#[derive(Clone, Debug)]
pub struct Identity {
    pub kind: &'static str,
    pub path: PathBuf,
    pub build: Option<&'static str>,
    pub how: How,
    /// What it declares, as written (`file:///ref/GRCh38.fa`).
    pub declared_text: Option<String>,
    /// Contig lengths it states (none for annotations, whose extents aren't lengths).
    pub lengths: Vec<(String, u64)>,
}

impl Identity {
    /// A file's identity from what it declares (a header value or its name) and the contig
    /// lengths it states.
    pub fn new(kind: &'static str, path: PathBuf, declared_text: Option<String>, lengths: Vec<(String, u64)>) -> Self {
        let declared = declared_text.as_deref().and_then(build_named);
        let inferred = crate::catalog::build_of(lengths.iter().map(|(n, l)| (n.as_str(), *l)));
        let (build, how) = match (declared, inferred) {
            (Some(b), _) => (Some(b), How::Declared),
            (None, Some(b)) => (Some(b), How::Inferred),
            (None, None) => (None, How::Unknown),
        };
        Self { kind, path, build, how, declared_text, lengths }
    }
}

/// The build a name or header value names (`hg38.ncbiRefSeq.gtf.gz`, `GRCh37`, `b37`), when it
/// names exactly one.
pub fn build_named(text: &str) -> Option<&'static str> {
    let name = text.to_ascii_lowercase();
    let words: Vec<&str> = name.split(|c: char| !c.is_ascii_alphanumeric()).filter(|w| !w.is_empty()).collect();
    let mut found = None;
    for g in crate::genomes::GENOMES {
        let names = std::iter::once(g.id).chain(g.aliases.iter().copied()).filter(|n| !matches!(*n, "human" | "mouse"));
        if names.into_iter().any(|n| words.contains(&n)) {
            if found.is_some_and(|f| f != g.id) {
                return None; // names two builds (a liftover file): don't guess
            }
            found = Some(g.id);
        }
    }
    found
}

/// `chr1`, `1` → `1`; `chrM`, `MT` → `mt`: names compared across files.
pub fn canonical_contig(name: &str) -> String {
    let n = name.strip_prefix("chr").unwrap_or(name).to_ascii_lowercase();
    if n == "m" { "mt".into() } else { n }
}

pub fn display_name(p: &Path) -> String {
    p.file_name().map_or_else(|| p.display().to_string(), |n| n.to_string_lossy().into_owned())
}

fn commas(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    for (i, ch) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(ch);
    }
    out
}

/// What two files' contig lengths say: `None` when they agree on every contig both name.
fn compare(a: &Identity, b: &Identity) -> Option<String> {
    let theirs: BTreeMap<String, u64> = b.lengths.iter().map(|(n, l)| (canonical_contig(n), *l)).collect();
    let mut differ = Vec::new();
    let mut shared = 0;
    for (n, l) in &a.lengths {
        if let Some(&m) = theirs.get(&canonical_contig(n)) {
            shared += 1;
            if m != *l && *l > 0 && m > 0 {
                differ.push((n.clone(), *l, m));
            }
        }
    }
    let (name, la, lb) = differ.first()?.clone();
    let more = if differ.len() > 1 { format!(" ({} of {shared} shared contigs differ)", differ.len()) } else { String::new() };
    Some(format!(
        "{} and {} are different assemblies: {name} is {} bp in the {} but {} bp in the {}{more}.",
        display_name(&a.path),
        display_name(&b.path),
        commas(la),
        a.kind,
        commas(lb),
        b.kind,
    ))
}

/// Every disagreement among `ids`, as sentences.
pub fn conflicts(ids: &[Identity]) -> Vec<String> {
    let mut out = Vec::new();
    for id in ids {
        if id.how == How::Declared
            && let (Some(declared), Some(inferred)) = (id.build, crate::catalog::build_of(id.lengths.iter().map(|(n, l)| (n.as_str(), *l))))
            && declared != inferred
        {
            out.push(format!("{} declares {declared}, but its chr1 length is {inferred}'s.", display_name(&id.path)));
        }
    }
    // Each file against the first (the reference when one is open), so one odd file is one
    // finding, not one per other file.
    if let Some((a, rest)) = ids.split_first() {
        for b in rest {
            // Lengths are the stronger evidence; builds are compared when there are none to compare.
            let by_length = (!a.lengths.is_empty() && !b.lengths.is_empty()).then(|| compare(a, b)).flatten();
            if let Some(c) = by_length {
                out.push(c);
            } else if let (Some(x), Some(y)) = (a.build, b.build)
                && x != y
                && (a.lengths.is_empty() || b.lengths.is_empty())
            {
                out.push(format!("{} is {y} ({}) but the {} {} is {x} ({}).", display_name(&b.path), b.how.text(), a.kind, display_name(&a.path), a.how.text()));
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn id(kind: &'static str, path: &str, build: Option<&'static str>, how: How, lengths: &[(&str, u64)]) -> Identity {
        Identity { kind, path: path.into(), build, how, declared_text: None, lengths: lengths.iter().map(|(n, l)| (n.to_string(), *l)).collect() }
    }

    #[test]
    fn hg19_and_hg38_are_told_apart_by_length() {
        let bam = id("alignments", "/d/hg19.bam", Some("hg19"), How::Inferred, &[("chr1", 249_250_621), ("chr2", 243_199_373)]);
        let vcf = id("variants", "/d/calls.vcf.gz", Some("hg38"), How::Inferred, &[("1", 248_956_422), ("2", 242_193_529)]);
        let c = conflicts(&[bam.clone(), vcf]);
        assert_eq!(c.len(), 1, "{c:?}");
        assert!(c[0].contains("chr1 is 249,250,621 bp") && c[0].contains("2 of 2 shared"), "{}", c[0]);
        // The same assembly under other contig names agrees.
        let same = id("variants", "/d/b.vcf.gz", Some("hg19"), How::Inferred, &[("1", 249_250_621), ("2", 243_199_373)]);
        assert!(conflicts(&[bam.clone(), same]).is_empty());
        // An hg38 annotation (by name) with an hg19 BAM.
        let genes = Identity::new("annotation", "/d/hg38.ncbiRefSeq.gtf.gz".into(), Some("hg38.ncbiRefSeq.gtf.gz".into()), Vec::new());
        assert_eq!((genes.build, genes.how), (Some("hg38"), How::Declared));
        assert_eq!(conflicts(&[bam.clone(), genes]).len(), 1);
        // A custom assembly is not a conflict.
        let custom = id("alignments", "/d/x.bam", None, How::Unknown, &[("ctg1", 5_000)]);
        assert!(conflicts(&[custom.clone(), id("annotation", "/d/x.gff", None, How::Unknown, &[])]).is_empty());
    }

    #[test]
    fn a_declared_build_that_its_lengths_contradict() {
        let vcf = id("variants", "/d/calls.vcf.gz", Some("hg19"), How::Declared, &[("chr1", 248_956_422)]);
        let c = conflicts(&[vcf]);
        assert!(c.len() == 1 && c[0].contains("declares hg19") && c[0].contains("hg38"), "{c:?}");
    }

    #[test]
    fn builds_named_in_text() {
        assert_eq!(build_named("hg38.ncbiRefSeq.gtf.gz"), Some("hg38"));
        assert_eq!(build_named("Homo_sapiens.GRCh37.87.gtf"), Some("hg19"));
        assert_eq!(build_named("hg19_to_hg38.chain"), None);
        assert_eq!(build_named("x.gff"), None);
    }
}
