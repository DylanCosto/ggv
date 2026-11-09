//! Core genomic types shared by every GGV crate. Nothing here does I/O, so it builds for WASM.

mod locus;

pub use locus::{Locus, LocusParseError};

/// A reference sequence (contig) as declared by a file's header.
#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct Contig {
    pub name: String,
    pub length: u64,
}

/// A locus resolved against a file: a known contig and a concrete, 0-based half-open range.
#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct Region {
    /// The contig name as the file spells it (which may differ from what the user typed,
    /// e.g. `chr7` → `7`).
    pub contig: String,
    /// 0-based inclusive start.
    pub start: u64,
    /// 0-based exclusive end.
    pub end: u64,
}

impl Region {
    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    pub fn is_empty(&self) -> bool {
        self.end == self.start
    }

    /// Whether a 0-based position falls inside the region.
    pub fn contains(&self, pos: u64) -> bool {
        (self.start..self.end).contains(&pos)
    }
}

impl std::fmt::Display for Region {
    /// Formats in 1-based inclusive coordinates, the way people and IGV write loci.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}-{}", self.contig, self.start + 1, self.end)
    }
}

/// Alternative spellings of a contig name to try when the exact name is not in a file:
/// `chr7` ⇄ `7`; the mitochondrion as `chrM`, `MT` or `M` (b37 references say `MT`, some
/// aligners' BAMs `M`).
pub fn contig_aliases(name: &str) -> Vec<String> {
    let mut out = Vec::new();
    match name {
        "chrM" | "chrMT" => return vec!["MT".to_owned(), "M".to_owned(), "chrM".to_owned()].into_iter().filter(|x| x != name).collect(),
        "MT" => return vec!["M".to_owned(), "chrM".to_owned()],
        "M" => return vec!["MT".to_owned(), "chrM".to_owned()],
        _ => {}
    }
    if let Some(bare) = name.strip_prefix("chr") {
        out.push(bare.to_owned());
    } else {
        out.push(format!("chr{name}"));
    }
    out
}

/// A junction's splice-site motif from its intron's first two and last two bases (on the
/// forward strand), and the strand that motif implies: GT-AG, GC-AG and AT-AC are read on +,
/// their reverse complements (CT-AC, CT-GC, GT-AT) on −; anything else is non-canonical (`.`).
pub fn splice_motif(first: [u8; 2], last: [u8; 2]) -> (&'static str, char) {
    match (first.map(|b| b.to_ascii_uppercase()), last.map(|b| b.to_ascii_uppercase())) {
        ([b'G', b'T'], [b'A', b'G']) => ("GT/AG", '+'),
        ([b'C', b'T'], [b'A', b'C']) => ("GT/AG", '-'),
        ([b'G', b'C'], [b'A', b'G']) => ("GC/AG", '+'),
        ([b'C', b'T'], [b'G', b'C']) => ("GC/AG", '-'),
        ([b'A', b'T'], [b'A', b'C']) => ("AT/AC", '+'),
        ([b'G', b'T'], [b'A', b'T']) => ("AT/AC", '-'),
        _ => ("non-canonical", '.'),
    }
}

#[cfg(test)]
mod tests {
    use super::contig_aliases;

    #[test]
    fn aliases_cover_chr_prefix_and_every_mitochondrion_spelling() {
        assert_eq!(contig_aliases("chr7"), ["7"]);
        assert_eq!(contig_aliases("7"), ["chr7"]);
        for (name, others) in [("M", ["MT", "chrM"]), ("MT", ["M", "chrM"])] {
            assert_eq!(contig_aliases(name), others);
        }
        assert_eq!(contig_aliases("chrM"), ["MT", "M"]);
    }
}
