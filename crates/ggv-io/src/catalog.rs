//! Built-in tracks: gene annotations, ClinVar and missense pathogenicity scores for each
//! built-in genome build, by name instead of finding the files.
//!
//! Like the genomes ([`crate::genomes`]), each is downloaded once into the genome cache folder
//! (`<cache>/<build>/tracks`), resuming partial files. Scores are UCSC's bigWigs, one per
//! alternate base (the score of a change to A, C, G or T), so they open as four signal tracks.
//! Sources too large to keep locally (CADD at 39 GB, phyloP at 10 GB) are left out until GGV
//! reads remote files.

use std::path::{Path, PathBuf};

use crate::genomes::{Progress, cache_dir, download};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Group {
    Genes,
    Variants,
    Scores,
}

impl Group {
    pub fn label(self) -> &'static str {
        match self {
            Self::Genes => "Gene annotations",
            Self::Variants => "Variants",
            Self::Scores => "Pathogenicity scores",
        }
    }
}

/// Where a track's files come from.
pub enum Files {
    /// These URLs, each saved under its name.
    Fixed(&'static [(&'static str, &'static str)]),
    /// The newest `<stem>_YYYYMMDD.vcf.gz` and its `.tbi` in a folder (ClinVar's dated weekly
    /// releases: a fixed name would change under a resumed download).
    Dated { dir: &'static str, stem: &'static str },
}

pub struct Track {
    /// What to type: `clinvar`.
    pub id: &'static str,
    /// The genome build ([`crate::genomes::GENOMES`] id).
    pub build: &'static str,
    pub group: Group,
    pub name: &'static str,
    pub about: &'static str,
    /// Approximate download size, in bytes (for the list).
    pub size: u64,
    pub files: Files,
}

/// One bigWig per alternate base from UCSC's `gbdb/<build>/<dir>`.
macro_rules! per_base {
    ($build:literal, $dir:literal, $name:literal) => {
        &[
            (concat!("https://hgdownload.soe.ucsc.edu/gbdb/", $build, "/", $dir, "/a.bw"), concat!($name, ".A.bw")),
            (concat!("https://hgdownload.soe.ucsc.edu/gbdb/", $build, "/", $dir, "/c.bw"), concat!($name, ".C.bw")),
            (concat!("https://hgdownload.soe.ucsc.edu/gbdb/", $build, "/", $dir, "/g.bw"), concat!($name, ".G.bw")),
            (concat!("https://hgdownload.soe.ucsc.edu/gbdb/", $build, "/", $dir, "/t.bw"), concat!($name, ".T.bw")),
        ]
    };
}

/// RefSeq genes as UCSC packages them for `build`.
macro_rules! refseq {
    ($build:literal) => {
        Files::Fixed(&[(
            concat!("https://hgdownload.soe.ucsc.edu/goldenPath/", $build, "/bigZips/genes/", $build, ".ncbiRefSeq.gtf.gz"),
            concat!($build, ".ncbiRefSeq.gtf.gz"),
        )])
    };
}

const ALPHAMISSENSE: &str = "AlphaMissense (Cheng et al. 2023): likelihood a missense change is pathogenic, 0-1; one track per alternate base";
const REVEL: &str = "REVEL (Ioannidis et al. 2016): ensemble missense pathogenicity score, 0-1; one track per alternate base";

pub const TRACKS: &[Track] = &[
    // GRCh38
    Track {
        id: "gencode",
        build: "hg38",
        group: Group::Genes,
        name: "GENCODE 50 (comprehensive)",
        about: "Every GENCODE/Ensembl transcript, with MANE and canonical tags",
        size: 125_000_000,
        files: Files::Fixed(&[("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz", "gencode.v50.annotation.gtf.gz")]),
    },
    Track {
        id: "gencode-basic",
        build: "hg38",
        group: Group::Genes,
        name: "GENCODE 50 (basic)",
        about: "Full-length protein-coding and lncRNA transcripts only",
        size: 80_000_000,
        files: Files::Fixed(&[("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.basic.annotation.gtf.gz", "gencode.v50.basic.annotation.gtf.gz")]),
    },
    Track {
        id: "mane",
        build: "hg38",
        group: Group::Genes,
        name: "MANE Select 1.5",
        about: "One matched Ensembl/RefSeq transcript per protein-coding gene (plus MANE Plus Clinical)",
        size: 9_000_000,
        files: Files::Fixed(&[("https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz", "MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz")]),
    },
    Track { id: "refseq", build: "hg38", group: Group::Genes, name: "NCBI RefSeq", about: "RefSeq transcripts (NM_/NR_/XM_), as UCSC packages them", size: 42_000_000, files: refseq!("hg38") },
    Track {
        id: "clinvar",
        build: "hg38",
        group: Group::Variants,
        name: "ClinVar (newest weekly release)",
        about: "Submitted clinical significance (CLNSIG), review status, conditions and gene for each variant",
        size: 202_000_000,
        files: Files::Dated { dir: "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38", stem: "clinvar" },
    },
    Track { id: "alphamissense", build: "hg38", group: Group::Scores, name: "AlphaMissense", about: ALPHAMISSENSE, size: 1_855_000_000, files: Files::Fixed(per_base!("hg38", "alphaMissense", "AlphaMissense")) },
    Track { id: "revel", build: "hg38", group: Group::Scores, name: "REVEL", about: REVEL, size: 1_875_000_000, files: Files::Fixed(per_base!("hg38", "revel", "REVEL")) },
    // GRCh37
    Track {
        id: "gencode",
        build: "hg19",
        group: Group::Genes,
        name: "GENCODE 50 (lifted to GRCh37)",
        about: "GENCODE's current release mapped back to GRCh37",
        size: 160_000_000,
        files: Files::Fixed(&[("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/GRCh37_mapping/gencode.v50lift37.annotation.gtf.gz", "gencode.v50lift37.annotation.gtf.gz")]),
    },
    Track { id: "refseq", build: "hg19", group: Group::Genes, name: "NCBI RefSeq", about: "RefSeq transcripts (NM_/NR_/XM_), as UCSC packages them", size: 20_000_000, files: refseq!("hg19") },
    Track {
        id: "clinvar",
        build: "hg19",
        group: Group::Variants,
        name: "ClinVar (newest weekly release)",
        about: "Submitted clinical significance (CLNSIG), review status, conditions and gene for each variant",
        size: 202_000_000,
        files: Files::Dated { dir: "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37", stem: "clinvar" },
    },
    Track { id: "alphamissense", build: "hg19", group: Group::Scores, name: "AlphaMissense", about: ALPHAMISSENSE, size: 1_810_000_000, files: Files::Fixed(per_base!("hg19", "alphaMissense", "AlphaMissense")) },
    Track { id: "revel", build: "hg19", group: Group::Scores, name: "REVEL", about: REVEL, size: 1_885_000_000, files: Files::Fixed(per_base!("hg19", "revel", "REVEL")) },
    // T2T-CHM13
    Track { id: "refseq", build: "hs1", group: Group::Genes, name: "NCBI RefSeq", about: "RefSeq annotation of T2T-CHM13, as UCSC packages it", size: 39_000_000, files: refseq!("hs1") },
    // Mouse
    Track {
        id: "gencode",
        build: "mm39",
        group: Group::Genes,
        name: "GENCODE M39",
        about: "Every GENCODE/Ensembl mouse transcript",
        size: 92_000_000,
        files: Files::Fixed(&[("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M39/gencode.vM39.annotation.gtf.gz", "gencode.vM39.annotation.gtf.gz")]),
    },
    Track { id: "refseq", build: "mm39", group: Group::Genes, name: "NCBI RefSeq", about: "RefSeq transcripts, as UCSC packages them", size: 31_000_000, files: refseq!("mm39") },
    Track {
        id: "gencode",
        build: "mm10",
        group: Group::Genes,
        name: "GENCODE M25",
        about: "The last GENCODE release on GRCm38",
        size: 29_000_000,
        files: Files::Fixed(&[("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz", "gencode.vM25.annotation.gtf.gz")]),
    },
    Track { id: "refseq", build: "mm10", group: Group::Genes, name: "NCBI RefSeq", about: "RefSeq transcripts, as UCSC packages them", size: 28_000_000, files: refseq!("mm10") },
    // Yeast
    Track { id: "refseq", build: "sacCer3", group: Group::Genes, name: "NCBI RefSeq", about: "RefSeq genes, as UCSC packages them", size: 400_000, files: refseq!("sacCer3") },
];

/// The built-in build these contigs are, from the first chromosome's length (`chr1`, `1` or
/// yeast's `chrI`); `None` for anything else.
pub fn build_of<'a>(contigs: impl IntoIterator<Item = (&'a str, u64)>) -> Option<&'static str> {
    const FIRST: [(u64, &str); 6] =
        [(248_956_422, "hg38"), (249_250_621, "hg19"), (248_387_328, "hs1"), (195_154_279, "mm39"), (195_471_971, "mm10"), (230_218, "sacCer3")];
    contigs
        .into_iter()
        .filter(|(name, _)| matches!(*name, "chr1" | "1" | "chrI" | "I"))
        .find_map(|(_, len)| FIRST.iter().find(|(l, _)| *l == len).map(|(_, b)| *b))
}

/// The tracks for `build`, in the catalogue's order (grouped).
pub fn for_build(build: &str) -> impl Iterator<Item = &'static Track> + '_ {
    TRACKS.iter().filter(move |t| t.build.eq_ignore_ascii_case(build))
}

/// A track by build and id (`hg38`, `clinvar`).
pub fn find(build: &str, id: &str) -> Option<&'static Track> {
    for_build(build).find(|t| t.id.eq_ignore_ascii_case(id))
}

/// Where a build's tracks are kept.
pub fn dir(build: &str) -> Option<PathBuf> {
    Some(cache_dir()?.join(build).join("tracks"))
}

fn is_index(name: &str) -> bool {
    name.ends_with(".tbi") || name.ends_with(".csi")
}

/// The newest `<stem>_YYYYMMDD.vcf.gz` among `names` that also has its `.tbi`.
fn newest_dated<'a>(names: impl IntoIterator<Item = &'a str>, stem: &str) -> Option<String> {
    let names: Vec<&str> = names.into_iter().collect();
    let dated = |n: &str| {
        n.strip_prefix(stem).and_then(|r| r.strip_prefix('_')).and_then(|r| r.strip_suffix(".vcf.gz")).is_some_and(|d| d.len() == 8 && d.bytes().all(|b| b.is_ascii_digit()))
    };
    names.iter().copied().filter(|n| dated(n) && names.contains(&format!("{n}.tbi").as_str())).max().map(str::to_owned)
}

/// The track's files to open (indexes left out: they're found beside their file) if they are
/// all downloaded.
pub fn cached(t: &Track) -> Option<Vec<PathBuf>> {
    let dir = dir(t.build)?;
    match &t.files {
        Files::Fixed(files) => {
            let all = files.iter().all(|(_, name)| dir.join(name).is_file());
            all.then(|| files.iter().filter(|(_, n)| !is_index(n)).map(|(_, n)| dir.join(n)).collect())
        }
        Files::Dated { stem, .. } => {
            let names: Vec<String> = std::fs::read_dir(&dir).ok()?.filter_map(|e| e.ok()?.file_name().into_string().ok()).collect();
            newest_dated(names.iter().map(String::as_str), stem).map(|n| vec![dir.join(n)])
        }
    }
}

/// File names linked from an HTTP folder listing.
fn listing(url: &str) -> std::io::Result<Vec<String>> {
    let out = std::process::Command::new("curl").args(["-sSfL", "--max-time", "60", url]).output()?;
    if !out.status.success() {
        return Err(std::io::Error::other(format!("listing {url} failed: {}", String::from_utf8_lossy(&out.stderr).trim())));
    }
    let html = String::from_utf8_lossy(&out.stdout);
    Ok(html.split("href=\"").skip(1).filter_map(|s| s.split('"').next()).map(|s| s.rsplit('/').next().unwrap_or(s).to_owned()).collect())
}

/// Downloads whatever of the track is missing (resuming partial files) and returns its files
/// to open.
pub fn fetch(t: &Track, progress: Progress) -> std::io::Result<Vec<PathBuf>> {
    let dir = dir(t.build).ok_or_else(|| std::io::Error::other("no cache folder (set GGV_GENOME_DIR)"))?;
    std::fs::create_dir_all(&dir)?;
    let get = |url: &str, dest: &Path, progress: Progress| if dest.is_file() { Ok(()) } else { download(url, dest, progress) };
    match &t.files {
        Files::Fixed(files) => {
            for (url, name) in *files {
                get(url, &dir.join(name), progress)?;
            }
        }
        Files::Dated { dir: remote, stem } => {
            if cached(t).is_none() {
                let names = listing(&format!("{remote}/"))?;
                let name = newest_dated(names.iter().map(String::as_str), stem)
                    .ok_or_else(|| std::io::Error::other(format!("no {stem}_YYYYMMDD.vcf.gz with an index in {remote}")))?;
                // The data first: an index without its file would look like a finished download.
                get(&format!("{remote}/{name}"), &dir.join(&name), progress)?;
                get(&format!("{remote}/{name}.tbi"), &dir.join(format!("{name}.tbi")), progress)?;
            }
        }
    }
    cached(t).ok_or_else(|| std::io::Error::other(format!("{}: files missing after download", t.name)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builds_from_first_chromosome() {
        assert_eq!(build_of([("chr1", 248_956_422), ("chr2", 242_193_529)]), Some("hg38"));
        assert_eq!(build_of([("MT", 16_569), ("1", 249_250_621)]), Some("hg19"));
        assert_eq!(build_of([("chrI", 230_218)]), Some("sacCer3"));
        assert_eq!(build_of([("chr1", 1000)]), None);
        assert_eq!(build_of([("chr20", 64_444_167)]), None, "a single other chromosome isn't enough");
    }

    #[test]
    fn catalogue_is_consistent() {
        for t in TRACKS {
            assert!(crate::genomes::find(t.build).is_some_and(|g| g.id == t.build), "{}: unknown build {}", t.id, t.build);
            assert_eq!(find(t.build, t.id).map(|x| x.name), Some(t.name), "{} {}: duplicate id", t.build, t.id);
            if let Files::Fixed(files) = &t.files {
                assert!(files.iter().all(|(url, name)| url.starts_with("https://") && !name.contains('/')), "{}", t.name);
            }
        }
        assert!(matches!(find("hg38", "ClinVar").unwrap().files, Files::Dated { .. }));
        let am = find("hg38", "alphamissense").unwrap();
        let Files::Fixed(files) = am.files else { panic!() };
        assert_eq!(files[3], ("https://hgdownload.soe.ucsc.edu/gbdb/hg38/alphaMissense/t.bw", "AlphaMissense.T.bw"));
    }

    #[test]
    fn newest_dated_release_with_index() {
        let names = ["clinvar.vcf.gz", "clinvar_20261025.vcf.gz", "clinvar_20261025.vcf.gz.tbi", "clinvar_20261101.vcf.gz", "clinvar_20261101.vcf.gz.md5", "clinvar_papu.vcf.gz"];
        assert_eq!(newest_dated(names, "clinvar").as_deref(), Some("clinvar_20261025.vcf.gz"), "the newest has no index yet");
        let names = ["clinvar_20261101.vcf.gz", "clinvar_20261101.vcf.gz.tbi", "clinvar_20261025.vcf.gz", "clinvar_20261025.vcf.gz.tbi"];
        assert_eq!(newest_dated(names, "clinvar").as_deref(), Some("clinvar_20261101.vcf.gz"));
    }
}
