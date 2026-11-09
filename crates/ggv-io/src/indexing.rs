//! Builds missing indexes, so a file without one can still be opened: BAM (`.bai`), CRAM
//! (`.crai`), bgzipped VCF (`.tbi`) and FASTA (`.fai`, plus `.gzi` when bgzipped). Each writes the
//! index next to the file, the way samtools/tabix name them.

use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write},
    path::{Path, PathBuf},
};

use pasta::{
    bam::{self, bai},
    bgzf, cram, fasta,
    variant::index::{IndexFormat, tabix::{TabixBuilder, TabixConfig}},
};

use crate::{Error, Result, error::IoContext};

/// What kind of index a file needs, from its name; `None` when GGV can't build one.
pub fn index_kind(path: &Path) -> Option<&'static str> {
    let name = path.file_name()?.to_string_lossy().to_ascii_lowercase();
    if name.ends_with(".bam") {
        Some("bai")
    } else if name.ends_with(".cram") {
        Some("crai")
    } else if name.ends_with(".vcf.gz") || name.ends_with(".vcf.bgz") {
        Some("tbi")
    } else if [".fa", ".fasta", ".fna"].iter().any(|e| name.ends_with(e) || name.ends_with(&format!("{e}.gz")) || name.ends_with(&format!("{e}.bgz"))) {
        Some("fai")
    } else {
        None
    }
}

/// Builds the index for `path` (see [`index_kind`]) and returns where it was written.
pub fn build_index(path: &Path) -> Result<PathBuf> {
    let dst = path.with_added_extension(index_kind(path).ok_or_else(|| Error::Io {
        path: path.to_owned(),
        cause: std::io::Error::other("GGV can index BAM, CRAM, bgzipped VCF and FASTA (plain or bgzipped) files"),
    })?);
    match index_kind(path) {
        Some("bai") => {
            match bam::fs::index(path).at(path)? {
                bam::Index::Bai(index) => bai::fs::write(&dst, &index).at(&dst)?,
                // Contigs too long for BAI: pasta builds a CSI, saved as such.
                bam::Index::Csi(index) => {
                    let csi = path.with_added_extension("csi");
                    pasta::csi::fs::write(&csi, &index).at(&csi)?;
                    return Ok(csi);
                }
            }
        }
        Some("crai") => {
            let index = cram::fs::index(path).at(path)?;
            cram::crai::fs::write(&dst, &index).at(&dst)?;
        }
        Some("tbi") => index_vcf(path, &dst)?,
        Some("fai") if is_gzip(path)? => {
            // samtools faidx: a .gzi of the BGZF blocks, and a .fai in uncompressed offsets.
            let gzi = path.with_added_extension("gzi");
            let blocks = bgzf_blocks(path)?;
            let mut w = BufWriter::new(File::create(&gzi).at(&gzi)?);
            // The first block (0, 0) is implied.
            let entries = &blocks[1.min(blocks.len())..];
            w.write_all(&(entries.len() as u64).to_le_bytes()).at(&gzi)?;
            for (c, u) in entries {
                w.write_all(&c.to_le_bytes()).at(&gzi)?;
                w.write_all(&u.to_le_bytes()).at(&gzi)?;
            }
            w.flush().at(&gzi)?;
            write_fai(BufReader::new(bgzf::io::Reader::new(File::open(path).at(path)?)), path, &dst)?;
        }
        Some("fai") => {
            write_fai(BufReader::new(File::open(path).at(path)?), path, &dst)?;
        }
        _ => unreachable!("index_kind checked above"),
    }
    Ok(dst)
}

/// Writes the `.fai` of the FASTA text `reader` gives (offsets into that text).
fn write_fai(reader: impl BufRead, path: &Path, dst: &Path) -> Result<()> {
    let mut indexer = fasta::io::Indexer::new(reader);
    let mut records = Vec::new();
    while let Some(r) = indexer.index_record().map_err(|e| Error::Io { path: path.to_owned(), cause: std::io::Error::other(e.to_string()) })? {
        records.push(r);
    }
    fasta::fai::fs::write(dst, &fasta::fai::Index::from(records)).at(dst)
}

/// Whether `path` is gzip-compressed (and if so, that it is BGZF, which can be read at a position).
pub(crate) fn is_gzip(path: &Path) -> Result<bool> {
    let mut header = [0u8; 18];
    let mut f = File::open(path).at(path)?;
    let n = f.read(&mut header).at(path)?;
    if n < 2 || header[0..2] != [0x1f, 0x8b] {
        return Ok(false);
    }
    if n < 18 || header[2..4] != [8, 4] || header[12..14] != *b"BC" {
        return Err(Error::Io { path: path.to_owned(), cause: std::io::Error::other("gzip but not bgzip: recompress it with `bgzip` to index it") });
    }
    Ok(true)
}

/// Each BGZF block's (compressed, uncompressed) start offset, the empty end-of-file block left out.
fn bgzf_blocks(path: &Path) -> Result<Vec<(u64, u64)>> {
    let mut f = BufReader::new(File::open(path).at(path)?);
    let len = f.get_ref().metadata().at(path)?.len();
    let (mut c, mut u) = (0u64, 0u64);
    let mut blocks = Vec::new();
    let bad = |why: &str, at: u64| Error::Io { path: path.to_owned(), cause: std::io::Error::other(format!("not valid bgzip ({why} at byte {at})")) };
    while c < len {
        let mut header = [0u8; 18];
        f.read_exact(&mut header).map_err(|_| bad("a cut-off block", c))?;
        if header[0..4] != [0x1f, 0x8b, 8, 4] || header[12..14] != *b"BC" {
            return Err(bad("a block that isn't BGZF", c));
        }
        let bsize = u64::from(u16::from_le_bytes([header[16], header[17]])) + 1;
        if bsize < 26 || c + bsize > len {
            return Err(bad("a cut-off block", c));
        }
        f.seek(SeekFrom::Start(c + bsize - 4)).at(path)?;
        let mut isize = [0u8; 4];
        f.read_exact(&mut isize).at(path)?;
        let isize = u64::from(u32::from_le_bytes(isize));
        if isize > 0 {
            blocks.push((c, u));
        }
        c += bsize;
        u += isize;
    }
    Ok(blocks)
}

/// Tabix-indexes a bgzipped VCF: each line with the BGZF virtual offset where it ends.
fn index_vcf(path: &Path, dst: &Path) -> Result<()> {
    let mut reader = bgzf::io::Reader::new(File::open(path).at(path)?);
    let mut builder = TabixBuilder::new(TabixConfig::VCF, 0);
    let mut line = Vec::new();
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line).at(path)? == 0 {
            break;
        }
        while line.last().is_some_and(|b| *b == b'\n' || *b == b'\r') {
            line.pop();
        }
        builder.push_line(&line, u64::from(reader.virtual_position())).at(path)?;
    }
    let tabix = builder.finish(u64::from(reader.virtual_position()));
    tabix.index.save(dst, IndexFormat::Tbi).at(dst)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kinds() {
        assert_eq!(index_kind(Path::new("a.BAM")), Some("bai"));
        assert_eq!(index_kind(Path::new("a.cram")), Some("crai"));
        assert_eq!(index_kind(Path::new("calls.vcf.gz")), Some("tbi"));
        assert_eq!(index_kind(Path::new("calls.vcf")), None, "plain VCF must be bgzipped first");
        assert_eq!(index_kind(Path::new("ref.fa")), Some("fai"));
        assert_eq!(index_kind(Path::new("ref.fa.gz")), Some("fai"), "and a .gzi");
        assert_eq!(index_kind(Path::new("ref.fasta.bgz")), Some("fai"));
    }

    /// A bgzipped FASTA indexed here gets the same .fai and .gzi samtools wrote for the fixture.
    #[test]
    fn bgzipped_fasta_indexes_match_samtools() {
        let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio");
        let src = data.join("GRCh38.chr20.slices.fa.gz");
        if !src.exists() {
            return;
        }
        let dir = std::env::temp_dir().join(format!("ggv-faidx-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("ref.fa.gz");
        std::fs::copy(&src, &fa).unwrap();
        assert_eq!(build_index(&fa).unwrap(), fa.with_added_extension("fai"));
        for ext in ["fai", "gzi"] {
            let (ours, theirs) = (std::fs::read(fa.with_added_extension(ext)).unwrap(), std::fs::read(src.with_added_extension(ext)).unwrap());
            assert!(ours == theirs, ".{ext} differs ({} vs {} bytes)", ours.len(), theirs.len());
        }
        let region = ggv_model::Region { contig: "chr20".into(), start: 58_880_000, end: 58_880_100 };
        assert_eq!(crate::ReferenceFasta::open(&fa).unwrap().fetch(&region).unwrap(), crate::ReferenceFasta::open(&src).unwrap().fetch(&region).unwrap());
        // Plain gzip can't be read at a position: refused, with what to do.
        let gz = dir.join("plain.fa.gz");
        let mut e = flate2::write::GzEncoder::new(File::create(&gz).unwrap(), flate2::Compression::fast());
        e.write_all(b">c
ACGT
").unwrap();
        e.finish().unwrap();
        assert!(build_index(&gz).unwrap_err().to_string().contains("bgzip"));
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn missing_file_is_not_a_missing_index() {
        // A BAM or CRAM that isn't there is reported as missing, not as lacking an index.
        let dir = std::env::temp_dir();
        for missing in [dir.join("ggv-no-such-file.bam"), dir.join("ggv-no-such-file.cram")] {
            match crate::BamFile::open(&missing) {
                Err(Error::Io { cause, .. }) => assert_eq!(cause.kind(), std::io::ErrorKind::NotFound),
                Err(e) => panic!("{e}"),
                Ok(_) => panic!("opened a missing file"),
            }
        }
    }

    #[test]
    fn rebuilt_indexes_answer_like_the_originals() {
        // GGV_TEST_INDEX_DIR: a scratch directory. Copies the HG002 fixtures (GGV_TEST_BAM,
        // GGV_TEST_VCF) there without their indexes, rebuilds them, and checks queries agree.
        let (Ok(dir), Ok(bam), Ok(vcf)) =
            (std::env::var("GGV_TEST_INDEX_DIR"), std::env::var("GGV_TEST_BAM"), std::env::var("GGV_TEST_VCF"))
        else {
            return;
        };
        let dir = PathBuf::from(dir);
        std::fs::create_dir_all(&dir).unwrap();
        let region = ggv_model::Region { contig: "chr20".into(), start: 1_890_000, end: 1_940_000 };

        let bam_copy = dir.join("copy.bam");
        std::fs::copy(&bam, &bam_copy).unwrap();
        let _ = std::fs::remove_file(bam_copy.with_added_extension("bai"));
        assert!(matches!(crate::BamFile::open(&bam_copy), Err(Error::MissingIndex { .. })));
        build_index(&bam_copy).unwrap();
        let count = |p: &Path| {
            let mut n = 0;
            crate::BamFile::open(p).unwrap().for_each_record(&region, |_| { n += 1; Ok(()) }).unwrap();
            n
        };
        assert_eq!(count(&bam_copy), count(Path::new(&bam)));

        let vcf_copy = dir.join("copy.vcf.gz");
        std::fs::copy(&vcf, &vcf_copy).unwrap();
        let _ = std::fs::remove_file(vcf_copy.with_added_extension("tbi"));
        build_index(&vcf_copy).unwrap();
        let sites = |p: &Path| crate::VariantFile::open(p).unwrap().query(&region, true).unwrap().len();
        assert_eq!(sites(&vcf_copy), sites(Path::new(&vcf)));
        eprintln!("rebuilt BAM and VCF indexes answer {} reads / {} sites like the originals", count(&bam_copy), sites(&vcf_copy));
    }
}
