use std::{
    fs::File,
    io,
    path::{Path, PathBuf},
};

use ggv_model::{Contig, Locus, Region};
use pasta::{
    bam::{self, bai},
    bgzf, cram, csi, fasta,
    core::{Position, Region as PastaRegion},
    sam,
    sam::header::record::value::map::read_group::tag as rg_tag,
};

use crate::{
    Result,
    error::{Error, IoContext},
    resolve_locus,
};

/// An indexed BAM or CRAM file opened for region queries. Both yield the same encoded
/// records, so everything downstream (pileup, layout, depth) is format-agnostic.
pub struct BamFile {
    path: PathBuf,
    /// The FASTA a CRAM decodes against (kept to reopen it).
    reference_path: Option<PathBuf>,
    index_path: PathBuf,
    inner: Inner,
    header: sam::Header,
    contigs: Vec<Contig>,
    /// The file's size when opened (a BAM's index must point inside it).
    len: u64,
    /// A BAM without the BGZF end-of-file marker: probably cut short.
    missing_eof: bool,
}

/// The empty BGZF block every complete BAM ends with.
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

fn ends_with_eof_marker(path: &Path, len: u64) -> bool {
    use std::io::{Read, Seek, SeekFrom};
    let mut tail = [0u8; 28];
    len >= 28
        && File::open(path).and_then(|mut f| {
            f.seek(SeekFrom::Start(len - 28))?;
            f.read_exact(&mut tail)
        }).is_ok()
        && tail == BGZF_EOF
}

enum Inner {
    /// The index is shared between handles on the same file ([`BamFile::reopen`]).
    Bam { reader: bam::io::Reader<bgzf::io::Reader<File>>, index: std::sync::Arc<bam::Index> },
    /// `has_reference`: records can't be decoded without it (the header can).
    Cram { reader: cram::io::Reader<File>, index: cram::crai::Index, has_reference: bool },
}

fn is_cram(path: &Path) -> bool {
    path.extension().is_some_and(|e| e.eq_ignore_ascii_case("cram"))
}

fn contigs_of(header: &sam::Header) -> Vec<Contig> {
    header.reference_sequences().iter().map(|(name, rs)| Contig { name: name.to_string(), length: rs.length().get() as u64 }).collect()
}

fn info_of(path: &Path, index: String, header: &sam::Header) -> BamInfo {
    let read_groups = header
        .read_groups()
        .iter()
        .map(|(id, rg)| ReadGroup {
            id: id.to_string(),
            sample: rg.other_fields().get(&rg_tag::SAMPLE).map(|v| v.to_string()),
            platform: rg.other_fields().get(&rg_tag::PLATFORM).map(|v| v.to_string()),
        })
        .collect();
    BamInfo { path: path.display().to_string(), index, contigs: contigs_of(header), read_groups }
}

/// What an agent or UI needs to know about a BAM before querying it.
#[derive(Clone, Debug, serde::Serialize)]
pub struct BamInfo {
    pub path: String,
    pub index: String,
    pub contigs: Vec<Contig>,
    pub read_groups: Vec<ReadGroup>,
}

#[derive(Clone, Debug, serde::Serialize)]
pub struct ReadGroup {
    pub id: String,
    pub sample: Option<String>,
    pub platform: Option<String>,
}

impl BamFile {
    /// Opens a BAM and its index. The index is looked for at `x.bam.bai`, `x.bam.csi` and
    /// `x.bai`, the same places samtools and IGV look.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        Self::open_with_reference(path, None::<&Path>)
    }

    /// Opens a BAM, or a CRAM with the indexed FASTA it was compressed against (`.crai` index
    /// next to it). The reference is ignored for BAM.
    pub fn open_with_reference(path: impl AsRef<Path>, reference: Option<impl AsRef<Path>>) -> Result<Self> {
        let path = path.as_ref().to_owned();
        // The file itself first: a missing BAM is not a missing index.
        let meta = std::fs::metadata(&path).at(&path)?;
        if meta.is_dir() {
            return Err(Error::Io { path, cause: std::io::Error::other("is a folder, not an alignment file") });
        }
        if is_cram(&path) {
            return Self::open_cram(path, reference.as_ref().map(AsRef::as_ref));
        }
        let (index, index_path) = read_index(&path)?;

        let mut reader = File::open(&path).map(bam::io::Reader::new).at(&path)?;
        let header = reader.read_header().at(&path)?;
        let mut file = Self::new(path, index_path, Inner::Bam { reader, index: std::sync::Arc::new(index) }, header);
        file.len = meta.len();
        file.missing_eof = !ends_with_eof_marker(&file.path, file.len);
        Ok(file)
    }

    fn open_cram(path: PathBuf, reference: Option<&Path>) -> Result<Self> {
        if let Some(r) = reference.filter(|r| r.to_string_lossy().to_ascii_lowercase().ends_with(".2bit")) {
            return Err(Error::Io { path: r.to_owned(), cause: std::io::Error::other("CRAM decoding needs the reference as FASTA (.fa with .fai); a .2bit works for BAM") });
        }
        let candidates = [path.with_added_extension("crai"), path.with_extension("crai")];
        let Some(index_path) = candidates.iter().find(|p| p.exists()).cloned() else {
            return Err(Error::MissingIndex { path: path.clone(), tried: "x.cram.crai, x.crai".into(), how: "create one with `samtools index`".into() });
        };
        let index = cram::crai::fs::read(&index_path).at(&index_path)?;
        let repository = match reference {
            Some(r) => {
                let fasta_reader = fasta::io::indexed_reader::Builder::default().build_from_path(r).at(r)?;
                fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(fasta_reader))
            }
            None => fasta::Repository::default(),
        };
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_path(&path)
            .at(&path)?;
        let header = reader.read_header().at(&path)?;
        let mut file = Self::new(path, index_path, Inner::Cram { reader, index, has_reference: reference.is_some() }, header);
        file.reference_path = reference.map(Path::to_owned);
        Ok(file)
    }

    fn new(path: PathBuf, index_path: PathBuf, inner: Inner, header: sam::Header) -> Self {
        let contigs = contigs_of(&header);
        Self { path, reference_path: None, index_path, inner, header, contigs, len: 0, missing_eof: false }
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    /// A warning when the file looks cut short (a BAM without its end-of-file marker). Reads
    /// the index places past the end are an error when queried; this flags the rest.
    pub fn warning(&self) -> Option<String> {
        self.missing_eof.then(|| format!("{}: no BGZF end-of-file marker; the file may be truncated (an incomplete download or copy?)", self.path.display()))
    }

    /// The same file opened again (with the same reference): another handle for another
    /// thread, since a reader keeps its own position.
    pub fn reopen(&self) -> Result<Self> {
        // A BAM's parsed index and header are shared, not read again: chunked passes reopen the
        // file for every thread on every wide load, and a whole-genome index is tens of MB.
        if let Inner::Bam { index, .. } = &self.inner {
            let mut reader = File::open(&self.path).map(bam::io::Reader::new).at(&self.path)?;
            reader.read_header().at(&self.path)?; // positions the reader past the header
            return Ok(Self {
                path: self.path.clone(),
                reference_path: self.reference_path.clone(),
                index_path: self.index_path.clone(),
                inner: Inner::Bam { reader, index: index.clone() },
                header: self.header.clone(),
                contigs: self.contigs.clone(),
                len: self.len,
                missing_eof: self.missing_eof,
            });
        }
        Self::open_with_reference(&self.path, self.reference_path.as_ref())
    }

    pub fn header(&self) -> &sam::Header {
        &self.header
    }

    pub fn contigs(&self) -> &[Contig] {
        &self.contigs
    }

    pub fn resolve(&self, locus: &Locus) -> Result<Region> {
        resolve_locus(&self.contigs, locus)
    }

    /// Parses and resolves locus text; see [`crate::resolve_locus_str`].
    pub fn resolve_str(&self, text: &str) -> Result<Region> {
        crate::resolve_locus_str(&self.contigs, text)
    }

    pub fn info(&self) -> BamInfo {
        info_of(&self.path, self.index_path.display().to_string(), &self.header)
    }

    /// [`Self::info`] from the header alone, with index `none`: for a file without an index
    /// (unaligned or unsorted reads can't have one), which [`Self::open`] refuses.
    pub fn header_info(path: impl AsRef<Path>) -> Result<BamInfo> {
        let path = path.as_ref();
        if std::fs::metadata(path).at(path)?.is_dir() {
            return Err(Error::Io { path: path.to_owned(), cause: std::io::Error::other("is a folder, not an alignment file") });
        }
        let header = if is_cram(path) {
            cram::io::reader::Builder::default().build_from_path(path).at(path)?.read_header().at(path)?
        } else {
            File::open(path).map(bam::io::Reader::new).at(path)?.read_header().at(path)?
        };
        Ok(info_of(path, "none".into(), &header))
    }

    /// Approximate bytes of alignment data per window of `contig`, read from the index alone
    /// (no records are decoded), as `(window length, bytes per window)`. BAI/CSI: the distance
    /// between the first record offsets of consecutive windows; CRAI: each slice's size spread
    /// over the span it covers. Proportional to read density, so it gives whole-chromosome
    /// coverage at a glance; calibrate it with a small real sample to get depth.
    pub fn index_density(&self, contig: &str) -> Option<(u64, Vec<f64>)> {
        let id = self.header.reference_sequences().get_index_of(contig.as_bytes())?;
        let len = self.contigs.get(id)?.length;
        match &self.inner {
            Inner::Bam { index, .. } => {
                use csi::BinningIndex;
                let (min_shift, depth) = (index.min_shift(), index.depth());
                let window = 1u64 << min_shift;
                let n = len.div_ceil(window) as usize;
                let virtual_offsets: Vec<bgzf::VirtualPosition> = match &**index {
                    bam::Index::Bai(ix) => {
                        let rs = ix.reference_sequences().get(id)?;
                        (0..=n).map(|w| rs.min_offset(min_shift, depth, Position::try_from(w * window as usize + 1).expect(">= 1"))).collect()
                    }
                    bam::Index::Csi(ix) => {
                        let rs = ix.reference_sequences().get(id)?;
                        (0..=n).map(|w| rs.min_offset(min_shift, depth, Position::try_from(w * window as usize + 1).expect(">= 1"))).collect()
                    }
                };
                // A virtual offset in compressed bytes: the block's start plus the offset within
                // the block scaled by this file's compression ratio. Windows inside one block
                // differ only in the second part, so the ratio must be the file's own.
                let ratio = compression_ratio(&self.path, &virtual_offsets).unwrap_or(0.3);
                let at = |v: &bgzf::VirtualPosition| v.compressed() as f64 + f64::from(v.uncompressed()) * ratio;
                let offsets: Vec<f64> = virtual_offsets.iter().map(at).collect();
                // Windows past the last record point back at offset 0 (or stay put): no data.
                Some((window, offsets.windows(2).map(|w| (w[1] - w[0]).max(0.0)).collect()))
            }
            Inner::Cram { index, .. } => {
                const WINDOW: u64 = 1 << 14;
                let mut bytes = vec![0.0; len.div_ceil(WINDOW) as usize];
                for r in index.iter().filter(|r| r.reference_sequence_id() == Some(id)) {
                    let Some(start) = r.alignment_start() else { continue };
                    let s = usize::from(start) as u64 - 1;
                    let e = (s + r.alignment_span() as u64).min(len).max(s + 1);
                    let per_base = r.slice_length() as f64 / (e - s) as f64;
                    let mut w = s / WINDOW;
                    while w * WINDOW < e && (w as usize) < bytes.len() {
                        let (a, b) = ((w * WINDOW).max(s), ((w + 1) * WINDOW).min(e));
                        bytes[w as usize] += per_base * (b - a) as f64;
                        w += 1;
                    }
                }
                Some((WINDOW, bytes))
            }
        }
    }

    /// Calls `f` with every record that overlaps `region`, in coordinate order. The record
    /// buffer is reused between calls; copy out what you need to keep.
    pub fn for_each_record(
        &mut self,
        region: &Region,
        mut f: impl FnMut(&bam::Record) -> io::Result<()>,
    ) -> Result<()> {
        let start = Position::try_from(region.start as usize + 1).expect("start is >= 1");
        let end = Position::try_from(region.end as usize).expect("end is >= 1");
        let pasta_region = PastaRegion::new(region.contig.as_str(), start..=end);

        let path = self.path.clone();
        let mut record = bam::Record::default();
        match &mut self.inner {
            Inner::Bam { reader, index } => {
                check_chunks(&path, self.len, &self.header, index, &pasta_region)?;
                let mut query = reader.query(&self.header, &**index, &pasta_region).at(&path)?;
                while query.read_record(&mut record).at(&path)? != 0 {
                    f(&record).at(&path)?;
                }
            }
            Inner::Cram { has_reference: false, .. } => {
                return Err(Error::Io {
                    path,
                    cause: io::Error::other("CRAM needs the reference FASTA it was made with to decode reads: open it (in the app, Open > Reference; on the command line, -f)"),
                });
            }
            Inner::Cram { reader, index, .. } => {
                let mut query = reader.query(&self.header, index, &pasta_region).at(&path)?;
                while query.read_record(&mut record).at(&path)? != 0 {
                    f(&record).at(&path)?;
                }
            }
        }
        Ok(())
    }
}

impl BamFile {
    /// One record as a SAM line (no trailing newline), with this file's header for names.
    pub fn sam_line(&self, record: &bam::Record) -> io::Result<String> {
        use sam::alignment::io::Write as _;
        let mut writer = sam::io::Writer::new(Vec::new());
        writer.write_alignment_record(&self.header, record)?;
        let line = String::from_utf8_lossy(&writer.into_inner()).into_owned();
        Ok(line.trim_end_matches(['\n', '\r']).to_owned())
    }

    /// SAM lines for records overlapping `region` that `keep` accepts, at most `limit`.
    pub fn sam_lines(&mut self, region: &Region, limit: usize, mut keep: impl FnMut(&bam::Record) -> bool) -> Result<Vec<String>> {
        let mut out = Vec::new();
        let mut lines = Vec::new();
        let result = self.for_each_record(region, |r| {
            if keep(r) {
                lines.push(r.clone());
                if lines.len() >= limit {
                    return Err(io::Error::other(STOP));
                }
            }
            Ok(())
        });
        stopped_ok(result)?;
        for r in &lines {
            out.push(self.sam_line(r).at(&self.path)?);
        }
        Ok(out)
    }

    /// The first record named `name` in `region`, as a SAM line. Checks `cancel` every 64k
    /// records so a long scan can be stopped.
    pub fn find_by_name(&mut self, region: &Region, name: &str, cancel: &std::sync::atomic::AtomicBool) -> Result<Option<String>> {
        use std::sync::atomic::Ordering;
        let mut found = None;
        let mut n = 0u32;
        let result = self.for_each_record(region, |r| {
            n = n.wrapping_add(1);
            if n & 0xffff == 0 && cancel.load(Ordering::Relaxed) {
                return Err(io::Error::other(STOP));
            }
            if r.name().is_some_and(|x| <[u8]>::eq(x.as_ref(), name.as_bytes())) {
                found = Some(r.clone());
                return Err(io::Error::other(STOP));
            }
            Ok(())
        });
        stopped_ok(result)?;
        found.map(|r| self.sam_line(&r).at(&self.path)).transpose()
    }
}

/// Fails when the index places `region`'s records past the end of the file: reading would
/// stop at the end and quietly return fewer reads (or none) instead of an error.
fn check_chunks(path: &Path, len: u64, header: &sam::Header, index: &bam::Index, region: &PastaRegion) -> Result<()> {
    use csi::BinningIndex;
    let Some(id) = header.reference_sequences().get_index_of(region.name()) else { return Ok(()) };
    let Ok(chunks) = index.query(id, region.interval()) else { return Ok(()) };
    // A chunk ends at the start of the block after its last record: the end of the file at most.
    let past = |v: bgzf::VirtualPosition| v.compressed() > len || (v.compressed() == len && v.uncompressed() > 0);
    if chunks.iter().any(|c| c.start().compressed() >= len || past(c.end())) {
        return Err(Error::Io {
            path: path.to_owned(),
            cause: io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!("the index places reads past the end of the file ({len} bytes): the file is truncated, or the index belongs to another file"),
            ),
        });
    }
    Ok(())
}

/// Marks an early stop from inside a record callback (not a real error).
const STOP: &str = "ggv: stop";

fn stopped_ok(result: Result<()>) -> Result<()> {
    match result {
        Err(Error::Io { cause, .. }) if cause.to_string() == STOP => Ok(()),
        other => other,
    }
}

/// Compressed bytes per uncompressed byte, from the headers and trailers of up to 32 BGZF
/// blocks at the given offsets (each block records both sizes).
fn compression_ratio(path: &Path, offsets: &[bgzf::VirtualPosition]) -> Option<f64> {
    use std::io::{Read, Seek, SeekFrom};
    let mut starts: Vec<u64> = offsets.iter().map(|v| v.compressed()).filter(|&c| c > 0).collect();
    starts.dedup();
    if starts.is_empty() {
        return None;
    }
    let step = starts.len().div_ceil(32);
    let mut file = File::open(path).ok()?;
    let (mut compressed, mut uncompressed) = (0u64, 0u64);
    for &c in starts.iter().step_by(step) {
        let mut header = [0u8; 18];
        file.seek(SeekFrom::Start(c)).ok()?;
        file.read_exact(&mut header).ok()?;
        // gzip magic, FEXTRA, and the BGZF "BC" subfield holding BSIZE (block size - 1).
        if header[0..4] != [0x1f, 0x8b, 8, 4] || header[12..14] != *b"BC" {
            continue;
        }
        let bsize = u64::from(u16::from_le_bytes([header[16], header[17]])) + 1;
        let mut isize = [0u8; 4];
        file.seek(SeekFrom::Start(c + bsize - 4)).ok()?;
        file.read_exact(&mut isize).ok()?;
        compressed += bsize;
        uncompressed += u64::from(u32::from_le_bytes(isize));
    }
    (uncompressed > 0).then(|| compressed as f64 / uncompressed as f64)
}

fn read_index(bam: &Path) -> Result<(bam::Index, PathBuf)> {
    let candidates = [
        (bam.with_added_extension("bai"), false),
        (bam.with_added_extension("csi"), true),
        (bam.with_extension("bai"), false),
    ];

    for (candidate, is_csi) in &candidates {
        let result = if *is_csi {
            csi::fs::read(candidate).map(bam::Index::Csi)
        } else {
            bai::fs::read(candidate).map(bam::Index::Bai)
        };
        match result {
            Ok(index) => return Ok((index, candidate.clone())),
            Err(e) if e.kind() == io::ErrorKind::NotFound => continue,
            Err(e) => return Err(Error::Io { path: candidate.clone(), cause: e }),
        }
    }

    Err(Error::MissingIndex {
        path: bam.to_owned(),
        tried: candidates
            .iter()
            .map(|(p, _)| p.file_name().unwrap_or_default().to_string_lossy().into_owned())
            .collect::<Vec<_>>()
            .join(", "),
        how: "create one with `samtools index`".into(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Copies of the bundled BAM: whole, without its end-of-file marker, and cut short (each
    /// with the original index).
    #[test]
    fn truncated_files_are_reported() {
        let src = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio/HG002.revio.slices.bam");
        let bytes = std::fs::read(&src).unwrap();
        let bai = std::fs::read(src.with_added_extension("bai")).unwrap();
        let dir = std::env::temp_dir().join(format!("ggv-truncated-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let copy = |name: &str, len: usize| {
            let p = dir.join(name);
            std::fs::write(&p, &bytes[..len]).unwrap();
            std::fs::write(p.with_added_extension("bai"), &bai).unwrap();
            BamFile::open(&p).unwrap()
        };
        let count = |b: &mut BamFile| -> Result<usize> {
            let region = b.resolve_str("chr20:58885000-58886000")?;
            let mut n = 0;
            b.for_each_record(&region, |_| {
                n += 1;
                Ok(())
            })?;
            Ok(n)
        };

        let mut whole = copy("whole.bam", bytes.len());
        assert!(whole.warning().is_none());
        assert_eq!(count(&mut whole).unwrap(), 34);

        // Only the marker missing: the reads are all there, with a warning.
        let mut no_eof = copy("no-eof.bam", bytes.len() - 28);
        assert!(no_eof.warning().is_some());
        assert_eq!(count(&mut no_eof).unwrap(), 34);

        // The index points past the end: an error, not zero reads.
        let mut cut = copy("cut.bam", bytes.len() - 2_000_000);
        assert!(cut.warning().is_some());
        assert!(count(&mut cut).unwrap_err().to_string().contains("past the end of the file"));

        // Cut inside the block holding the region's reads: an error too.
        let region = whole.resolve_str("chr20:58885000-58886000").unwrap();
        let mut first = None;
        let mut offsets = Vec::new();
        if let Inner::Bam { index, .. } = &whole.inner {
            use csi::BinningIndex;
            let id = whole.header.reference_sequences().get_index_of(region.contig.as_bytes()).unwrap();
            let pasta = PastaRegion::new(region.contig.as_str(), Position::try_from(region.start as usize + 1).unwrap()..=Position::try_from(region.end as usize).unwrap());
            offsets = index.query(id, pasta.interval()).unwrap().iter().map(|c| c.start().compressed()).collect();
            first = offsets.iter().min().copied();
        }
        let mut mid = copy("mid.bam", first.unwrap() as usize + 100);
        assert!(count(&mut mid).is_err(), "{offsets:?}");
        std::fs::remove_dir_all(&dir).ok();
    }
}
