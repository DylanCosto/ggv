//! A binary copy of a parsed annotation, so a large GTF or GFF3 is parsed once rather than on
//! every open: GENCODE takes seconds and hundreds of MB to parse, and its cache opens in
//! milliseconds holding a few MB.
//!
//! The cache keeps each contig's transcripts as a separately compressed block, read the first
//! time something asks about that contig. Gene spans and gene names stay in memory; transcript
//! ids are looked up in a sorted table in the file. It lives in the user's cache directory
//! (`GGV_ANNOTATION_CACHE` names another directory, or `off`), never next to the data, and is
//! keyed by the source's path, size and modification time, and by the parsing code: a changed
//! file, or an GGV whose parser changed, parses it again.
//!
//! Layout (little-endian): magic, version, the header's offset; the contig blocks; the
//! transcript-id table and its key bytes; then the header (source, contigs, genes, gene keys).

use std::{
    fs::File,
    io::{Read, Seek, SeekFrom, Write},
    path::{Path, PathBuf},
    sync::{Mutex, OnceLock},
    time::UNIX_EPOCH,
};

use crate::{Annotation, GeneSpan, Keys, Loaded, Part, Strand, Transcript, Translation};

const MAGIC: &[u8; 8] = b"GVANNOT\0";
const VERSION: u32 = 3;
/// The parsing code, hashed when GGV is built: a cache written by an GGV that parsed differently
/// isn't used (nobody has to remember to bump `VERSION` when parsing changes).
const PARSER: u64 = fnv(include_str!("parse.rs").as_bytes()) ^ fnv(include_str!("lib.rs").as_bytes()).rotate_left(1) ^ fnv(include_str!("cache.rs").as_bytes()).rotate_left(2)
    ^ fnv(include_str!("mane.rs").as_bytes()).rotate_left(3);

/// FNV-1a: stable across builds and platforms.
const fn fnv(bytes: &[u8]) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    let mut i = 0;
    while i < bytes.len() {
        h = (h ^ bytes[i] as u64).wrapping_mul(0x0100_0000_01b3);
        i += 1;
    }
    h
}
/// Smaller files parse faster than a cache would open.
const MIN_SOURCE_BYTES: u64 = 4 << 20;

/// Where `source`'s cache goes, or `None` when caching is off or the file is small.
pub(crate) fn path_for(source: &Path) -> Option<PathBuf> {
    // `GV_ANNOTATION_CACHE`: its name before the rename to GGV.
    let dir = match std::env::var_os("GGV_ANNOTATION_CACHE").or_else(|| std::env::var_os("GV_ANNOTATION_CACHE")) {
        Some(d) if d == "off" => return None,
        Some(d) => PathBuf::from(d),
        None => {
            // As `ggv_io::folders::cache_dir`: the cache folder from before the rename (`gv`)
            // while there is no `ggv` one.
            let base = cache_dir()?;
            let (dir, old) = (base.join("ggv"), base.join("gv"));
            if !dir.exists() && old.is_dir() { old } else { dir }.join("annotations")
        }
    };
    let meta = std::fs::metadata(source).ok()?;
    if meta.len() < MIN_SOURCE_BYTES {
        return None;
    }
    let abs = std::fs::canonicalize(source).ok()?;
    // A hash of the absolute path names the file (the header holds the full key).
    let hash = fnv(abs.to_string_lossy().as_bytes());
    let name = abs.file_name().unwrap_or_default().to_string_lossy();
    Some(dir.join(format!("{name}.{hash:016x}.ggvannot")))
}

fn cache_dir() -> Option<PathBuf> {
    if cfg!(windows) {
        std::env::var_os("LOCALAPPDATA").map(PathBuf::from)
    } else if cfg!(target_os = "macos") {
        std::env::var_os("HOME").map(|h| PathBuf::from(h).join("Library/Caches"))
    } else {
        std::env::var_os("XDG_CACHE_HOME").map(PathBuf::from).or_else(|| std::env::var_os("HOME").map(|h| PathBuf::from(h).join(".cache")))
    }
}

/// What the cache must match: the source's absolute path, size and modification time, and the
/// parser that read it.
#[derive(PartialEq, Eq)]
pub(crate) struct SourceKey {
    path: String,
    size: u64,
    mtime_ns: u128,
    parser: u64,
}

impl SourceKey {
    /// Taken before parsing, so a file that changes while it is parsed doesn't get a cache
    /// claiming its new state.
    pub(crate) fn of(source: &Path) -> std::io::Result<Self> {
        let meta = std::fs::metadata(source)?;
        let mtime_ns = meta.modified()?.duration_since(UNIX_EPOCH).map_or(0, |d| d.as_nanos());
        Ok(Self { path: std::fs::canonicalize(source)?.to_string_lossy().into_owned(), size: meta.len(), mtime_ns, parser: PARSER })
    }
}

// ---- writing ----

#[derive(Default)]
struct Buf(Vec<u8>);

impl Buf {
    fn u8(&mut self, v: u8) {
        self.0.push(v);
    }
    fn u16(&mut self, v: u16) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn u32(&mut self, v: u32) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn u64(&mut self, v: u64) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    fn u128(&mut self, v: u128) {
        self.0.extend_from_slice(&v.to_le_bytes());
    }
    /// Strings longer than 64 KB are cut (no annotation field comes close).
    fn str(&mut self, s: &str) {
        let b = &s.as_bytes()[..s.len().min(u16::MAX as usize)];
        self.u16(b.len() as u16);
        self.0.extend_from_slice(b);
    }
}

fn encode_transcripts(ts: &[Transcript]) -> std::io::Result<Vec<u8>> {
    let mut b = Buf::default();
    for t in ts {
        b.u64(t.start);
        b.u64(t.end);
        b.u8(match t.strand {
            Strand::Forward => 0,
            Strand::Reverse => 1,
            Strand::Unknown => 2,
        });
        b.u8(t.tags);
        b.u8(u8::from(t.representative));
        match t.cds {
            Some((s, e)) => {
                b.u8(1);
                b.u64(s);
                b.u64(e);
            }
            None => b.u8(0),
        }
        for s in [&t.id, &t.gene_id, &t.gene_name, &t.biotype] {
            b.str(s);
        }
        b.u32(t.exons.len() as u32);
        for &(s, e) in &t.exons {
            // Exons sit inside their transcript, so offsets from its start fit in 32 bits; a
            // file where one doesn't isn't cached (it is still read as parsed).
            let offset = |p: u64| {
                p.checked_sub(t.start)
                    .and_then(|d| u32::try_from(d).ok())
                    .ok_or_else(|| std::io::Error::other(format!("transcript {}: an exon outside it at {p}", t.id)))
            };
            b.u32(offset(s)?);
            b.u32(offset(e)?);
        }
        let tr = &t.translation;
        b.u8(tr.phase);
        b.u8(tr.table);
        b.u8(u8::from(t.exons_from_span));
        b.u16(tr.selenocysteine.len().min(u16::MAX as usize) as u16);
        for &(s, e) in tr.selenocysteine.iter().take(u16::MAX as usize) {
            b.u64(s);
            b.u64(e);
        }
    }
    let mut z = flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::fast());
    z.write_all(&b.0)?;
    z.finish()
}

/// Writes `ann` (fully in memory, parsed from the source `key` describes) as its cache,
/// atomically.
pub(crate) fn write(ann: &Annotation, key: &SourceKey, cache: &Path) -> std::io::Result<()> {
    let mut body = Buf::default();
    body.0.extend_from_slice(MAGIC);
    body.u32(VERSION);
    body.u64(0); // the header's offset, filled in below

    let mut blocks = Vec::with_capacity(ann.parts.len());
    for part in &ann.parts {
        let ts = part.loaded.get().map_or(&[][..], |l| &l.transcripts[..]);
        let start = body.0.len() as u64;
        body.0.extend_from_slice(&encode_transcripts(ts)?);
        blocks.push((start, body.0.len() as u64 - start));
    }

    // Transcript ids: fixed-size entries (key offset, key length, gene), then the key bytes.
    let Keys::Memory(tx) = &ann.tx_keys else { return Err(std::io::Error::other("annotation already cached")) };
    let table = body.0.len() as u64;
    let blob = table + 14 * tx.len() as u64;
    let mut at = blob;
    for (k, g) in tx {
        body.u64(at);
        body.u16(k.len() as u16);
        body.u32(*g);
        at += k.len() as u64;
    }
    for (k, _) in tx {
        body.0.extend_from_slice(k.as_bytes());
    }

    let header = body.0.len() as u64;
    body.0[12..20].copy_from_slice(&header.to_le_bytes());
    body.str(&key.path);
    body.u64(key.size);
    body.u128(key.mtime_ns);
    body.u64(key.parser);
    body.u32(ann.contigs.len() as u32);
    for ((name, part), (off, len)) in ann.contigs.iter().zip(&ann.parts).zip(&blocks) {
        body.str(name);
        body.u32(part.count as u32);
        body.u64(part.length);
        body.u64(*off);
        body.u64(*len);
    }
    body.u32(ann.genes.len() as u32);
    for g in &ann.genes {
        body.str(&g.name);
        body.u32(ann.contigs.iter().position(|c| *c == g.contig).unwrap_or(0) as u32);
        body.u64(g.start);
        body.u64(g.end);
    }
    body.u32(ann.gene_keys.len() as u32);
    for (k, g) in &ann.gene_keys {
        body.str(k);
        body.u32(*g);
    }
    body.u32(tx.len() as u32);
    body.u64(table);

    std::fs::create_dir_all(cache.parent().unwrap_or(Path::new(".")))?;
    let tmp = cache.with_extension(format!("tmp{}", std::process::id()));
    std::fs::write(&tmp, &body.0)?;
    std::fs::rename(&tmp, cache).inspect_err(|_| {
        let _ = std::fs::remove_file(&tmp);
    })
}

// ---- reading ----

struct Cursor<'a>(&'a [u8]);

impl Cursor<'_> {
    fn take<const N: usize>(&mut self) -> Option<[u8; N]> {
        let (head, rest) = self.0.split_first_chunk::<N>()?;
        self.0 = rest;
        Some(*head)
    }
    fn u8(&mut self) -> Option<u8> {
        self.take::<1>().map(|b| b[0])
    }
    fn u16(&mut self) -> Option<u16> {
        self.take().map(u16::from_le_bytes)
    }
    fn u32(&mut self) -> Option<u32> {
        self.take().map(u32::from_le_bytes)
    }
    fn u64(&mut self) -> Option<u64> {
        self.take().map(u64::from_le_bytes)
    }
    fn u128(&mut self) -> Option<u128> {
        self.take().map(u128::from_le_bytes)
    }
    fn bytes(&mut self, n: usize) -> Option<&[u8]> {
        let (head, rest) = (self.0.get(..n)?, self.0.get(n..)?);
        self.0 = rest;
        Some(head)
    }
    fn str(&mut self) -> Option<Box<str>> {
        let n = self.u16()? as usize;
        Some(String::from_utf8_lossy(self.bytes(n)?).into())
    }
}

/// The open cache file: contig blocks and the transcript-id table are read from it on demand.
pub(crate) struct Reader {
    file: Mutex<File>,
    path: PathBuf,
}

impl Reader {
    fn read_at(&self, offset: u64, len: usize) -> std::io::Result<Vec<u8>> {
        let mut f = self.file.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
        f.seek(SeekFrom::Start(offset))?;
        let mut out = vec![0; len];
        f.read_exact(&mut out)?;
        Ok(out)
    }

    /// Decodes one contig's transcripts (empty, with a note, if the cache can't be read).
    pub(crate) fn load(&self, contig: u32, offset: u64, len: u64) -> Loaded {
        let decoded = self.read_at(offset, len as usize).ok().and_then(|z| {
            let mut raw = Vec::new();
            flate2::read::ZlibDecoder::new(&z[..]).read_to_end(&mut raw).ok()?;
            decode_transcripts(&raw, contig)
        });
        let transcripts = decoded.unwrap_or_else(|| {
            eprintln!("ggv: {}: unreadable; delete it to rebuild", self.path.display());
            Vec::new()
        });
        Loaded::new(transcripts)
    }

    /// Binary search of the transcript-id table (`count` entries at `table`).
    pub(crate) fn find(&self, key: &str, table: u64, count: u32) -> Option<u32> {
        let (mut lo, mut hi) = (0u64, u64::from(count));
        while lo < hi {
            let mid = (lo + hi) / 2;
            let e = self.read_at(table + 14 * mid, 14).ok()?;
            let mut c = Cursor(&e);
            let (at, len, gene) = (c.u64()?, c.u16()?, c.u32()?);
            let k = self.read_at(at, len as usize).ok()?;
            match k.as_slice().cmp(key.as_bytes()) {
                std::cmp::Ordering::Equal => return Some(gene),
                std::cmp::Ordering::Less => lo = mid + 1,
                std::cmp::Ordering::Greater => hi = mid,
            }
        }
        None
    }
}

fn decode_transcripts(raw: &[u8], contig: u32) -> Option<Vec<Transcript>> {
    let mut c = Cursor(raw);
    let mut out = Vec::new();
    while !c.0.is_empty() {
        let (start, end) = (c.u64()?, c.u64()?);
        let strand = match c.u8()? {
            0 => Strand::Forward,
            1 => Strand::Reverse,
            _ => Strand::Unknown,
        };
        let (tags, representative) = (c.u8()?, c.u8()? != 0);
        let cds = if c.u8()? != 0 { Some((c.u64()?, c.u64()?)) } else { None };
        let (id, gene_id, gene_name, biotype) = (c.str()?, c.str()?, c.str()?, c.str()?);
        let n = c.u32()? as usize;
        let mut exons = Vec::with_capacity(n);
        for _ in 0..n {
            exons.push((start + u64::from(c.u32()?), start + u64::from(c.u32()?)));
        }
        let (phase, table, exons_from_span) = (c.u8()?, c.u8()?, c.u8()? != 0);
        let n = c.u16()? as usize;
        let mut selenocysteine = Vec::with_capacity(n);
        for _ in 0..n {
            selenocysteine.push((c.u64()?, c.u64()?));
        }
        let translation = Translation { phase, table, selenocysteine };
        out.push(Transcript { contig, start, end, strand, id, gene_id, gene_name, biotype, exons, cds, tags, representative, translation, exons_from_span });
    }
    Some(out)
}

/// Opens `cache` if it is a current cache of `source`.
pub(crate) fn open(source: &Path, cache: &Path) -> Option<Annotation> {
    let key = SourceKey::of(source).ok()?;
    let mut file = File::open(cache).ok()?;
    let mut pre = [0u8; 20];
    file.read_exact(&mut pre).ok()?;
    let mut c = Cursor(&pre);
    if c.take::<8>()? != *MAGIC || c.u32()? != VERSION {
        return None;
    }
    let header_at = c.u64()?;
    file.seek(SeekFrom::Start(header_at)).ok()?;
    let mut header = Vec::new();
    file.read_to_end(&mut header).ok()?;
    let mut c = Cursor(&header);
    let found = SourceKey { path: c.str()?.into(), size: c.u64()?, mtime_ns: c.u128()?, parser: c.u64()? };
    if found != key {
        return None;
    }
    let n = c.u32()? as usize;
    let mut contigs = Vec::with_capacity(n);
    let mut parts = Vec::with_capacity(n);
    for _ in 0..n {
        contigs.push(c.str()?.into_string());
        let (count, length, offset, len) = (c.u32()? as usize, c.u64()?, c.u64()?, c.u64()?);
        parts.push(Part { count, length, block: Some((offset, len)), loaded: OnceLock::new() });
    }
    let n = c.u32()? as usize;
    let mut genes = Vec::with_capacity(n);
    for _ in 0..n {
        let name = c.str()?.into_string();
        let contig = contigs.get(c.u32()? as usize)?.clone();
        genes.push(GeneSpan { name, contig, start: c.u64()?, end: c.u64()? });
    }
    let n = c.u32()? as usize;
    let mut gene_keys = Vec::with_capacity(n);
    for _ in 0..n {
        gene_keys.push((c.str()?, c.u32()?));
    }
    let (count, table) = (c.u32()?, c.u64()?);
    Some(Annotation {
        contigs,
        parts,
        genes,
        gene_keys,
        tx_keys: Keys::Cached { count, table },
        reader: Some(Reader { file: Mutex::new(file), path: cache.to_owned() }),
    })
}
