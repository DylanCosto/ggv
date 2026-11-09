//! Gene annotations: GTF, GFF3 and BED (plain or gzipped) parsed once into transcripts, with an
//! interval index per contig and a name index for gene search.
//!
//! Queries are a binary search plus a short backwards scan over a running maximum of ends, so
//! finding what overlaps a view costs microseconds even for a full GENCODE annotation. Large
//! files are cached after their first parse (see `cache`), so later opens take milliseconds
//! and only the contigs in use are held in memory.

mod cache;
pub mod coding;
mod mane;
mod parse;

use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, Read},
    path::Path,
    sync::OnceLock,
};

use ggv_model::{Contig, contig_aliases};

pub use mane::without_copy_suffix;
pub use parse::Format;

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("{path}: {cause}")]
    Io { path: String, cause: std::io::Error },
    #[error("{path}: no transcripts or features found (is it GTF, GFF3 or BED?)")]
    Empty { path: String },
}

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
}

/// Tags that pick the transcript to show for a gene when only one is shown.
pub const TAG_MANE_SELECT: u8 = 1;
pub const TAG_CANONICAL: u8 = 2;
pub const TAG_BASIC: u8 = 4;
pub const TAG_MANE_PLUS_CLINICAL: u8 = 8;

#[derive(Clone, Debug, PartialEq)]
pub struct Transcript {
    pub contig: u32,
    /// 0-based, half-open.
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub id: Box<str>,
    pub gene_id: Box<str>,
    pub gene_name: Box<str>,
    pub biotype: Box<str>,
    /// Sorted, non-overlapping, 0-based half-open.
    pub exons: Vec<(u64, u64)>,
    /// Coding span (first CDS start to last CDS end), if coding.
    pub cds: Option<(u64, u64)>,
    pub tags: u8,
    /// The one transcript shown for its gene in the collapsed view.
    pub representative: bool,
    /// How the coding sequence is read: frame, genetic code, selenocysteines.
    pub translation: Translation,
    /// The file gave no exon or CDS rows: the exon is the transcript's span, not known
    /// structure.
    pub exons_from_span: bool,
}

/// What reading a coding sequence needs beyond its exons.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Translation {
    /// Bases before the first complete codon (the phase of the 5'-most CDS row): non-zero
    /// when the coding sequence is incomplete at its start.
    pub phase: u8,
    /// NCBI genetic code: 1 (or 0) standard, 2 vertebrate mitochondrial.
    pub table: u8,
    /// Codons read as selenocysteine rather than stop (genomic, 0-based half-open).
    pub selenocysteine: Vec<(u64, u64)>,
}

/// Whether `contig` names the mitochondrial genome (`chrM`, `MT`, `M`, `chrMT`).
pub fn is_mitochondrial(contig: &str) -> bool {
    let c = contig.strip_prefix("chr").unwrap_or(contig);
    c.eq_ignore_ascii_case("M") || c.eq_ignore_ascii_case("MT")
}

impl Transcript {
    pub fn cds_len(&self) -> u64 {
        let Some((cs, ce)) = self.cds else { return 0 };
        self.exons.iter().map(|(s, e)| e.min(&ce).saturating_sub(*s.max(&cs))).sum()
    }
}

/// An annotated intron (0-based, half-open).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Intron {
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub gene: Box<str>,
}

/// A gene's extent, for search.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GeneSpan {
    pub name: String,
    pub contig: String,
    pub start: u64,
    pub end: u64,
}

pub struct Annotation {
    contigs: Vec<String>,
    /// Per contig: its transcripts, decoded from the cache the first time they're needed.
    parts: Vec<Part>,
    genes: Vec<GeneSpan>,
    /// Lowercased gene names and gene ids (with and without version) → gene, sorted.
    gene_keys: Vec<(Box<str>, u32)>,
    /// Lowercased transcript ids (with and without version) → gene.
    tx_keys: Keys,
    /// The cache file the parts and transcript ids are read from, if cached.
    reader: Option<cache::Reader>,
}

struct Part {
    count: usize,
    /// The end of the last feature: the contig's length as far as genes go.
    length: u64,
    /// Where its compressed transcripts are in the cache file.
    block: Option<(u64, u64)>,
    loaded: OnceLock<Loaded>,
}

struct Loaded {
    /// Sorted by start.
    transcripts: Vec<Transcript>,
    /// Running maximum of `end` (the interval index).
    max_end: Vec<u64>,
}

impl Loaded {
    fn new(transcripts: Vec<Transcript>) -> Self {
        let max_end = transcripts.iter().scan(0, |m, t| {
            *m = t.end.max(*m);
            Some(*m)
        });
        Self { max_end: max_end.collect(), transcripts }
    }
}

enum Keys {
    Memory(Vec<(Box<str>, u32)>),
    /// `count` entries of a sorted table at `table` in the cache file.
    Cached { count: u32, table: u64 },
}

impl Annotation {
    /// Reads a file, picking the format from its extension (`.gtf`, `.gff`, `.gff3`, `.bed`,
    /// optionally followed by `.gz`/`.bgz`).
    pub fn read(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let err = |cause| Error::Io { path: path.display().to_string(), cause };
        let file = File::open(path).map_err(err)?;
        let name = path.file_name().unwrap_or_default().to_string_lossy().to_ascii_lowercase();
        let gz = name.ends_with(".gz") || name.ends_with(".bgz");
        let format = Format::from_name(name.trim_end_matches(".gz").trim_end_matches(".bgz"));
        // Decompress on a second thread so it overlaps with parsing.
        let reader: Box<dyn Read> =
            if gz { Box::new(Pipe::spawn(flate2::read::MultiGzDecoder::new(file))) } else { Box::new(file) };
        let ann = Self::parse(BufReader::with_capacity(1 << 20, reader), format).map_err(err)?;
        if ann.is_empty() {
            return Err(Error::Empty { path: path.display().to_string() });
        }
        Ok(ann)
    }

    /// Like [`read`](Self::read), but for a large file opens its cache when it has one, and
    /// otherwise parses it and writes one (falling back to the parsed copy if that fails).
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let Some(cache) = cache::path_for(path) else { return Self::read(path) };
        if let Some(ann) = cache::open(path, &cache) {
            return Ok(ann);
        }
        let key = cache::SourceKey::of(path);
        let ann = Self::read(path)?;
        match key.and_then(|key| cache::write(&ann, &key, &cache)) {
            // Reopen from the cache, which holds only what's in use.
            Ok(()) => Ok(cache::open(path, &cache).unwrap_or(ann)),
            Err(e) => {
                eprintln!("ggv: couldn't cache {} at {}: {e}", path.display(), cache.display());
                Ok(ann)
            }
        }
    }

    /// Whether this annotation reads from a cache.
    pub fn is_cached(&self) -> bool {
        self.reader.is_some()
    }

    pub fn parse(reader: impl BufRead, format: Format) -> std::io::Result<Self> {
        let parsed = parse::parse(reader, format)?;
        Ok(Self::index(parsed.contigs, parsed.transcripts))
    }

    fn index(contigs: Vec<String>, mut transcripts: Vec<Transcript>) -> Self {
        // MANE from the file's tags or, for files without them (UCSC's ncbiRefSeq), from the
        // shipped list: tagged only when the version matches too.
        let hits: Vec<Option<mane::Hit>> = transcripts.iter().map(|t| mane::lookup(&t.id)).collect();
        for (t, hit) in transcripts.iter_mut().zip(&hits) {
            match hit {
                Some(mane::Hit { kind: mane::Kind::Select, exact: true }) => t.tags |= TAG_MANE_SELECT,
                Some(mane::Hit { kind: mane::Kind::PlusClinical, exact: true }) => t.tags |= TAG_MANE_PLUS_CLINICAL,
                _ => {}
            }
        }
        // Pick one transcript per gene: MANE Select (the listed version over another), then
        // MANE Plus Clinical, then Ensembl canonical, then a curated RefSeq transcript over a
        // predicted model (XM_/XR_), coding first, then the longest coding sequence (for genes
        // MANE doesn't list; the oldest accession, tried before, picks short isoforms such as
        // TTN's N2B), then the longest span, then the oldest accession.
        let level = |t: &Transcript, hit: Option<mane::Hit>, tag: u8, kind: mane::Kind| -> u8 {
            match hit {
                _ if t.tags & tag != 0 => 2,
                Some(h) if h.kind == kind => 1,
                _ => 0,
            }
        };
        let rank = |i: usize| {
            let t = &transcripts[i];
            (
                level(t, hits[i], TAG_MANE_SELECT, mane::Kind::Select),
                level(t, hits[i], TAG_MANE_PLUS_CLINICAL, mane::Kind::PlusClinical),
                t.tags & TAG_CANONICAL != 0,
                !is_predicted(&t.id),
                t.cds.is_some(),
                t.cds_len(),
                t.end - t.start,
                std::cmp::Reverse(curated_number(&t.id).unwrap_or(u64::MAX)),
            )
        };
        let mut best: HashMap<(u32, Box<str>), usize> = HashMap::new();
        for (i, t) in transcripts.iter().enumerate() {
            let key = (t.contig, t.gene_id.clone());
            match best.get(&key) {
                Some(&j) if rank(j) >= rank(i) => {}
                _ => {
                    best.insert(key, i);
                }
            }
        }
        for i in best.into_values() {
            transcripts[i].representative = true;
        }

        transcripts.sort_by(|a, b| (a.contig, a.start, a.end).cmp(&(b.contig, b.start, b.end)));

        // Gene spans, for search.
        let mut gene_index: HashMap<(u32, Box<str>), usize> = HashMap::new();
        let mut genes: Vec<GeneSpan> = Vec::new();
        // Key → (gene, whether a gene name or id gave it).
        let mut names: HashMap<String, (usize, bool)> = HashMap::new();
        for t in &transcripts {
            let g = *gene_index.entry((t.contig, t.gene_id.clone())).or_insert_with(|| {
                genes.push(GeneSpan {
                    name: t.gene_name.to_string(),
                    contig: contigs[t.contig as usize].clone(),
                    start: t.start,
                    end: t.end,
                });
                genes.len() - 1
            });
            genes[g].start = genes[g].start.min(t.start);
            genes[g].end = genes[g].end.max(t.end);
            for (key, of_gene) in [(&*t.gene_name, true), (&*t.gene_id, true), (&*t.id, false)] {
                // `NM_000314.8_2` is also found as `NM_000314.8` and `NM_000314` (a UCSC copy).
                let base = if of_gene { key } else { without_copy_suffix(key) };
                for k in [key, base, strip_version(base)] {
                    let k = k.to_ascii_lowercase();
                    // Prefer primary contigs (no `_` as in `chr6_GL000250v2_alt`) for names on
                    // several contigs, and chrX for a pseudo-autosomal gene, whichever row
                    // comes first: analysis sets mask chrY's PARs, so reads and calls are on X.
                    match names.get_mut(&k) {
                        Some((existing, gene_key)) if copy_rank(&genes[*existing].contig) <= copy_rank(&genes[g].contig) => {
                            *gene_key |= of_gene;
                        }
                        Some(slot) => *slot = (g, slot.1 | of_gene),
                        None => {
                            names.insert(k, (g, of_gene));
                        }
                    }
                }
            }
        }
        let (mut gene_keys, mut tx_keys): (Vec<_>, Vec<_>) = names.into_iter().partition(|(_, (_, of_gene))| *of_gene);
        let keys = |v: &mut Vec<(String, (usize, bool))>| -> Vec<(Box<str>, u32)> {
            v.sort_unstable_by(|a, b| a.0.cmp(&b.0));
            v.drain(..).map(|(k, (g, _))| (k.into_boxed_str(), g as u32)).collect()
        };
        let (gene_keys, tx_keys) = (keys(&mut gene_keys), keys(&mut tx_keys));

        let mut parts: Vec<Vec<Transcript>> = (0..contigs.len()).map(|_| Vec::new()).collect();
        for t in transcripts {
            parts[t.contig as usize].push(t);
        }
        let parts = parts
            .into_iter()
            .map(|ts| {
                let loaded = Loaded::new(ts);
                let length = loaded.max_end.last().copied().unwrap_or(0);
                Part { count: loaded.transcripts.len(), length, block: None, loaded: OnceLock::from(loaded) }
            })
            .collect();
        Self { contigs, parts, genes, gene_keys, tx_keys: Keys::Memory(tx_keys), reader: None }
    }

    /// A contig's transcripts, decoding them from the cache on first use.
    fn part(&self, c: usize) -> &Loaded {
        let part = &self.parts[c];
        part.loaded.get_or_init(|| match (&self.reader, part.block) {
            (Some(r), Some((offset, len))) => r.load(c as u32, offset, len),
            _ => Loaded::new(Vec::new()),
        })
    }

    pub fn len(&self) -> usize {
        self.parts.iter().map(|p| p.count).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn gene_count(&self) -> usize {
        self.genes.len()
    }

    pub fn contig_names(&self) -> &[String] {
        &self.contigs
    }

    /// Contigs with a length taken from the last annotated feature, for browsing genes without
    /// a BAM.
    pub fn contigs(&self) -> Vec<Contig> {
        self.contigs
            .iter()
            .zip(&self.parts)
            .filter(|(_, p)| p.count > 0)
            .map(|(name, p)| Contig { name: name.clone(), length: p.length })
            .collect()
    }

    /// This annotation's index for `contig`, trying `chr` aliases.
    pub fn contig_index(&self, contig: &str) -> Option<usize> {
        self.contigs.iter().position(|c| c == contig).or_else(|| {
            contig_aliases(contig).into_iter().find_map(|a| self.contigs.iter().position(|c| *c == a))
        })
    }

    /// Transcripts overlapping `[start, end)` on `contig` (any naming), in start order.
    pub fn query(&self, contig: &str, start: u64, end: u64) -> impl Iterator<Item = &Transcript> {
        let (ts, lo, hi) = match self.contig_index(contig).filter(|&c| self.parts[c].count > 0) {
            Some(c) => {
                let part = self.part(c);
                let ts = &part.transcripts[..];
                let hi = ts.partition_point(|t| t.start < end);
                // Everything before `lo` ends at or before `start` (the running max says so);
                // between `lo` and `hi`, short transcripts may still end early, so filter.
                let lo = part.max_end[..hi].partition_point(|&m| m <= start);
                (ts, lo, hi)
            }
            None => (&[][..], 0, 0),
        };
        ts[lo..hi].iter().filter(move |t| t.end > start)
    }

    /// The annotated introns of transcripts overlapping `[start, end)` on `contig`, each once per
    /// gene and strand (genes sharing an intron each keep theirs): 0-based half-open intron, the
    /// gene, and its strand.
    pub fn introns(&self, contig: &str, start: u64, end: u64) -> Vec<Intron> {
        let mut out: Vec<Intron> = Vec::new();
        for t in self.query(contig, start, end) {
            for w in t.exons.windows(2) {
                let (s, e) = (w[0].1, w[1].0);
                if s < e && !out.iter().any(|i| i.start == s && i.end == e && i.strand == t.strand && i.gene == t.gene_name) {
                    out.push(Intron { start: s, end: e, strand: t.strand, gene: t.gene_name.clone() });
                }
            }
        }
        out.sort_by_key(|i| (i.start, i.end));
        out
    }

    /// Every gene, in file order (for name search).
    pub fn genes(&self) -> &[GeneSpan] {
        &self.genes
    }

    /// Finds a gene by name, gene id or transcript id (case-insensitive, version optional).
    pub fn find(&self, name: &str) -> Option<&GeneSpan> {
        let k = name.trim().to_ascii_lowercase();
        let g = self.key(&k).or_else(|| self.key(strip_version(&k)))?;
        self.genes.get(g as usize)
    }

    /// Every copy of the gene `name` finds: the same gene on other primary contigs too (a
    /// pseudo-autosomal gene such as SHOX is on chrX and chrY), the one found first.
    pub fn find_all(&self, name: &str) -> Vec<&GeneSpan> {
        let Some(g) = self.find(name) else { return Vec::new() };
        let mut out = vec![g];
        out.extend(self.genes.iter().filter(|o| !std::ptr::eq(*o, g) && o.name == g.name && o.contig != g.contig && !o.contig.contains('_')));
        out
    }

    /// The transcript with this id (case-insensitive, version and UCSC copy suffix optional),
    /// when `id` is a transcript id rather than a gene name or gene id. A transcript on several
    /// sequences is found on a primary contig.
    pub fn find_transcript(&self, id: &str) -> Option<&Transcript> {
        let k = id.trim().to_ascii_lowercase();
        let gene_key = |k: &str| self.gene_keys.binary_search_by(|(x, _)| (**x).cmp(k)).is_ok();
        if gene_key(&k) || gene_key(strip_version(&k)) {
            return None;
        }
        let g = self.find(&k)?;
        let matches = |t: &&Transcript| {
            let full = t.id.to_ascii_lowercase();
            let base = without_copy_suffix(&full);
            full == k || base == k || strip_version(base) == k
        };
        self.query(&g.contig, g.start, g.end).find(matches)
    }

    fn key(&self, k: &str) -> Option<u32> {
        let sorted = |v: &[(Box<str>, u32)]| v.binary_search_by(|(x, _)| (**x).cmp(k)).ok().map(|i| v[i].1);
        sorted(&self.gene_keys).or_else(|| match (&self.tx_keys, &self.reader) {
            (Keys::Memory(v), _) => sorted(v),
            (Keys::Cached { count, table }, Some(r)) => r.find(k, *table, *count),
            (Keys::Cached { .. }, None) => None,
        })
    }
}

/// A reader fed by another thread in 4 MB blocks, so producing (e.g. decompressing) and
/// consuming run in parallel.
struct Pipe {
    rx: std::sync::mpsc::Receiver<std::io::Result<Vec<u8>>>,
    block: Vec<u8>,
    pos: usize,
}

impl Pipe {
    fn spawn(mut inner: impl Read + Send + 'static) -> Self {
        const BLOCK: usize = 4 << 20;
        let (tx, rx) = std::sync::mpsc::sync_channel(4);
        std::thread::spawn(move || {
            loop {
                let mut block = vec![0u8; BLOCK];
                let mut filled = 0;
                let result = loop {
                    match inner.read(&mut block[filled..]) {
                        Ok(0) => break Ok(()),
                        Ok(n) => {
                            filled += n;
                            if filled == BLOCK {
                                break Ok(());
                            }
                        }
                        Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
                        Err(e) => break Err(e),
                    }
                };
                block.truncate(filled);
                let done = filled < BLOCK;
                let send = match result {
                    Ok(()) => tx.send(Ok(block)),
                    Err(e) => tx.send(Err(e)),
                };
                if send.is_err() || done {
                    return;
                }
            }
        });
        Self { rx, block: Vec::new(), pos: 0 }
    }
}

impl Read for Pipe {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        while self.pos == self.block.len() {
            match self.rx.recv() {
                Ok(block) => {
                    self.block = block?;
                    self.pos = 0;
                    if self.block.is_empty() {
                        return Ok(0);
                    }
                }
                Err(_) => return Ok(0), // producer finished
            }
        }
        let n = buf.len().min(self.block.len() - self.pos);
        buf[..n].copy_from_slice(&self.block[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

/// Ensembl and RefSeq ids, whose version may be left out (`ENST…`, `NM_004333`).
fn versioned(base: &str) -> bool {
    let has = |p: &str| base.get(..p.len()).is_some_and(|b| b.eq_ignore_ascii_case(p));
    has("ENS") || ["NM_", "NR_", "XM_", "XR_"].into_iter().any(has)
}

/// A RefSeq model predicted by the genome pipeline rather than curated (`XM_`, `XR_`).
pub fn is_predicted(id: &str) -> bool {
    id.starts_with("XM_") || id.starts_with("XR_")
}

/// The number of a curated RefSeq accession (`NM_004333.6` → 4333; also `NR_`), whatever
/// version or UCSC copy suffix (`NM_000451.4_2`) follows.
fn curated_number(id: &str) -> Option<u64> {
    let rest = id.strip_prefix("NM_").or_else(|| id.strip_prefix("NR_"))?;
    let digits = &rest[..rest.bytes().position(|b| !b.is_ascii_digit()).unwrap_or(rest.len())];
    digits.parse().ok()
}

/// How a contig ranks for a name found on several (lower first): a primary contig, then chrY
/// (its PARs copy chrX's), then the rest (`_alt`, `_fix` and the like).
fn copy_rank(contig: &str) -> u8 {
    if contig.contains('_') {
        2
    } else {
        u8::from(contig.strip_prefix("chr").unwrap_or(contig).eq_ignore_ascii_case("Y"))
    }
}

/// `ENSG00000012048.23` → `ENSG00000012048`.
fn strip_version(id: &str) -> &str {
    match id.rsplit_once('.') {
        Some((base, v)) if !v.is_empty() && v.bytes().all(|b| b.is_ascii_digit()) && versioned(base) => base,
        _ => id,
    }
}

#[cfg(test)]
mod tests;
