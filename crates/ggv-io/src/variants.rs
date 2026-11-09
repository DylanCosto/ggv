//! Indexed VCF/BCF access through pasta's HTSlib-compatible variant engine: bgzipped VCF with a
//! TBI or CSI index, or BCF with a CSI index.

use std::path::{Path, PathBuf};

use ggv_model::{Region, contig_aliases};
use pasta::variant::{
    Header, Record, Unpack,
    encoding::{INT32_MISSING, INT32_VECTOR_END},
    index::{Index, IndexFormat, locate, tabix::Tabix},
    io::{InputFormat, Query, Reader},
    vcf,
};

use crate::{Error, Result};

/// A sample's genotype, simplified for display.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize)]
pub enum GenotypeClass {
    HomRef,
    Het,
    HomAlt,
    /// Any allele missing (`./.`, `0/.`).
    NoCall,
    /// The sample has no GT at this record (or the file has no GT).
    Unknown,
}

#[derive(Clone, Debug, PartialEq, serde::Serialize)]
pub struct SampleCall {
    pub class: GenotypeClass,
    /// e.g. `0/1`, `1|0`, `./.`.
    pub gt: String,
    pub dp: Option<i32>,
    pub gq: Option<i32>,
    /// Allele depths, REF first; a missing value (`.`) is `None`, kept in its place so the rest
    /// still line up with their alleles. Empty when AD is absent or wholly missing.
    pub ad: Vec<Option<i32>>,
}

impl SampleCall {
    /// AD as VCF writes it: `5,.`; `None` when there is none.
    pub fn ad_text(&self) -> Option<String> {
        (!self.ad.is_empty()).then(|| self.ad.iter().map(|v| v.map_or_else(|| ".".to_owned(), |v| v.to_string())).collect::<Vec<_>>().join(","))
    }
}

#[derive(Clone, Debug, Default, PartialEq, serde::Serialize)]
pub struct Variant {
    /// 0-based start and exclusive end (start + REF length, or INFO/END).
    pub start: u64,
    pub end: u64,
    pub id: String,
    pub reference: String,
    pub alts: Vec<String>,
    pub qual: Option<f32>,
    /// `PASS`, filter names, or empty for `.`.
    pub filters: Vec<String>,
    /// The INFO column as text (truncated past 8,000 characters).
    pub info: String,
    /// Frequency of any non-reference allele at the site: INFO/AF summed over the ALTs, else
    /// ALT alleles over called alleles in the genotypes. Per-ALT values are in `alt_af`.
    pub af: Option<f32>,
    /// Each ALT's frequency, in ALT order: INFO/AF, else from the genotypes.
    pub alt_af: Vec<Option<f32>>,
    /// INFO/SVTYPE, read from the whole record (not the shortened `info`).
    pub svtype: Option<String>,
    /// INFO/SVLEN as written, read from the whole record: one per ALT when well formed, `None`
    /// for `.`; empty when absent.
    pub svlen: Vec<Option<i64>>,
    pub samples: Vec<SampleCall>,
}

/// INFO kept per record, in characters: whole for ClinVar's longest (condition database lists),
/// bounded for files with huge annotations.
const INFO_MAX: usize = 8000;

/// An `##INFO` header line: its ID and Description.
#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize)]
pub struct InfoField {
    pub id: String,
    pub description: String,
}

impl Variant {
    /// INFO's value for `key`: `Some("")` for a flag, `None` when absent.
    pub fn info_value(&self, key: &str) -> Option<&str> {
        self.info.split(';').find_map(|kv| match kv.split_once('=') {
            Some((k, v)) if k == key => Some(v),
            None if kv == key => Some(""),
            _ => None,
        })
    }

    /// INFO as key/value pairs, in file order (flags with an empty value).
    pub fn info_pairs(&self) -> impl Iterator<Item = (&str, &str)> {
        self.info.split(';').filter(|kv| !kv.is_empty() && *kv != ".").map(|kv| kv.split_once('=').unwrap_or((kv, "")))
    }

    /// AF for display: the site's, then each ALT's when there are several
    /// (`0.910 (C 0.010, G 0.900)`); `.` when unknown.
    pub fn af_text(&self) -> String {
        let f = |a: Option<f32>| a.map_or_else(|| ".".to_owned(), |a| format!("{a:.3}"));
        let mut s = f(self.af);
        if self.alts.len() > 1 && self.alt_af.iter().any(Option::is_some) {
            let each: Vec<String> = self.alts.iter().zip(&self.alt_af).map(|(alt, a)| format!("{alt} {}", f(*a))).collect();
            s += &format!(" ({})", each.join(", "));
        }
        s
    }

    pub fn is_pass(&self) -> bool {
        self.filters.is_empty() || self.filters.iter().all(|f| f == "PASS")
    }

    /// A structural variant's type and length, when this is one: the first ALT that is an
    /// insertion or deletion of at least `min_len` bases (from the alleles), or a symbolic
    /// `<DEL>`, `<INS>`, `<DUP>`, `<INV>`, `<CNV>` (length from that ALT's SVLEN, else INFO/END).
    /// Breakends are `BND` with no length.
    pub fn structural(&self, min_len: u64) -> Option<(&'static str, Option<u64>)> {
        self.alts.iter().enumerate().filter(|(_, a)| !is_reference_block(a)).find_map(|(i, a)| self.structural_alt(i, a, min_len))
    }

    fn structural_alt(&self, i: usize, alt: &str, min_len: u64) -> Option<(&'static str, Option<u64>)> {
        if is_breakend(alt) {
            return Some(("BND", None));
        }
        if let Some(t) = alt.strip_prefix('<') {
            let t = t.trim_end_matches('>').split(':').next().unwrap_or("");
            let name = match t {
                "DEL" => "DEL",
                "INS" => "INS",
                "DUP" => "DUP",
                "INV" => "INV",
                "CNV" => "CNV",
                "BND" => return Some(("BND", None)),
                _ => return None,
            };
            let span = self.end.saturating_sub(self.start + 1);
            return Some((name, self.svlen_of(i).or((name != "INS").then_some(span))));
        }
        let (r, a) = (self.reference.len() as u64, alt.len() as u64);
        match r.abs_diff(a) {
            d if d >= min_len && a > r => Some(("INS", Some(d))),
            d if d >= min_len => Some(("DEL", Some(d))),
            _ => None,
        }
    }

    /// ALT `i`'s SVLEN as a length; `None` when absent, `.`, or SVLEN doesn't give one value
    /// per ALT (see `svlen_mismatch`).
    pub fn svlen_of(&self, i: usize) -> Option<u64> {
        if self.svlen.len() != self.alts.len() {
            return None;
        }
        self.svlen.get(i).copied().flatten().map(i64::unsigned_abs)
    }

    /// SVLEN is present but doesn't give one value per ALT.
    pub fn svlen_mismatch(&self) -> bool {
        !self.svlen.is_empty() && self.svlen.len() != self.alts.len()
    }

    /// One ALT's type: SNV, MNV, INS, DEL, BND (breakend), REF (gVCF block), or the symbolic
    /// ALT's name.
    pub fn alt_kind(&self, alt: &str) -> String {
        if is_reference_block(alt) {
            return "REF".into();
        }
        if is_breakend(alt) {
            return "BND".into();
        }
        if alt.starts_with('<') {
            return alt.trim_matches(|c| c == '<' || c == '>').to_owned();
        }
        match (self.reference.len(), alt.len()) {
            (1, 1) => "SNV".into(),
            (a, b) if a == b => "MNV".into(),
            (a, b) if a < b => "INS".into(),
            _ => "DEL".into(),
        }
    }

    /// A short type label over all ALTs: their shared `alt_kind`, MIXED when they differ, REF
    /// for a reference block. A spanning deletion (`*`) is left out.
    pub fn kind(&self) -> String {
        let mut kinds = self.alts.iter().filter(|a| !is_reference_block(a) && a.as_str() != "*").map(|a| self.alt_kind(a));
        let Some(first) = kinds.next() else { return "REF".into() };
        if kinds.all(|k| k == first) { first } else { "MIXED".into() }
    }

    /// The bases the variant changes, 0-based half-open, for drawing it over reads. VCF indels
    /// and symbolic SVs start with an unchanged padding base (`GC`→`G` at POS deletes the base
    /// after POS), so it is skipped as IGV does when every allele shares it. An insertion
    /// becomes an empty span at the boundary where the bases go.
    pub fn shown(&self) -> (u64, u64) {
        let first = self.reference.as_bytes().first().map(u8::to_ascii_uppercase);
        let alts: Vec<&String> = self.alts.iter().filter(|a| !is_reference_block(a) && a.as_str() != "*").collect();
        let padded = first.is_some()
            && !alts.is_empty()
            && alts.iter().any(|a| a.starts_with('<') || a.len() != self.reference.len())
            && alts.iter().all(|a| a.starts_with('<') || a.as_bytes().first().map(u8::to_ascii_uppercase) == first);
        if padded {
            let start = self.start + 1;
            (start, self.end.max(start))
        } else {
            (self.start, self.end)
        }
    }
}

/// A breakend ALT: `G]chr2:100]`, `[chr2:100[G`, or a single breakend `G.` / `.G`.
pub fn is_breakend(alt: &str) -> bool {
    alt.contains('[') || alt.contains(']') || (alt.len() > 1 && (alt.starts_with('.') || alt.ends_with('.')))
}

/// A breakend ALT's mate: its contig and 1-based position. `None` for a single breakend.
pub fn breakend_mate(alt: &str) -> Option<(&str, u64)> {
    let open = alt.find(['[', ']'])?;
    let bracket = alt.as_bytes()[open] as char;
    let inner = &alt[open + 1..];
    let inner = &inner[..inner.find(bracket)?];
    let (contig, pos) = inner.rsplit_once(':')?;
    Some((contig, pos.parse().ok()?))
}

/// How to index a VCF: only a bgzipped one can be (a plain-text one is compressed first).
fn vcf_index_advice(path: &Path) -> String {
    let mut magic = [0u8; 2];
    let gzipped = std::fs::File::open(path).and_then(|mut f| std::io::Read::read_exact(&mut f, &mut magic)).is_ok() && magic == [0x1f, 0x8b];
    if gzipped {
        "create one with `tabix -p vcf` (or `bcftools index -t`); it must be bgzip-compressed, not plain gzip".into()
    } else {
        "a plain-text VCF can't be indexed: compress it with `bgzip` (or `bcftools view -Oz -o x.vcf.gz`), then index that with `tabix -p vcf`".into()
    }
}

/// INFO's value for `key` in INFO text.
fn info_text_value<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    info.split(';').find_map(|kv| kv.split_once('=').filter(|(k, _)| *k == key).map(|(_, v)| v))
}

/// gVCF reference-block ALTs.
fn is_reference_block(alt: &str) -> bool {
    alt == "<*>" || alt == "<NON_REF>" || alt == "."
}

enum IndexKind {
    Tabix(Box<Tabix>),
    Csi(Box<Index>),
}

pub struct VariantFile {
    path: PathBuf,
    reader: Reader,
    header: Header,
    index: IndexKind,
    samples: Vec<String>,
    /// Contig names as the index knows them, for alias matching.
    contigs: Vec<String>,
    line: Vec<u8>,
    record: Record,
    gts: Vec<i32>,
    ints: Vec<i32>,
    text: Vec<u8>,
}

impl VariantFile {
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref().to_owned();
        let p = path.to_string_lossy().into_owned();
        let io_err = |e: std::io::Error| Error::Io { path: path.clone(), cause: e };
        let mut reader = Reader::open(&p, 0).map_err(|e| io_err(std::io::Error::other(e.to_string())))?;
        let header = reader.read_header().map_err(|e| io_err(std::io::Error::other(e.to_string())))?;

        let index = match reader.format() {
            InputFormat::Vcf => {
                let idx = locate(&p, IndexFormat::Tbi).or_else(|| locate(&p, IndexFormat::Csi)).ok_or_else(|| Error::MissingIndex {
                    path: path.clone(),
                    tried: ".tbi, .csi".into(),
                    how: vcf_index_advice(&path),
                })?;
                IndexKind::Tabix(Box::new(Tabix::load(&idx).map_err(io_err)?))
            }
            InputFormat::Bcf => {
                let idx = locate(&p, IndexFormat::Csi)
                    .ok_or_else(|| Error::MissingIndex { path: path.clone(), tried: ".csi".into(), how: "create one with `bcftools index`".into() })?;
                IndexKind::Csi(Box::new(Index::read(&idx).map_err(io_err)?))
            }
        };

        let samples = header.samples().iter().map(|s| String::from_utf8_lossy(s).into_owned()).collect();
        let contigs = match &index {
            IndexKind::Tabix(t) => t.names().iter().map(|n| String::from_utf8_lossy(n).into_owned()).collect(),
            IndexKind::Csi(_) => header.contig_names().iter().map(|n| String::from_utf8_lossy(n).into_owned()).collect(),
        };

        Ok(Self {
            path,
            reader,
            header,
            index,
            samples,
            contigs,
            line: Vec::new(),
            record: Record::new(),
            gts: Vec::new(),
            ints: Vec::new(),
            text: Vec::new(),
        })
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    pub fn samples(&self) -> &[String] {
        &self.samples
    }

    /// The `##INFO` fields the header declares, in header order.
    pub fn info_fields(&self) -> Vec<InfoField> {
        use pasta::variant::header::HeaderRecordType;
        let text = |b: &[u8]| String::from_utf8_lossy(b).trim_matches('"').to_owned();
        self.header
            .records()
            .filter(|r| r.kind() == HeaderRecordType::Info)
            .filter_map(|r| Some(InfoField { id: text(r.get(b"ID")?), description: r.get(b"Description").map(text).unwrap_or_default() }))
            .collect()
    }

    /// The `##contig` lines that give a length, in header order.
    pub fn header_contigs(&self) -> Vec<ggv_model::Contig> {
        use pasta::variant::header::HeaderRecordType;
        self.header
            .records()
            .filter(|r| r.kind() == HeaderRecordType::Contig)
            .filter_map(|r| {
                let name = String::from_utf8_lossy(r.get(b"ID")?).into_owned();
                let length = std::str::from_utf8(r.get(b"length")?).ok()?.trim().parse().ok()?;
                Some(ggv_model::Contig { name, length })
            })
            .collect()
    }

    /// What the header says the assembly is: `##reference`, `##assembly`, or a contig's
    /// `assembly=` (as written, e.g. `file:///ref/GRCh38.fa` or `b37`).
    pub fn declared_assembly(&self) -> Option<String> {
        use pasta::variant::header::HeaderRecordType;
        let text = |b: &[u8]| String::from_utf8_lossy(b).trim_matches('"').to_owned();
        let mut contig_assembly = None;
        for r in self.header.records() {
            match r.kind() {
                HeaderRecordType::Generic if matches!(r.key(), b"reference" | b"assembly") => {
                    if let Some(v) = r.value().filter(|v| !v.is_empty()) {
                        return Some(text(v));
                    }
                }
                HeaderRecordType::Contig if contig_assembly.is_none() => contig_assembly = r.get(b"assembly").map(text),
                _ => {}
            }
        }
        contig_assembly
    }

    /// The contigs the index (or a BCF's header) names.
    pub fn contig_names(&self) -> &[String] {
        &self.contigs
    }

    /// This file's name for `contig`, trying `chr` aliases.
    pub fn contig_name(&self, contig: &str) -> Option<&str> {
        self.contigs.iter().find(|c| *c == contig).map(String::as_str).or_else(|| {
            contig_aliases(contig).into_iter().find_map(|a| self.contigs.iter().find(|c| **c == a).map(String::as_str))
        })
    }

    /// Records overlapping `region` (any contig naming), in file order. gVCF reference blocks
    /// (ALT only `<*>`/`<NON_REF>`) are skipped unless `reference_blocks` is set.
    pub fn query(&mut self, region: &Region, reference_blocks: bool) -> Result<Vec<Variant>> {
        let Some(contig) = self.contig_name(&region.contig).map(str::to_owned) else { return Ok(Vec::new()) };
        let (beg, end) = (region.start as i64, region.end as i64);
        let err = |e: std::io::Error| Error::Io { path: self.path.clone(), cause: e };
        let mut out = Vec::new();

        match &mut self.index {
            IndexKind::Tabix(tbx) => {
                let Some(tid) = tbx.name_to_id(contig.as_bytes()) else { return Ok(out) };
                let Some(plan) = tbx.index.query(tid as usize, beg, end) else { return Ok(out) };
                let mut q = Query::new(plan);
                while self.reader.query_next_line(&mut q, tbx, &mut self.line).map_err(err)? {
                    let mut line = std::mem::take(&mut self.line);
                    let parsed = vcf::parse(&mut self.header, &mut line, &mut self.record);
                    self.line = line;
                    if parsed.is_ok() {
                        push(&self.header, &mut self.record, &mut self.gts, &mut self.ints, &mut self.text, reference_blocks, &mut out);
                    }
                }
            }
            IndexKind::Csi(index) => {
                let Some(tid) = self.header.contig_id(contig.as_bytes()) else { return Ok(out) };
                let Some(plan) = index.query(tid as usize, beg, end) else { return Ok(out) };
                let mut q = Query::new(plan);
                while self.reader.query_next_bcf(&mut q, &self.header, &mut self.record).map_err(err)? {
                    push(&self.header, &mut self.record, &mut self.gts, &mut self.ints, &mut self.text, reference_blocks, &mut out);
                }
            }
        }
        Ok(out)
    }
}

#[allow(clippy::too_many_arguments)]
fn push(
    header: &Header,
    rec: &mut Record,
    gts: &mut Vec<i32>,
    ints: &mut Vec<i32>,
    text: &mut Vec<u8>,
    reference_blocks: bool,
    out: &mut Vec<Variant>,
) {
    // Records decode lazily (as in HTSlib): unpack the site fields before reading them.
    rec.unpack(Unpack::SHR);
    let alleles: Vec<String> = rec.alleles().map(|a| String::from_utf8_lossy(a).into_owned()).collect();
    let (reference, alts) = match alleles.split_first() {
        Some((r, a)) => (r.clone(), a.to_vec()),
        None => (String::new(), Vec::new()),
    };
    if !reference_blocks && alts.iter().all(|a| is_reference_block(a)) {
        return;
    }

    let filters = rec
        .filters()
        .iter()
        .map(|&id| header.id_name(id as u32).map_or_else(|| id.to_string(), |n| String::from_utf8_lossy(n).into_owned()))
        .collect();

    // INFO as text: format the record and take column 8. Fields that give the variant's
    // meaning are read from it whole, before it is shortened for display.
    text.clear();
    let (mut svtype, mut svlen) = (None, Vec::new());
    let info = if vcf::format_record(header, rec, text).is_ok() {
        let s = String::from_utf8_lossy(text);
        let full = s.split('\t').nth(7).unwrap_or(".").trim_end();
        svtype = info_text_value(full, "SVTYPE").map(str::to_owned);
        if let Some(v) = info_text_value(full, "SVLEN") {
            svlen = v.split(',').map(|x| x.trim().parse::<i64>().ok()).collect();
        }
        let mut info: String = full.chars().take(INFO_MAX).collect();
        if info.chars().count() == INFO_MAX {
            info.push('…');
        }
        info
    } else {
        String::new()
    };

    let n = header.sample_count();
    let mut samples = vec![
        SampleCall { class: GenotypeClass::Unknown, gt: ".".into(), dp: None, gq: None, ad: Vec::new() };
        n
    ];
    if n > 0 && let Ok(per) = rec.get_genotypes(header, gts).map(|total| total / n.max(1)) {
        for (s, call) in samples.iter_mut().enumerate() {
            (call.class, call.gt) = classify(&gts[s * per..(s + 1) * per]);
        }
    }
    for (tag, set) in [(&b"DP"[..], 0), (b"GQ", 1)] {
        if n > 0 && let Ok(total) = rec.get_format_i32(header, tag, ints) {
            let per = total / n;
            for (s, call) in samples.iter_mut().enumerate() {
                let v = ints[s * per];
                let v = (v != INT32_MISSING && v != INT32_VECTOR_END).then_some(v);
                if set == 0 { call.dp = v } else { call.gq = v }
            }
        }
    }
    if n > 0 && let Ok(total) = rec.get_format_i32(header, b"AD", ints) {
        let per = total / n;
        for (s, call) in samples.iter_mut().enumerate() {
            // htslib's missing value is i32::MIN: never a depth.
            let ad: Vec<Option<i32>> = ints[s * per..(s + 1) * per].iter().copied().take_while(|&v| v != INT32_VECTOR_END).map(|v| (v != INT32_MISSING).then_some(v)).collect();
            call.ad = if ad.iter().all(Option::is_none) { Vec::new() } else { ad };
        }
    }

    let mut info_af = Vec::new();
    let info_af: Vec<Option<f32>> = match rec.get_info_f32(header, b"AF", &mut info_af) {
        Ok(k) if k > 0 => info_af.iter().map(|&x| (f32::to_bits(x) != pasta::variant::encoding::FLOAT_MISSING).then_some(x)).collect(),
        _ => Vec::new(),
    };
    let (af, alt_af) = allele_frequencies(&info_af, alts.len(), &samples);

    let start = rec.pos().max(0) as u64;
    out.push(Variant {
        start,
        end: (rec.pos() + rec.rlen().max(1)).max(0) as u64,
        id: String::from_utf8_lossy(rec.id()).into_owned(),
        reference,
        alts,
        qual: rec.qual(),
        filters,
        info,
        af,
        alt_af,
        svtype,
        svlen,
        samples,
    });
}

/// The site's non-reference frequency and each ALT's, from INFO/AF values (one per ALT; a
/// single value for several ALTs is taken as the site's) or else from the genotypes.
fn allele_frequencies(info_af: &[Option<f32>], n_alts: usize, samples: &[SampleCall]) -> (Option<f32>, Vec<Option<f32>>) {
    if info_af.len() == n_alts && info_af.iter().any(Option::is_some) {
        let sum: f32 = info_af.iter().flatten().sum();
        return (Some(sum.min(1.0)), info_af.to_vec());
    }
    if info_af.len() == 1 && info_af[0].is_some() {
        return (info_af[0], vec![None; n_alts]);
    }
    // From genotypes: each ALT's alleles over called alleles.
    let mut counts = vec![0u32; n_alts];
    let mut called = 0u32;
    for c in samples {
        for a in c.gt.split(['/', '|']) {
            if let Ok(i) = a.parse::<usize>() {
                called += 1;
                if let Some(n) = i.checked_sub(1).and_then(|k| counts.get_mut(k)) {
                    *n += 1;
                }
            }
        }
    }
    if called == 0 {
        return (None, vec![None; n_alts]);
    }
    let f = |n: u32| n as f32 / called as f32;
    (Some(f(counts.iter().sum())), counts.into_iter().map(|n| Some(f(n))).collect())
}

/// HTSlib genotype encoding: each allele is `(index + 1) << 1 | phased`; 0 means missing;
/// `INT32_VECTOR_END` pads shorter ploidies.
fn classify(gt: &[i32]) -> (GenotypeClass, String) {
    let mut text = String::new();
    let mut alleles = Vec::new();
    for (k, &v) in gt.iter().enumerate() {
        if v == INT32_VECTOR_END {
            break;
        }
        if k > 0 {
            text.push(if v & 1 == 1 { '|' } else { '/' });
        }
        if v == INT32_MISSING || v >> 1 == 0 {
            text.push('.');
            alleles.push(None);
        } else {
            let a = (v >> 1) - 1;
            text.push_str(&a.to_string());
            alleles.push(Some(a));
        }
    }
    let class = if alleles.is_empty() {
        GenotypeClass::Unknown
    } else if alleles.iter().any(Option::is_none) {
        GenotypeClass::NoCall
    } else if alleles.iter().all(|a| *a == Some(0)) {
        GenotypeClass::HomRef
    } else if alleles.windows(2).all(|w| w[0] == w[1]) {
        GenotypeClass::HomAlt
    } else {
        GenotypeClass::Het
    };
    (class, if text.is_empty() { ".".into() } else { text })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn enc(a: i32, phased: bool) -> i32 {
        ((a + 1) << 1) | i32::from(phased)
    }

    #[test]
    fn genotype_classes() {
        assert_eq!(classify(&[enc(0, false), enc(1, false)]), (GenotypeClass::Het, "0/1".into()));
        assert_eq!(classify(&[enc(1, false), enc(0, true)]), (GenotypeClass::Het, "1|0".into()));
        assert_eq!(classify(&[enc(1, false), enc(1, false)]), (GenotypeClass::HomAlt, "1/1".into()));
        assert_eq!(classify(&[enc(0, false), enc(0, false)]), (GenotypeClass::HomRef, "0/0".into()));
        assert_eq!(classify(&[0, 0]), (GenotypeClass::NoCall, "./.".into()));
        assert_eq!(classify(&[enc(1, false), enc(2, false)]), (GenotypeClass::Het, "1/2".into()));
        // Haploid, padded with a vector end.
        assert_eq!(classify(&[enc(1, false), INT32_VECTOR_END]), (GenotypeClass::HomAlt, "1".into()));
        assert_eq!(classify(&[]), (GenotypeClass::Unknown, ".".into()));
    }

    #[test]
    fn kinds() {
        let v = |r: &str, a: &str| Variant {
            start: 0,
            end: 1,
            id: ".".into(),
            reference: r.into(),
            alts: vec![a.into()],
            qual: None,
            filters: vec![],
            info: String::new(),
            af: None,
            samples: vec![],
            ..Default::default()
        };
        assert_eq!(v("A", "G").kind(), "SNV");
        assert_eq!(v("AC", "A").kind(), "DEL");
        assert_eq!(v("A", "ACT").kind(), "INS");
        assert_eq!(v("AC", "GT").kind(), "MNV");
        assert_eq!(v("A", "<DEL>").kind(), "DEL");
        assert_eq!(v("A", "<*>").kind(), "REF");
        // Breakends are BND, not insertions by text length.
        assert_eq!(v("A", "A]chr1:700]").kind(), "BND");
        assert_eq!(v("A", "[chr1:700[A").kind(), "BND");
        assert_eq!(v("A", "A.").kind(), "BND");
        assert_eq!(v("A", "A]chr1:700]").structural(50), Some(("BND", None)));
    }

    #[test]
    fn kind_covers_every_alt() {
        let v = |alts: &[&str]| Variant { reference: "A".into(), alts: alts.iter().map(|a| (*a).to_owned()).collect(), ..Default::default() };
        assert_eq!(v(&["C", "<DEL>"]).kind(), "MIXED");
        assert_eq!(v(&["<DEL>", "C"]).kind(), "MIXED");
        assert_eq!(v(&["C", "G"]).kind(), "SNV");
        assert_eq!(v(&["C", "*"]).kind(), "SNV");
        assert_eq!(v(&["C", "<DEL>"]).alt_kind("<DEL>"), "DEL");
    }

    #[test]
    fn svlen_is_per_alt() {
        let v = Variant {
            start: 119,
            end: 170,
            reference: "A".into(),
            alts: vec!["<DEL>".into(), "<INS>".into()],
            svlen: vec![Some(-50), Some(100)],
            ..Default::default()
        };
        assert_eq!(v.svlen_of(0), Some(50));
        assert_eq!(v.svlen_of(1), Some(100));
        assert!(!v.svlen_mismatch());
        // One value for two ALTs: no length is made up for either.
        let short = Variant { svlen: vec![Some(-50)], ..v.clone() };
        assert!(short.svlen_mismatch());
        assert_eq!(short.svlen_of(1), None);
        // An SNV first doesn't hide the SV after it.
        let mixed = Variant { alts: vec!["C".into(), "<DEL>".into()], svlen: vec![None, Some(-50)], ..v };
        assert_eq!(mixed.structural(50), Some(("DEL", Some(50))));
        assert_eq!(breakend_mate("A]chrQA:700]"), Some(("chrQA", 700)));
        assert_eq!(breakend_mate("[chr2:5[T"), Some(("chr2", 5)));
        assert_eq!(breakend_mate("A."), None);
    }

    #[test]
    fn allele_frequency_per_alt() {
        let call = |gt: &str| SampleCall { class: GenotypeClass::Het, gt: gt.into(), dp: None, gq: None, ad: vec![] };
        // INFO/AF=0.01,0.9: both kept, the site's is their sum.
        let (af, alt) = allele_frequencies(&[Some(0.01), Some(0.9)], 2, &[]);
        assert_eq!(alt, vec![Some(0.01), Some(0.9)]);
        assert!((af.unwrap() - 0.91).abs() < 1e-6);
        // From genotypes the same way: 1/2 is half C, half G, all non-reference.
        let (af, alt) = allele_frequencies(&[], 2, &[call("1/2")]);
        assert_eq!((af, alt), (Some(1.0), vec![Some(0.5), Some(0.5)]));
        let (af, alt) = allele_frequencies(&[], 1, &[call("0/1"), call("0/0")]);
        assert_eq!((af, alt), (Some(0.25), vec![Some(0.25)]));
        let v = Variant { alts: vec!["C".into(), "G".into()], af: Some(0.91), alt_af: vec![Some(0.01), Some(0.9)], ..Default::default() };
        assert_eq!(v.af_text(), "0.910 (C 0.010, G 0.900)");
    }

    /// Writes `body` under a minimal header as a bgzipped, tabix-indexed VCF.
    fn write_vcf(name: &str, body: &str) -> PathBuf {
        use std::io::Write;
        let dir = std::env::temp_dir().join(format!("ggv-io-variants-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(name);
        let text = format!(
            "##fileformat=VCFv4.2\n##contig=<ID=chrQA,length=2000>\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"t\">\n##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"l\">\n##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">\n##INFO=<ID=NOTE,Number=1,Type=String,Description=\"n\">\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"f\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{body}"
        );
        let mut w = pasta::bgzf::io::Writer::new(std::fs::File::create(&path).unwrap());
        w.write_all(text.as_bytes()).unwrap();
        w.finish().unwrap();
        crate::indexing::build_index(&path).unwrap();
        path
    }

    #[test]
    fn a_missing_index_says_how_to_make_one() {
        let indexed = write_vcf("noindex.vcf.gz", "chrQA	120	.	A	C	.	PASS	.
");
        let tbi = indexed.with_added_extension("tbi");
        std::fs::remove_file(&tbi).unwrap();
        let e = VariantFile::open(&indexed).err().unwrap().to_string();
        assert!(e.contains("tabix -p vcf") && !e.contains("samtools"), "{e}");
        // A plain-text VCF must be bgzipped before it can be indexed.
        let plain = indexed.with_file_name("plain.vcf");
        std::fs::write(&plain, "##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
").unwrap();
        let e = VariantFile::open(&plain).err().unwrap().to_string();
        assert!(e.contains("bgzip") && !e.contains("samtools"), "{e}");
    }

    #[test]
    fn semantic_info_survives_long_info() {
        // An INFO longer than the display limit: SVLEN after it is still read.
        let long = "x".repeat(INFO_MAX + 100);
        let body = format!("chrQA\t120\tlate\tA\t<INS>\t.\tPASS\tSVTYPE=INS;NOTE={long};SVLEN=100;AF=0.01\nchrQA\t300\tmulti\tA\tC,G\t.\tPASS\tAF=0.01,0.9\n");
        let path = write_vcf("late.vcf.gz", &body);
        let mut f = VariantFile::open(&path).unwrap();
        let v = f.query(&Region { contig: "chrQA".into(), start: 0, end: 2000 }, false).unwrap();
        assert!(v[0].info.ends_with('…'));
        assert_eq!(v[0].svtype.as_deref(), Some("INS"));
        assert_eq!(v[0].svlen_of(0), Some(100));
        assert_eq!(v[0].structural(50), Some(("INS", Some(100))));
        assert_eq!(v[1].alt_af, vec![Some(0.01), Some(0.9)]);
    }

    #[test]
    fn reads_fixture_vcfs() {
        // GGV_TEST_VCF=file.vcf.gz GGV_TEST_REGIONS=chr20:1-200000 checks real files.
        let (Ok(path), Ok(regions)) = (std::env::var("GGV_TEST_VCF"), std::env::var("GGV_TEST_REGIONS")) else { return };
        let mut f = VariantFile::open(&path).unwrap();
        for r in regions.split(',') {
            let locus: ggv_model::Locus = r.parse().unwrap();
            let region = Region { contig: locus.contig, start: locus.start.unwrap() - 1, end: locus.end.unwrap() };
            let v = f.query(&region, false).unwrap();
            eprintln!("{r}: {} variants, samples {:?}", v.len(), f.samples());
            for x in v.iter().take(3) {
                eprintln!("  {}:{} {}>{:?} {:?} af={:?} {:?}", region.contig, x.start + 1, x.reference, x.alts, x.filters, x.af, x.samples.iter().map(|s| &s.gt).collect::<Vec<_>>());
            }
        }
    }
}

#[cfg(test)]
mod shown_tests {
    use super::Variant;

    fn v(pos1: u64, reference: &str, alts: &[&str], end: Option<u64>) -> Variant {
        let start = pos1 - 1;
        Variant {
            start,
            end: end.unwrap_or(start + reference.len() as u64),
            id: ".".into(),
            reference: reference.into(),
            alts: alts.iter().map(|a| (*a).to_owned()).collect(),
            qual: None,
            filters: vec![],
            info: String::new(),
            af: None,
            samples: vec![],
            ..Default::default()
        }
    }

    #[test]
    fn padding_base_is_skipped() {
        // chr20:58,854,245 GC>G deletes the C at 58,854,246 (0-based 58,854,245).
        assert_eq!(v(58_854_245, "GC", &["G"], None).shown(), (58_854_245, 58_854_246));
        // An insertion is the boundary after the padding base.
        assert_eq!(v(100, "A", &["ACGT"], None).shown(), (100, 100));
        // Symbolic SVs: POS is the padding base, END the last deleted base.
        assert_eq!(v(1000, "N", &["<DEL>"], Some(1500)).shown(), (1000, 1500));
        // SNVs, MNVs and indels without a shared first base stay as written.
        assert_eq!(v(10, "A", &["G"], None).shown(), (9, 10));
        assert_eq!(v(10, "AC", &["GT"], None).shown(), (9, 11));
        assert_eq!(v(10, "AC", &["T"], None).shown(), (9, 11));
        // Reference-only (RefCall ".") and gVCF blocks have nothing to skip.
        assert_eq!(v(10, "A", &["."], None).shown(), (9, 10));
    }
}
