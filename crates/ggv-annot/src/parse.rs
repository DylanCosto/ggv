use std::{collections::HashMap, io::BufRead};

use crate::{Strand, TAG_BASIC, TAG_CANONICAL, TAG_MANE_PLUS_CLINICAL, TAG_MANE_SELECT, Transcript, Translation, is_mitochondrial};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    Gtf,
    Gff3,
    Bed,
}

impl Format {
    /// From a (lowercased, decompressed-extension) file name; GTF if unknown.
    pub fn from_name(name: &str) -> Self {
        if name.ends_with(".gff3") || name.ends_with(".gff") {
            Self::Gff3
        } else if name.ends_with(".bed") || name.ends_with(".bed12") {
            Self::Bed
        } else {
            Self::Gtf
        }
    }
}

pub(crate) struct Parsed {
    pub contigs: Vec<String>,
    pub transcripts: Vec<Transcript>,
}

#[derive(Default)]
struct Builder {
    contig: u32,
    span: Option<(u64, u64)>,
    strand: Option<Strand>,
    id: String,
    gene_id: String,
    gene_name: String,
    biotype: String,
    exons: Vec<(u64, u64)>,
    cds: Vec<(u64, u64)>,
    /// CDS rows with their phase, to find the frame at the 5' end.
    phases: Vec<(u64, u64, u8)>,
    selenocysteine: Vec<(u64, u64)>,
    tags: u8,
}

impl Builder {
    /// Adds an exon, CDS, start/stop codon or selenocysteine row. `phase` is the row's phase
    /// column, `attrs` its attributes (for `transl_except`).
    fn add(&mut self, kind: &str, range: (u64, u64), phase: &str, attrs: &str) {
        match kind {
            "exon" => self.exons.push(range),
            "Selenocysteine" => self.selenocysteine.push(range),
            _ => {
                self.cds.push(range);
                if kind == "CDS" {
                    self.phases.push((range.0, range.1, phase.parse().ok().filter(|&p| p < 3).unwrap_or(0)));
                    if attrs.contains("aa:Sec") {
                        self.selenocysteine.extend(sec_exceptions(attrs));
                    }
                }
            }
        }
    }

    fn finish(mut self) -> Option<Transcript> {
        let exons_from_span = self.exons.is_empty() && self.cds.is_empty();
        if self.exons.is_empty() {
            // Some files give only CDS for coding transcripts.
            self.exons = self.cds.clone();
        }
        if self.exons.is_empty() && self.span.is_none() {
            return None;
        }
        self.exons.sort_unstable();
        let mut exons: Vec<(u64, u64)> = Vec::with_capacity(self.exons.len());
        for (s, e) in self.exons {
            match exons.last_mut() {
                Some(last) if s <= last.1 => last.1 = last.1.max(e),
                _ => exons.push((s, e)),
            }
        }
        let cds = self.cds.iter().fold(None, |acc: Option<(u64, u64)>, &(s, e)| {
            Some(acc.map_or((s, e), |(a, b)| (a.min(s), b.max(e))))
        });
        let (mut start, mut end) = self.span.unwrap_or((u64::MAX, 0));
        if let (Some(f), Some(l)) = (exons.first(), exons.last()) {
            start = start.min(f.0);
            end = end.max(l.1);
        }
        if exons.is_empty() {
            exons.push((start, end));
        }
        // The 5'-most CDS row's phase says where the first complete codon starts.
        let first = if self.strand == Some(Strand::Reverse) {
            self.phases.iter().max_by_key(|p| p.1)
        } else {
            self.phases.iter().min_by_key(|p| p.0)
        };
        self.selenocysteine.sort_unstable();
        self.selenocysteine.dedup();
        let translation = Translation { phase: first.map_or(0, |p| p.2), table: 1, selenocysteine: self.selenocysteine };
        let gene_name = if self.gene_name.is_empty() { self.gene_id.clone() } else { self.gene_name };
        let gene_id = if self.gene_id.is_empty() { gene_name.clone() } else { self.gene_id };
        Some(Transcript {
            contig: self.contig,
            start,
            end,
            strand: self.strand.unwrap_or(Strand::Unknown),
            id: self.id.into(),
            gene_id: gene_id.into(),
            gene_name: gene_name.into(),
            biotype: self.biotype.into(),
            exons,
            cds,
            tags: self.tags,
            representative: false,
            translation,
            exons_from_span,
        })
    }
}

#[derive(Default)]
struct Contigs {
    names: Vec<String>,
    index: HashMap<String, u32>,
}

impl Contigs {
    fn get(&mut self, name: &str) -> u32 {
        if let Some(&i) = self.index.get(name) {
            return i;
        }
        let i = self.names.len() as u32;
        self.names.push(name.to_owned());
        self.index.insert(name.to_owned(), i);
        i
    }
}

fn strand(s: &str) -> Strand {
    match s {
        "+" => Strand::Forward,
        "-" => Strand::Reverse,
        _ => Strand::Unknown,
    }
}

fn tag_bits(tag: &str) -> u8 {
    match tag.trim() {
        "MANE_Select" => TAG_MANE_SELECT,
        "MANE_Plus_Clinical" => TAG_MANE_PLUS_CLINICAL,
        "Ensembl_canonical" => TAG_CANONICAL,
        "basic" => TAG_BASIC,
        _ => 0,
    }
}

/// 1-based inclusive text coordinates → 0-based half-open.
fn span(start: &str, end: &str) -> Option<(u64, u64)> {
    let s: u64 = start.parse().ok()?;
    let e: u64 = end.parse().ok()?;
    (s >= 1 && e >= s).then(|| (s - 1, e))
}

pub(crate) fn parse(mut reader: impl BufRead, format: Format) -> std::io::Result<Parsed> {
    let mut contigs = Contigs::default();
    let mut line = Vec::with_capacity(4096);
    let mut transcripts = match format {
        Format::Gtf => {
            let mut gtf = Gtf::default();
            while next_line(&mut reader, &mut line)? {
                gtf.line(&String::from_utf8_lossy(&line), &mut contigs);
            }
            gtf.finish()
        }
        Format::Gff3 => {
            let mut gff = Gff3::default();
            while next_line(&mut reader, &mut line)? {
                let text = String::from_utf8_lossy(&line);
                if text.starts_with("##FASTA") {
                    break;
                }
                gff.line(&text, &mut contigs);
            }
            gff.finish()
        }
        Format::Bed => {
            let mut out = Vec::new();
            while next_line(&mut reader, &mut line)? {
                if let Some(t) = bed_line(&String::from_utf8_lossy(&line), &mut contigs) {
                    out.push(t);
                }
            }
            out
        }
    };
    for t in &mut transcripts {
        if is_mitochondrial(&contigs.names[t.contig as usize]) {
            t.translation.table = 2;
        }
    }
    Ok(Parsed { contigs: contigs.names, transcripts })
}

/// Reads one line without its line ending; false at end of input.
fn next_line(reader: &mut impl BufRead, buf: &mut Vec<u8>) -> std::io::Result<bool> {
    buf.clear();
    if reader.read_until(b'\n', buf)? == 0 {
        return Ok(false);
    }
    while buf.last().is_some_and(|b| *b == b'\n' || *b == b'\r') {
        buf.pop();
    }
    Ok(true)
}

/// The 9 tab-separated columns of a GTF/GFF line, or `None` for headers and blank lines.
fn columns(line: &str) -> Option<[&str; 9]> {
    if line.is_empty() || line.starts_with('#') {
        return None;
    }
    let mut it = line.splitn(9, '\t');
    let mut out = [""; 9];
    for slot in &mut out {
        *slot = it.next()?;
    }
    Some(out)
}

#[derive(Default)]
struct Gtf {
    builders: Vec<Builder>,
    by_id: HashMap<String, usize>,
}

impl Gtf {
    fn line(&mut self, line: &str, contigs: &mut Contigs) {
        let Some([seq, _, kind, start, end, _, strand_s, phase, attrs]) = columns(line) else { return };
        if !matches!(kind, "exon" | "CDS" | "start_codon" | "stop_codon" | "Selenocysteine" | "transcript") {
            return;
        }
        let Some(range) = span(start, end) else { return };

        // Fast path: most lines are exons/CDS of a transcript already seen, which need only the
        // transcript id. Finding it directly skips splitting the other ~20 attributes.
        if kind != "transcript"
            && let Some(tid) = quoted_attr(attrs, "transcript_id")
            && let Some(&i) = self.by_id.get(tid)
        {
            let b = &mut self.builders[i];
            if b.gene_name.is_empty() || b.gene_id.is_empty() || b.biotype.is_empty() {
                // Rare: the first line lacked something; take the slow path to fill it in.
            } else {
                b.add(kind, range, phase, attrs);
                return;
            }
        }

        let mut tid = None;
        let (mut gene_id, mut gene_name, mut biotype, mut tags) = (None, None, None, 0u8);
        for field in attrs.split(';') {
            let Some((key, value)) = field.trim().split_once(' ') else { continue };
            let value = value.trim().trim_matches('"');
            match key {
                "transcript_id" => tid = Some(value),
                "gene_id" => gene_id = Some(value),
                "gene_name" => gene_name = Some(value),
                "transcript_type" | "transcript_biotype" => biotype = Some(value),
                "gene_type" | "gene_biotype" if biotype.is_none() => biotype = Some(value),
                "tag" => tags |= tag_bits(value),
                _ => {}
            }
        }
        let Some(tid) = tid else { return };

        let i = match self.by_id.get(tid) {
            Some(&i) => i,
            None => {
                self.builders.push(Builder {
                    contig: contigs.get(seq),
                    id: tid.to_owned(),
                    gene_id: gene_id.unwrap_or_default().to_owned(),
                    gene_name: gene_name.unwrap_or_default().to_owned(),
                    biotype: biotype.unwrap_or_default().to_owned(),
                    ..Default::default()
                });
                self.by_id.insert(tid.to_owned(), self.builders.len() - 1);
                self.builders.len() - 1
            }
        };
        let b = &mut self.builders[i];
        // Fill in what the transcript's first line didn't say.
        for (slot, value) in [(&mut b.gene_id, gene_id), (&mut b.gene_name, gene_name), (&mut b.biotype, biotype)] {
            if slot.is_empty()
                && let Some(v) = value
            {
                *slot = v.to_owned();
            }
        }
        b.strand.get_or_insert(strand(strand_s));
        b.tags |= tags;
        match kind {
            "transcript" => b.span = Some(range),
            _ => b.add(kind, range, phase, attrs),
        }
    }

    fn finish(self) -> Vec<Transcript> {
        self.builders.into_iter().filter_map(Builder::finish).collect()
    }
}

/// The value of `key "value"` in GTF attributes, matching `key` only as a whole word.
fn quoted_attr<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    let bytes = attrs.as_bytes();
    let mut from = 0;
    while let Some(off) = memchr::memmem::find(&bytes[from..], key.as_bytes()) {
        let at = from + off;
        let before_ok = at == 0 || matches!(bytes[at - 1], b' ' | b';' | b'\t');
        let rest = &attrs[at + key.len()..];
        if before_ok && let Some(rest) = rest.strip_prefix(" \"") {
            return rest.split_once('"').map(|(v, _)| v);
        }
        from = at + key.len();
    }
    None
}

/// A GFF3 feature with an `ID`, kept to resolve `Parent` links at the end.
struct Node {
    contig: u32,
    span: (u64, u64),
    strand: Strand,
    parent: Option<String>,
    name: Option<String>,
    gene_id: Option<String>,
    transcript_id: Option<String>,
    biotype: Option<String>,
    tags: u8,
}

#[derive(Default)]
struct Gff3 {
    nodes: HashMap<String, Node>,
    /// Transcript (parent) ID → its exon / CDS parts.
    parts: HashMap<String, Builder>,
    order: Vec<String>,
    /// Node IDs in file order.
    node_order: Vec<String>,
}

impl Gff3 {
    fn line(&mut self, line: &str, contigs: &mut Contigs) {
        let Some([seq, _, kind, start, end, _, strand_s, phase, attrs]) = columns(line) else { return };
        let Some(range) = span(start, end) else { return };
        let contig = contigs.get(seq);

        let (mut id, mut parents, mut name, mut gene_id, mut transcript_id, mut biotype, mut tags) =
            (None, None, None, None, None, None, 0u8);
        for field in attrs.split(';') {
            let Some((key, value)) = field.split_once('=') else { continue };
            match key.trim() {
                "ID" => id = Some(value),
                "Parent" => parents = Some(value),
                "Name" | "gene_name" => name = Some(percent_decode(value)),
                "gene_id" => gene_id = Some(value.to_owned()),
                "transcript_id" => transcript_id = Some(value.to_owned()),
                "biotype" | "transcript_type" | "gene_type" | "gene_biotype" | "transcript_biotype" => {
                    biotype = Some(value.to_owned());
                }
                "tag" | "tags" => tags |= value.split(',').map(tag_bits).fold(0, |a, b| a | b),
                _ => {}
            }
        }

        if matches!(kind, "exon" | "CDS" | "start_codon" | "stop_codon" | "Selenocysteine") {
            for parent in parents.into_iter().flat_map(|p| p.split(',')) {
                let b = self.parts.entry(parent.to_owned()).or_insert_with(|| {
                    self.order.push(parent.to_owned());
                    Builder { contig, strand: Some(strand(strand_s)), ..Default::default() }
                });
                b.add(kind, range, phase, attrs);
            }
            return;
        }

        if let Some(id) = id {
            if !self.nodes.contains_key(id) {
                self.node_order.push(id.to_owned());
            }
            self.nodes.insert(
                id.to_owned(),
                Node {
                    contig,
                    span: range,
                    strand: strand(strand_s),
                    parent: parents.and_then(|p| p.split(',').next()).map(str::to_owned),
                    name,
                    gene_id,
                    transcript_id,
                    biotype,
                    tags,
                },
            );
        }
    }

    fn finish(mut self) -> Vec<Transcript> {
        // Transcripts (children of a top-level gene) with no exon or CDS rows: kept as their
        // span, flagged as such, rather than dropped.
        for id in std::mem::take(&mut self.node_order) {
            let Some(n) = self.nodes.get(&id) else { continue };
            let under_gene = n.parent.as_ref().and_then(|g| self.nodes.get(g)).is_some_and(|g| g.parent.is_none());
            if under_gene && !self.parts.contains_key(&id) {
                self.parts.insert(id.clone(), Builder { contig: n.contig, strand: Some(n.strand), ..Default::default() });
                self.order.push(id);
            }
        }
        let mut out = Vec::with_capacity(self.order.len());
        for parent_id in self.order {
            let Some(mut b) = self.parts.remove(&parent_id) else { continue };
            let node = self.nodes.get(&parent_id);
            let gene = node.and_then(|n| n.parent.as_ref()).and_then(|g| self.nodes.get(g));
            let short = |id: &str| id.split_once(':').map_or(id, |(_, s)| s).to_owned();

            b.id = node.and_then(|n| n.transcript_id.clone()).unwrap_or_else(|| short(&parent_id));
            if let Some(n) = node {
                b.contig = n.contig;
                b.span = Some(n.span);
                b.strand = Some(n.strand);
                b.tags = n.tags;
                b.biotype = n.biotype.clone().unwrap_or_default();
            }
            // Exons hanging straight off a gene: the gene is its own transcript.
            let gene = gene.or(node.filter(|n| n.parent.is_none()));
            let gene_key = node.and_then(|n| n.parent.clone()).unwrap_or_else(|| parent_id.clone());
            b.gene_id = gene.and_then(|g| g.gene_id.clone()).unwrap_or_else(|| short(&gene_key));
            b.gene_name = gene.and_then(|g| g.name.clone()).unwrap_or_else(|| b.gene_id.clone());
            if b.biotype.is_empty() {
                b.biotype = gene.and_then(|g| g.biotype.clone()).unwrap_or_default();
            }
            out.extend(b.finish());
        }
        out
    }
}

fn bed_line(line: &str, contigs: &mut Contigs) -> Option<Transcript> {
    if line.is_empty() || line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
        return None;
    }
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 3 {
        return None;
    }
    let start: u64 = f[1].parse().ok()?;
    let end: u64 = f[2].parse().ok()?;
    if end <= start {
        return None;
    }
    let name = f.get(3).map_or_else(|| format!("{}:{}-{}", f[0], start + 1, end), |s| (*s).to_owned());
    let mut b = Builder {
        contig: contigs.get(f[0]),
        span: Some((start, end)),
        strand: Some(strand(f.get(5).copied().unwrap_or("."))),
        id: name.clone(),
        gene_id: name.clone(),
        gene_name: name,
        ..Default::default()
    };
    if let (Some(ts), Some(te)) = (f.get(6).and_then(|s| s.parse::<u64>().ok()), f.get(7).and_then(|s| s.parse::<u64>().ok()))
        && te > ts
    {
        b.cds.push((ts, te));
    }
    if let (Some(sizes), Some(starts)) = (f.get(10), f.get(11)) {
        for (size, offset) in sizes.split(',').zip(starts.split(',')) {
            if let (Ok(size), Ok(offset)) = (size.parse::<u64>(), offset.parse::<u64>()) {
                b.exons.push((start + offset, start + offset + size));
            }
        }
    }
    if b.exons.is_empty() {
        b.exons.push((start, end));
    }
    // A thick range is a coding span, not a list of exons.
    let thick = b.cds.pop();
    let mut t = b.finish()?;
    t.cds = thick;
    Some(t)
}

/// Selenocysteine codons from a `transl_except` value such as `(pos:1234..1236,aa:Sec)` or
/// `(pos:complement(1234..1236),aa:Sec)`, possibly percent-encoded.
fn sec_exceptions(attrs: &str) -> Vec<(u64, u64)> {
    let text = percent_decode(attrs);
    text.split("(pos:")
        .skip(1)
        .filter_map(|e| {
            let aa = e.find("aa:")?;
            if !e[aa..].starts_with("aa:Sec") {
                return None;
            }
            let body = e[..aa].trim_start_matches("complement(");
            let (a, rest) = body.split_once("..")?;
            let b: String = rest.chars().take_while(char::is_ascii_digit).collect();
            span(a.trim(), &b)
        })
        .collect()
}

/// Decodes `%XX` escapes (GFF3 encodes `;`, `=`, `,` and others in values).
fn percent_decode(s: &str) -> String {
    if !s.contains('%') {
        return s.to_owned();
    }
    let hex = |b: u8| (b as char).to_digit(16).map(|d| d as u8);
    let bytes = s.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'%'
            && i + 2 < bytes.len()
            && let (Some(hi), Some(lo)) = (hex(bytes[i + 1]), hex(bytes[i + 2]))
        {
            out.push(hi << 4 | lo);
            i += 3;
        } else {
            out.push(bytes[i]);
            i += 1;
        }
    }
    String::from_utf8_lossy(&out).into_owned()
}
