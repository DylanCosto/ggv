//! The tools `ggv mcp` offers, and what each returns.
//!
//! Every answer is compact JSON built for a model to read: summaries first (counts, fractions,
//! medians, warnings), raw reads only when asked for and capped, and a `query` block saying
//! which files, filters and region produced it. Numbers come from the same engines as the
//! desktop app and `ggv query`, so a person checking the app sees the same thing.
//!
//! Open files are kept between calls (keyed by path), and an annotation is parsed once, so
//! repeated questions about one sample cost milliseconds.

use std::{
    collections::HashMap,
    path::{Path, PathBuf},
    sync::Arc,
};

use anyhow::{Context, Result, anyhow, bail};
use ggv_annot::{Annotation, coding};
use ggv_io::{BamFile, ReferenceFasta, VariantFile};
use ggv_layout::{AlignedRead, ViewOptions, load_view};
use ggv_model::Region;
use ggv_query::{PileupOptions, PositionStats, SummaryOptions};
use serde_json::{Value, json};

/// Files used when a call names none (from `ggv mcp -b … -f … -v … -g …`, or `open`).
#[derive(Clone, Debug, Default)]
pub struct Defaults {
    pub bams: Vec<PathBuf>,
    pub reference: Option<PathBuf>,
    pub vcfs: Vec<PathBuf>,
    pub genes: Vec<PathBuf>,
    /// Signal tracks (bigWig, bedGraph, SEG).
    pub signals: Vec<PathBuf>,
}

pub enum Output {
    Json(Value),
    Image { png: Vec<u8>, note: String },
}

pub struct State {
    defaults: Defaults,
    /// `ggv mcp --attach`: the desktop window whose state fills in what calls leave out.
    pub desktop: Option<super::desktop::Desktop>,
    bams: HashMap<(PathBuf, Option<PathBuf>), BamFile>,
    /// Behind a `RefCell`: consequences look bases up from inside a shared-borrow closure.
    fastas: std::cell::RefCell<HashMap<PathBuf, ReferenceFasta>>,
    vcfs: HashMap<PathBuf, VariantFile>,
    annotations: HashMap<PathBuf, Arc<Annotation>>,
    /// Each cached file's size and modification time when opened: a file replaced on disk (a
    /// pipeline re-run) is opened again rather than read through a stale handle.
    stamps: std::cell::RefCell<HashMap<PathBuf, Stamp>>,
}

/// One pass over `region`: the exact counts at the positions `want` picks (0-based) and the
/// reads that pass the filters, decoded as `load_view` decodes them (`seq` is the reference
/// for exactly `region`).
fn counts_and_reads(
    bam: &mut BamFile,
    region: &Region,
    seq: Option<&[u8]>,
    opts: &PileupOptions,
    want: impl Fn(u64) -> bool,
) -> Result<(HashMap<u64, PositionStats>, Vec<AlignedRead>)> {
    let mut stats = HashMap::new();
    let mut reads = Vec::new();
    let refslice = seq.map(|r| (r, region.start));
    ggv_query::pileup_with_reads(
        bam,
        region,
        seq,
        opts,
        |read| {
            if read.position() >= 0 {
                let mut r = ggv_layout::decode_with(read, refslice, opts.min_base_quality);
                if refslice.is_none() {
                    r.bases = Some(ggv_layout::read_bases(read));
                }
                reads.push(r);
            }
        },
        |s| {
            if want(s.pos - 1) {
                stats.insert(s.pos - 1, s);
            }
            Ok(())
        },
    )?;
    Ok((stats, reads))
}

/// Counts and reads at bases of a several-base `evidence` call, read once per nearby group:
/// keyed by file index, contig and 0-based position.
type Prefetch = HashMap<(usize, String, u64), (Option<PositionStats>, Arc<Vec<AlignedRead>>)>;

/// Regions wider than this get no reads (as in the app).
const READ_WINDOW: u64 = 100_000;
const MAX_LOCI: usize = 50;
/// Heterozygous sites `phase` lists unless asked for more: a satellite region has hundreds, at
/// about 100 characters each.
const DEFAULT_PHASE_SITES: u64 = 50;
const MAX_PHASE_SITES: u64 = 5000;
/// `repeats` warns that its region may be too narrow once this many reads were set aside as
/// partial and they outnumber the reads measured.
const MIN_PARTIAL_FOR_WARNING: u32 = 3;
/// Bases of a REF or ALT allele an answer spells out. A longer one (an SV written out in full)
/// is cut to its first [`ALLELE_HEAD`] bases with its length, which `ref_length` / `alt_lengths`
/// also give.
const MAX_ALLELE_CHARS: usize = 60;
const ALLELE_HEAD: usize = 20;

// ---- Tool descriptions ----

fn files_schema() -> Value {
    json!({
        "bams": {"type": "array", "items": {"type": "string"}, "description": "Alignment files (BAM/CRAM). Default: the server's."},
        "reference": {"type": "string", "description": "Indexed reference FASTA. Default: the server's."},
    })
}

fn filter_schema() -> Value {
    json!({
        "filter": {"type": "string", "description": "Keep only reads matching this expression, e.g. \"mapq >= 20 && !dup\", \"hp == 1\", or \"!hp\" for reads with no HP tag (a missing tag equals nothing, so \"hp == 0\" matches none). Fields: mapq, flag, start, end, len, tlen, insert, strand, reverse, dup, secondary, supplementary, qcfail, paired, proper, read1, read2, mate_unmapped, mate_elsewhere, hp, rg, name, mismatches, deletions, insertions, clips, mods, split; operators && || ! == != < <= > >= ~ (contains) & (bits)."},
        "min_mapq": {"type": "integer", "minimum": 0, "maximum": 255, "description": "Minimum mapping quality (default 0)."},
        "include_duplicates": {"type": "boolean", "description": "Count duplicate reads (default false: dropped, like IGV and samtools)."},
        "exclude_flags": {"type": ["string", "integer"], "description": "Drop reads with any of these SAM flag bits, e.g. \"0x704\" (the default: unmapped, secondary, QC-fail, duplicate) or \"0xF04\" to drop supplementary too. Overrides include_duplicates."},
        "min_base_quality": {"type": "integer", "minimum": 0, "maximum": 255, "description": "Bases below this quality are not counted (default 0)."}
    })
}

/// `compare`'s read filters for one of its two files, in place of the call's own.
fn per_file_filters(which: &str) -> Value {
    json!({
        "type": "object",
        "description": format!("Read filters for `{which}` alone, in place of the call's own (e.g. {{\"filter\": \"hp == 1\"}} against {{\"filter\": \"hp == 2\"}} compares one file's haplotypes)."),
        "properties": {"filter": {"type": "string"}, "min_mapq": {"type": "integer", "minimum": 0, "maximum": 255}, "include_duplicates": {"type": "boolean"}, "exclude_flags": {"type": ["string", "integer"]}, "min_base_quality": {"type": "integer", "minimum": 0, "maximum": 255}},
    })
}

fn schema(required: &[&str], props: &[Value]) -> Value {
    let mut p = serde_json::Map::new();
    for v in props {
        if let Value::Object(m) = v {
            p.extend(m.clone());
        }
    }
    json!({"type": "object", "properties": p, "required": required})
}

fn tool(name: &str, title: &str, description: &str, input: Value) -> Value {
    json!({"name": name, "title": title, "description": description, "inputSchema": input, "annotations": {"readOnlyHint": true, "openWorldHint": false}})
}

/// The arguments every data tool reads to find its files and place a locus (contig names and
/// lengths come from the alignment file or the reference; a region given as a gene name, from
/// the annotations). Published for each tool, so that naming one is never refused as unknown.
fn add_file_arguments(tool: &mut Value) {
    let common = json!({
        "bams": {"type": "array", "items": {"type": "string"}, "description": "Alignment files (BAM/CRAM). Default: the server's."},
        "bam": {"type": "string", "description": "One alignment file: the same as `bams` with one item."},
        "reference": {"type": "string", "description": "Indexed reference FASTA (or .2bit). Default: the server's."},
        "genes": {"type": "array", "items": {"type": "string"}, "description": "Annotations (GTF/GFF3/BED), for a region given as a gene name. Default: the server's."},
    });
    let (Some(props), Value::Object(common)) = (tool["inputSchema"]["properties"].as_object_mut(), common) else { return };
    for (name, schema) in common {
        props.entry(name).or_insert(schema);
    }
}

pub fn list() -> Vec<Value> {
    let mut tools = data_tools();
    // `info` and `open` take a list of files of any kind instead.
    for tool in tools.iter_mut().filter(|t| t["name"] != "info" && t["name"] != "open") {
        add_file_arguments(tool);
    }
    tools
}

fn data_tools() -> Vec<Value> {
    let locus = json!({"locus": {"type": "string", "description": "A base, e.g. chr17:43,092,919 (1-based)."}});
    let region = json!({"region": {"type": "string", "description": "e.g. chr17:43,090,000-43,100,000 (1-based, inclusive)."}});
    vec![
        tool(
            "info",
            "Files and contigs",
            "The files open by default and what is in them: samples and read groups of each alignment file, contigs (the main ones, and how many in all), VCF samples, annotation size. Call this first.",
            schema(&[], &[json!({"files": {"type": "array", "items": {"type": "string"}, "description": "Other files to describe (BAM/CRAM, VCF/BCF, GTF/GFF3/BED, FASTA, bigWig/bedGraph/SEG)."}})]),
        ),
        {
            // Changes what later calls use and, attached, the window's tracks.
            let mut open = tool(
                "open",
                "Set default files",
                "Adds files to the defaults used when a call names none. Files are sorted by extension: .bam/.cram, .vcf.gz/.bcf, .gtf/.gff3/.bed, .fa/.fasta/.2bit (the reference), .bw/.bigwig/.bedgraph/.seg (signal tracks).",
                schema(&["files"], &[json!({"files": {"type": "array", "items": {"type": "string"}}, "replace": {"type": "boolean", "description": "Replace the defaults instead of adding (default false). Attached: close the window's files not given (the reference stays unless a FASTA is given)."}})]),
            );
            open["annotations"] = json!({"readOnlyHint": false, "destructiveHint": false, "idempotentHint": true, "openWorldHint": false});
            open
        },
        tool(
            "evidence",
            "Evidence at a base",
            "Everything about one base (or up to 50): for each alignment file, depth and reads per allele (A/C/G/T/N/DEL) by strand, and per allele the median MAPQ, median distance of the base from the read end, and soft-clipped share; mean MAPQ and base quality; warnings for common artifacts (one-strand support, alt only near read ends, alt reads with lower MAPQ, low depth, many soft clips); the VCF records there with genotypes; and each change's consequence for the genes (HGVS c./p., missense, splice, frameshift...).",
            schema(&[], &[locus.clone(), json!({"loci": {"type": "array", "items": {"type": "string"}, "maxItems": MAX_LOCI, "description": "Several bases at once."}, "vcfs": {"type": "array", "items": {"type": "string"}}, "genes": {"type": "array", "items": {"type": "string"}, "description": "Annotations for consequences (default: the server's or window's; [] for none)."}, "transcripts": {"type": "string", "enum": ["main", "all"], "description": "main (default): one transcript per gene, its MANE Select (from the file's tags or NCBI/EMBL-EBI's MANE list; else MANE Plus Clinical, else canonical, else the curated one with the longest coding sequence), or where that one doesn't reach, the gene's best one there (curated before predicted XM_/XR_ models); all: every transcript there. Answers say how many were left out."}}), files_schema(), filter_schema()]),
        ),
        tool(
            "reads",
            "Reads, grouped",
            "Reads at a base or in a region (up to 100 kb), grouped by allele at the base (default for a base), strand, HP tag (haplotype), allele-based phase, or not at all. Each group: count, strands, median MAPQ, insert size, soft clips, and (at a base) median distance from the read end. Raw reads are listed only when `list` is given (at most 50 per group).",
            schema(&[], &[locus, region.clone(), json!({"bam": {"type": "string"}, "reference": {"type": "string"}, "group_by": {"type": "string", "enum": ["allele", "strand", "haplotype", "phase", "none"]}, "list": {"type": "integer", "minimum": 0, "maximum": 50, "description": "Raw reads to list per group (default 0)."}}), filter_schema()]),
        ),
        tool(
            "summarize",
            "Region summary",
            "For a region (any size, though large ones take time): reads used and why others were dropped, mean/median/min/max depth, breadth at 1/10/20/30x, and the most variable positions with their alt allele, fraction and strand split. Same as `ggv query summary`.",
            schema(&["region"], &[region.clone(), json!({"max_sites": {"type": "integer", "minimum": 0, "maximum": 100, "description": "Variable positions to list (default 10)."}}), files_schema(), filter_schema()]),
        ),
        tool(
            "coverage",
            "Depth profile",
            "Depth across a region in bins: counted from the reads up to 5 Mb, estimated from the index beyond (whole chromosomes in milliseconds). Returns mean (and max) depth per bin.",
            schema(&["region"], &[region.clone(), json!({"bins": {"type": "integer", "minimum": 1, "maximum": 1000, "description": "Default 50."}}), files_schema(), filter_schema()]),
        ),
        tool(
            "variants",
            "VCF records",
            "VCF/BCF records in a region: position, ID, REF/ALT, type (MIXED when the ALTs differ in type, with alt_types per ALT; BND for breakends), QUAL, FILTER, AF (the site's frequency of any non-reference allele: INFO/AF summed over the ALTs, else from the genotypes; alt_af per ALT when there are several) and each sample's GT/DP/GQ/AD.",
            schema(&["region"], &[region.clone(), json!({"vcfs": {"type": "array", "items": {"type": "string"}}, "max": {"type": "integer", "minimum": 1, "maximum": 500, "description": "Default 50."}, "pass_only": {"type": "boolean"}})]),
        ),
        tool(
            "genes",
            "Genes and transcripts",
            "Transcripts in a region, or a gene by name or ID: strand, span, exons, coding span, biotype, MANE/canonical tags. One transcript per gene unless `all` is set. Each file's answer gives the `total` found and whether the list was `truncated` at `max`.",
            schema(&[], &[region, json!({"name": {"type": "string", "description": "Gene name, gene ID or transcript ID."}, "genes": {"type": "array", "items": {"type": "string"}}, "all": {"type": "boolean"}, "max": {"type": "integer", "minimum": 1, "maximum": 10000, "description": "Transcripts to list per file (default 100)."}})]),
        ),
        tool(
            "consequence",
            "What a change does",
            "The consequence of a change (VCF-style REF and ALT at a 1-based position) for every transcript there: HGVS c. and p., effect (missense, synonymous, stop gained, frameshift, splice donor...), exon or intron, and `warnings` when the result is uncertain. Indels are shifted 3' as HGVS requires, within their exon or intron.",
            schema(&["locus", "ref", "alt"], &[json!({"locus": {"type": "string"}, "ref": {"type": "string"}, "alt": {"type": "string"}, "genes": {"type": "array", "items": {"type": "string"}}, "reference": {"type": "string"}, "transcripts": {"type": "string", "enum": ["main", "all"], "description": "main (default): one transcript per gene, its MANE Select (from the file's tags or NCBI/EMBL-EBI's MANE list; else MANE Plus Clinical, else canonical, else the curated one with the longest coding sequence), or where that one doesn't reach, the gene's best one there (curated before predicted XM_/XR_ models); all: every transcript there. Answers say how many were left out."}})]),
        ),
        tool(
            "junctions",
            "Splice junctions",
            "RNA-seq splice junctions in a region (up to 10 Mb), for every alignment file side by side: each junction's intron (1-based, inclusive), length, strand and motif (GT/AG, GC/AG, AT/AC or non-canonical, from the reference), whether the gene annotation knows it (known, novel (known sites) — e.g. exon skipping, novel donor, novel acceptor, novel) and its genes; and per file its reads, uniquely mapped reads, and its share of the reads leaving its donor and entering its acceptor. Shifts in those shares between samples are aberrant splicing.",
            schema(&[], &[json!({"region": {"type": "string", "description": "e.g. chr17:7,668,421-7,687,490 or a gene."}, "genes": {"type": "array", "items": {"type": "string"}}, "min_reads": {"type": "integer", "minimum": 1, "description": "Keep junctions with at least this many reads in some file (default 2)."}, "min_anchor": {"type": "integer", "minimum": 1, "description": "Count a read only if it aligns this many bases on each side (default 1)."}, "max": {"type": "integer", "minimum": 1, "maximum": 2000, "description": "At most this many junctions, the most-read (default 200)."}}), files_schema(), filter_schema()]),
        ),
        tool(
            "svs",
            "Structural variants against the reads",
            "Structural variants (≥ `min_len` bp, default 50) in a region up to 2 Mb: each VCF record (sequence-resolved, or symbolic <DEL>/<INS>/<DUP>/<INV>, each ALT on its own with its own SVLEN) with its type, length and the VCF's genotype per sample (`genotypes`, by sample name), and for every alignment file the reads that support it (a CIGAR deletion or insertion of about its size near the breakpoint, or a split alignment across it; for a tandem <DUP>, an insertion of about its size in or beside it, or a split read from its end back to its start), the reads that span it showing the reference, reads with a different allele there (another-sized event), clipped reads at its breakpoints, the allele fraction and a genotype from the reads, and both split by HP tag. With `discover` (default true), clusters of large deletions and insertions in the reads that no record in the region explains (all of them, not only those tested): what the call set may have missed. Breakends (BND) have no length, so `min_len` doesn't apply: a pair on the same contig is tested once (its `mate` listed with it), as a deletion or tandem duplication when its brackets say so (`t[p[` with p after t is a deletion), else as a rearrangement between the two breakpoints; others are listed under `breakends_not_tested`. An <INS> without SVLEN can't be tested and is named in `warnings`.",
            schema(&["region"], &[json!({"region": {"type": "string"}, "vcfs": {"type": "array", "items": {"type": "string"}}, "min_len": {"type": "integer", "minimum": 20, "description": "Default 50."}, "max": {"type": "integer", "minimum": 1, "maximum": 200, "description": "Records to test, the longest first (default 30)."}, "discover": {"type": "boolean"}, "min_reads": {"type": "integer", "minimum": 1, "description": "Reads a discovered cluster needs (default 3)."}}), files_schema(), filter_schema()]),
        ),
        tool(
            "repeats",
            "Tandem repeat lengths",
            "How long a tandem repeat is in the reads (e.g. HTT CAG, FMR1 CGG, a VNTR): for every alignment file, each read that spans the repeat `region` with flank on both sides gives its length on that read's haplotype (read bases between the repeat's ends). Returns the reference length, the alleles (unphased reads: one per cluster of read lengths, each of two reads or more, so a diploid sample's two alleles and any mosaic third; with phased reads one per HP tag, untagged reads joining the one they match and forming their own allele when two or more match neither (a lone one is counted as `unassigned`); every spanning read is in one or in `unassigned`: median length, range, reads, HP split, change from the reference, and copies of `motif` if given), the median length per HP tag, and reads that end inside the repeat (too short to measure it), with a `warning` when they outnumber the reads measured (the region is probably narrower than the repeat tract). Long reads span most repeats whole, so expansions are measured, not inferred.",
            schema(&["region"], &[json!({"region": {"type": "string", "description": "The repeat itself, e.g. chr4:3,074,877-3,074,939 (1-based, inclusive)."}, "motif": {"type": "string", "description": "The repeat unit, e.g. CAG: copies are reported as length / motif length."}, "list": {"type": "boolean", "description": "Also list every read's length (default false)."}}), files_schema(), filter_schema()]),
        ),
        tool(
            "phase",
            "Haplotypes from the reads",
            "Splits reads into two haplotypes by the alleles they carry at heterozygous sites found in the region (up to 100 kb; needs a reference). Returns the phase blocks (per block its span, sites and reads per haplotype: A and B only match up within a block), reads per haplotype, and for the first `max_sites` sites their block and which allele each haplotype carries (`het_sites` is the count of all of them). Agrees with trio phasing on the HG002 fixtures.",
            schema(&["region"], &[json!({"region": {"type": "string"}, "bam": {"type": "string"}, "reference": {"type": "string"}, "max_sites": {"type": "integer", "minimum": 0, "maximum": MAX_PHASE_SITES, "description": "Sites to list, leftmost first (default 50; a repeat-rich region has hundreds)."}}), filter_schema()]),
        ),
        tool(
            "compare",
            "Differences between two files",
            "Positions in a region (up to 1 Mb) where an allele's fraction differs between two alignment files (e.g. tumor and normal) by more than a threshold, with enough depth in both: the allele (A, C, G, T, N, DEL, or INS for reads with an insertion right after the base), both fractions and depths. The region is read against `a` and `b`.",
            schema(&["region"], &[json!({"region": {"type": "string"}, "a": {"type": "string", "description": "Default: the first default BAM."}, "b": {"type": "string", "description": "Default: the second."}, "threshold": {"type": "number", "minimum": 0.01, "maximum": 1.0, "description": "Default 0.2."}, "min_depth": {"type": "integer", "minimum": 1, "maximum": 4294967295u64, "description": "Default 10."}, "max": {"type": "integer", "minimum": 1, "maximum": 500, "description": "Default 50."}, "reference": {"type": "string"},
                "a_filters": per_file_filters("a"), "b_filters": per_file_filters("b")}), filter_schema()]),
        ),
        tool(
            "signal",
            "Signal values",
            "Values of signal tracks (bigWig, bedGraph, SEG: ChIP/ATAC/RNA coverage, conservation, GC, methylation summaries, copy number) across a region in equal bins: mean, min and max over the bases with data, per file (and per sample for SEG). Wide regions use bigWig zoom levels, so whole chromosomes are quick.",
            schema(&["region"], &[region_schema_for_signal(), json!({"signals": {"type": "array", "items": {"type": "string"}, "description": "Signal files. Default: the server's."}, "bins": {"type": "integer", "minimum": 1, "maximum": 1000, "description": "Default 20."}})]),
        ),
        tool(
            "render",
            "Picture",
            "A PNG of a base (with 50 bp each side, outlined) or a region, every alignment file stacked, as IGV would show it. For clients that can see images, or to show a person; the numbers are in `evidence`.",
            schema(&[], &[json!({"locus": {"type": "string"}, "region": {"type": "string"}, "width": {"type": "integer", "minimum": 300, "maximum": 2400, "description": "Default 1000."}, "color_by": {"type": "string", "description": "Read colours: none (grey, default), strand, haplotype (HP tag), insert-size (pair orientation and insert size), base-mods (methylation, MM/ML), or tag:XX (one colour per value of any tag, e.g. tag:CB)."}, "group_by": {"type": "string", "description": "Group reads into labelled bands: none (default), strand, haplotype (HP tag), read group, alleles (heterozygous sites), base (for a single base: what each read shows there), or tag:XX (by any tag's value, e.g. tag:CB, tag:PS)."}, "link": {"type": "boolean", "description": "Put each split read's parts (primary and supplementary) on one row, joined. Default false."}, "sashimi": {"type": "boolean", "description": "Draw splice junctions (RNA-seq) under each file's coverage: an arc per intron with its read count, + strand above and − below. Default false."}, "min_junction_reads": {"type": "integer", "minimum": 1, "description": "With sashimi: leave out junctions with fewer reads (default 2)."}}), files_schema(), filter_schema()]),
        ),
    ]
}

// ---- Helpers ----

fn region_schema_for_signal() -> Value {
    json!({"region": {"type": "string", "description": "e.g. chr17:43,090,000-43,100,000 or chr17 (1-based, inclusive)."}})
}

/// A pane's state as the reply to a command that changed it (`navigate`, `display`, `open`):
/// without the per-file assembly table, which `view` gives and which only matters here when
/// the files disagree.
fn pane_reply(mut pane: Value) -> Value {
    let agree = pane["assembly"]["conflicts"].as_array().is_none_or(Vec::is_empty);
    if agree && let Some(p) = pane.as_object_mut() {
        p.remove("assembly");
    }
    pane
}

fn pane_check(state: &Value, pane: u64) -> Result<()> {
    let n = state["panes"].as_array().map_or(0, Vec::len);
    if pane == 0 || pane as usize > n {
        bail!("no pane {pane}: the window has {n}");
    }
    Ok(())
}

/// Adds what came from the desktop window to an answer.
fn annotate(out: Output, info: Value) -> Output {
    match out {
        Output::Json(mut v) => {
            if v.is_object() {
                v["desktop"] = info;
            } else {
                v = json!({"result": v, "desktop": info});
            }
            Output::Json(v)
        }
        Output::Image { png, note } => Output::Image { png, note: format!("{note}\nfrom the GGV window: {info}") },
    }
}

/// A list of paths: a list of strings, or one string. Anything else is refused, never dropped.
fn strs(v: &Value, key: &str) -> Result<Vec<PathBuf>> {
    match v.get(key) {
        None | Some(Value::Null) => Ok(Vec::new()),
        Some(Value::String(s)) => Ok(vec![PathBuf::from(s)]),
        Some(Value::Array(a)) => a.iter().enumerate().map(|(i, x)| x.as_str().map(PathBuf::from).ok_or_else(|| anyhow!("`{key}[{i}]` must be a string, not {x}"))).collect(),
        Some(other) => bail!("`{key}` must be a list of strings, not {other}"),
    }
}

/// An optional non-negative integer argument; a supplied value of another type is an error.
fn uint_arg(v: &Value, key: &str) -> Result<Option<u64>> {
    match v.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(x) => x.as_u64().map(Some).ok_or_else(|| anyhow!("`{key}` must be a non-negative integer, not {x}")),
    }
}

fn uint_or(v: &Value, key: &str, default: u64) -> Result<u64> {
    Ok(uint_arg(v, key)?.unwrap_or(default))
}

/// As [`uint_or`], checked to fit 32 bits.
fn u32_or(v: &Value, key: &str, default: u32) -> Result<u32> {
    let n = uint_or(v, key, u64::from(default))?;
    u32::try_from(n).map_err(|_| anyhow!("`{key}` must be at most {}, not {n}", u32::MAX))
}

fn bool_arg(v: &Value, key: &str) -> Result<Option<bool>> {
    match v.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(x) => x.as_bool().map(Some).ok_or_else(|| anyhow!("`{key}` must be true or false, not {x}")),
    }
}

fn str_arg<'a>(v: &'a Value, key: &str) -> Option<&'a str> {
    v.get(key).and_then(Value::as_str).filter(|s| !s.trim().is_empty())
}

/// An error opening `path`, with the path in front unless its message names the file already
/// (most do: said twice, the path reads as two failures).
fn named(e: impl Into<anyhow::Error>, path: &Path) -> anyhow::Error {
    let e = e.into();
    let shown = path.display().to_string();
    if format!("{e:#}").contains(&shown) { e } else { e.context(shown) }
}

fn file_name(p: &Path) -> String {
    p.file_name().map_or_else(|| p.display().to_string(), |n| n.to_string_lossy().into_owned())
}

fn round(x: f64, places: i32) -> f64 {
    let m = 10f64.powi(places);
    (x * m).round() / m
}

fn median(mut v: Vec<f64>) -> Option<f64> {
    if v.is_empty() {
        return None;
    }
    v.sort_by(f64::total_cmp);
    let n = v.len();
    Some(if n % 2 == 1 { v[n / 2] } else { (v[n / 2 - 1] + v[n / 2]) / 2.0 })
}

/// Bin `i` of `n` over `len` bases covers offsets `[i * len / n, (i + 1) * len / n)`.
fn bin_bounds(i: usize, n: usize, len: u64) -> (u64, u64) {
    (i as u64 * len / n as u64, (i as u64 + 1) * len / n as u64)
}

/// The bin of offset `o`, with the same boundaries as [`bin_bounds`].
fn bin_of(o: u64, n: usize, len: u64) -> usize {
    (((o + 1) * n as u64 - 1) / len) as usize
}

/// When a file's contig ends before the region asked for: says what was used.
fn clipped(asked: &Region, local: &Region) -> Option<String> {
    (local.end < asked.end).then(|| {
        format!(
            "{} ends at {} in this file: only {}:{}-{} ({} bp of the {} asked for) is reported",
            local.contig,
            local.end,
            local.contig,
            local.start + 1,
            local.end,
            local.len(),
            asked.len()
        )
    })
}

/// The read filters a call gives. A supplied filter of the wrong type is an error: it is never
/// treated as left out.
fn filters(args: &Value) -> Result<PileupOptions> {
    let mut o = PileupOptions::default();
    if bool_arg(args, "include_duplicates")? == Some(true) {
        o.exclude_flags &= !0x400;
    }
    match args.get("exclude_flags") {
        None | Some(Value::Null) => {}
        Some(Value::String(s)) if s.trim().is_empty() => {}
        Some(Value::String(s)) => o.exclude_flags = crate::parse_flags(s.trim()).map_err(|e| anyhow!("exclude_flags: {e}"))? | 0x4,
        Some(Value::Number(n)) => o.exclude_flags = u16::try_from(n.as_u64().unwrap_or(u64::MAX)).map_err(|_| anyhow!("exclude_flags: not a 16-bit flag value"))? | 0x4,
        Some(other) => bail!("`exclude_flags` must be a string such as \"0x704\", not {other}"),
    }
    let quality = |key: &str| -> Result<Option<u8>> {
        uint_arg(args, key)?.map(|q| u8::try_from(q).map_err(|_| anyhow!("`{key}` must be at most 255, not {q}"))).transpose()
    };
    if let Some(q) = quality("min_mapq")? {
        o.min_mapping_quality = q;
    }
    if let Some(q) = quality("min_base_quality")? {
        o.min_base_quality = q;
    }
    match args.get("filter") {
        None | Some(Value::Null) => {}
        Some(Value::String(t)) if t.trim().is_empty() => {}
        Some(Value::String(text)) => {
            let expr = ggv_layout::expr::Expr::parse(text).map_err(|e| anyhow!("filter: {e}"))?;
            o.read_filter = Some(ggv_layout::expr::predicate(Arc::new(expr), text));
        }
        Some(other) => bail!("`filter` must be an expression string such as \"mapq >= 20\", not {other}"),
    }
    Ok(o)
}

fn filters_json(o: &PileupOptions) -> Value {
    json!({
        "exclude_flags": format!("0x{:x}", o.exclude_flags),
        "min_mapq": o.min_mapping_quality,
        "min_base_quality": o.min_base_quality,
        "expression": o.read_filter.as_ref().map(|p| p.text.clone()),
    })
}

/// What one read shows at a base.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
enum Allele {
    Base(u8),
    Deletion,
    Skip,
    /// Covered, but below the minimum base quality.
    LowQuality,
    /// Covered, but the base can't be told: no reference to compare with and no stored base.
    Unknown,
}

impl Allele {
    fn name(self) -> String {
        match self {
            Self::Base(b) => (b as char).to_string(),
            Self::Deletion => "DEL".into(),
            Self::Skip => "N-skip".into(),
            Self::LowQuality => "below min_base_quality".into(),
            Self::Unknown => "base unknown".into(),
        }
    }
}

/// The allele a read carries at `pos` and how far the base is from the read's nearer end
/// (in read bases, soft clips included), or `None` if the read doesn't cover it.
fn read_at(r: &AlignedRead, pos: u64, ref_base: Option<u8>) -> Option<(Allele, u32)> {
    if pos < r.start || pos >= r.end {
        return None;
    }
    if r.ref_skips.iter().any(|d| (d.pos..d.pos + u64::from(d.len)).contains(&pos)) {
        return Some((Allele::Skip, 0));
    }
    if r.is_low_quality(pos) {
        // Covered, but the base is below the minimum quality: not counted as any allele.
        return Some((Allele::LowQuality, 0));
    }
    // Every base the read has: clips, aligned blocks and insertions (a terminal one too).
    let qlen = r.soft_clip_left + r.blocks.iter().map(|b| b.len).sum::<u32>() + r.insertions.iter().map(|i| i.len).sum::<u32>() + r.soft_clip_right;
    if let Some(b) = r.blocks.iter().find(|b| b.ref_start <= pos && pos < b.ref_start + u64::from(b.len)) {
        let off = b.query_start + (pos - b.ref_start) as u32;
        let end = off.min(qlen.saturating_sub(1 + off));
        // A mismatch, else the reference base; loaded without a reference, the read's own base.
        let base = match r.mismatches.iter().find(|m| m.pos == pos) {
            Some(m) => Some(m.base),
            None if r.bases.is_some() => r.base_at(pos),
            None => ref_base,
        };
        return Some((base.map_or(Allele::Unknown, Allele::Base), end));
    }
    if r.deletions.iter().any(|d| (d.pos..d.pos + u64::from(d.len)).contains(&pos)) {
        // Distance of the deletion from the end: the bases before it.
        let before = r.blocks.iter().filter(|b| b.ref_start < pos).map(|b| b.query_start + b.len).max().unwrap_or(0);
        return Some((Allele::Deletion, before.min(qlen.saturating_sub(before))));
    }
    None
}

/// Statistics for a group of reads.
fn group_json(reads: &[(&AlignedRead, u32)], at_base: bool) -> Value {
    let n = reads.len();
    let fwd = reads.iter().filter(|(r, _)| !r.reverse).count();
    let mapq = median(reads.iter().map(|(r, _)| f64::from(r.mapq)).collect());
    let insert = median(reads.iter().filter(|(r, _)| r.mate.is_some_and(|m| m.same_contig) && r.template_length != 0).map(|(r, _)| f64::from(r.template_length.unsigned_abs())).collect());
    let clipped = reads.iter().filter(|(r, _)| r.soft_clip_left + r.soft_clip_right > 0).count();
    let mut g = json!({
        "reads": n,
        "forward": fwd,
        "reverse": n - fwd,
        "mapq_median": mapq,
        "mapq0": reads.iter().filter(|(r, _)| r.mapq == 0).count(),
        "insert_median": insert,
        "soft_clipped_share": if n == 0 { 0.0 } else { round(clipped as f64 / n as f64, 3) },
    });
    if at_base {
        g["end_distance_median"] = json!(median(reads.iter().map(|(_, d)| f64::from(*d)).collect()));
    }
    g
}

fn read_json(r: &AlignedRead, at: Option<(Allele, u32)>) -> Value {
    let mut v = json!({
        "name": r.name,
        "start": r.start + 1,
        "end": r.end,
        "strand": if r.reverse { "-" } else { "+" },
        "mapq": r.mapq,
        "flags": r.flags,
        "mismatches": r.mismatches.len(),
        "deletions": r.deletions.len(),
        "insertions": r.insertions.len(),
        "soft_clips": [r.soft_clip_left, r.soft_clip_right],
    });
    if let Some(h) = r.haplotype {
        v["hp"] = json!(h);
    }
    if let Some(m) = r.mate {
        v["mate"] = if m.same_contig { json!({"pos": m.pos + 1, "strand": if m.reverse { "-" } else { "+" }, "insert": r.template_length}) } else { json!({"other_contig": true}) };
    }
    if !r.supplementary.is_empty() {
        v["split"] = json!(r.supplementary.iter().map(|s| format!("{}:{}", s.contig, s.pos + 1)).collect::<Vec<_>>());
    }
    if let Some((a, d)) = at {
        v["allele"] = json!(a.name());
        v["end_distance"] = json!(d);
    }
    v
}

impl State {
    pub fn new(defaults: Defaults) -> Self {
        Self { defaults, desktop: None, bams: HashMap::new(), fastas: Default::default(), vcfs: HashMap::new(), annotations: HashMap::new(), stamps: Default::default() }
    }

    fn bam_paths(&self, args: &Value) -> Result<Vec<PathBuf>> {
        let mut v = strs(args, "bams")?;
        v.extend(strs(args, "bam")?);
        if v.is_empty() {
            v = self.defaults.bams.clone();
        }
        if v.is_empty() {
            bail!("no alignment files: pass `bams`, or start the server with -b FILE (or call `open`)");
        }
        Ok(v)
    }

    /// The reference to use: `reference`, else the default. One that is named (or set as the
    /// default) but can't be opened is an error, never a quiet switch to reference-free
    /// answers; with none at all, answers are reference-free.
    fn reference_path(&self, args: &Value) -> Result<Option<PathBuf>> {
        let path = match args.get("reference") {
            None | Some(Value::Null) => self.defaults.reference.clone(),
            Some(Value::String(s)) if s.trim().is_empty() => self.defaults.reference.clone(),
            Some(Value::String(s)) => Some(PathBuf::from(s)),
            Some(other) => bail!("`reference` must be a path, not {other}"),
        };
        if let Some(p) = &path {
            self.with_fasta(p, |_| ()).with_context(|| format!("reference {}", p.display()))?;
        }
        Ok(path)
    }

    /// The read filters in `args`. A filter testing mismatches gets the reference over each
    /// whole read, and is refused without a reference (every read would have none).
    fn read_filters(&self, args: &Value) -> Result<PileupOptions> {
        let mut o = filters(args)?;
        if o.read_filter.as_ref().is_some_and(|p| p.needs_reference) {
            let r = self.reference_path(args)?.ok_or_else(|| anyhow!("filter: mismatches need a reference; pass `reference` or start the server with -f"))?;
            o.reference = Some(ggv_query::ReadReference::from_path(r));
        }
        Ok(o)
    }

    /// Whether `path` changed on disk since it was cached (and records it as it is now).
    fn changed(&self, path: &Path) -> bool {
        let now = stamp(path);
        let mut stamps = self.stamps.borrow_mut();
        let was = stamps.insert(path.to_owned(), now);
        was.is_some_and(|w| w != now)
    }

    fn bam(&mut self, path: &Path, reference: Option<&PathBuf>) -> Result<&mut BamFile> {
        let key = (path.to_owned(), reference.cloned());
        if self.changed(path) {
            self.bams.retain(|(p, _), _| p != path);
        }
        if !self.bams.contains_key(&key) {
            let b = BamFile::open_with_reference(path, reference).map_err(|e| named(e, path))?;
            self.bams.insert(key.clone(), b);
        }
        Ok(self.bams.get_mut(&key).expect("inserted"))
    }

    /// Runs `f` with the (cached) reference at `path`.
    fn with_fasta<T>(&self, path: &Path, f: impl FnOnce(&mut ReferenceFasta) -> T) -> Result<T> {
        let changed = self.changed(path);
        let mut map = self.fastas.borrow_mut();
        if changed {
            map.remove(path);
        }
        if !map.contains_key(path) {
            // Its errors name the file already.
            let fa = ReferenceFasta::open(path)?;
            map.insert(path.to_owned(), fa);
        }
        Ok(f(map.get_mut(path).expect("inserted")))
    }

    fn vcf(&mut self, path: &Path) -> Result<&mut VariantFile> {
        if self.changed(path) {
            self.vcfs.remove(path);
        }
        if !self.vcfs.contains_key(path) {
            let f = VariantFile::open(path).map_err(|e| named(e, path))?;
            self.vcfs.insert(path.to_owned(), f);
        }
        Ok(self.vcfs.get_mut(path).expect("inserted"))
    }

    fn annotation(&mut self, path: &Path) -> Result<Arc<Annotation>> {
        if self.changed(path) {
            self.annotations.remove(path);
        }
        if let Some(a) = self.annotations.get(path) {
            return Ok(a.clone());
        }
        let a = Arc::new(Annotation::open(path).map_err(|e| named(e, path))?);
        self.annotations.insert(path.to_owned(), a.clone());
        Ok(a)
    }

    fn vcf_paths(&self, args: &Value) -> Result<Vec<PathBuf>> {
        let v = strs(args, "vcfs")?;
        Ok(if v.is_empty() { self.defaults.vcfs.clone() } else { v })
    }

    /// The annotations to use: `genes` when given (an empty list means none), else the defaults.
    fn gene_paths(&self, args: &Value) -> Result<Vec<PathBuf>> {
        match args.get("genes") {
            None | Some(Value::Null) => Ok(self.defaults.genes.clone()),
            _ => strs(args, "genes"),
        }
    }

    /// Resolves a locus against the first alignment file's contigs, or the reference's; a gene
    /// name against the annotations (`genes`, the defaults, else the built-in RefSeq for the
    /// files' build when it's downloaded), as the tools' descriptions promise.
    fn resolve(&mut self, text: &str, args: &Value) -> Result<Region> {
        self.resolve_also(text, args, Vec::new())
    }

    /// [`resolve`](Self::resolve), also knowing `more` contigs (a VCF's or an annotation's, for
    /// the tools that read only those).
    fn resolve_also(&mut self, text: &str, args: &Value, more: Vec<ggv_model::Contig>) -> Result<Region> {
        let reference = self.reference_path(args)?;
        let mut contigs: Vec<ggv_model::Contig> = if let Ok(bams) = self.bam_paths(args) {
            self.bam(&bams[0], reference.as_ref())?.contigs().to_vec()
        } else if let Some(r) = &reference {
            self.with_fasta(r, |f| f.contigs().to_vec())?
        } else if !more.is_empty() {
            Vec::new()
        } else {
            bail!("no file to resolve {text} against: pass `bams` or `reference`")
        };
        // The alignment file's or reference's own lengths first.
        for c in more {
            if ggv_io::resolve_locus_str(&contigs, &c.name).is_err() {
                contigs.push(c);
            }
        }
        let as_locus = ggv_io::resolve_locus_str(&contigs, text);
        if as_locus.is_ok() || text.contains(':') {
            return Ok(as_locus?);
        }
        let mut annotations = self.gene_paths(args)?;
        // `genes: []` means no annotations, as documented: only a call that names none falls
        // back to the built-in ones.
        let explicit = matches!(args.get("genes"), Some(Value::Array(_) | Value::String(_)));
        if annotations.is_empty()
            && !explicit
            && let Some(build) = ggv_io::catalog::build_of(contigs.iter().map(|c| (c.name.as_str(), c.length)))
        {
            // The built-in genome's own genes (downloaded with it), then the catalogue's.
            annotations.extend(ggv_io::genomes::find(build).and_then(ggv_io::genomes::local).map(|l| l.genes).filter(|p| p.is_file()));
            annotations.extend(ggv_io::catalog::for_build(build).filter(|t| t.group == ggv_io::catalog::Group::Genes).filter_map(ggv_io::catalog::cached).flatten());
        }
        for path in &annotations {
            let Ok(ann) = self.annotation(path) else { continue };
            if let Some(g) = ann.find(text) {
                let at = format!("{}:{}-{}", g.contig, g.start + 1, g.end);
                if let Ok(r) = ggv_io::resolve_locus_str(&contigs, &at) {
                    return Ok(r);
                }
            }
        }
        if annotations.is_empty() {
            bail!("`{text}` is not a contig in these files; to use a gene name, pass `genes` (an annotation)")
        }
        bail!("`{text}` is not a contig in these files or a gene in {}", annotations.iter().map(|p| file_name(p)).collect::<Vec<_>>().join(", "))
    }

    /// The reference bases for `region`: `None` with no reference, an error when the one given
    /// can't be read there (naming the file and region).
    fn fetch(&self, args: &Value, region: &Region) -> Result<Option<Vec<u8>>> {
        let Some(r) = self.reference_path(args)? else { return Ok(None) };
        let seq = self.with_fasta(&r, |f| f.fetch(region))?.with_context(|| format!("reference {} at {region}", r.display()))?;
        Ok(Some(seq))
    }

    /// The tools, as listed to the client: attached, places and files are optional (they come
    /// from the window), and `view` / `navigate` are added.
    pub fn list(&self) -> Vec<Value> {
        let mut tools = list();
        if self.desktop.is_none() {
            return tools;
        }
        let place = json!({
            "pane": {"type": "integer", "minimum": 1, "description": "Which pane of the window (default: the active one)."},
            "window": {"type": "integer", "description": "Which GGV window, by process id (see `view`; default: the focused or last used one)."},
        });
        for t in &mut tools {
            if let Some(req) = t["inputSchema"]["required"].as_array_mut() {
                req.retain(|r| r != "region" && r != "locus");
            }
            if let (Some(p), Some(extra)) = (t["inputSchema"]["properties"].as_object_mut(), place.as_object()) {
                p.extend(extra.clone());
            }
            let d = t["description"].as_str().unwrap_or("").to_owned();
            t["description"] = json!(format!("{d} Attached to the GGV app: what you leave out (files, base, region, read filters) comes from what the person is looking at, and the answer's `desktop` block says what was used."));
        }
        tools.push(tool(
            "view",
            "What the person is looking at",
            "The GGV desktop window's state: for each pane its view (region), outlined base, selection, files, each file's read filter and display options, and which pane is active; plus the other GGV windows open. Call this first when attached.",
            schema(&[], &[place.clone()]),
        ));
        let mut nav = tool(
            "navigate",
            "Move the view",
            "Moves a pane of the GGV window to a locus (a single base is centred and outlined, as when typed), a region, or a gene; with `select`, selects the region instead (its statistics open in the Region window). Returns the pane's new state.",
            schema(&["target"], &[json!({"target": {"type": "string", "description": "e.g. chr17:43,092,919, chr17:43,090,000-43,100,000, or BRCA1."}, "select": {"type": "boolean", "description": "Select the region rather than just going there."}}), place.clone()]),
        );
        nav["annotations"] = json!({"readOnlyHint": false, "idempotentHint": true, "openWorldHint": false});
        tools.push(nav);
        let mut display = tool(
            "display",
            "Change how the window shows reads",
            "Changes the GGV window's read display, as its toolbar would (every pane follows): colour, grouping, sort, soft clips, squished rows, all transcripts, read filters for every file or one file's filter expression, and what colours variant sites. Give only what should change; a bad value changes nothing. Returns the pane's new state with `before`: the display settings it replaced, which display accepts as they are (pass them back to restore every setting, including sort, its position and per-file variant colours). A file name shared by two open files is refused: give the path.",
            schema(&[], &[json!({
                "color": {"type": "string", "anyOf": [{"enum": ["none", "strand", "haplotype", "insert size", "base modifications"]}, {"pattern": "^tag:[A-Za-z][A-Za-z0-9]$"}], "description": "Or `tag:XX`: colour by the value of any two-letter read tag, e.g. tag:HP or tag:RG."},
                "variant_color": {"type": "string", "description": "Colour variant sites by `af` (allele frequency) or an INFO key, e.g. CLNSIG (ClinVar significance: pathogenic red, uncertain dark blue, benign green, conflicting light blue, other grey) or SVTYPE. ClinVar files start on CLNSIG."},
                "variant_file": {"type": "string", "description": "With `variant_color`: only this VCF (path or file name); else every VCF declaring the key."},
                "group": {"type": "string", "anyOf": [{"enum": ["none", "strand", "haplotype", "read group", "alleles", "base"]}, {"pattern": "^tag:[A-Za-z][A-Za-z0-9]$"}], "description": "`base` groups by the base each read shows at `pos` (else the outlined base, else the view centre); `tag:XX` groups by any two-letter read tag's value, e.g. tag:HP."},
                "sort": {"anyOf": [{"type": "string", "enum": ["none", "base", "strand", "MAPQ", "start", "insert size"]}, {"type": "object", "properties": {"by": {"type": "string"}, "pos": {"type": "integer", "minimum": 1}}, "required": ["by", "pos"]}, {"type": "null"}], "description": "Sort at `pos` (else the outlined base, else the view centre); or {by, pos} as the state reports it; null for none."},
                "group_pos": {"type": "integer", "minimum": 1, "description": "With group `base`: the 1-based position (as the state reports it)."},
                "variant_colors": {"type": "object", "additionalProperties": {"type": "string"}, "description": "Per VCF (path or file name): `af` or an INFO key, as the state reports them."},
                "pos": {"type": "integer", "minimum": 1, "description": "1-based position for grouping by base or sorting."},
                "soft_clips": {"type": "boolean"},
                "squished": {"type": "boolean"},
                "junctions": {"type": "boolean", "description": "Splice junctions (RNA-seq) as arcs over each track's coverage, with read counts."},
                "link_splits": {"type": "boolean", "description": "Split reads' parts on one row, joined."},
                "read_window": {"type": "integer", "minimum": 10000, "maximum": 50000000, "description": "Show reads for views up to this many bases (the app's default is 1000000; wider views show coverage). Remembered by the app."},
                "all_transcripts": {"type": "boolean"},
                "filter": {"type": "object", "description": "For every file: min_mapq, hide_duplicates, hide_secondary, hide_supplementary, hide_qc_fail, expression (\"\" clears).", "properties": {
                    "min_mapq": {"type": "integer", "minimum": 0, "maximum": 255},
                    "hide_duplicates": {"type": "boolean"}, "hide_secondary": {"type": "boolean"}, "hide_supplementary": {"type": "boolean"}, "hide_qc_fail": {"type": "boolean"},
                    "expression": {"type": "string"}
                }},
                "file": {"type": "string", "description": "With `expression`: one file's own filter (path or file name)."},
                "expression": {"type": "string", "description": "With `file`: that file's filter expression (\"\" clears it)."}
            }), place.clone()]),
        );
        display["annotations"] = json!({"readOnlyHint": false, "idempotentHint": true, "openWorldHint": false});
        tools.push(display);
        let changes = |mut t: Value, idempotent: bool| {
            t["annotations"] = json!({"readOnlyHint": false, "idempotentHint": idempotent, "openWorldHint": false});
            t
        };
        tools.push(changes(
            tool(
                "select",
                "Select or clear a region",
                "Selects a region in a pane of the GGV window (its statistics open in the Region window, and region tools default to it), or clears the selection with `clear`.",
                schema(&[], &[json!({"region": {"type": "string", "description": "e.g. chr17:43,090,000-43,100,000 or a gene."}, "clear": {"type": "boolean"}}), place.clone()]),
            ),
            true,
        ));
        tools.push(changes(
            tool(
                "bookmarks",
                "Bookmarks",
                "The GGV window's bookmarks (shared by its panes, saved with the session): `list` (default); `add` one at `locus` (a base, region or gene; else the selection, the outlined base or the view) with a `note`; `remove` one by `index` (1-based) or `locus`; `clear` them; or `export` them as BED to `path` (ending .bed; a file open in GGV, or its index, is refused). Returns the bookmarks.",
                schema(&[], &[json!({"action": {"type": "string", "enum": ["list", "add", "remove", "clear", "export"]}, "locus": {"type": "string"}, "note": {"type": "string"}, "index": {"type": "integer", "minimum": 1}, "path": {"type": "string"}}), place.clone()]),
            ),
            false,
        ));
        tools.push(changes(
            tool(
                "session",
                "Save or load a session",
                "`save` writes the window's session (files, view, outlined base, display, filters, bookmarks, track layout) to `path` (a .json file, e.g. case.ggv.json); `load` opens one, replacing what is open.",
                schema(&["action", "path"], &[json!({"action": {"type": "string", "enum": ["save", "load"]}, "path": {"type": "string"}}), place.clone()]),
            ),
            false,
        ));
        tools.push(changes(
            tool(
                "tracks",
                "Tracks",
                "The pane's tracks by kind (alignments, variants, genes, signals) in order, and the reference; `close` one, or `move` it `up` or `down` among its kind. Name a track by its path, or its file name when that is unique.",
                schema(&[], &[json!({"action": {"type": "string", "enum": ["list", "close", "move"]}, "track": {"type": "string"}, "direction": {"type": "string", "enum": ["up", "down"]}}), place.clone()]),
            ),
            false,
        ));
        tools.push(changes(
            tool(
                "overlays",
                "Mark regions on the view",
                "Draws labelled, tinted bands across every track of the GGV window to show the person what you mean (e.g. \"deletion at the run start, insertion at the site\"): `add` one over `region` (a base, region or gene; else the selection, the outlined base or the view) with a short `label` and a `color`; `list` (default); `remove` one by `index`; or `clear` them. Saved with the session; the person can clear them from the toolbar.",
                schema(&[], &[json!({"action": {"type": "string", "enum": ["list", "add", "remove", "clear"]}, "region": {"type": "string"}, "label": {"type": "string", "description": "Up to 120 characters."}, "color": {"type": "string", "enum": ["yellow", "red", "green", "blue", "purple", "orange", "grey"], "description": "Default yellow."}, "index": {"type": "integer", "minimum": 1}}), place.clone()]),
            ),
            false,
        ));
        tools.push(tool(
            "events",
            "Wait for the person",
            "What the person does in the GGV window, as it happens: outlining a base (`outline`), selecting or clearing a region (`select`, `unselect`), adding a bookmark (`bookmark`) and triage decisions (`review`). Waits up to `wait_s` seconds (default 30, at most 300) for the next ones, or returns at once with `wait_s: 0`. Pass the returned `cursor` as `since` to pick up where you left off, and `kinds` to wait for some only. Use it to answer as the person points: wait for `outline`, then call `evidence` (which uses the outlined base).",
            schema(&[], &[json!({"wait_s": {"type": "number", "minimum": 0, "maximum": 300}, "since": {"type": "integer", "minimum": 0, "description": "The `cursor` of a previous answer: events after it. Default: only events from now on."}, "kinds": {"type": "array", "items": {"type": "string", "enum": ["outline", "select", "unselect", "bookmark", "review"]}}}), place.clone()]),
        ));
        tools.push(changes(
            tool(
                "review",
                "Variant review",
                "Variant triage in the GGV window, as its Variant review section does. With no `action` (or `status`): for each open VCF the sites marked real, artifact or unsure (with notes), their counts and the file they are saved in, and `current`: the site being reviewed. `start` reviews a VCF (`vcf`: path or file name; default the first) from the outlined site or the next one; `decide` marks the current site `real`, `artifact` or `unsure` (with an optional `note`), saves it and moves to the next site; `note` sets the current site's note; `next` / `prev` step without deciding; `stop` ends the review. Moving to a site takes a moment: when `current.site` is null and `finding_next` is true, call again. Decisions are saved as they are made, next to the VCF or, where that can't be written, in GGV's folder (`file` says where).",
                schema(&[], &[json!({
                    "action": {"type": "string", "enum": ["status", "start", "decide", "note", "next", "prev", "stop"]},
                    "vcf": {"type": "string"},
                    "decision": {"type": "string", "enum": ["real", "artifact", "unsure"]},
                    "note": {"type": "string"}
                }), place.clone()]),
            ),
            false,
        ));
        tools.push(changes(
            tool(
                "queue",
                "Review queue",
                "A list of places to step through and decide pass, fail or unsure, as the window's Review queue does. `open` a list file at `path` (BED, VCF, or one locus per line) or `paste` `text`; `decide` the current item (`decision`, optional `note`; saves and moves to the next undecided item); `note`; `next`, `prev`, `next_undecided`, or `go` to item `index` (1-based); `export` the decisions to `path` (.bed or .tsv; a file open in GGV, or its index, is refused); `stop`. Returns the queue: position, counts, the file decisions are saved in, every item with its decision, and how many lines of the list were not places (`lines_not_understood`). Reopening the same list resumes its decisions.",
                schema(&[], &[json!({
                    "action": {"type": "string", "enum": ["status", "open", "paste", "decide", "note", "next", "prev", "next_undecided", "go", "export", "stop"]},
                    "path": {"type": "string"}, "text": {"type": "string"},
                    "decision": {"type": "string", "enum": ["pass", "fail", "unsure"]},
                    "note": {"type": "string"}, "index": {"type": "integer", "minimum": 1}
                }), place.clone()]),
            ),
            false,
        ));
        tools.push(changes(
            tool(
                "panes",
                "Split view",
                "The GGV window's panes (side-by-side views of the same files): `list` (default); `split` adds a pane at `locus` (else where the target pane is; at most 4); `close` closes the target `pane`; `activate` makes it the one the person's keys and the toolbar act on; `link` sets whether panes scroll and zoom together (`linked`). Returns the pane count, each pane's region, the active pane and whether linked; a new pane reaches its locus as it loads: until then its `region` is null and `going_to` says where it was sent (`view` reports it `ready` once there).",
                schema(&[], &[json!({
                    "action": {"type": "string", "enum": ["list", "split", "close", "activate", "link"]},
                    "locus": {"type": "string"}, "linked": {"type": "boolean"}
                }), place.clone()]),
            ),
            false,
        ));
        tools.push(changes(
            tool(
                "action",
                "Run a window action",
                "Runs one of the GGV window's own actions by `name`, exactly as its menu, key or command palette would: history (`back`, `forward`), `zoom_in`/`zoom_out`, stepping (`next_variant`, `prev_variant`, `next_diff`, `prev_diff`), panels and toggles (`toggle_evidence`, `toggle_squish`, `toggle_dark`, `toggle_auto_reload`, `reload_all`…), `clear_selection` (clears the selection, outlined base, read and search marks), clipboard copies (`copy_link`, `copy_fasta`, `copy_sam`, `copy_image`), `add_bookmark`, `color:strand`, `group:haplotype`, `scale:log` and more. With no `name` (or `list`), every action's name, description and key. Actions that open a file dialog are refused: use `open`, `session`, `queue` or `bookmarks` with a path.",
                schema(&[], &[json!({"name": {"type": "string"}}), place.clone()]),
            ),
            false,
        ));
        tools.push(changes(
            tool(
                "rows",
                "Scroll read rows",
                "Scrolls an alignment track's stacked read rows in the GGV window, as the person does with the wheel or the track's scrollbar, or by dragging: to `row` (1-based, shown at the top), `by` a number of rows (negative is up), `to` the `top` or `bottom`, to the row holding the read named `read` (also outlined), or to the first row of `group` (when reads are grouped). `track`: 1-based, a path or a file name (default: every track for `to`, else the first). Without any, each track's row count, the rows on screen and its groups.",
                schema(&[], &[json!({
                    "track": {"type": ["integer", "string"]}, "row": {"type": "integer", "minimum": 1}, "by": {"type": "integer"},
                    "to": {"type": "string", "enum": ["top", "bottom"]}, "read": {"type": "string"}, "group": {"type": "string"}
                }), place.clone()]),
            ),
            false,
        ));
        tools.push(changes(
            tool(
                "search",
                "Find a sequence in the window",
                "The window's sequence search (primers, oligos; both strands, IUPAC codes): `start` one for `sequence` with up to `mismatches` (0-3) in `scope` `genome` (default), `contig` or `view`; `status` (default) gives progress and hits; `go` to hit `index` (selects it); `stop` cancels and closes it. Hits are marked in the window as the person would see them.",
                schema(&[], &[json!({
                    "action": {"type": "string", "enum": ["start", "status", "go", "stop"]},
                    "sequence": {"type": "string"}, "mismatches": {"type": "integer", "minimum": 0, "maximum": 3},
                    "scope": {"type": "string", "enum": ["genome", "contig", "view"]}, "index": {"type": "integer", "minimum": 1}
                }), place.clone()]),
            ),
            false,
        ));
        tools.push(changes(
            tool(
                "show_compare",
                "Show a comparison in the window",
                "The window's visible tumor/normal comparison of two alignment tracks, as the person would set it up: `open` (default) with tracks `a` and `b` (1-based, path or file name; default the first two), `threshold` (allele-fraction difference, 0-1) and `min_depth`; `next` / `prev` moves to the next difference (scanning the files past what is loaded); `close`. For the numbers without changing the window, use the `compare` tool.",
                schema(&[], &[json!({
                    "action": {"type": "string", "enum": ["open", "next", "prev", "close"]},
                    "a": {"type": ["integer", "string"]}, "b": {"type": ["integer", "string"]},
                    "threshold": {"type": "number", "minimum": 0.01, "maximum": 1}, "min_depth": {"type": "integer", "minimum": 1, "maximum": 1000}
                }), place]),
            ),
            false,
        ));
        tools
    }

    /// Sends a window command and returns what it gave back, with the window's process id.
    fn window_command(&mut self, window: Option<u32>, pane: Option<u64>, method: &str, mut params: Value) -> Result<Output> {
        let desktop = self.desktop.as_mut().expect("attached");
        let w = desktop.window(window).map_err(|e| anyhow!(e))?;
        if let Some(p) = pane {
            pane_check(&w.state, p)?;
            params["pane"] = json!(p);
        }
        let state = desktop.command(w.pid, method, params).map_err(|e| anyhow!(e))?;
        let mut out = match state.get("result") {
            Some(r) => r.clone(),
            None => {
                let after = super::desktop::Window { pid: w.pid, state };
                after.pane(pane).map_err(|e| anyhow!(e))?.clone()
            }
        };
        out["window"] = json!(w.pid);
        Ok(Output::Json(out))
    }

    pub fn call(&mut self, name: &str, args: &Value) -> Result<Output> {
        if self.desktop.is_none() {
            return self.call_with(name, args);
        }
        let window = args.get("window").and_then(Value::as_u64).map(|w| w as u32);
        let pane = args.get("pane").and_then(Value::as_u64);
        let mut args = args.clone();
        if let Some(m) = args.as_object_mut() {
            m.remove("window");
            m.remove("pane");
        }
        let desktop = self.desktop.as_mut().expect("attached");
        match name {
            "view" => {
                let others: Vec<Value> = desktop
                    .windows()
                    .into_iter()
                    .map(|w| json!({"window": w.pid, "focused": w.state["focused"], "idle_ms": w.state["idle_ms"], "panes": w.state["panes"].as_array().map(|p| p.iter().map(|x| x["view"]["region"].clone()).collect::<Vec<_>>())}))
                    .collect();
                let w = desktop.window(window).map_err(|e| anyhow!(e))?;
                let mut state = w.state;
                if let Some(p) = pane {
                    pane_check(&state, p)?;
                }
                state["window"] = json!(w.pid);
                state["windows"] = json!(others);
                return Ok(Output::Json(state));
            }
            "display" => {
                let w = desktop.window(window).map_err(|e| anyhow!(e))?;
                let mut params = args.clone();
                if let Some(p) = pane {
                    params["pane"] = json!(p);
                }
                let state = desktop.command(w.pid, "display", params).map_err(|e| anyhow!(e))?;
                let before = state["before"].clone();
                let after = super::desktop::Window { pid: w.pid, state };
                let mut p = pane_reply(after.pane(pane).map_err(|e| anyhow!(e))?.clone());
                p["window"] = json!(w.pid);
                p["before"] = before;
                return Ok(Output::Json(p));
            }
            "select" => {
                if args.get("clear").and_then(Value::as_bool) == Some(true) {
                    return self.window_command(window, pane, "unselect", json!({}));
                }
                let region = str_arg(&args, "region").ok_or_else(|| anyhow!("give a `region` to select, or `clear`"))?.to_owned();
                return self.window_command(window, pane, "select", json!({"target": region}));
            }
            "bookmarks" => return self.window_command(window, pane, "bookmarks", args),
            "overlays" => return self.window_command(window, pane, "overlays", args),
            "events" => {
                let out = desktop.events_tool(window, &args, &super::desktop::Cancel::default()).map_err(|e| anyhow!(e))?;
                return Ok(Output::Json(out));
            }
            "tracks" => return self.window_command(window, pane, "tracks", args),
            "review" => return self.window_command(window, pane, "review", args),
            "queue" => return self.window_command(window, pane, "queue", args),
            "panes" => return self.window_command(window, pane, "panes", args),
            "action" => return self.window_command(window, pane, "action", args),
            "search" => return self.window_command(window, pane, "search", args),
            "show_compare" => return self.window_command(window, pane, "compare", args),
            "rows" => return self.window_command(window, pane, "rows", args),
            "session" => {
                let path = str_arg(&args, "path").ok_or_else(|| anyhow!("give a `path`"))?.to_owned();
                return match str_arg(&args, "action") {
                    Some("save") => self.window_command(window, pane, "save_session", json!({"path": path})),
                    Some("load") => {
                        if !Path::new(&path).is_file() {
                            bail!("{path}: no such file");
                        }
                        // Only a session: `open` would add any other file to the workspace.
                        if !path.to_ascii_lowercase().ends_with(".json") {
                            bail!("{path}: not a session file (a .json, e.g. case.ggv.json); use `open` for data files");
                        }
                        self.window_command(window, pane, "open", json!({"files": [path]}))
                    }
                    _ => bail!("give `action`: save or load"),
                };
            }
            "navigate" => {
                let target = str_arg(&args, "target").ok_or_else(|| anyhow!("give a `target`: a locus, region or gene"))?;
                let w = desktop.window(window).map_err(|e| anyhow!(e))?;
                let method = if args.get("select").and_then(Value::as_bool) == Some(true) { "select" } else { "goto" };
                let mut params = json!({"target": target});
                if let Some(p) = pane {
                    params["pane"] = json!(p);
                }
                let state = desktop.command(w.pid, method, params).map_err(|e| anyhow!(e))?;
                let after = super::desktop::Window { pid: w.pid, state };
                let mut p = pane_reply(after.pane(None).map_err(|e| anyhow!(e))?.clone());
                p["window"] = json!(w.pid);
                return Ok(Output::Json(p));
            }
            _ => {}
        }
        let w = match desktop.window(window) {
            Ok(w) => w,
            // No window: the server's own defaults, if it has any.
            Err(e) if window.is_none() && !self.defaults.bams.is_empty() => {
                let out = self.call_with(name, &args)?;
                return Ok(annotate(out, json!({"attached": false, "note": format!("{e}; used the server's files")})));
            }
            Err(e) => bail!(e),
        };
        let p = w.pane(pane).map_err(|e| anyhow!(e))?.clone();
        if name == "open" {
            let files = strs(&args, "files")?;
            if files.is_empty() {
                bail!("give `files`");
            }
            let mut params = json!({"files": files});
            if args.get("replace").and_then(Value::as_bool) == Some(true) {
                params["replace"] = json!(true);
            }
            if let Some(n) = pane {
                params["pane"] = json!(n);
            }
            let state = self.desktop.as_mut().expect("attached").command(w.pid, "open", params).map_err(|e| anyhow!(e))?;
            let after = super::desktop::Window { pid: w.pid, state };
            let mut p = pane_reply(after.pane(None).map_err(|e| anyhow!(e))?.clone());
            p["window"] = json!(w.pid);
            p["note"] = json!("opened in the GGV window (files load in the background; `view` shows when they are ready)");
            return Ok(Output::Json(p));
        }
        let (filled, applied, notes) = super::desktop::fill(name, &args, &p);
        let server_defaults = std::mem::replace(&mut self.defaults, super::desktop::defaults_of(&p));
        let out = self.call_with(name, &filled);
        self.defaults = server_defaults;
        let mut info = json!({"window": w.pid, "pane": p["pane"], "applied": applied});
        if !notes.is_empty() {
            info["notes"] = json!(notes);
        }
        Ok(annotate(out?, info))
    }

    fn call_with(&mut self, name: &str, args: &Value) -> Result<Output> {
        let v = match name {
            "info" => self.info(args)?,
            "open" => self.open(args)?,
            "evidence" => self.evidence(args)?,
            "reads" => self.reads(args)?,
            "summarize" => self.summarize(args)?,
            "coverage" => self.coverage(args)?,
            "variants" => self.variants(args)?,
            "genes" => self.genes(args)?,
            "consequence" => self.consequence(args)?,
            "phase" => self.phase(args)?,
            "junctions" => self.junctions(args)?,
            "svs" => self.svs(args)?,
            "repeats" => self.repeats(args)?,
            "compare" => self.compare(args)?,
            "signal" => self.signal(args)?,
            "render" => {
                let out = self.render(args)?;
                return Ok(match (out, self.file_warnings(name, args)) {
                    (Output::Image { png, note }, w) if !w.is_empty() => Output::Image { png, note: format!("{note}
warning: {}", w.join("; ")) },
                    (out, _) => out,
                });
            }
            _ => bail!("unknown tool {name}"),
        };
        let warnings = self.file_warnings(name, args);
        Ok(Output::Json(match v {
            Value::Object(mut m) if !warnings.is_empty() => {
                m.insert("file_warnings".into(), json!(warnings));
                Value::Object(m)
            }
            v => v,
        }))
    }

    /// Warnings about the files this call used: an alignment file that looks cut short, and
    /// files that aren't the same assembly (an hg19 BAM against an hg38 reference or genes).
    fn file_warnings(&mut self, name: &str, args: &Value) -> Vec<String> {
        let paths = self.used_bams(name, args);
        let mut w: Vec<String> = self.bams.iter().filter(|((p, _), _)| paths.contains(p)).filter_map(|(_, b)| b.warning()).collect();
        w.sort();
        w.dedup();
        w.extend(self.assembly_warning(&paths, args));
        w
    }

    /// The alignment files a call reads: `compare`'s pair, else `bams`/`bam` or the defaults.
    fn used_bams(&self, name: &str, args: &Value) -> Vec<PathBuf> {
        let defaults = self.bam_paths(args).unwrap_or_default();
        if name != "compare" {
            return defaults;
        }
        let a = str_arg(args, "a").map(PathBuf::from).or_else(|| defaults.first().cloned());
        let b = str_arg(args, "b").map(PathBuf::from).or_else(|| defaults.get(1).cloned());
        a.into_iter().chain(b).collect()
    }

    /// Where the reference, the alignment files read, the VCFs and the annotations disagree on
    /// the assembly (the same check as the desktop's), as one warning.
    fn assembly_warning(&mut self, bams: &[PathBuf], args: &Value) -> Option<String> {
        use ggv_io::assembly::{Identity, display_name};
        let lengths = |c: &[ggv_model::Contig]| c.iter().map(|c| (c.name.clone(), c.length)).collect::<Vec<_>>();
        let mut ids = Vec::new();
        if let Some(r) = self.reference_path(args).ok().flatten()
            && let Ok(l) = self.with_fasta(&r, |f| lengths(f.contigs()))
        {
            ids.push(Identity::new("reference", r.clone(), Some(display_name(&r)), l));
        }
        for p in bams {
            // Only files this call opened (one that failed has its own error).
            if let Some((_, b)) = self.bams.iter().find(|((q, _), _)| q == p) {
                ids.push(Identity::new("alignments", p.clone(), None, lengths(b.contigs())));
            }
        }
        for p in self.vcf_paths(args).unwrap_or_default() {
            if let Ok(v) = self.vcf(&p) {
                let declared = v.declared_assembly().or_else(|| Some(display_name(&p)));
                ids.push(Identity::new("variants", p.clone(), declared, lengths(&v.header_contigs())));
            }
        }
        for p in self.gene_paths(args).unwrap_or_default() {
            ids.push(Identity::new("annotation", p.clone(), Some(display_name(&p)), Vec::new()));
        }
        let conflicts = ggv_io::assembly::conflicts(&ids);
        (!conflicts.is_empty()).then(|| {
            format!("assembly mismatch: {} Until the files match, positions, reference bases, variant calls and consequences may be wrong.", conflicts.join(" "))
        })
    }

    // ---- info / open ----

    fn describe(&mut self, path: &Path) -> Value {
        let name = path.to_string_lossy().to_ascii_lowercase();
        let kind = crate::kind_of(&name);
        let result: Result<Value> = (|| {
            Ok(match kind {
                "alignments" => {
                    let reference = self.defaults.reference.clone();
                    let b = self.bam(path, reference.as_ref())?;
                    let info = b.info();
                    let contigs = b.contigs();
                    let main: Vec<Value> = contigs.iter().filter(|c| c.length >= 1_000_000 && !c.name.contains('_')).take(40).map(|c| json!([c.name, c.length])).collect();
                    json!({
                        "path": path, "kind": kind,
                        "samples": info.read_groups.iter().filter_map(|r| r.sample.clone()).collect::<std::collections::BTreeSet<_>>(),
                        "read_groups": info.read_groups.iter().map(|r| json!({"id": r.id, "sample": r.sample, "platform": r.platform})).collect::<Vec<_>>(),
                        "contigs": contigs.len(), "main_contigs": main, "index": info.index,
                    })
                }
                "variants" => {
                    let v = self.vcf(path)?;
                    json!({"path": path, "kind": kind, "samples": v.samples()})
                }
                "genes" => {
                    let a = self.annotation(path)?;
                    json!({"path": path, "kind": kind, "genes": a.gene_count(), "transcripts": a.len(), "contigs": a.contig_names().len()})
                }
                "reference" => self.with_fasta(path, |f| {
                    json!({"path": path, "kind": kind, "contigs": f.contigs().len(), "main_contigs": f.contigs().iter().filter(|c| c.length >= 1_000_000 && !c.name.contains('_')).take(40).map(|c| json!([c.name, c.length])).collect::<Vec<_>>()})
                })?,
                "signal" => {
                    let s = ggv_io::signal::Signal::open(path)?;
                    let contigs = s.contigs();
                    let main: Vec<Value> = contigs.iter().filter(|c| !c.name.contains('_')).take(40).map(|c| json!([c.name, c.length])).collect();
                    // One row for bigWig/bedGraph; SEG has one per sample.
                    let rows = s.rows();
                    let samples = if rows.len() == 1 && rows[0].is_empty() { Value::Null } else { json!(rows) };
                    json!({"path": path, "kind": kind, "samples": samples, "contigs": contigs.len(), "main_contigs": main})
                }
                _ => json!({"path": path, "kind": "unknown"}),
            })
        })();
        result.unwrap_or_else(|e| json!({"path": path, "error": format!("{e:#}")}))
    }

    fn info(&mut self, args: &Value) -> Result<Value> {
        let mut files: Vec<PathBuf> = self.defaults.bams.clone();
        files.extend(self.defaults.vcfs.clone());
        files.extend(self.defaults.genes.clone());
        files.extend(self.defaults.reference.clone());
        files.extend(self.defaults.signals.clone());
        files.extend(strs(args, "files")?);
        let described: Vec<Value> = files.iter().map(|p| self.describe(p)).collect();
        Ok(json!({
            "defaults": {"bams": self.defaults.bams, "reference": self.defaults.reference, "vcfs": self.defaults.vcfs, "genes": self.defaults.genes, "signals": self.defaults.signals},
            "files": described,
        }))
    }

    fn open(&mut self, args: &Value) -> Result<Value> {
        let files = strs(args, "files")?;
        if files.is_empty() {
            bail!("give `files`");
        }
        // All or nothing: check every file (it exists, is a known kind, and opens) before
        // changing the defaults, so a failed call leaves them as they were.
        let mut next = if bool_arg(args, "replace")? == Some(true) { Defaults::default() } else { self.defaults.clone() };
        for p in &files {
            let kind = crate::kind_of(&p.to_string_lossy().to_ascii_lowercase());
            if kind == "reference" {
                next.reference = Some(p.clone());
            }
        }
        for p in files {
            if !p.exists() {
                bail!("{} does not exist; the defaults are unchanged", p.display());
            }
            // Idempotent, as advertised: a file already open (by any spelling of its path) isn't
            // added again.
            let add = |list: &mut Vec<PathBuf>| {
                if !list.iter().any(|q| same_file(q, &p)) {
                    list.push(p.clone());
                }
            };
            let opened = match crate::kind_of(&p.to_string_lossy().to_ascii_lowercase()) {
                "alignments" => self.bam(&p, next.reference.as_ref()).map(|_| add(&mut next.bams)),
                "variants" => self.vcf(&p).map(|_| add(&mut next.vcfs)),
                "genes" => self.annotation(&p).map(|_| add(&mut next.genes)),
                "signal" => ggv_io::signal::Signal::open(&p).map(|_| add(&mut next.signals)).map_err(anyhow::Error::from),
                "reference" => self.with_fasta(&p, |_| ()),
                _ => Err(anyhow!("not a BAM/CRAM, VCF/BCF, GTF/GFF3/BED, FASTA or bigWig/bedGraph/SEG")),
            };
            if let Err(e) = opened {
                let e = format!("{e:#}");
                let named = e.contains(&*p.to_string_lossy());
                bail!("{}{e}; the defaults are unchanged", if named { String::new() } else { format!("{}: ", p.display()) });
            }
        }
        self.defaults = next;
        self.info(&json!({}))
    }

    // ---- evidence ----

    fn evidence(&mut self, args: &Value) -> Result<Value> {
        let mut loci: Vec<String> = strs(args, "loci")?.into_iter().map(|p| p.to_string_lossy().into_owned()).collect();
        loci.extend(str_arg(args, "locus").map(String::from));
        if loci.is_empty() {
            bail!("give a `locus` (or `loci`)");
        }
        if loci.len() > MAX_LOCI {
            bail!("at most {MAX_LOCI} loci per call");
        }
        let opts = self.read_filters(args)?;
        let bams = self.bam_paths(args)?;
        let pre = if loci.len() > 1 { self.prefetch(&loci, args, &bams, &opts) } else { HashMap::new() };
        let mut out = Vec::new();
        for l in &loci {
            out.push(self.evidence_at(l, args, &bams, &opts, &pre).unwrap_or_else(|e| json!({"locus": l, "error": format!("{e:#}")})));
        }
        Ok(json!({
            "query": {"bams": bams, "filters": filters_json(&opts)},
            "results": out,
        }))
    }

    /// For several bases at once: per file, the bases within 50 kb of each other are read in one
    /// pass over their span (counts at each base, and the reads), instead of one pass per base.
    /// Counts at a base don't depend on the region the pileup ran over, so the answers are the
    /// same; only the decoding is shared.
    fn prefetch(&mut self, loci: &[String], args: &Value, bams: &[PathBuf], opts: &PileupOptions) -> Prefetch {
        const SPAN: u64 = 50_000;
        let mut by_contig: std::collections::BTreeMap<String, Vec<u64>> = Default::default();
        for l in loci {
            if let Ok(r) = self.resolve(l, args)
                && r.len() == 1
            {
                by_contig.entry(r.contig).or_default().push(r.start);
            }
        }
        // A reference that can't be read is reported by each base's own call.
        let Ok(reference) = self.reference_path(args) else { return Prefetch::new() };
        let mut out = Prefetch::new();
        for (contig, mut positions) in by_contig {
            positions.sort_unstable();
            positions.dedup();
            let mut groups: Vec<Vec<u64>> = Vec::new();
            for p in positions {
                match groups.last_mut() {
                    Some(g) if p - g[0] <= SPAN => g.push(p),
                    _ => groups.push(vec![p]),
                }
            }
            for group in groups.into_iter().filter(|g| g.len() > 1) {
                let span = Region { contig: contig.clone(), start: group[0], end: group[group.len() - 1] + 1 };
                let seq = self.fetch(args, &span).ok().flatten();
                for (bi, path) in bams.iter().enumerate() {
                    let Ok(bam) = self.bam(path, reference.as_ref()) else { continue };
                    let Ok(local) = bam.resolve_str(&format!("{}:{}-{}", span.contig, span.start + 1, span.end)) else { continue };
                    let wanted: std::collections::HashSet<u64> = group.iter().copied().collect();
                    let Ok((mut stats, reads)) = counts_and_reads(bam, &local, seq.as_deref(), opts, |p| wanted.contains(&p)) else { continue };
                    let reads = Arc::new(reads);
                    for &p in &group {
                        out.insert((bi, span.contig.clone(), p), (stats.remove(&p), reads.clone()));
                    }
                }
            }
        }
        out
    }

    fn evidence_at(&mut self, locus: &str, args: &Value, bams: &[PathBuf], opts: &PileupOptions, pre: &Prefetch) -> Result<Value> {
        let region = self.resolve(locus, args)?;
        if region.len() != 1 {
            bail!("{locus} is not a single base (use `summarize` for regions)");
        }
        let pos = region.start;
        let reference = self.reference_path(args)?;
        let ref_base = self.fetch(args, &region)?.and_then(|s| s.first().copied());
        let mut tracks = Vec::new();
        let mut changes: Vec<(u64, String, String, String)> = Vec::new();
        for (bi, path) in bams.iter().enumerate() {
            let bam = self.bam(path, reference.as_ref())?;
            let local = match bam.resolve_str(&format!("{}:{}", region.contig, pos + 1)) {
                Ok(r) => r,
                Err(_) => {
                    tracks.push(json!({"file": file_name(path), "path": path, "error": format!("{} is not in this file", region.contig)}));
                    continue;
                }
            };
            // Exact counts (samtools-compatible), and the reads for read-level detail per allele
            // with the same filters: shared with nearby bases of the same call when it named
            // several (see `prefetch`), else read here.
            let fetched;
            let (stats, reads): (Option<PositionStats>, &[AlignedRead]) = match pre.get(&(bi, region.contig.clone(), pos)) {
                Some((s, r)) => (s.clone(), r.as_slice()),
                None => {
                    let refslice = ref_base.map(|b| vec![b]);
                    let (mut stats, reads) = counts_and_reads(bam, &local, refslice.as_deref(), opts, |_| true)?;
                    fetched = reads;
                    (stats.remove(&pos), fetched.as_slice())
                }
            };
            let mut by: std::collections::BTreeMap<Allele, Vec<(&AlignedRead, u32)>> = Default::default();
            for r in reads {
                if let Some((a, d)) = read_at(r, pos, ref_base) {
                    by.entry(a).or_default().push((r, d));
                }
            }
            let Some(s) = stats else {
                tracks.push(json!({"file": file_name(path), "path": path, "depth": 0}));
                continue;
            };
            let depth = f64::from(s.depth.max(1));
            let count = |name: &str| -> (u32, u32, u32) {
                let b = &s.bases;
                let c = match name {
                    "A" => b.a,
                    "C" => b.c,
                    "G" => b.g,
                    "T" => b.t,
                    "N" => b.n,
                    _ => s.deletions,
                };
                (c.total, c.fwd, c.rev)
            };
            let mut alleles = Vec::new();
            for name in ["A", "C", "G", "T", "N", "DEL"] {
                let (n, f, r) = count(name);
                if n == 0 {
                    continue;
                }
                let key = if name == "DEL" { Allele::Deletion } else { Allele::Base(name.as_bytes()[0]) };
                let mut a = json!({"allele": name, "reads": n, "fraction": round(f64::from(n) / depth, 3), "forward": f, "reverse": r});
                if s.ref_base.is_some_and(|rb| rb.to_string() == name) {
                    a["reference"] = json!(true);
                }
                if let Some(g) = by.get(&key) {
                    let gj = group_json(g, true);
                    for k in ["mapq_median", "end_distance_median", "soft_clipped_share"] {
                        a[k] = gj[k].clone();
                    }
                }
                alleles.push(a);
            }
            // Warnings a person would see by eye.
            let mut warnings = Vec::new();
            if s.depth < 8 {
                warnings.push(format!("low depth ({})", s.depth));
            }
            if s.mapping_quality_zero as f64 > 0.1 * depth {
                warnings.push(format!("{} of {} reads have MAPQ 0 (repetitive or ambiguous)", s.mapping_quality_zero, s.depth));
            }
            if let Some(alt) = &s.top_alt {
                if alt.count.total >= 4 && (alt.count.fwd == 0 || alt.count.rev == 0) {
                    // One-strand support only means something when the other reads here come
                    // from the strand the alt lacks: at the edge of a captured region, or with
                    // amplicons, every read can be on one strand, alt or not.
                    let b = &s.bases;
                    let (all_f, all_r) = [b.a, b.c, b.g, b.t, b.n, s.deletions].iter().fold((0u32, 0u32), |(f, r), c| (f + c.fwd, r + c.rev));
                    let (other_f, other_r) = (all_f.saturating_sub(alt.count.fwd), all_r.saturating_sub(alt.count.rev));
                    let lacking = if alt.count.fwd == 0 { other_f } else { other_r };
                    let other = other_f + other_r;
                    if other >= 4 && f64::from(lacking) >= 0.2 * f64::from(other) {
                        warnings.push(format!(
                            "alt {} on one strand only ({} forward / {} reverse) while the other reads here are on both ({other_f} / {other_r})",
                            alt.allele, alt.count.fwd, alt.count.rev
                        ));
                    } else if lacking == 0 {
                        warnings.push(format!(
                            "strand bias can't be judged: every read here is on the {} strand",
                            if alt.count.fwd == 0 { "reverse" } else { "forward" }
                        ));
                    }
                }
                let key = if alt.allele == "DEL" { Allele::Deletion } else { Allele::Base(alt.allele.as_bytes()[0]) };
                let reference_key = s.ref_base.map(|b| Allele::Base(b as u8));
                let med = |k: Option<Allele>, f: fn(&(&AlignedRead, u32)) -> f64| k.and_then(|k| by.get(&k)).and_then(|g| median(g.iter().map(f).collect()));
                let (alt_end, ref_end) = (med(Some(key), |x| f64::from(x.1)), med(reference_key, |x| f64::from(x.1)));
                if let (Some(a), Some(r)) = (alt_end, ref_end)
                    && alt.count.total >= 4
                    && a <= 10.0
                    && r >= 25.0
                {
                    warnings.push(format!("alt {} sits near read ends (median {a} bases from the end vs {r} for the reference)", alt.allele));
                }
                let (alt_q, ref_q) = (med(Some(key), |x| f64::from(x.0.mapq)), med(reference_key, |x| f64::from(x.0.mapq)));
                if let (Some(a), Some(r)) = (alt_q, ref_q)
                    && alt.count.total >= 4
                    && a + 20.0 <= r
                {
                    warnings.push(format!("alt {} reads map worse (median MAPQ {a} vs {r})", alt.allele));
                }
                if let Some(r) = s.ref_base.filter(|c| "ACGT".contains(*c)) {
                    if alt.allele.len() == 1 {
                        changes.push((pos, r.to_string(), alt.allele.to_owned(), format!("{} reads", file_name(path))));
                    }
                }
            }
            if (s.soft_clipped_left + s.soft_clipped_right) as f64 > 0.2 * depth {
                warnings.push(format!("{} reads are soft-clipped right here (a breakpoint or misalignment?)", s.soft_clipped_left + s.soft_clipped_right));
            }
            tracks.push(json!({
                "file": file_name(path),
                "path": path,
                "depth": s.depth,
                "alleles": alleles,
                "insertions_after": s.insertions_after.total,
                "mapq_mean": round(s.mean_mapping_quality, 1),
                "base_quality_mean": round(s.mean_base_quality, 1),
                "low_base_quality": s.low_base_quality,
                "warnings": warnings,
            }));
        }

        // VCF records at the base.
        let mut records = Vec::new();
        let mut not_annotated = Vec::new();
        let mut warnings = Vec::new();
        for path in self.vcf_paths(args)? {
            let samples = self.vcf(&path)?.samples().to_vec();
            for v in self.vcf(&path)?.query(&region, false)?.into_iter().filter(|v| v.start <= pos && pos < v.end) {
                // A REF the reference disagrees with makes every consequence of it suspect.
                let span = Region { contig: region.contig.clone(), start: v.start, end: v.start + v.reference.len() as u64 };
                if let Some(found) = self.fetch(args, &span)?
                    && !found.eq_ignore_ascii_case(v.reference.as_bytes())
                {
                    warnings.push(format!(
                        "{}: REF {} at {}:{} does not match the reference ({}); its consequences assume the VCF's REF",
                        file_name(&path),
                        allele_text(&v.reference),
                        region.contig,
                        v.start + 1,
                        allele_text(&String::from_utf8_lossy(&found))
                    ));
                }
                for a in v.alts.iter().filter(|a| *a != "*" && *a != "." && *a != "<*>" && *a != "<NON_REF>") {
                    match unannotatable(&v.reference, a) {
                        None => changes.push((v.start, v.reference.clone(), a.clone(), file_name(&path))),
                        Some(what) => not_annotated.push(json!({"change": format!("{}:{} {}>{}", region.contig, v.start + 1, allele_text(&v.reference), allele_text(a)), "from": file_name(&path),
                            "reason": format!("{what} doesn't say which bases change, so its effect on a transcript isn't computed; the genes here are listed")})),
                    }
                }
                records.push(variant_json(&v, &samples, &path));
            }
        }

        // Consequences of each change, for the genes there.
        let mut consequences = Vec::new();
        let mut seen = std::collections::HashSet::new();
        let fetch_args = args.clone();
        let all = all_transcripts(args)?;
        let mut hidden = 0;
        for gpath in self.gene_paths(args)? {
            let ann = self.annotation(&gpath)?;
            let before = consequences.len();
            for (p, r, a, source) in &changes {
                if !seen.insert((gpath.clone(), *p, r.clone(), a.clone())) {
                    continue;
                }
                let (shown, more) = pick_transcripts(transcripts_over(&ann, &region.contig, *p, *p + r.len() as u64), all);
                hidden += more;
                for t in shown {
                    let c = coding::consequence(t, *p, r, a, |q| self.base_at(&fetch_args, &region.contig, q));
                    if let Some(c) = c {
                        let mut v = json!({"gene": c.gene, "transcript": c.transcript, "hgvs_c": c.hgvs_c, "hgvs_p": c.hgvs_p, "effect": c.effect, "where": c.location.feature, "change": format!("{}:{} {}>{}", region.contig, p + 1, allele_text(r), allele_text(a)), "from": source});
                        if !c.warnings.is_empty() {
                            v["warnings"] = json!(c.warnings);
                        }
                        consequences.push(v);
                    }
                }
            }
            // Nothing annotated (no change, or none that could be): the genes at the base.
            if consequences.len() == before {
                for t in transcripts_at(&ann, &region.contig, pos).into_iter().take(3) {
                    if let Some(l) = coding::locate(t, pos) {
                        consequences.push(json!({"gene": t.gene_name, "transcript": ggv_annot::without_copy_suffix(&t.id), "where": l.feature, "region": l.kind.name(), "c": format!("{}{}", if t.cds.is_some() { "c." } else { "n." }, l.c)}));
                    }
                }
            }
        }

        let mut out = json!({"locus": format!("{}:{}", region.contig, pos + 1), "ref": ref_base.map(|b| (b as char).to_string()), "tracks": tracks});
        if !records.is_empty() {
            out["variants"] = json!(records);
        }
        if !warnings.is_empty() {
            out["warnings"] = json!(warnings);
        }
        if !consequences.is_empty() {
            out["genes"] = json!(consequences);
        }
        if hidden > 0 {
            out["other_transcripts"] = json!(hidden);
        }
        if !not_annotated.is_empty() {
            out["not_annotated"] = json!(not_annotated);
        }
        Ok(out)
    }

    fn base_at(&self, args: &Value, contig: &str, p: u64) -> Option<u8> {
        // Lenient: callers have already fetched from the same reference (so a bad one has been
        // reported), and transcripts may reach past a contig's end.
        self.fetch(args, &Region { contig: contig.to_owned(), start: p, end: p + 1 }).ok().flatten().and_then(|s| s.first().copied())
    }

    // ---- reads ----

    fn reads(&mut self, args: &Value) -> Result<Value> {
        let opts = self.read_filters(args)?;
        if str_arg(args, "locus").is_some() && str_arg(args, "region").is_some() {
            bail!("give either `locus` (one base) or `region`, not both");
        }
        let at = str_arg(args, "locus").or_else(|| str_arg(args, "region")).ok_or_else(|| anyhow!("give a `locus` or a `region`"))?.to_owned();
        let region = self.resolve(&at, args)?;
        if region.len() > READ_WINDOW {
            bail!("{at} is {} bp; reads are listed for up to {READ_WINDOW} bp", region.len());
        }
        let single = region.len() == 1;
        let path = self.bam_paths(args)?.remove(0);
        let reference = self.reference_path(args)?;
        let seq = self.fetch(args, &region)?;
        let ref_base = if single { seq.as_ref().and_then(|s| s.first().copied()) } else { None };
        let group_by = str_arg(args, "group_by").unwrap_or(if single { "allele" } else { "none" }).to_owned();
        if group_by == "allele" && !single {
            bail!("grouping by allele needs one base: pass `locus` (e.g. {}:{}) instead of a {}-bp region", region.contig, region.start + 1, region.len());
        }
        let list = uint_or(args, "list", 0)?.min(50) as usize;
        let bam = self.bam(&path, reference.as_ref())?;
        let local = bam.resolve_str(&format!("{}:{}-{}", region.contig, region.start + 1, region.end))?;
        let view = load_view(bam, &local, seq.clone(), String::new(), &ViewOptions { filters: opts.clone(), max_rows: 0, ..Default::default() })?;
        let reads: Vec<&AlignedRead> = view.reads.as_ref().map(|(r, _)| r.iter().collect()).unwrap_or_default();
        let phasing = (group_by == "phase").then(|| {
            let all: Vec<AlignedRead> = reads.iter().map(|r| (*r).clone()).collect();
            let p = ggv_layout::phase::phase(&all, &(0..all.len()).collect::<Vec<_>>());
            (p, all.len())
        });
        let mut groups: std::collections::BTreeMap<String, Vec<(&AlignedRead, u32, Option<(Allele, u32)>)>> = Default::default();
        for (i, r) in reads.iter().enumerate() {
            let here = if single { read_at(r, region.start, ref_base) } else { None };
            if single && here.is_none() {
                continue;
            }
            let key = match group_by.as_str() {
                "allele" => here.map_or("not covering".into(), |(a, _)| {
                    let n = a.name();
                    if Some(n.as_bytes()[0]) == ref_base && n.len() == 1 { format!("{n} (ref)") } else { n }
                }),
                "strand" => if r.reverse { "-".into() } else { "+".into() },
                "haplotype" => r.haplotype.map_or("no HP".into(), |h| format!("HP {h}")),
                "phase" => match phasing.as_ref().and_then(|(p, _)| p.hap[i]) {
                    Some(true) => "phase A".into(),
                    Some(false) => "phase B".into(),
                    None => "unphased".into(),
                },
                _ => "all".into(),
            };
            groups.entry(key).or_default().push((r, here.map_or(0, |h| h.1), here));
        }
        let out: Vec<Value> = groups
            .iter()
            .map(|(k, g)| {
                let pairs: Vec<(&AlignedRead, u32)> = g.iter().map(|(r, d, _)| (*r, *d)).collect();
                let mut j = group_json(&pairs, single);
                j["group"] = json!(k);
                if list > 0 {
                    j["list"] = json!(g.iter().take(list).map(|(r, _, h)| read_json(r, *h)).collect::<Vec<_>>());
                }
                j
            })
            .collect();
        let mut result = json!({
            "query": {"bam": path, "at": at, "group_by": group_by, "filters": filters_json(&opts)},
            "reads": reads.len(),
            "read_counts": view.read_counts,
            "groups": out,
        });
        if let Some((p, _)) = &phasing {
            result["phase"] = json!({"het_sites": p.sites.len(), "blocks": p.blocks});
        }
        Ok(result)
    }

    // ---- summarize / coverage ----

    fn summarize(&mut self, args: &Value) -> Result<Value> {
        let opts = self.read_filters(args)?;
        let text = str_arg(args, "region").ok_or_else(|| anyhow!("give a `region`"))?.to_owned();
        let region = self.resolve(&text, args)?;
        let max_sites = uint_or(args, "max_sites", 10)?.min(100) as usize;
        // The reference for the whole region, at any size (as `ggv query summary`): without it,
        // alleles are judged against the commonest base instead, which changes the answer.
        let seq = self.fetch(args, &region)?;
        let reference = self.reference_path(args)?;
        let mut out = Vec::new();
        for path in self.bam_paths(args)? {
            let bam = self.bam(&path, reference.as_ref())?;
            let Ok(local) = bam.resolve_str(&format!("{}:{}-{}", region.contig, region.start + 1, region.end)) else {
                out.push(json!({"file": file_name(&path), "path": path, "error": "contig not in this file"}));
                continue;
            };
            let so = SummaryOptions { pileup: opts.clone(), max_variable_sites: max_sites, ..Default::default() };
            let s = ggv_query::summarize(bam, &local, seq.as_deref(), &so)?;
            let mut v = serde_json::to_value(&s)?;
            v["file"] = json!(file_name(&path));
            v["path"] = json!(path);
            if let Some(m) = v.get_mut("mean_depth") {
                *m = json!(round(s.mean_depth, 2));
            }
            if let Some(note) = clipped(&region, &local) {
                v["note"] = json!(note);
            }
            out.push(v);
        }
        Ok(json!({"query": {"region": text, "reference": reference, "filters": filters_json(&opts)}, "results": out}))
    }

    fn coverage(&mut self, args: &Value) -> Result<Value> {
        let opts = self.read_filters(args)?;
        let text = str_arg(args, "region").ok_or_else(|| anyhow!("give a `region`"))?.to_owned();
        let region = self.resolve(&text, args)?;
        let bins = uint_or(args, "bins", 50)?.clamp(1, 1000) as usize;
        let reference = self.reference_path(args)?;
        let mut out = Vec::new();
        for path in self.bam_paths(args)? {
            let bam = self.bam(&path, reference.as_ref())?;
            let Ok(local) = bam.resolve_str(&format!("{}:{}-{}", region.contig, region.start + 1, region.end)) else {
                out.push(json!({"file": file_name(&path), "path": path, "error": "contig not in this file"}));
                continue;
            };
            // Bins cover what this file has of the region: its contig may end sooner.
            let span = json!(format!("{}:{}-{}", local.contig, local.start + 1, local.end));
            let note = clipped(&region, &local);
            // A base-quality minimum needs each base's quality: count with the pileup engine
            // (as `evidence` and `summarize` do). The index estimate can't apply it.
            if opts.min_base_quality > 0 && local.len() <= 5_000_000 {
                let nb = bins.min(local.len() as usize).max(1);
                let (mut sum, mut max) = (vec![0u64; nb], vec![0u32; nb]);
                let len = local.len();
                ggv_query::pileup(bam, &local, None, &opts, |s| {
                    let i = bin_of(s.pos - 1 - local.start, nb, len);
                    sum[i] += u64::from(s.depth);
                    max[i] = max[i].max(s.depth);
                    Ok(())
                })?;
                let values: Vec<Value> = (0..nb)
                    .map(|i| {
                        let (a, b) = bin_bounds(i, nb, len);
                        json!([round(sum[i] as f64 / (b - a).max(1) as f64, 1), max[i]])
                    })
                    .collect();
                let mut r = json!({"file": file_name(&path), "path": path, "region": span, "kind": "counted", "bins": nb, "bin_size": len as f64 / nb as f64, "depth": values});
                if let Some(n) = note {
                    r["clipped"] = json!(n);
                }
                out.push(r);
                continue;
            }
            let view = load_view(bam, &local, None, String::new(), &ViewOptions { filters: opts.clone(), bins, max_read_window: 0, ..Default::default() })?;
            let (kind, values) = match &view.coverage {
                ggv_layout::Coverage::Binned(b) => ("counted", b.iter().map(|x| json!([round(x.mean, 1), x.max])).collect::<Vec<_>>()),
                ggv_layout::Coverage::Estimated(b) => ("estimated from the index", b.iter().map(|x| json!(round(x.mean, 1))).collect()),
                ggv_layout::Coverage::PerBase(c) => ("counted", c.iter().map(|x| json!(x.depth)).collect()),
                ggv_layout::Coverage::Unavailable => ("unavailable", Vec::new()),
            };
            let mut r = json!({"file": file_name(&path), "path": path, "region": span, "kind": kind, "bins": values.len(), "bin_size": local.len() as f64 / values.len().max(1) as f64, "depth": values});
            if let Some(n) = note {
                r["clipped"] = json!(n);
            }
            if kind.starts_with("estimated") {
                r["note"] = json!(format!(
                    "An estimate: read density from the index, calibrated by counting samples in each bin (by density stratum) with the read filters{}. Typically within a few percent of counted depth, but in repeats, above all with a MAPQ filter, single bins can be off by a fifth to a half (the index can't see MAPQ). Use `summarize`, or a region under 5 Mb, for counted depth.",
                    if opts.min_base_quality > 0 { "; the base-quality minimum can't be applied to it" } else { "" }
                ));
            }
            out.push(r);
        }
        Ok(json!({"query": {"region": text, "filters": filters_json(&opts), "format": "per bin: [mean, max] when counted, mean when estimated"}, "results": out}))
    }

    // ---- variants / genes / consequence ----

    fn variants(&mut self, args: &Value) -> Result<Value> {
        let text = str_arg(args, "region").ok_or_else(|| anyhow!("give a `region`"))?.to_owned();
        let max = uint_or(args, "max", 50)?.clamp(1, 500) as usize;
        let pass_only = bool_arg(args, "pass_only")? == Some(true);
        let paths = self.vcf_paths(args)?;
        if paths.is_empty() {
            bail!("no VCF: pass `vcfs`, or start the server with -v FILE");
        }
        // A gene name or a misspelt contig is resolved or refused like in the other region tools,
        // never searched as an empty contig.
        let mut known = Vec::new();
        for path in &paths {
            let v = self.vcf(path)?;
            let mut contigs = v.header_contigs();
            let unstated: Vec<_> = v.contig_names().iter().filter(|n| !contigs.iter().any(|c| c.name == **n)).cloned().collect();
            contigs.extend(unstated.into_iter().map(|name| ggv_model::Contig { name, length: u64::MAX / 2 }));
            known.extend(contigs);
        }
        let region = self.resolve_also(&text, args, known)?;
        let mut out = Vec::new();
        for path in paths {
            let samples = self.vcf(&path)?.samples().to_vec();
            let all = self.vcf(&path)?.query(&region, false)?;
            let kept: Vec<_> = all.iter().filter(|v| !pass_only || v.is_pass()).collect();
            out.push(json!({
                "file": file_name(&path),
                "records": kept.len(),
                "shown": kept.iter().take(max).map(|v| variant_json(v, &samples, &path)).collect::<Vec<_>>(),
            }));
        }
        Ok(json!({"query": {"region": text, "searched": region.to_string(), "pass_only": pass_only}, "results": out}))
    }

    fn genes(&mut self, args: &Value) -> Result<Value> {
        let paths = self.gene_paths(args)?;
        if paths.is_empty() {
            bail!("no annotation: pass `genes`, or start the server with -g FILE");
        }
        let all = bool_arg(args, "all")? == Some(true);
        let max = uint_or(args, "max", 100)?.clamp(1, 10_000) as usize;
        // A region is resolved like in the other region tools (a gene name works, an unknown
        // contig is refused), against the annotations' contigs too.
        let region = match (str_arg(args, "name"), str_arg(args, "region")) {
            (Some(_), _) => None,
            (None, Some(text)) => {
                let mut known = Vec::new();
                for path in &paths {
                    // Extents aren't lengths: a region past the last gene is just empty.
                    known.extend(self.annotation(path)?.contigs().into_iter().map(|c| ggv_model::Contig { length: u64::MAX / 2, ..c }));
                }
                Some(self.resolve_also(text, args, known)?)
            }
            (None, None) => bail!("give a `region` or a `name`"),
        };
        let mut out = Vec::new();
        for path in paths {
            let ann = self.annotation(&path)?;
            let mut note = None;
            let (contig, found): (String, Vec<&ggv_annot::Transcript>) = if let Some(name) = str_arg(args, "name") {
                // A version the file doesn't have: the accession's version it has, said so.
                let other_version = || {
                    let (base, v) = name.trim().rsplit_once('.')?;
                    let t = ann.find_transcript(base).filter(|_| !v.is_empty() && v.bytes().all(|b| b.is_ascii_digit()))?;
                    let found = ggv_annot::without_copy_suffix(&t.id);
                    Some((t, format!("{} is not in {}: it has {found}, listed instead", name.trim(), file_name(&path))))
                };
                if let Some(t) = ann.find_transcript(name) {
                    // A transcript id asks for that transcript, whether or not it's the main one.
                    (ann.contig_names()[t.contig as usize].clone(), vec![t])
                } else if let Some((t, n)) = other_version() {
                    note = Some(n);
                    (ann.contig_names()[t.contig as usize].clone(), vec![t])
                } else {
                    // A gene name asks for that gene only, not everything overlapping it, and
                    // for every copy of it (pseudo-autosomal genes are on chrX and chrY).
                    let spans = ann.find_all(name);
                    let first = spans.first().ok_or_else(|| anyhow!("{name} is not in {}", file_name(&path)))?;
                    if spans.len() > 1 {
                        let at: Vec<String> = spans.iter().map(|g| format!("{}:{}-{}", g.contig, g.start + 1, g.end)).collect();
                        note = Some(format!("{} is annotated on {} contigs ({}): each copy is listed", first.name, spans.len(), at.join(", ")));
                    }
                    let found = spans
                        .iter()
                        .flat_map(|g| ann.query(&g.contig, g.start, g.end).filter(|t| t.gene_name.eq_ignore_ascii_case(&g.name)))
                        .filter(|t| all || t.representative)
                        .collect();
                    (first.contig.clone(), found)
                }
            } else {
                let r = region.as_ref().expect("a region");
                let contig = std::iter::once(r.contig.clone()).chain(ggv_model::contig_aliases(&r.contig)).find(|c| ann.contig_index(c).is_some()).unwrap_or(r.contig.clone());
                let found = ann.query(&contig, r.start, r.end).filter(|t| all || t.representative).collect();
                (contig, found)
            };
            let total = found.len();
            let ts: Vec<Value> = found
                .into_iter()
                .take(max)
                .map(|t| {
                    let mut v = json!({
                        "gene": t.gene_name, "gene_id": t.gene_id, "transcript": ggv_annot::without_copy_suffix(&t.id), "biotype": (!t.biotype.is_empty()).then_some(&t.biotype), "contig": ann.contig_names()[t.contig as usize],
                        "strand": match t.strand { ggv_annot::Strand::Forward => "+", ggv_annot::Strand::Reverse => "-", ggv_annot::Strand::Unknown => "." },
                        "start": t.start + 1, "end": t.end, "exons": t.exons.len(),
                        "coding": t.cds.map(|(s, e)| json!({"start": s + 1, "end": e, "length": t.cds_len()})),
                        "mane_select": t.tags & ggv_annot::TAG_MANE_SELECT != 0,
                        "canonical": t.tags & ggv_annot::TAG_CANONICAL != 0,
                    });
                    if t.tags & ggv_annot::TAG_MANE_PLUS_CLINICAL != 0 {
                        v["mane_plus_clinical"] = json!(true);
                    }
                    // No exon rows in the file: the one exon is the transcript's span.
                    if t.exons_from_span {
                        v["exons_from_span"] = json!(true);
                    }
                    v
                })
                .collect();
            let mut r = json!({"file": file_name(&path), "contig": contig, "total": total, "shown": ts.len(), "truncated": total > ts.len(), "transcripts": ts});
            if total > ts.len() {
                note = Some(format!("{} of {total} transcripts shown, in position order: raise `max` (up to 10000) or ask for a narrower region", ts.len()));
            }
            if let Some(n) = note {
                r["note"] = json!(n);
            }
            out.push(r);
        }
        Ok(json!({"results": out}))
    }

    fn consequence(&mut self, args: &Value) -> Result<Value> {
        let text = str_arg(args, "locus").ok_or_else(|| anyhow!("give a `locus`"))?.to_owned();
        let (r, a) = (str_arg(args, "ref").unwrap_or("").to_owned(), str_arg(args, "alt").unwrap_or("").to_owned());
        if r.is_empty() || a.is_empty() {
            bail!("give `ref` and `alt` (VCF style: for indels include the base before)");
        }
        if a.contains(',') {
            bail!("ALT {a} lists several alleles (a multi-allelic record): give one ALT at a time");
        }
        if !ggv_annot::coding::is_sequence(&r) || !ggv_annot::coding::is_sequence(&a) {
            bail!("{r}>{a} is not a sequence change: symbolic alleles (<DEL>, <DUP>, breakends) don't say which bases change, so their effect on a transcript isn't computed here. Pass the REF and ALT sequences instead");
        }
        if r.eq_ignore_ascii_case(&a) {
            bail!("REF and ALT are both {r}: that is not a change, so it has no consequence");
        }
        let l: ggv_model::Locus = text.parse().map_err(|e| anyhow!("{text}: {e}"))?;
        let pos = l.start.ok_or_else(|| anyhow!("{text} needs a position"))? - 1;
        // A range must be REF's own span: anything else names two different places.
        if let Some(end) = l.end
            && end > pos + 1
            && end != pos + r.len() as u64
        {
            bail!("{text} covers {} bases but REF {r} has {}: give REF's first base (or the range it covers)", end - pos, r.len());
        }
        let paths = self.gene_paths(args)?;
        if paths.is_empty() {
            bail!("no annotation: pass `genes`, or start the server with -g FILE");
        }
        let mut out = Vec::new();
        let fetch_args = args.clone();
        // Check REF against the reference, when there is one.
        let mut warning = None;
        if let Some(found) = self.fetch(args, &Region { contig: l.contig.clone(), start: pos, end: pos + r.len() as u64 })? {
            if !found.eq_ignore_ascii_case(r.as_bytes()) {
                warning = Some(format!("REF {r} does not match the reference here ({})", String::from_utf8_lossy(&found)));
            }
        }
        let all = all_transcripts(args)?;
        let mut hidden = 0;
        for path in paths {
            let ann = self.annotation(&path)?;
            // Every transcript the REF touches, not only the one at its first base.
            let (shown, more) = pick_transcripts(transcripts_over(&ann, &l.contig, pos, pos + r.len() as u64), all);
            hidden += more;
            for t in shown {
                if let Some(c) = coding::consequence(t, pos, &r, &a, |q| self.base_at(&fetch_args, &l.contig, q)) {
                    let mut v = json!({"gene": c.gene, "transcript": c.transcript, "hgvs_c": c.hgvs_c, "hgvs_p": c.hgvs_p, "effect": c.effect, "where": c.location.feature, "codons": c.codons});
                    if !c.warnings.is_empty() {
                        v["warnings"] = json!(c.warnings);
                    }
                    out.push(v);
                }
            }
        }
        let mut v = json!({"change": format!("{}:{} {r}>{a}", l.contig, pos + 1), "warning": warning, "consequences": out});
        if hidden > 0 {
            v["other_transcripts"] = json!(hidden);
        }
        Ok(v)
    }

    // ---- tandem repeats ----

    fn repeats(&mut self, args: &Value) -> Result<Value> {
        let opts = self.read_filters(args)?;
        let text = str_arg(args, "region").ok_or_else(|| anyhow!("give the repeat's `region`"))?.to_owned();
        let region = self.resolve(&text, args)?;
        if region.len() > 50_000 {
            bail!("{text} is {} bp: give the repeat itself (up to 50 kb)", region.len());
        }
        let motif_len = str_arg(args, "motif").map(|m| m.trim().len() as f64).filter(|&n| n > 0.0);
        let copies = |len: u64| motif_len.map(|m| (len as f64 / m * 10.0).round() / 10.0);
        let list = bool_arg(args, "list")? == Some(true);
        let reference = self.reference_path(args)?;
        let mut files = Vec::new();
        for path in self.bam_paths(args)? {
            let bam = self.bam(&path, reference.as_ref())?;
            let local = bam.resolve_str(&format!("{}:{}-{}", region.contig, region.start + 1, region.end))?;
            let r = ggv_query::repeat_lengths(bam, &local.contig, local.start, local.end, &opts)?;
            let alleles: Vec<Value> = r
                .alleles
                .iter()
                .map(|a| {
                    json!({"length": a.length, "copies": copies(a.length), "change": a.length as i64 - r.reference_length as i64,
                        "range": [a.min, a.max], "reads": a.reads,
                        "by_haplotype": {"none": a.by_haplotype[0], "1": a.by_haplotype[1], "2": a.by_haplotype[2]}})
                })
                .collect();
            let mut f = json!({
                "file": file_name(&path), "spanning": r.spanning, "partial": r.partial, "skipping": r.skipping, "alleles": alleles,
                "haplotype_lengths": {"none": r.haplotype_lengths[0], "1": r.haplotype_lengths[1], "2": r.haplotype_lengths[2]},
            });
            if r.unassigned > 0 {
                f["unassigned"] = json!(r.unassigned);
            }
            // Few reads measured and many set aside is what a region narrower than the repeat
            // tract looks like: an indel placed at the tract's other end falls on a flank.
            if r.partial >= MIN_PARTIAL_FOR_WARNING && r.partial > r.spanning {
                f["warning"] = json!(format!(
                    "only {} of {} reads here span the region with {} aligned bases either side; the other {} end inside it or have an indel on a flank. If the region is part of a longer repeat tract, give the whole tract: aligners may place a length change anywhere in it",
                    r.spanning,
                    r.spanning + r.partial,
                    ggv_query::repeat_flank(),
                    r.partial
                ));
            }
            if list {
                f["read_lengths"] = json!(r.reads.iter().map(|x| json!([x.length, x.haplotype])).collect::<Vec<_>>());
            }
            files.push(f);
        }
        let ref_len = region.len();
        Ok(json!({
            "query": {"region": region.to_string(), "filters": filters_json(&opts), "flank": ggv_query::repeat_flank()},
            "reference_length": ref_len, "reference_copies": copies(ref_len),
            "files": files,
        }))
    }

    // ---- structural variants ----

    fn svs(&mut self, args: &Value) -> Result<Value> {
        let opts = self.read_filters(args)?;
        let text = str_arg(args, "region").ok_or_else(|| anyhow!("give a `region`"))?.to_owned();
        let region = self.resolve(&text, args)?;
        if region.len() > 2_000_000 {
            bail!("{text} is {} bp; SVs are tested over up to 2 Mb", region.len());
        }
        let min_len = uint_or(args, "min_len", 50)?.max(20);
        let max = usize::try_from(uint_or(args, "max", 30)?).unwrap_or(usize::MAX);
        let min_reads = u32_or(args, "min_reads", 3)?;
        let min_event = u32::try_from(min_len).map_err(|_| anyhow!("`min_len` must be at most {}, not {min_len}", u32::MAX))?;
        let reference = self.reference_path(args)?;
        let bams = self.bam_paths(args)?;

        // Every SV record in the region; the longest `max` are tested against the reads, and
        // all of them explain discovered events.
        let mut records: Vec<SvRecord> = Vec::new();
        let mut breakends = Vec::new();
        let mut warnings = Vec::new();
        // Records not listed on their own: <INS> without a length, a breakend pair's mates.
        let (mut unsized_ins, mut mates) = (0, 0);
        for path in self.vcf_paths(args)? {
            let samples = self.vcf(&path)?.samples().to_vec();
            for v in self.vcf(&path)?.query(&region, false)? {
                let id = if v.id.is_empty() || v.id == "." { Value::Null } else { json!(v.id) };
                if v.svlen_mismatch() {
                    warnings.push(format!(
                        "{}:{} in {}: SVLEN has {} values for {} ALTs, so it is not used (lengths from END; insertions skipped)",
                        region.contig,
                        v.start + 1,
                        file_name(&path),
                        v.svlen.len(),
                        v.alts.len()
                    ));
                }
                let genotypes = sample_genotypes(&v, &samples);
                let filter = if v.is_pass() { "PASS".to_owned() } else { v.filters.join(";") };
                for (i, alt) in v.alts.iter().enumerate() {
                    if alt.trim_start_matches('<').split([':', '>']).next() == Some("INS") && alt.starts_with('<') && v.svlen_of(i).is_none() && !v.svlen_mismatch() {
                        unsized_ins += 1;
                        warnings.push(format!(
                            "{}:{} in {}: {alt} has no SVLEN, so its length is unknown and it isn't tested (reads that carry it are listed in `unexplained`)",
                            region.contig,
                            v.start + 1,
                            file_name(&path)
                        ));
                        continue;
                    }
                    let found = if ggv_io::is_breakend(alt) {
                        match breakend_site(&v, alt, &region.contig) {
                            Ok(site) => Some((site, "BND".to_owned())),
                            Err(reason) => {
                                breakends.push(json!({"pos": v.start + 1, "alt": alt, "id": id, "genotypes": genotypes, "filter": filter, "vcf": file_name(&path), "reason": reason}));
                                None
                            }
                        }
                    } else {
                        sv_site(&v, i, alt, min_len)
                    };
                    let Some((site, kind)) = found else { continue };
                    let vcf = file_name(&path);
                    // A breakend pair is one event: its second mate joins the first.
                    if kind == "BND"
                        && let Some(first) = records.iter_mut().find(|o| o.kind == "BND" && o.mate.is_none() && o.vcf == vcf && same_breakends(&o.site, &site))
                    {
                        first.mate = Some(json!({"pos": v.start + 1, "id": id, "alt": alt}));
                        mates += 1;
                        continue;
                    }
                    let mate_alt = (kind == "BND").then(|| alt.clone());
                    records.push(SvRecord { site, kind, id: id.clone(), genotypes: genotypes.clone(), filter: filter.clone(), vcf, alt: mate_alt, vcf_pos: v.start + 1, mate: None });
                }
            }
        }
        let total = records.len();
        let tested = pick_tested(&records, max);

        // Each record's read support, per file, tested in parallel (a handle on the file per
        // thread; the records are independent).
        let mut support: Vec<Vec<ggv_query::SvSupport>> = vec![Vec::new(); tested.len()];
        let sites: Vec<ggv_query::SvSite> = tested.iter().map(|&i| records[i].site).collect();
        for path in &bams {
            let bam = self.bam(path, reference.as_ref())?;
            let contig = bam.resolve_str(&format!("{}:1", region.contig))?.contig;
            let threads = std::thread::available_parallelism().map_or(1, |n| n.get()).min(8).min(sites.len()).max(1);
            let mut handles: Vec<ggv_io::BamFile> = (1..threads).filter_map(|_| bam.reopen().ok()).collect();
            let results: Vec<Result<ggv_query::SvSupport>> = {
                let mut all: Vec<&mut ggv_io::BamFile> = vec![bam];
                all.extend(handles.iter_mut());
                let n = all.len();
                std::thread::scope(|scope| {
                    let jobs: Vec<_> = all
                        .into_iter()
                        .enumerate()
                        .map(|(t, b)| {
                            let (sites, contig, opts) = (&sites, &contig, &opts);
                            scope.spawn(move || {
                                (t..sites.len()).step_by(n).map(|i| (i, ggv_query::sv_support(b, contig, &sites[i], opts).map_err(anyhow::Error::from))).collect::<Vec<_>>()
                            })
                        })
                        .collect();
                    let mut flat: Vec<(usize, Result<ggv_query::SvSupport>)> = jobs.into_iter().flat_map(|j| j.join().expect("sv worker")).collect();
                    flat.sort_by_key(|(i, _)| *i);
                    flat.into_iter().map(|(_, r)| r).collect()
                })
            };
            for (i, r) in results.into_iter().enumerate() {
                support[i].push(r?);
            }
        }

        let mut out_records = Vec::new();
        for (i, &ri) in tested.iter().enumerate() {
            let r = &records[ri];
            let site = &r.site;
            // Other records for the same event (a call set may list one SV per haplotype).
            let same_site: Vec<Value> = records
                .iter()
                .enumerate()
                .filter(|(j, o)| *j != ri && o.site.kind == site.kind && o.site.start.abs_diff(site.start) <= site.slack() && o.site.len.abs_diff(site.len) * 10 <= o.site.len.max(site.len) * 3)
                .map(|(_, o)| json!({"pos": o.pos(), "type": o.kind, "length": o.length(), "genotypes": o.genotypes}))
                .collect();
            let mut per_file = Vec::new();
            for (path, s) in bams.iter().zip(&support[i]) {
                per_file.push(json!({
                    "file": file_name(path), "supporting": s.supporting, "reference": s.reference, "other_allele": s.other_allele, "clipped": s.clipped,
                    "fraction": s.fraction(), "genotype": s.genotype(), "median_event_len": s.median_event_len,
                    "in_cigar": s.in_cigar, "split": s.split, "strands": [s.forward, s.reverse],
                    "by_haplotype": {"support": {"none": s.support_by_haplotype[0], "1": s.support_by_haplotype[1], "2": s.support_by_haplotype[2]},
                                     "reference": {"none": s.reference_by_haplotype[0], "1": s.reference_by_haplotype[1], "2": s.reference_by_haplotype[2]}},
                }));
            }
            let mut rec = json!({
                "pos": r.pos(), "end": r.end(), "type": r.kind, "length": r.length(),
                "id": r.id, "genotypes": r.genotypes, "filter": r.filter, "vcf": r.vcf,
                "reads": per_file,
                "same_site": if same_site.is_empty() { Value::Null } else { json!({"records": same_site, "note": "other records for the same event here (e.g. one per haplotype): the reads count them together"}) },
            });
            if let Some(alt) = &r.alt {
                rec["alt"] = json!(alt);
                rec["tested_as"] = json!(site.kind);
            }
            if let Some(mate) = &r.mate {
                rec["mate"] = mate.clone();
            }
            out_records.push(rec);
        }

        let mut unexplained = Vec::new();
        if bool_arg(args, "discover")? != Some(false) {
            for path in &bams {
                let bam = self.bam(path, reference.as_ref())?;
                let local = bam.resolve_str(&format!("{}:{}-{}", region.contig, region.start + 1, region.end))?;
                let (events, _) = ggv_query::sv_events(bam, &local, min_event, min_reads, &opts)?;
                for e in events {
                    let site = ggv_query::SvSite { kind: e.kind, start: e.pos, end: e.pos, len: e.len };
                    let known = explained(&records, &site);
                    if !known {
                        unexplained.push(json!({"file": file_name(path), "type": e.kind, "pos": e.pos + 1, "length": e.len, "reads": e.reads,
                            "by_haplotype": {"none": e.by_haplotype[0], "1": e.by_haplotype[1], "2": e.by_haplotype[2]}}));
                    }
                }
            }
        }
        let mut out = json!({
            "query": {"region": region.to_string(), "filters": filters_json(&opts), "min_len": min_len},
            "records": out_records, "tested": out_records.len(), "records_in_region": total + breakends.len() + unsized_ins + mates,
            "unexplained": unexplained,
        });
        if total > tested.len() {
            out["not_tested"] = json!({"records": total - tested.len(), "note": format!("only the {} longest records were tested against the reads (`max`); `unexplained` is still checked against all {total}", tested.len())});
        }
        if !breakends.is_empty() {
            out["breakends_not_tested"] = json!(breakends);
        }
        if !warnings.is_empty() {
            out["warnings"] = json!(warnings);
        }
        Ok(out)
    }

    // ---- junctions ----

    fn junctions(&mut self, args: &Value) -> Result<Value> {
        let reads = self.read_filters(args)?;
        let text = str_arg(args, "region").ok_or_else(|| anyhow!("give a `region`"))?.to_owned();
        let region = match self.resolve(&text, args) {
            Ok(r) => r,
            Err(e) => {
                // A gene name: its span from the annotation.
                let paths = self.gene_paths(args)?;
                let found = paths.first().map(|p| self.annotation(p)).transpose()?.and_then(|a| a.find(&text).cloned());
                let g = found.ok_or(e)?;
                Region { contig: g.contig, start: g.start, end: g.end }
            }
        };
        if region.len() > 10_000_000 {
            bail!("{text} is {} bp; junctions are counted over up to 10 Mb", region.len());
        }
        let opts = ggv_query::JunctionOptions { reads: reads.clone(), min_anchor: u32_or(args, "min_anchor", 1)? };
        let reference = self.reference_path(args)?;
        let mut files = Vec::new();
        let mut names = Vec::new();
        for path in self.bam_paths(args)? {
            let bam = self.bam(&path, reference.as_ref())?;
            let local = bam.resolve_str(&format!("{}:{}-{}", region.contig, region.start + 1, region.end))?;
            files.push(ggv_query::junctions(bam, &local, &opts)?);
            names.push(path);
        }
        let ann = self.gene_paths(args)?.first().map(|p| self.annotation(p)).transpose()?;
        let min_reads = u32_or(args, "min_reads", 2)?;
        let max = uint_or(args, "max", 200)?.clamp(1, 2000) as usize;
        let fetch_args = args.clone();
        let rep = crate::junctions::report(&region.contig, &files, |q| self.base_at(&fetch_args, &region.contig, q), ann.as_deref(), min_reads, max);
        let file_info: Vec<Value> = names
            .iter()
            .zip(&files)
            .map(|(p, f)| json!({"file": file_name(p), "path": p, "reads_used": f.reads.used, "spliced_reads": f.spliced_reads}))
            .collect();
        let junctions: Vec<Value> = rep
            .junctions
            .iter()
            .map(|j| {
                json!({
                    "intron": format!("{}:{}-{}", region.contig, j.start, j.end), "length": j.length, "strand": j.strand.to_string(),
                    "motif": j.motif, "annotation": j.annotation, "genes": j.genes,
                    "samples": j.samples.iter().map(|c| json!({"reads": c.reads, "unique": c.unique, "donor_share": c.donor_share, "acceptor_share": c.acceptor_share})).collect::<Vec<_>>(),
                })
            })
            .collect();
        Ok(json!({
            "query": {"region": region.to_string(), "filters": filters_json(&reads), "min_reads": min_reads, "min_anchor": opts.min_anchor},
            "files": file_info,
            "shown": junctions.len(),
            "total": rep.total,
            "junctions": junctions,
        }))
    }

    // ---- phase / compare ----

    fn phase(&mut self, args: &Value) -> Result<Value> {
        let opts = self.read_filters(args)?;
        let text = str_arg(args, "region").ok_or_else(|| anyhow!("give a `region`"))?.to_owned();
        let region = self.resolve(&text, args)?;
        if region.len() > READ_WINDOW {
            bail!("phase works on up to {READ_WINDOW} bp at a time");
        }
        let seq = self.fetch(args, &region)?.ok_or_else(|| anyhow!("phase needs a reference (sites come from mismatches)"))?;
        let path = self.bam_paths(args)?.remove(0);
        let reference = self.reference_path(args)?;
        let bam = self.bam(&path, reference.as_ref())?;
        let local = bam.resolve_str(&format!("{}:{}-{}", region.contig, region.start + 1, region.end))?;
        let view = load_view(bam, &local, Some(seq.clone()), String::new(), &ViewOptions { filters: opts.clone(), max_rows: 0, ..Default::default() })?;
        let reads = view.reads.as_ref().map(|(r, _)| r.as_slice()).unwrap_or(&[]);
        let p = ggv_layout::phase::phase(reads, &(0..reads.len()).collect::<Vec<_>>());
        let max_sites = uint_or(args, "max_sites", DEFAULT_PHASE_SITES)?.min(MAX_PHASE_SITES) as usize;
        let sites: Vec<Value> = p
            .sites
            .iter()
            .zip(&p.site_blocks)
            .take(max_sites)
            .map(|(&(pos, alt), &block)| {
                let r = seq.get((pos - region.start) as usize).copied();
                let mut counts = [[0u32; 2]; 2]; // [hap A, hap B] x [ref, alt]
                for (i, read) in reads.iter().enumerate() {
                    let Some(h) = p.hap[i] else { continue };
                    if let Some((Allele::Base(b), _)) = read_at(read, pos, r) {
                        let k = usize::from(b == alt);
                        if b == alt || Some(b) == r {
                            counts[usize::from(!h)][k] += 1;
                        }
                    }
                }
                let carries = |c: [u32; 2]| if c[1] > c[0] { (alt as char).to_string() } else { r.map_or("?".into(), |b| (b as char).to_string()) };
                json!({"pos": pos + 1, "block": block + 1, "ref": r.map(|b| (b as char).to_string()), "alt": (alt as char).to_string(), "A": carries(counts[0]), "B": carries(counts[1]), "A_reads": counts[0], "B_reads": counts[1]})
            })
            .collect();
        let a = p.hap.iter().filter(|h| **h == Some(true)).count();
        let b = p.hap.iter().filter(|h| **h == Some(false)).count();
        // Each block: its sites' span and its own reads per haplotype.
        let blocks: Vec<Value> = (0..p.blocks)
            .map(|k| {
                let at: Vec<u64> = p.sites.iter().zip(&p.site_blocks).filter(|(_, b)| **b == k).map(|(s, _)| s.0 + 1).collect();
                let of = |h: Option<bool>| p.read_blocks.iter().zip(&p.hap).filter(|(rb, hh)| **rb == Some(k) && **hh == h).count();
                json!({"block": k + 1, "start": at.first(), "end": at.last(), "sites": at.len(), "phase_A": of(Some(true)), "phase_B": of(Some(false)), "unphased": of(None)})
            })
            .collect();
        Ok(json!({
            "query": {"bam": path, "region": text, "filters": filters_json(&opts)},
            "reads": reads.len(), "phase_A": a, "phase_B": b, "unphased": reads.len() - a - b,
            "het_sites": p.sites.len(), "blocks": p.blocks,
            "note": "A and B are consistent within a block, not across blocks: compare sites (and reads) only within one `block`. phase_A/phase_B at the top add up every block; each block lists its own. A_reads/B_reads are [ref, alt] counts. A block's start and end are its first and last sites (1-based).",
            "phase_blocks": blocks,
            "sites_shown": sites.len(),
            "sites_truncated": (sites.len() < p.sites.len()).then(|| format!("the leftmost {} of {} sites: raise `max_sites` (up to {MAX_PHASE_SITES}) or narrow the region for the rest", sites.len(), p.sites.len())),
            "sites": sites,
        }))
    }

    fn compare(&mut self, args: &Value) -> Result<Value> {
        // Each file's filters: its own (`a_filters` / `b_filters`) where given, else the call's.
        let own_filters = |key: &str| -> Result<Option<Value>> {
            match args.get(key) {
                None | Some(Value::Null) => Ok(None),
                Some(Value::Object(own)) => {
                    let mut merged = args.clone();
                    if let Some(m) = merged.as_object_mut() {
                        for k in ["filter", "min_mapq", "include_duplicates", "exclude_flags", "min_base_quality"] {
                            m.remove(k);
                        }
                        m.extend(own.clone());
                    }
                    Ok(Some(merged))
                }
                Some(other) => bail!("`{key}` must be an object of read filters, not {other}"),
            }
        };
        let opts = self.read_filters(args)?;
        let per_file = [own_filters("a_filters")?, own_filters("b_filters")?];
        let file_opts = [
            per_file[0].as_ref().map(|f| self.read_filters(f)).transpose()?.unwrap_or_else(|| opts.clone()),
            per_file[1].as_ref().map(|f| self.read_filters(f)).transpose()?.unwrap_or_else(|| opts.clone()),
        ];
        let same_filters = per_file.iter().all(Option::is_none) || filters_json(&file_opts[0]) == filters_json(&file_opts[1]);
        let text = str_arg(args, "region").ok_or_else(|| anyhow!("give a `region`"))?.to_owned();
        // The pair first, then the region against them: an unrelated default file must neither
        // refuse the region nor shorten it.
        let (given_a, given_b) = (str_arg(args, "a").map(PathBuf::from), str_arg(args, "b").map(PathBuf::from));
        let defaults = if given_a.is_some() && given_b.is_some() { Vec::new() } else { self.bam_paths(args)? };
        let a = given_a.or_else(|| defaults.first().cloned()).ok_or_else(|| anyhow!("give `a` and `b`"))?;
        let b = given_b.or_else(|| defaults.get(1).cloned()).ok_or_else(|| anyhow!("give `b` (or start the server with two BAMs)"))?;
        let mut pair = args.clone();
        if let Some(m) = pair.as_object_mut() {
            m.remove("bam");
            m.insert("bams".into(), json!([a, b]));
        }
        let region = self.resolve(&text, &pair)?;
        if region.len() > 1_000_000 {
            bail!("compare works on up to 1 Mb at a time");
        }
        let threshold = match args.get("threshold") {
            None | Some(Value::Null) => 0.2,
            Some(x) => x.as_f64().filter(|t| (0.01..=1.0).contains(t)).ok_or_else(|| anyhow!("`threshold` must be a number from 0.01 to 1, not {x}"))?,
        };
        let min_depth = u32_or(args, "min_depth", 10)?.max(1);
        let max = uint_or(args, "max", 50)?.clamp(1, 500) as usize;
        let reference = self.reference_path(args)?;
        let seq = self.fetch(args, &region)?;
        // Per position: reads per allele (A, C, G, T, N, DEL), reads with an insertion after
        // the base, and depth.
        let mut counts: Vec<Vec<Option<([u32; 6], u32, u32)>>> = Vec::new();
        for (path, opts) in [&a, &b].into_iter().zip(&file_opts) {
            let bam = self.bam(path, reference.as_ref())?;
            let local = bam.resolve_str(&format!("{}:{}-{}", region.contig, region.start + 1, region.end)).with_context(|| path.display().to_string())?;
            let mut v = vec![None; region.len() as usize];
            ggv_query::pileup(bam, &local, seq.as_deref(), opts, |s| {
                if let Some(slot) = v.get_mut((s.pos - 1 - local.start) as usize) {
                    let x = &s.bases;
                    *slot = Some(([x.a.total, x.c.total, x.g.total, x.t.total, x.n.total, s.deletions.total], s.insertions_after.total, s.depth));
                }
                Ok(())
            })?;
            counts.push(v);
        }
        const NAMES: [&str; 6] = ["A", "C", "G", "T", "N", "DEL"];
        // Strictly more than the threshold; the tolerance keeps exact ratios (1/5 vs 0.2) from
        // flipping on floating-point noise.
        let differs = |d: f64| d > threshold + 1e-9;
        let mut sites = Vec::new();
        let mut total = 0;
        for i in 0..region.len() as usize {
            let (Some((ca, ia, da)), Some((cb, ib, db))) = (counts[0][i], counts[1][i]) else { continue };
            if da < min_depth || db < min_depth {
                continue;
            }
            let pos = region.start + i as u64 + 1;
            let ref_base = seq.as_ref().and_then(|s| s.get(i)).map(|b| (*b as char).to_string());
            // The allele with the biggest change; of a gained/lost pair, the rarer one.
            let best = (0..6)
                .map(|k| {
                    let (fa, fb) = (f64::from(ca[k]) / f64::from(da), f64::from(cb[k]) / f64::from(db));
                    (k, fa, fb, (fa - fb).abs())
                })
                .max_by(|x, y| x.3.total_cmp(&y.3).then((y.1 + y.2).total_cmp(&(x.1 + x.2))))
                .expect("six alleles");
            let mut found = Vec::new();
            if differs(best.3) {
                found.push((NAMES[best.0], best.1, best.2));
            }
            // Insertions sit between bases, beside whatever base the read shows: compared on
            // their own, as reads with an insertion right after this base.
            let (fa, fb) = (f64::from(ia) / f64::from(da), f64::from(ib) / f64::from(db));
            if differs((fa - fb).abs()) {
                found.push(("INS", fa, fb));
            }
            for (allele, fa, fb) in found {
                total += 1;
                if sites.len() < max {
                    let mut site = json!({
                        "pos": pos, "ref": ref_base, "allele": allele,
                        "a_fraction": round(fa, 3), "b_fraction": round(fb, 3),
                        "a_depth": da, "b_depth": db,
                    });
                    if allele == "INS" {
                        site["note"] = json!("reads with an insertion (of any sequence) right after this base");
                    }
                    sites.push(site);
                }
            }
        }
        Ok(json!({
            "query": {"a": a, "b": b, "region": text, "compared": region.to_string(), "threshold": threshold, "min_depth": min_depth,
                      "filters": if same_filters { filters_json(&file_opts[0]) } else { json!({"a": filters_json(&file_opts[0]), "b": filters_json(&file_opts[1])}) },
                      "alleles": "A, C, G, T, N, DEL (per base) and INS (an insertion after the base)"},
            "differences": total,
            "shown": sites,
            "warning": (same_file(&a, &b) && same_filters).then_some("`a` and `b` are the same file under the same filters: no differences is expected, not a finding"),
        }))
    }

    // ---- render ----

    fn signal(&mut self, args: &Value) -> Result<Value> {
        let text = str_arg(args, "region").ok_or_else(|| anyhow!("give a `region`"))?;
        let bins = uint_or(args, "bins", 20)?.clamp(1, 1000) as usize;
        let mut files = strs(args, "signals")?;
        if files.is_empty() {
            files = self.defaults.signals.clone();
        }
        if files.is_empty() {
            bail!("no signal files: pass `signals`, start the server with -s, or (attached) open one in the window");
        }
        // A text signal file knows only where its data ends: the reference's (or alignment
        // file's) contig lengths cut a region past a contig's end, as in `coverage`.
        let reference = self.reference_path(args)?;
        let known: Vec<ggv_model::Contig> = if let Ok(bams) = self.bam_paths(args) {
            self.bam(&bams[0], reference.as_ref()).map(|b| b.contigs().to_vec()).unwrap_or_default()
        } else if let Some(r) = &reference {
            self.with_fasta(r, |f| f.contigs().to_vec()).unwrap_or_default()
        } else {
            Vec::new()
        };
        let mut results = Vec::new();
        for path in &files {
            let mut sig = ggv_io::signal::Signal::open(path).with_context(|| path.display().to_string())?;
            let mut region = match sig.resolve(text) {
                Ok(r) => r,
                Err(e) => {
                    results.push(json!({"file": file_name(path), "path": path, "error": e.to_string()}));
                    continue;
                }
            };
            let asked = region.clone();
            if let Ok(whole) = ggv_io::resolve_locus_str(&known, &region.contig) {
                if region.start >= whole.end {
                    results.push(json!({"file": file_name(path), "path": path, "error": format!("{text} starts past the end of {} (length {})", region.contig, whole.end)}));
                    continue;
                }
                region.end = region.end.min(whole.end);
            }
            let note = clipped(&asked, &region);
            // At most one bin per base (as `ggv query signal`): finer bins fall between bases
            // and come back empty.
            let bins = ggv_io::signal::bin_edges(region.start, region.end, bins).len();
            let mut rows = Vec::new();
            for (k, row) in sig.rows().iter().enumerate() {
                let bs = sig.bins(k, &region.contig, region.start, region.end, bins)?;
                let r3 = |x: f64| round(x, 3);
                let with_data: Vec<&ggv_io::signal::Bin> = bs.iter().flatten().collect();
                rows.push(json!({
                    "row": row,
                    "mean": bs.iter().map(|b| b.map(|b| r3(b.mean))).collect::<Vec<_>>(),
                    "max": bs.iter().map(|b| b.map(|b| r3(b.max))).collect::<Vec<_>>(),
                    "overall": if with_data.is_empty() { Value::Null } else { json!({
                        "mean": r3(with_data.iter().map(|b| b.mean * b.covered as f64).sum::<f64>() / with_data.iter().map(|b| b.covered as f64).sum::<f64>()),
                        "min": r3(with_data.iter().map(|b| b.min).fold(f64::INFINITY, f64::min)),
                        "max": r3(with_data.iter().map(|b| b.max).fold(f64::NEG_INFINITY, f64::max)),
                    }) },
                }));
            }
            let mut r = json!({"file": file_name(path), "path": path, "region": format!("{}:{}-{}", region.contig, region.start + 1, region.end), "bins": bins, "bin_size": round(region.len() as f64 / bins as f64, 1), "rows": rows});
            if let Some(n) = note {
                r["note"] = json!(n);
            }
            results.push(r);
        }
        Ok(json!({"query": {"signals": files, "region": text, "bins": bins}, "results": results}))
    }

    fn render(&mut self, args: &Value) -> Result<Output> {
        let opts = self.read_filters(args)?;
        let text = str_arg(args, "locus").or_else(|| str_arg(args, "region")).ok_or_else(|| anyhow!("give a `locus` or a `region`"))?.to_owned();
        let target = self.resolve(&text, args)?;
        let width = uint_or(args, "width", 1000)?.clamp(300, 2400) as u32;
        let link = bool_arg(args, "link")?.unwrap_or(false);
        let sashimi = bool_arg(args, "sashimi")?.unwrap_or(false);
        let min_junction_reads = u32_or(args, "min_junction_reads", 2)?;
        let tag_of = |key: &str| -> Result<Option<[u8; 2]>> {
            match str_arg(args, key).and_then(|v| v.strip_prefix("tag:")) {
                Some(t) => ggv_layout::parse_tag_name(t).map(Some).ok_or_else(|| anyhow!("{key}: `tag:{t}` needs a two-character tag name")),
                None => Ok(None),
            }
        };
        let (color_tag, group_tag) = (tag_of("color_by")?, tag_of("group_by")?);
        if let (Some(a), Some(b)) = (color_tag, group_tag)
            && a != b
        {
            bail!("color_by and group_by can use one tag at a time");
        }
        let tag = group_tag.or(color_tag);
        let color_by = match str_arg(args, "color_by").unwrap_or("none") {
            c if c.starts_with("tag:") => ggv_render::ColorBy::Tag,
            "none" => ggv_render::ColorBy::None,
            "strand" => ggv_render::ColorBy::Strand,
            "haplotype" => ggv_render::ColorBy::Haplotype,
            "insert-size" => ggv_render::ColorBy::InsertSize,
            "base-mods" => ggv_render::ColorBy::BaseMods,
            other => bail!("color_by: `{other}` is not one of none, strand, haplotype, insert-size, base-mods"),
        };
        let group = match str_arg(args, "group_by").unwrap_or("none") {
            g if g.starts_with("tag:") => ggv_layout::GroupBy::Tag(tag.expect("parsed above")),
            "none" => ggv_layout::GroupBy::None,
            "strand" => ggv_layout::GroupBy::Strand,
            "haplotype" => ggv_layout::GroupBy::Haplotype,
            "read group" | "read-group" => ggv_layout::GroupBy::ReadGroup,
            "alleles" => ggv_layout::GroupBy::Alleles,
            "base" if target.len() == 1 => ggv_layout::GroupBy::BaseAt(target.start),
            "base" => bail!("group_by base needs a `locus` (a single base)"),
            other => bail!("group_by: `{other}` is not one of none, strand, haplotype, read group, alleles, base"),
        };
        let region = if target.len() == 1 {
            // Clamped to the contig, so the picture and its note stop at its end.
            self.resolve(&format!("{}:{}-{}", target.contig, target.start.saturating_sub(50) + 1, target.end + 50), args)?
        } else {
            target.clone()
        };
        let seq = if region.len() <= READ_WINDOW { self.fetch(args, &region)? } else { None };
        let reference = self.reference_path(args)?;
        let mut scenes = Vec::new();
        for path in self.bam_paths(args)? {
            let bam = self.bam(&path, reference.as_ref())?;
            let Ok(local) = bam.resolve_str(&format!("{}:{}-{}", region.contig, region.start + 1, region.end)) else { continue };
            let label = file_name(&path);
            let mut view = load_view(bam, &local, seq.clone(), label, &ViewOptions { filters: opts.clone(), bins: width as usize, max_rows: 40, tag, ..Default::default() })?;
            if group != ggv_layout::GroupBy::None || link {
                ggv_layout::repack(&mut view, 40, group, None, link);
            }
            scenes.push(ggv_render::build_scene(&view, &ggv_render::RenderOptions { width, color_by, link, highlight: Some((target.start, target.end)), sashimi, min_junction_reads, ..Default::default() }));
        }
        let png = crate::batch::stack(scenes).to_png(1.0).map_err(anyhow::Error::msg)?;
        Ok(Output::Image { png, note: format!("{}:{}-{} ({} bp), {}", region.contig, region.start + 1, region.end, region.len(), if target.len() == 1 { "the base outlined" } else { "as asked" }) })
    }
}

/// An allele as answers spell it: whole when short, else its first bases and its length.
/// Symbolic alleles and breakends are never cut.
fn allele_text(allele: &str) -> String {
    if allele.len() <= MAX_ALLELE_CHARS || !coding::is_sequence(allele) {
        return allele.to_owned();
    }
    format!("{}… ({} bp)", &allele[..ALLELE_HEAD], allele.len())
}

fn variant_json(v: &ggv_io::Variant, samples: &[String], path: &Path) -> Value {
    let gts: serde_json::Map<String, Value> = samples
        .iter()
        .zip(&v.samples)
        .map(|(s, c)| (s.clone(), json!({"gt": c.gt, "dp": c.dp, "gq": c.gq, "ad": if c.ad.is_empty() { Value::Null } else { json!(c.ad) }})))
        .collect();
    let info: String = v.info.chars().take(200).collect();
    let alts: Vec<String> = v.alts.iter().map(|a| allele_text(a)).collect();
    let mut out = json!({
        "file": file_name(path),
        "path": path,
        "pos": v.start + 1, "id": if v.id.is_empty() || v.id == "." { Value::Null } else { json!(v.id) },
        "ref": allele_text(&v.reference), "alt": alts, "type": v.kind(),
        // QUAL is a 32-bit float in the file: printed as one, not widened into noise digits.
        "qual": v.qual.map(|q| round(f64::from(q), 2)), "filter": if v.filters.is_empty() { json!(".") } else { json!(v.filters.join(",")) },
        "af": v.af, "genotypes": gts, "info": info,
    });
    if v.reference.len() > MAX_ALLELE_CHARS || v.alts.iter().any(|a| a.len() > MAX_ALLELE_CHARS) {
        out["ref_length"] = json!(v.reference.len());
        out["alt_lengths"] = json!(v.alts.iter().map(|a| coding::is_sequence(a).then_some(a.len())).collect::<Vec<_>>());
    }
    if v.alts.len() > 1 {
        out["alt_types"] = json!(v.alts.iter().map(|a| v.alt_kind(a)).collect::<Vec<_>>());
        out["alt_af"] = json!(v.alt_af);
    }
    out
}

/// Why a VCF change can't be annotated as a sequence change (a breakend, a symbolic allele);
/// `None` when it can.
fn unannotatable(reference: &str, alt: &str) -> Option<&'static str> {
    if ggv_io::is_breakend(alt) {
        Some("a breakend")
    } else if !coding::is_sequence(alt) || !coding::is_sequence(reference) {
        Some("a symbolic allele")
    } else {
        None
    }
}

/// Each VCF sample's GT, by sample name; `null` for a file without samples.
fn sample_genotypes(v: &ggv_io::Variant, samples: &[String]) -> Value {
    if samples.is_empty() {
        return Value::Null;
    }
    Value::Object(samples.iter().zip(&v.samples).map(|(s, c)| (s.clone(), json!(c.gt))).collect())
}

/// An SV record from a VCF, as `svs` tests it.
struct SvRecord {
    site: ggv_query::SvSite,
    kind: String,
    id: Value,
    /// The VCF's genotypes, by sample (not the reads').
    genotypes: Value,
    filter: String,
    vcf: String,
    /// A breakend's ALT, which names its mate.
    alt: Option<String>,
    /// The record's POS.
    vcf_pos: u64,
    /// A breakend's mate record, tested with it.
    mate: Option<Value>,
}

impl SvRecord {
    /// Where it is, 1-based: a breakend's own POS; else the SV's first base.
    fn pos(&self) -> u64 {
        if self.kind == "BND" { self.vcf_pos } else { self.site.start + 1 }
    }

    /// Its last base, 1-based: a breakend's other breakpoint; none for an insertion.
    fn end(&self) -> Option<u64> {
        let s = &self.site;
        if s.kind == ggv_query::SvKind::Rearrangement && self.kind == "BND" { Some(s.end + 1) } else { (s.end > s.start).then_some(s.end) }
    }

    /// Its length; none for a breakend that isn't a deletion or duplication.
    fn length(&self) -> Value {
        if self.site.kind == ggv_query::SvKind::Rearrangement && self.kind == "BND" { Value::Null } else { json!(self.site.len) }
    }
}

/// Whether two breakends name the same junction (the two mates of a pair).
fn same_breakends(a: &ggv_query::SvSite, b: &ggv_query::SvSite) -> bool {
    a.kind == b.kind && a.start == b.start && a.end == b.end
}

/// Which records to test against the reads: the `max` longest, in position order.
fn pick_tested(records: &[SvRecord], max: usize) -> Vec<usize> {
    let mut order: Vec<usize> = (0..records.len()).collect();
    order.sort_by_key(|&i| std::cmp::Reverse(records[i].site.len));
    order.truncate(max);
    order.sort_by_key(|&i| records[i].site.start);
    order
}

/// Whether some record (tested or not) accounts for an event found in the reads: one of its
/// kind and size near it, or a tandem duplication the reads carry as an insertion in or beside
/// it.
fn explained(records: &[SvRecord], e: &ggv_query::SvSite) -> bool {
    use ggv_query::SvKind;
    records.iter().map(|r| &r.site).any(|r| {
        let w = e.slack().max(r.slack());
        let near = match (r.kind, e.kind) {
            (SvKind::Duplication, SvKind::Insertion) => e.start + w >= r.start && e.start <= r.end + w,
            (a, b) => a == b && r.start.abs_diff(e.start) <= w,
        };
        near && r.len.abs_diff(e.len) * 10 <= r.len.max(e.len) * 3
    })
}

/// A breakend whose mate is on the same contig, as the event its brackets describe (VCF 4.2,
/// 5.4): `t[p[` joins t to the bases from p on, so with p after t it skips those between (a
/// deletion), and before t it repeats them (a tandem duplication); `]p]t` is the mirror image;
/// `t]p]` and `[p[t` join opposite strands, tested as a rearrangement between the two
/// breakpoints. Both mates of a pair give the same site. Else why it isn't tested.
fn breakend_site(v: &ggv_io::Variant, alt: &str, contig: &str) -> std::result::Result<ggv_query::SvSite, String> {
    use ggv_query::{SvKind, SvSite};
    let Some((mate_contig, mate_pos)) = ggv_io::breakend_mate(alt) else {
        return Err("single breakend (no mate): read support is not tested".into());
    };
    if mate_contig != contig && !ggv_model::contig_aliases(contig).iter().any(|a| a == mate_contig) {
        return Err(format!("mate on another contig ({mate_contig}:{mate_pos}): read support is not tested"));
    }
    let (a, b) = (v.start, mate_pos.saturating_sub(1));
    let (lo, hi) = (a.min(b), a.max(b));
    // Where the mate's sequence goes: after this base (`t[p[`, `t]p]`) or before it.
    let mate_after = !alt.starts_with(['[', ']']);
    let same_strand = if mate_after { alt.ends_with('[') } else { alt.starts_with(']') };
    let skips = mate_after == (b > a);
    Ok(match (same_strand, skips) {
        _ if a == b => SvSite { kind: SvKind::Rearrangement, start: lo, end: hi, len: 0 },
        (true, true) => SvSite { kind: SvKind::Deletion, start: lo + 1, end: hi, len: hi - lo - 1 },
        (true, false) => SvSite { kind: SvKind::Duplication, start: lo, end: hi + 1, len: hi - lo + 1 },
        (false, _) => SvSite { kind: SvKind::Rearrangement, start: lo, end: hi, len: 0 },
    })
}

/// The SV ALT `i` describes, if it is one of at least `min_len`: its site and type name.
/// Sequence alleles give a deletion or insertion from the length difference (after the
/// shared first base); symbolic ones use INFO/END and that ALT's SVLEN.
fn sv_site(v: &ggv_io::Variant, i: usize, alt: &str, min_len: u64) -> Option<(ggv_query::SvSite, String)> {
    use ggv_query::{SvKind, SvSite};
    let svlen = v.svlen_of(i);
    let (kind, start, end, len, name) = if alt.starts_with('<') {
        let t = alt.trim_matches(['<', '>']).split(':').next().unwrap_or("");
        let (start, end) = (v.start + 1, v.end.max(v.start + 1));
        match t {
            "DEL" => (SvKind::Deletion, start, end, svlen.unwrap_or(end - start), "DEL"),
            "INS" => (SvKind::Insertion, start, start, svlen?, "INS"),
            "DUP" => (SvKind::Duplication, start, end, svlen.unwrap_or(end - start), "DUP"),
            "INV" => (SvKind::Rearrangement, start, end, svlen.unwrap_or(end - start), "INV"),
            _ => return None,
        }
    } else if ggv_annot::coding::is_sequence(alt) {
        let (r, a) = (v.reference.len() as u64, alt.len() as u64);
        if a > r {
            (SvKind::Insertion, v.start + 1, v.start + 1, a - r, "INS")
        } else if r > a {
            (SvKind::Deletion, v.start + a, v.start + r, r - a, "DEL")
        } else {
            return None;
        }
    } else {
        return None;
    };
    (len >= min_len).then(|| (SvSite { kind, start, end, len }, name.to_owned()))
}

/// `transcripts`: "main" (the default) or "all".
fn all_transcripts(args: &Value) -> Result<bool> {
    match str_arg(args, "transcripts") {
        None | Some("main") => Ok(false),
        Some("all") => Ok(true),
        Some(other) => bail!("transcripts: `{other}` is not main or all"),
    }
}

/// Whether two paths name the same file (compared canonically, so `a/../x.bam` and `x.bam` are
/// one file).
fn same_file(a: &Path, b: &Path) -> bool {
    a == b || std::fs::canonicalize(a).is_ok_and(|ca| std::fs::canonicalize(b).is_ok_and(|cb| ca == cb))
}

/// One transcript for each gene here (its representative: MANE Select, else canonical, else the
/// longest curated coding one; where that one doesn't reach, the gene's best one here: curated
/// over a predicted XM_/XR_ model, coding, the longest coding sequence) unless `all`, and how
/// many were left out.
fn pick_transcripts(ts: Vec<&ggv_annot::Transcript>, all: bool) -> (Vec<&ggv_annot::Transcript>, usize) {
    if all {
        return (ts, 0);
    }
    let out = ggv_annot::one_per_gene(ts.iter().copied());
    let hidden = ts.len() - out.len();
    (out, hidden)
}

/// Transcripts over `pos`, trying the annotation's own contig naming; coding and
/// representative ones first.
fn transcripts_at<'a>(ann: &'a Annotation, contig: &str, pos: u64) -> Vec<&'a ggv_annot::Transcript> {
    crate::batch::transcripts(ann, contig, pos)
}

/// Like [`transcripts_at`], for every transcript overlapping `[start, end)`.
fn transcripts_over<'a>(ann: &'a Annotation, contig: &str, start: u64, end: u64) -> Vec<&'a ggv_annot::Transcript> {
    crate::batch::transcripts_over(ann, contig, start, end)
}

/// A file's size and modification time (none when it can't be read).
type Stamp = Option<(u64, std::time::SystemTime)>;

fn stamp(path: &Path) -> Stamp {
    let m = std::fs::metadata(path).ok()?;
    Some((m.len(), m.modified().ok()?))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ggv_query::{SvKind, SvSite};

    fn variant(start: u64, end: u64, reference: &str, alts: &[&str], svlen: &[Option<i64>]) -> ggv_io::Variant {
        ggv_io::Variant {
            start,
            end,
            reference: reference.into(),
            alts: alts.iter().map(|a| (*a).to_owned()).collect(),
            svlen: svlen.to_vec(),
            ..Default::default()
        }
    }

    fn record(kind: SvKind, start: u64, len: u64) -> SvRecord {
        SvRecord { site: SvSite { kind, start, end: start + len, len }, kind: "DEL".into(), id: Value::Null, genotypes: Value::Null, filter: "PASS".into(), vcf: "x".into(), alt: None, vcf_pos: start, mate: None }
    }

    #[test]
    fn each_alt_has_its_own_svlen() {
        // ALT <DEL>,<INS> with SVLEN=-50,100.
        let v = variant(119, 170, "A", &["<DEL>", "<INS>"], &[Some(-50), Some(100)]);
        let (del, kind) = sv_site(&v, 0, "<DEL>", 20).unwrap();
        assert_eq!((del.len, kind.as_str()), (50, "DEL"));
        let (ins, kind) = sv_site(&v, 1, "<INS>", 70).unwrap();
        assert_eq!((ins.len, kind.as_str()), (100, "INS"));
        assert!(sv_site(&v, 0, "<DEL>", 70).is_none());
        // SVLEN with too few values: the insertion has no length, the deletion uses END.
        let v = variant(119, 170, "A", &["<DEL>", "<INS>"], &[Some(-50)]);
        assert!(sv_site(&v, 1, "<INS>", 20).is_none());
        assert_eq!(sv_site(&v, 0, "<DEL>", 20).unwrap().0.len, 50);
        // A <DUP> is tested as a tandem duplication, not a rearrangement.
        let v = variant(119, 705, "A", &["<DUP>"], &[Some(586)]);
        assert_eq!(sv_site(&v, 0, "<DUP>", 20).unwrap().0.kind, SvKind::Duplication);
    }

    #[test]
    fn discovery_matches_untested_records() {
        // max 1 tests the 100-base deletion; the 50-base one still explains reads.
        let records = vec![record(SvKind::Deletion, 120, 50), record(SvKind::Deletion, 600, 100)];
        assert_eq!(pick_tested(&records, 1), vec![1]);
        assert_eq!(pick_tested(&records, 5), vec![0, 1]);
        assert!(explained(&records, &SvSite { kind: SvKind::Deletion, start: 120, end: 120, len: 50 }));
        assert!(!explained(&records, &SvSite { kind: SvKind::Insertion, start: 120, end: 120, len: 50 }));
        // A tandem duplication explains the insertion of its size the reads carry for it.
        let dup = vec![record(SvKind::Duplication, 1000, 586)];
        assert!(explained(&dup, &SvSite { kind: SvKind::Insertion, start: 1001, end: 1001, len: 586 }));
        assert!(explained(&dup, &SvSite { kind: SvKind::Insertion, start: 1590, end: 1590, len: 570 }));
        assert!(!explained(&dup, &SvSite { kind: SvKind::Insertion, start: 1001, end: 1001, len: 100 }));
    }

    #[test]
    fn genotypes_are_keyed_by_sample() {
        let call = |gt: &str| ggv_io::SampleCall { class: ggv_io::GenotypeClass::Unknown, gt: gt.into(), dp: None, gq: None, ad: vec![] };
        let v = ggv_io::Variant { samples: vec![call("1/1"), call("0/0")], ..Default::default() };
        assert_eq!(sample_genotypes(&v, &["TUMOR".into(), "NORMAL".into()]), json!({"NORMAL": "0/0", "TUMOR": "1/1"}));
        assert_eq!(sample_genotypes(&ggv_io::Variant::default(), &[]), Value::Null);
    }

    #[test]
    fn breakends_are_not_sequence_changes() {
        assert_eq!(unannotatable("C", "C]chrQA:700]"), Some("a breakend"));
        assert_eq!(unannotatable("C", "<DEL>"), Some("a symbolic allele"));
        assert_eq!(unannotatable("C", "CA"), None);
    }

    #[test]
    fn breakends_are_tested_or_explained() {
        let v = variant(119, 120, "A", &["A]chrQA:700]"], &[]);
        let site = breakend_site(&v, "A]chrQA:700]", "chrQA").unwrap();
        assert_eq!((site.kind, site.start, site.end), (SvKind::Rearrangement, 119, 699));
        assert!(breakend_site(&v, "A]chr2:700]", "chrQA").unwrap_err().contains("another contig"));
        assert!(breakend_site(&v, "A.", "chrQA").unwrap_err().contains("single breakend"));
        // A deletion written as a breakend pair (Manta/GRIDSS): both mates give the deletion
        // of 121..699 (1-based), the same site as `<DEL>` with END=699.
        let at = |pos: u64, alt: &str| breakend_site(&variant(pos - 1, pos, "A", &[alt], &[]), alt, "chrQA").unwrap();
        let (one, two) = (at(120, "A[chrQA:700["), at(700, "]chrQA:120]A"));
        assert_eq!((one.kind, one.start, one.end, one.len), (SvKind::Deletion, 120, 699, 579));
        assert!(same_breakends(&one, &two));
        let del = sv_site(&variant(119, 699, "A", &["<DEL>"], &[Some(-579)]), 0, "<DEL>", 20).unwrap().0;
        assert_eq!((del.start, del.end, del.len), (one.start, one.end, one.len));
        // A tandem duplication of 120..700: the mate's sequence comes before (or after) the
        // other way round.
        let (one, two) = (at(120, "]chrQA:700]A"), at(700, "A[chrQA:120["));
        assert_eq!((one.kind, one.start, one.end, one.len), (SvKind::Duplication, 119, 700, 581));
        assert!(same_breakends(&one, &two));
        // Opposite strands: an inversion's junction, tested as a rearrangement.
        assert_eq!(at(120, "[chrQA:700[A").kind, SvKind::Rearrangement);
    }

    #[test]
    fn transcripts_left_out_favour_curated_ones_for_each_gene() {
        let tx = |id: &str, gene: &str, cds: u64, representative: bool| ggv_annot::Transcript {
            contig: 0,
            start: 0,
            end: 1000,
            strand: ggv_annot::Strand::Forward,
            id: id.into(),
            gene_id: gene.into(),
            gene_name: gene.into(),
            biotype: "".into(),
            exons: vec![(0, 1000)],
            cds: (cds > 0).then_some((0, cds)),
            tags: 0,
            representative,
            translation: Default::default(),
            exons_from_span: false,
        };
        // GNAS: the main transcript doesn't reach the base; a predicted model comes first.
        let (xm, nesp, short, other) = (tx("XM_017027815.2", "GNAS", 900, false), tx("NM_016592.5", "GNAS", 600, false), tx("NM_001309861.2", "GNAS", 300, false), tx("NM_1.1", "OTHER", 0, true));
        let (shown, hidden) = pick_transcripts(vec![&xm, &short, &nesp, &other], false);
        assert_eq!(shown.iter().map(|t| &*t.id).collect::<Vec<_>>(), ["NM_016592.5", "NM_1.1"]);
        assert_eq!(hidden, 2);
    }

    /// A raw BAM record at 0-based `pos` with this CIGAR (`(len, op)`) and bases.
    fn raw_record(pos: i32, cigar: &[(u32, u32)], seq: &[u8]) -> ggv_io::pasta::pileup::Read {
        let code = |b: u8| match b {
            b'A' => 1u8,
            b'C' => 2,
            b'G' => 4,
            b'T' => 8,
            _ => 15,
        };
        let mut v = Vec::new();
        v.extend(0i32.to_le_bytes());
        v.extend(pos.to_le_bytes());
        v.extend([2u8, 60]);
        v.extend(0u16.to_le_bytes());
        v.extend((cigar.len() as u16).to_le_bytes());
        v.extend(0u16.to_le_bytes());
        v.extend((seq.len() as u32).to_le_bytes());
        v.extend([(-1i32).to_le_bytes(), (-1i32).to_le_bytes(), 0i32.to_le_bytes()].concat());
        v.extend(b"r\0");
        for (len, op) in cigar {
            v.extend((len << 4 | op).to_le_bytes());
        }
        for pair in seq.chunks(2) {
            v.push(code(pair[0]) << 4 | pair.get(1).map_or(0, |b| code(*b)));
        }
        v.extend(std::iter::repeat_n(30u8, seq.len()));
        ggv_io::pasta::pileup::Read::new(v).unwrap()
    }

    #[test]
    fn reads_without_a_reference_show_their_own_base() {
        let mut seq = vec![b'A'; 20];
        seq[9] = b'C';
        seq[10] = b'N';
        let raw = raw_record(100, &[(20, 0)], &seq);
        let mut r = ggv_layout::decode(&raw, None);
        // Nothing to compare with, and no bases kept: unknown, not N.
        assert!(matches!(read_at(&r, 109, None), Some((Allele::Unknown, _))));
        r.bases = Some(ggv_layout::read_bases(&raw));
        assert!(matches!(read_at(&r, 109, None), Some((Allele::Base(b'C'), 9))));
        assert!(matches!(read_at(&r, 110, None), Some((Allele::Base(b'N'), _))));
        assert!(matches!(read_at(&r, 101, None), Some((Allele::Base(b'A'), 1))));
    }

    #[test]
    fn read_end_distance_counts_terminal_insertions() {
        let seq = [b'A'; 60];
        // 60 read bases each time; reference 100.. and the base at offset 29.
        for (cigar, want) in [(vec![(60, 0)], 29), (vec![(30, 0), (30, 1)], 29), (vec![(30, 0), (20, 1), (10, 4)], 29), (vec![(10, 4), (20, 1), (30, 0)], 1)] {
            let r = ggv_layout::decode(&raw_record(100, &cigar, &seq), Some((&[b'A'; 100], 100)));
            let at = if cigar[0].1 == 4 { 128 } else { 129 };
            let (_, d) = read_at(&r, at, Some(b'A')).unwrap();
            assert_eq!(d, want, "{cigar:?}");
        }
    }
}
