mod batch;
mod junctions;
mod mcp;
#[cfg(test)]
mod tests_cli;
mod transcripts;

use std::{
    io::{self, BufWriter, Write},
    path::PathBuf,
};

use anyhow::{Result, bail};
use clap::{Args, Parser, Subcommand, ValueEnum};
use ggv_io::{BamFile, ReferenceFasta, signal::bin_edges};
use ggv_model::Region;
use ggv_query::{PileupOptions, PositionStats, SummaryOptions};

/// `0.9.0 (commit abc123def456, x86_64-pc-windows-msvc)`: the exact build, for bug reports.
const VERSION: &str = concat!(env!("CARGO_PKG_VERSION"), " (commit ", env!("GGV_COMMIT"), ", ", env!("GGV_TARGET"), ")");

/// GGV — a genome viewer built for people and agents. Today: structured queries.
#[derive(Parser)]
#[command(name = "ggv", version = VERSION, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Built-in genomes: list them (and which are downloaded), or download one ahead of time.
    Genomes {
        /// Download this genome (hg38, GRCh38, T2T, mouse...) now.
        #[arg(long, value_name = "NAME")]
        fetch: Option<String>,
    },
    /// The demo (real PacBio HiFi reads of GIAB HG002, 4.3 MB): downloads it once to the cache and
    /// prints where its files are, and commands to try.
    Demo,
    /// Show a BAM's contigs, read groups and index.
    Info {
        bam: PathBuf,
        #[arg(long, value_enum, default_value_t = Format::Tsv)]
        format: Format,
    },
    /// Build the index a file lacks, next to it: BAM (.bai, or .csi when a contig is too long
    /// for one), CRAM (.crai), bgzipped VCF (.tbi), FASTA (.fai, and .gzi when bgzipped). The
    /// same indexes samtools and tabix write, and the desktop app offers to build.
    Index {
        #[arg(required = true, value_name = "FILE")]
        files: Vec<PathBuf>,
        /// Build the index again when the file already has one.
        #[arg(long)]
        force: bool,
    },
    /// Ask questions about alignments; answers are data, not pictures.
    #[command(subcommand)]
    Query(Query),
    /// Draw a region to PNG or SVG, headless (no GPU or display needed).
    Render(RenderArgs),
    /// A picture and an evidence report for every place in a list (BED, VCF, or loci), plus a
    /// contact sheet (index.html) and summary.tsv. Headless and parallel.
    ///
    /// Exit status: 0 when every place was written, 3 when some failed or a line of the list is
    /// not a place (summary.tsv says which and why), 1 when none was written or the batch
    /// couldn't run.
    Batch(BatchArgs),
    /// Serve genome data to AI agents over the Model Context Protocol (stdio). Files given here
    /// are the defaults for every tool call.
    ///
    /// With --attach, calls are about what the person sees in a running GGV desktop window: its
    /// files, outlined base, selection, view and read filters fill in whatever a call leaves
    /// out, and `view` / `navigate` read and move the window.
    Mcp {
        /// Attach to the GGV desktop window (the focused or last used one, or the one with this
        /// process id). Files given here are used when no window is open.
        #[arg(long, value_name = "PID", num_args = 0..=1, default_missing_value = "0")]
        attach: Option<u32>,
        /// Alignment file (BAM/CRAM); repeatable.
        #[arg(short = 'b', long = "bam", value_name = "FILE")]
        bams: Vec<PathBuf>,
        /// Indexed reference FASTA.
        #[arg(short = 'f', long)]
        reference: Option<PathBuf>,
        /// Indexed VCF/BCF; repeatable.
        #[arg(short = 'v', long = "vcf", value_name = "FILE")]
        vcfs: Vec<PathBuf>,
        /// A built-in genome (see `ggv genomes`) for the reference and genes, downloaded once if
        /// needed; `-f` / `-g` win.
        #[arg(long, value_name = "NAME")]
        genome: Option<String>,
        /// Signal track (bigWig, bedGraph, SEG); repeatable.
        #[arg(short = 's', long = "signal", value_name = "FILE")]
        signals: Vec<PathBuf>,
        /// Gene annotation (GTF, GFF3, BED); repeatable.
        #[arg(short = 'g', long = "genes", value_name = "FILE")]
        genes: Vec<PathBuf>,
    },
}

/// What a file is, from its (lowercased) name.
pub(crate) fn kind_of(name: &str) -> &'static str {
    let base = name.trim_end_matches(".gz").trim_end_matches(".bgz");
    if base.ends_with(".bam") || base.ends_with(".cram") {
        "alignments"
    } else if base.ends_with(".vcf") || base.ends_with(".bcf") {
        "variants"
    } else if [".gtf", ".gff", ".gff3", ".bed"].iter().any(|e| base.ends_with(e)) {
        "genes"
    } else if [".fa", ".fasta", ".fna", ".2bit"].iter().any(|e| base.ends_with(e)) {
        "reference"
    } else if ggv_io::signal::is_signal(&base) {
        "signal"
    } else {
        "unknown"
    }
}

#[derive(Args)]
struct BatchArgs {
    /// Places: BED, VCF, or one locus per line (`chr1:100`, `chr1:100-200`, text after it is a
    /// label); gzipped is fine. A zero-length BED interval (an insertion site) is the base
    /// before it. Lines that are not places are failed rows of summary.tsv, with their line
    /// number.
    list: PathBuf,
    /// Alignment file (BAM/CRAM), one panel each; repeat for more.
    #[arg(short = 'b', long = "bam", required = true, value_name = "FILE")]
    bams: Vec<PathBuf>,
    /// Indexed reference FASTA (mismatches, allele counts against the reference, codons).
    #[arg(short = 'f', long)]
    reference: Option<PathBuf>,
    /// Indexed VCF/BCF whose records are listed (and whose changes are annotated); repeatable.
    #[arg(short = 'v', long = "vcf", value_name = "FILE")]
    vcfs: Vec<PathBuf>,
    /// Gene annotation (GTF, GFF3 or BED) for consequences (`c.` / `p.` changes).
    #[arg(short = 'g', long)]
    genes: Option<PathBuf>,
    /// A built-in genome (see `ggv genomes`) for the reference and genes, downloaded once if
    /// needed; `-f` / `-g` win.
    #[arg(long, value_name = "NAME")]
    genome: Option<String>,
    /// Output folder (created if needed).
    #[arg(short, long)]
    out: PathBuf,
    #[command(flatten)]
    filters: Filters,
    /// Image width in pixels.
    #[arg(long, default_value_t = 1200)]
    width: u32,
    /// Pixel ratio for the PNGs.
    #[arg(long, default_value_t = 1.0)]
    scale: f32,
    /// For single bases: bases to show on each side.
    #[arg(long, default_value_t = 50)]
    flank: u64,
    /// Maximum read rows per panel.
    #[arg(long, default_value_t = 60)]
    max_rows: usize,
    /// Colour reads by strand, HP tag, pair orientation and insert size, or methylation.
    #[arg(long, value_enum, default_value_t = ColorByArg::None)]
    color_by: ColorByArg,
    /// Group reads by strand, HP tag, read group, alleles at heterozygous sites, or (at
    /// single-base places) the base each read shows.
    #[arg(long, value_enum, default_value_t = GroupArg::None)]
    group: GroupArg,
    /// Put each split read's parts (primary and supplementary) on one row, joined.
    #[arg(long)]
    link: bool,
    /// The tag for `--color-by tag` / `--group tag` (CB, PS, UB, MI...).
    #[arg(long, value_name = "XX", value_parser = parse_tag_arg)]
    tag: Option<[u8; 2]>,
    /// Places worked on at once (default: the number of cores).
    #[arg(short, long)]
    jobs: Option<usize>,
}

#[derive(Args)]
struct RenderArgs {
    #[command(flatten)]
    input: Input,
    #[command(flatten)]
    filters: Filters,
    /// Output file: `.png` or `.svg`.
    #[arg(short, long)]
    output: PathBuf,
    /// Image width in pixels.
    #[arg(long, default_value_t = 1200)]
    width: u32,
    /// Pixel ratio for PNG (2 = sharper, twice the size).
    #[arg(long, default_value_t = 1.0)]
    scale: f32,
    /// For a single-base locus: bases to show on each side; the base itself is outlined.
    #[arg(long, default_value_t = 50)]
    flank: u64,
    #[arg(long, value_enum, default_value_t = ColorByArg::None)]
    color_by: ColorByArg,
    /// Group reads into bands, as the app does: by strand, HP tag, read group, the alleles
    /// reads carry at heterozygous sites, or (for a single-base locus) the base each read shows.
    #[arg(long, value_enum, default_value_t = GroupArg::None)]
    group: GroupArg,
    /// Put each split read's parts (primary and supplementary) on one row, joined.
    #[arg(long)]
    link: bool,
    /// Draw splice junctions (RNA-seq) as a sashimi band under the coverage: an arc per intron
    /// with its read count, + strand above the line and − below.
    #[arg(long)]
    sashimi: bool,
    /// With --sashimi: leave out junctions with fewer reads.
    #[arg(long, default_value_t = 2)]
    min_junction_reads: u32,
    /// The tag for `--color-by tag` / `--group tag` (CB, PS, UB, MI...).
    #[arg(long, value_name = "XX", value_parser = parse_tag_arg)]
    tag: Option<[u8; 2]>,
    /// For a single-base locus: sort reads by what they show there.
    #[arg(long, value_enum)]
    sort: Option<SortArg>,
    /// Maximum read rows; reads past this are counted but not drawn.
    #[arg(long, default_value_t = 100)]
    max_rows: usize,
    /// Show reads only for regions up to this many bases; wider regions show coverage only (the
    /// app's default too).
    #[arg(long, default_value_t = 1_000_000)]
    read_window: u64,
}

#[derive(Clone, Copy, PartialEq, Eq, ValueEnum)]
enum GroupArg {
    None,
    Strand,
    Haplotype,
    ReadGroup,
    Alleles,
    Base,
    /// By the value of `--tag`.
    Tag,
}

#[derive(Clone, Copy, ValueEnum)]
enum SortArg {
    Base,
    Strand,
    Mapq,
    Start,
    InsertSize,
}

fn parse_tag_arg(s: &str) -> Result<[u8; 2], String> {
    ggv_layout::parse_tag_name(s).ok_or_else(|| format!("`{s}` is not a two-character tag name (CB, PS...)"))
}

/// The tag the reads must carry for these options (an error when one needs `--tag`).
fn tag_needed(color: ColorByArg, group: GroupArg, tag: Option<[u8; 2]>) -> Result<Option<[u8; 2]>> {
    if matches!(color, ColorByArg::Tag) || group == GroupArg::Tag {
        tag.map(Some).ok_or_else(|| anyhow::anyhow!("--color-by tag and --group tag need --tag (e.g. --tag CB)"))
    } else {
        Ok(None)
    }
}

fn color_of(c: ColorByArg) -> ggv_render::ColorBy {
    match c {
        ColorByArg::Tag => ggv_render::ColorBy::Tag,
        ColorByArg::None => ggv_render::ColorBy::None,
        ColorByArg::Strand => ggv_render::ColorBy::Strand,
        ColorByArg::Haplotype => ggv_render::ColorBy::Haplotype,
        ColorByArg::InsertSize => ggv_render::ColorBy::InsertSize,
        ColorByArg::BaseMods => ggv_render::ColorBy::BaseMods,
    }
}

/// The layout's grouping for `--group` (`base` needs the single base, `tag` the `--tag`).
fn group_of(g: GroupArg, base: Option<u64>, tag: Option<[u8; 2]>) -> Result<ggv_layout::GroupBy> {
    use ggv_layout::GroupBy;
    Ok(match g {
        GroupArg::Tag => GroupBy::Tag(tag.ok_or_else(|| anyhow::anyhow!("--group tag needs --tag (e.g. --tag CB)"))?),
        GroupArg::None => GroupBy::None,
        GroupArg::Strand => GroupBy::Strand,
        GroupArg::Haplotype => GroupBy::Haplotype,
        GroupArg::ReadGroup => GroupBy::ReadGroup,
        GroupArg::Alleles => GroupBy::Alleles,
        GroupArg::Base => GroupBy::BaseAt(base.ok_or_else(|| anyhow::anyhow!("--group base needs a single-base locus"))?),
    })
}

fn sort_of(s: SortArg) -> ggv_layout::SortBy {
    use ggv_layout::SortBy;
    match s {
        SortArg::Base => SortBy::Base,
        SortArg::Strand => SortBy::Strand,
        SortArg::Mapq => SortBy::MappingQuality,
        SortArg::Start => SortBy::Start,
        SortArg::InsertSize => SortBy::InsertSize,
    }
}

#[derive(Clone, Copy, ValueEnum)]
enum ColorByArg {
    None,
    Strand,
    /// By HP tag (haplotype 1 red, 2 blue).
    Haplotype,
    /// Pair orientation and insert size, as IGV.
    InsertSize,
    /// Methylation and other base modifications (MM/ML tags): red likely, blue not.
    BaseMods,
    /// By the value of `--tag`: one colour per value.
    Tag,
}

#[derive(Subcommand)]
enum Query {
    /// Per-position depth, alleles by strand, MAPQ, clipping.
    Pileup {
        #[command(flatten)]
        input: Input,
        #[command(flatten)]
        filters: Filters,
        /// Also print positions no read covers (a single base is always printed).
        #[arg(short = 'a', long)]
        all_positions: bool,
        /// Only print positions whose alt fraction is at least this.
        #[arg(long, value_name = "FRACTION", value_parser = parse_fraction)]
        min_alt_fraction: Option<f64>,
        #[arg(long, value_enum, default_value_t = Format::Tsv)]
        format: Format,
    },
    /// Coverage, breadth, read filtering and the most variable sites across a region.
    Summary {
        #[command(flatten)]
        input: Input,
        #[command(flatten)]
        filters: Filters,
        /// Minimum alt fraction for a variable site.
        #[arg(long, default_value_t = 0.2, value_parser = parse_fraction)]
        min_alt_fraction: f64,
        /// Minimum depth for a variable site.
        #[arg(long, default_value_t = 8)]
        min_depth: u32,
        /// The most variable sites to list.
        #[arg(long, default_value_t = 25)]
        max_sites: usize,
        #[arg(long, value_enum, default_value_t = Format::Tsv)]
        format: Format,
    },
    /// Splice junctions (RNA-seq) in a region, for one or more files side by side: each
    /// junction's motif and strand, whether the annotation knows it, and per file its reads and
    /// its share of the reads from its donor and into its acceptor.
    Junctions {
        /// Indexed BAM/CRAM files (one set of columns each), then the locus.
        #[arg(required = true, num_args = 1..)]
        bams: Vec<PathBuf>,
        /// Locus: `chr17:7,668,421-7,687,490` (1-based, inclusive), or a contig.
        locus: String,
        /// Indexed reference FASTA: for splice-site motifs and strands (and CRAM).
        #[arg(short = 'f', long)]
        reference: Option<PathBuf>,
        /// Gene annotation (GTF, GFF3 or BED): marks junctions known or novel.
        #[arg(short = 'g', long)]
        genes: Option<PathBuf>,
        /// Keep junctions with at least this many reads in some file.
        #[arg(long, default_value_t = 1)]
        min_reads: u32,
        /// Count a read only if it aligns at least this many bases on each side of the junction.
        #[arg(long, default_value_t = 1)]
        min_anchor: u32,
        /// List at most this many junctions (the most-read).
        #[arg(long, default_value_t = 10_000)]
        max: usize,
        #[command(flatten)]
        filters: Filters,
        #[arg(long, value_enum, default_value_t = Format::Tsv)]
        format: Format,
    },
    /// Values of a signal track (bigWig, bedGraph, SEG) across a region, in equal bins: mean,
    /// min and max over the bases with data.
    Signal {
        /// A .bw / .bigwig, .bedGraph / .bg (gzipped too) or .seg file.
        file: PathBuf,
        /// Locus: `chr7` or `chr7:55,191,000-55,192,000` (1-based, inclusive).
        locus: String,
        /// Number of bins.
        #[arg(long, default_value_t = 100)]
        bins: usize,
        #[arg(long, value_enum, default_value_t = Format::Tsv)]
        format: Format,
    },
}

#[derive(Args)]
struct Input {
    /// Indexed BAM (.bai or .csi next to it).
    bam: PathBuf,
    /// Locus: `chr7`, `chr7:55,191,822` or `chr7:55191800-55191900` (1-based, inclusive).
    locus: String,
    /// Indexed reference FASTA, for reference bases and alt fractions against the reference.
    #[arg(short = 'f', long)]
    reference: Option<PathBuf>,
    /// A built-in genome (hg38, hg19, hs1, mm39, mm10, sacCer3; see `ggv genomes`) for the
    /// reference, downloaded once if needed. `-f` wins.
    #[arg(long, value_name = "NAME")]
    genome: Option<String>,
}

#[derive(Args)]
struct Filters {
    /// Skip reads with any of these flags (decimal or 0x hex).
    #[arg(long, value_parser = parse_flags, default_value = "0x704")]
    exclude_flags: u16,
    /// Skip reads without all of these flags.
    #[arg(long, value_parser = parse_flags, default_value = "0")]
    require_flags: u16,
    #[arg(short = 'q', long, default_value_t = 0)]
    min_mapq: u8,
    #[arg(short = 'Q', long, default_value_t = 0)]
    min_base_quality: u8,
    /// Keep only reads matching this expression, e.g. "mapq >= 20 && hp == 1 && !dup" (the
    /// desktop app's filter language; fields are listed in its filter window).
    #[arg(long, value_name = "EXPR")]
    filter: Option<String>,
}

impl Filters {
    /// The options, with `reference` for filters that test mismatches (refused without one:
    /// with no reference every read would have none).
    fn options(&self, include_zero_depth: bool, reference: Option<&std::path::Path>) -> PileupOptions {
        PileupOptions {
            exclude_flags: self.exclude_flags,
            require_flags: self.require_flags,
            min_mapping_quality: self.min_mapq,
            min_base_quality: self.min_base_quality,
            include_zero_depth,
            read_filter: self.filter.as_deref().map(|text| {
                let expr = std::sync::Arc::new(ggv_layout::expr::Expr::parse(text).unwrap_or_else(|e| {
                    eprintln!("ggv: --filter: {e}");
                    std::process::exit(2);
                }));
                if expr.uses_mismatches() && reference.is_none() {
                    eprintln!("ggv: --filter: mismatches need a reference; give one with -f");
                    std::process::exit(2);
                }
                ggv_layout::expr::predicate(expr, text)
            }),
            reference: reference.map(ggv_query::ReadReference::from_path),
        }
    }
}

#[derive(Clone, Copy, ValueEnum)]
enum Format {
    /// Tab-separated, for people and shell pipelines.
    Tsv,
    /// One JSON document.
    Json,
    /// One JSON object per line, streamed.
    Jsonl,
}

/// A fraction from 0 to 1 (NaN and anything outside would silently match nothing).
fn parse_fraction(s: &str) -> Result<f64, String> {
    let f: f64 = s.trim().parse().map_err(|e| format!("{e}"))?;
    if !(0.0..=1.0).contains(&f) {
        return Err(format!("{s} is not a fraction from 0 to 1"));
    }
    Ok(f)
}

fn parse_flags(s: &str) -> Result<u16, String> {
    let r = match s.strip_prefix("0x").or_else(|| s.strip_prefix("0X")) {
        Some(hex) => u16::from_str_radix(hex, 16),
        None => s.parse(),
    };
    r.map_err(|e| format!("invalid flags `{s}`: {e}"))
}

fn main() {
    if let Err(e) = run(Cli::parse()) {
        // A closed pipe (`ggv ... | head`) is not an error.
        if e.chain().any(|c| {
            c.downcast_ref::<io::Error>()
                .is_some_and(|e| e.kind() == io::ErrorKind::BrokenPipe)
        }) {
            return;
        }
        eprintln!("ggv: {e:#}");
        std::process::exit(1);
    }
}

/// A built-in genome's files, downloading what is missing (progress on stderr).
fn genome_files(name: &str) -> Result<ggv_io::genomes::Local> {
    let g = ggv_io::genomes::find(name).ok_or_else(|| {
        let ids: Vec<&str> = ggv_io::genomes::GENOMES.iter().map(|g| g.id).collect();
        anyhow::anyhow!("unknown genome `{name}`; built in: {}", ids.join(", "))
    })?;
    if let Some(l) = ggv_io::genomes::cached(g) {
        return Ok(l);
    }
    eprintln!("ggv: downloading {} (about {} MB, once) to {}", g.name, g.size / 1_000_000, ggv_io::genomes::cache_dir().map_or_else(String::new, |d| d.display().to_string()));
    let mut last = std::time::Instant::now() - std::time::Duration::from_secs(5);
    let mut finished = String::new();
    let mb = |b: u64| format!("{:.1}", b as f64 / 1e6);
    let l = ggv_io::genomes::fetch(g, &mut |file, done, total| {
        let complete = total > 0 && done >= total;
        if complete && finished == file {
            return;
        }
        if last.elapsed().as_secs_f32() >= 2.0 || complete {
            last = std::time::Instant::now();
            if complete {
                finished = file.to_owned();
            }
            if total > 0 {
                eprintln!("ggv:   {file}: {} / {} MB", mb(done), mb(total));
            } else {
                eprintln!("ggv:   {file}: {} MB", mb(done));
            }
        }
    })?;
    Ok(l)
}

/// `--genome` for `ggv batch`: the reference and genes `-f` / `-g` don't give, fetched (a
/// download the first time) only when one of them is missing.
fn batch_genome(a: &mut BatchArgs, files: impl FnOnce(&str) -> Result<ggv_io::genomes::Local>) -> Result<()> {
    if let Some(name) = &a.genome
        && (a.reference.is_none() || a.genes.is_none())
    {
        let l = files(name)?;
        a.reference.get_or_insert(l.reference);
        a.genes.get_or_insert(l.genes);
    }
    Ok(())
}

fn run(mut cli: Cli) -> Result<()> {
    // A built-in genome stands in for -f (and, where genes are used, -g).
    let input_genome = |input: &mut Input| -> Result<()> {
        if let (Some(name), None) = (&input.genome, &input.reference) {
            input.reference = Some(genome_files(name)?.reference);
        }
        Ok(())
    };
    match &mut cli.command {
        Command::Render(a) => input_genome(&mut a.input)?,
        Command::Query(Query::Pileup { input, .. } | Query::Summary { input, .. }) => input_genome(input)?,
        Command::Batch(a) => batch_genome(a, genome_files)?,
        _ => {}
    }
    // Not locked: `ggv mcp` writes replies from several threads.
    let stdout = io::stdout();
    let mut out = BufWriter::new(stdout);

    match cli.command {
        Command::Genomes { fetch } => {
            if let Some(name) = fetch {
                let l = genome_files(&name)?;
                writeln!(out, "{}
{}", l.reference.display(), l.genes.display())?;
            } else {
                writeln!(out, "id	name	download	status")?;
                for g in ggv_io::genomes::GENOMES {
                    let status = match ggv_io::genomes::cached(g) {
                        Some(l) => format!("ready: {}", l.reference.parent().map_or_else(String::new, |d| d.display().to_string())),
                        None => "not downloaded (fetched on first use)".to_owned(),
                    };
                    writeln!(out, "{}	{}	{} MB	{status}", g.id, g.name, g.size / 1_000_000)?;
                }
            }
        }
        Command::Index { files, force } => {
            let mut failed = 0;
            for file in &files {
                match build_missing_index(file, force) {
                    Ok(Indexed::Written(index)) => writeln!(out, "wrote {}", index.display())?,
                    Ok(Indexed::AlreadyThere(index)) => writeln!(out, "{} is already there (--force builds it again)", index.display())?,
                    Err(e) => {
                        failed += 1;
                        eprintln!("ggv: {e:#}");
                    }
                }
            }
            if failed > 0 {
                anyhow::bail!("{failed} of {} files could not be indexed", files.len());
            }
        }
        Command::Demo => {
            let d = match ggv_io::demo::cached() {
                Some(d) => d,
                None => {
                    eprintln!("downloading the demo ({:.1} MB, once)", ggv_io::demo::SIZE as f64 / 1e6);
                    ggv_io::demo::fetch(&mut |_, _, _| {}).map_err(|e| anyhow::anyhow!("downloading the demo: {e}"))?
                }
            };
            writeln!(out, "demo files in {}", d.dir.display())?;
            for f in d.files() {
                writeln!(out, "  {}", f.display())?;
            }
            writeln!(out, "\ntry:")?;
            for c in demo_commands(&d) {
                writeln!(out, "  {c}")?;
            }
        }
        Command::Render(args) => render(&args)?,
        Command::Mcp { attach, bams, mut reference, vcfs, mut genes, genome, signals } => {
            if let Some(name) = &genome {
                let l = genome_files(name)?;
                reference.get_or_insert(l.reference);
                if genes.is_empty() {
                    genes.push(l.genes);
                }
            }
            let mut state = mcp::tools::State::new(mcp::tools::Defaults { bams, reference, vcfs, genes, signals });
            state.desktop = attach.map(|pid| mcp::desktop::Desktop::new((pid != 0).then_some(pid)));
            mcp::run(state)?
        }
        Command::Batch(args) => {
            let outcome = batch::run(&batch::Batch {
            list: {
                check_image_size(args.width, args.scale)?;
                args.list
            },
            bams: args.bams,
            reference: args.reference.clone(),
            vcfs: args.vcfs,
            genes: args.genes,
            out: args.out,
            width: args.width,
            scale: args.scale,
            flank: args.flank,
            max_rows: args.max_rows,
            jobs: args.jobs.unwrap_or_else(|| std::thread::available_parallelism().map_or(4, |n| n.get())),
            pileup: args.filters.options(false, args.reference.as_deref()),
            color_by: color_of(args.color_by),
            group: group_of(args.group, Some(u64::MAX), args.tag)?,
            tag: tag_needed(args.color_by, args.group, args.tag)?,
            link: args.link,
            })?;
            if outcome.written == 0 {
                anyhow::bail!("no place was written ({} failed; see summary.tsv)", outcome.failed);
            }
            if outcome.failed > 0 {
                out.flush()?;
                std::process::exit(3);
            }
        }

        Command::Info { bam, format } => {
            // The header needs no index: unaligned or unsorted reads (`hifi_reads.bam`) can't
            // have one. Region queries still require it.
            let info = match open_bam(&bam, None) {
                Ok(b) => b.info(),
                Err(e) if matches!(e.downcast_ref(), Some(ggv_io::Error::MissingIndex { .. })) => BamFile::header_info(&bam)?,
                Err(e) => return Err(e),
            };
            match format {
                Format::Json => {
                    serde_json::to_writer_pretty(&mut out, &info)?;
                    writeln!(out)?;
                }
                Format::Jsonl => {
                    serde_json::to_writer(&mut out, &info)?;
                    writeln!(out)?;
                }
                Format::Tsv => {
                    writeln!(out, "path\t{}\nindex\t{}", info.path, info.index)?;
                    for rg in &info.read_groups {
                        writeln!(
                            out,
                            "read_group\t{}\tsample={}\tplatform={}",
                            rg.id,
                            rg.sample.as_deref().unwrap_or("-"),
                            rg.platform.as_deref().unwrap_or("-")
                        )?;
                    }
                    for c in &info.contigs {
                        writeln!(out, "contig\t{}\t{}", c.name, c.length)?;
                    }
                }
            }
        }

        Command::Query(Query::Pileup { input, filters, all_positions, min_alt_fraction, format }) => {
            let (mut bam, region, reference) = open(&input)?;
            // One base asked for by name is answered even when no read covers it: a row of
            // zeros says so, where no row at all reads as a failed query.
            let zero_depth = (all_positions || region.len() == 1) && min_alt_fraction.is_none();
            let opts = filters.options(zero_depth, input.reference.as_deref());
            let keep = |s: &PositionStats| min_alt_fraction.is_none_or(|f| s.exact_alt_fraction() >= f);

            match format {
                Format::Tsv => {
                    writeln!(out, "{}", TSV_HEADER)?;
                    ggv_query::pileup(&mut bam, &region, reference.as_deref(), &opts, |s| {
                        if keep(&s) { write_tsv(&mut out, &s) } else { Ok(()) }
                    })?;
                }
                Format::Jsonl => {
                    ggv_query::pileup(&mut bam, &region, reference.as_deref(), &opts, |s| {
                        if keep(&s) {
                            serde_json::to_writer(&mut out, &s)?;
                            writeln!(out)
                        } else {
                            Ok(())
                        }
                    })?;
                }
                Format::Json => {
                    let mut positions = Vec::new();
                    let reads =
                        ggv_query::pileup(&mut bam, &region, reference.as_deref(), &opts, |s| {
                            if keep(&s) {
                                positions.push(s);
                            }
                            Ok(())
                        })?;
                    #[derive(serde::Serialize)]
                    struct Doc {
                        region: String,
                        reads: ggv_query::ReadCounts,
                        positions: Vec<PositionStats>,
                    }
                    let doc = Doc { region: region.to_string(), reads, positions };
                    serde_json::to_writer_pretty(&mut out, &doc)?;
                    writeln!(out)?;
                }
            }
        }

        Command::Query(Query::Signal { file, locus, bins, format }) => {
            anyhow::ensure!((1..=100_000).contains(&bins), "--bins must be 1-100000");
            let mut sig = ggv_io::signal::Signal::open(&file)?;
            if let Some(w) = sig.warning() {
                eprintln!("ggv: warning: {}: {w}", file.display());
            }
            let region = sig.resolve(&locus)?;
            let edges = bin_edges(region.start, region.end, bins);
            let bins = edges.len();
            let rows = sig.rows();
            let mut json_rows = Vec::new();
            if matches!(format, Format::Tsv) {
                writeln!(out, "#row\tcontig\tstart\tend\tmean\tmin\tmax\tcovered")?;
            }
            for (k, row) in rows.iter().enumerate() {
                let bs = sig.bins(k, &region.contig, region.start, region.end, bins)?;
                for (b, &(s, e)) in bs.iter().zip(&edges) {
                    match format {
                        Format::Tsv => match b {
                            Some(b) => writeln!(out, "{row}\t{}\t{}\t{e}\t{:.4}\t{:.4}\t{:.4}\t{}", region.contig, s + 1, b.mean, b.min, b.max, b.covered)?,
                            None => writeln!(out, "{row}\t{}\t{}\t{e}\t.\t.\t.\t0", region.contig, s + 1)?,
                        },
                        Format::Jsonl => writeln!(out, "{}", serde_json::json!({"row": row, "contig": region.contig, "start": s + 1, "end": e, "mean": b.map(|b| b.mean), "min": b.map(|b| b.min), "max": b.map(|b| b.max), "covered": b.map_or(0, |b| b.covered)}))?,
                        Format::Json => {}
                    }
                }
                json_rows.push(serde_json::json!({"row": row, "mean": bs.iter().map(|b| b.map(|b| (b.mean * 1e4).round() / 1e4)).collect::<Vec<_>>()}));
            }
            if matches!(format, Format::Json) {
                serde_json::to_writer_pretty(&mut out, &serde_json::json!({"region": format!("{}:{}-{}", region.contig, region.start + 1, region.end), "bins": bins, "rows": json_rows}))?;
                writeln!(out)?;
            }
        }

        Command::Query(Query::Junctions { bams, locus, reference, genes, min_reads, min_anchor, max, filters, format }) => {
            let opts = ggv_query::JunctionOptions { reads: filters.options(false, reference.as_deref()), min_anchor };
            let mut fasta = reference.as_ref().map(ReferenceFasta::open).transpose()?;
            let mut files = Vec::new();
            let mut region = None;
            for path in &bams {
                let mut bam = open_bam(path, reference.as_deref())?;
                let r = bam.resolve_str(&locus)?;
                // A reference without the contig is refused, as by pileup and render: else
                // motifs were all "unknown" and `mismatches` filters counted none, silently.
                if let Some(f) = &mut fasta {
                    f.fetch(&Region { contig: r.contig.clone(), start: r.start, end: r.start + 1 })?;
                }
                files.push(ggv_query::junctions(&mut bam, &r, &opts)?);
                region.get_or_insert(r);
            }
            let region = region.expect("at least one file");
            let fasta = fasta.map(std::cell::RefCell::new);
            let base = |p: u64| {
                let r = Region { contig: region.contig.clone(), start: p, end: p + 1 };
                fasta.as_ref()?.borrow_mut().fetch(&r).ok()?.first().copied()
            };
            let ann = genes.as_ref().map(ggv_annot::Annotation::open).transpose()?;
            let rep = junctions::report(&region.contig, &files, base, ann.as_ref(), min_reads, max);
            let names: Vec<String> = bams.iter().map(|p| p.file_name().map_or_else(|| p.display().to_string(), |n| n.to_string_lossy().into_owned())).collect();
            match format {
                Format::Json | Format::Jsonl => {
                    let v = serde_json::json!({
                        "region": region.to_string(),
                        "files": names.iter().zip(&files).map(|(n, f)| serde_json::json!({"file": n, "reads": f.reads, "spliced_reads": f.spliced_reads})).collect::<Vec<_>>(),
                        "total": rep.total,
                        "junctions": rep.junctions,
                    });
                    if matches!(format, Format::Json) { serde_json::to_writer_pretty(&mut out, &v)? } else { serde_json::to_writer(&mut out, &v)? }
                    writeln!(out)?;
                }
                Format::Tsv => {
                    write!(out, "#contig\tstart\tend\tlength\tstrand\tmotif\tannotation\tgenes")?;
                    for n in &names {
                        write!(out, "\t{n}:reads\t{n}:donor_share\t{n}:acceptor_share")?;
                    }
                    writeln!(out)?;
                    let share = |x: Option<f64>| x.map_or_else(|| ".".to_owned(), |v| v.to_string());
                    for j in &rep.junctions {
                        let genes = if j.genes.is_empty() { ".".to_owned() } else { j.genes.join(",") };
                        write!(out, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", region.contig, j.start, j.end, j.length, j.strand, j.motif, j.annotation, genes)?;
                        for c in &j.samples {
                            write!(out, "\t{}\t{}\t{}", c.reads, share(c.donor_share), share(c.acceptor_share))?;
                        }
                        writeln!(out)?;
                    }
                }
            }
        }
        Command::Query(Query::Summary {
            input,
            filters,
            min_alt_fraction,
            min_depth,
            max_sites,
            format,
        }) => {
            let (mut bam, region, reference) = open(&input)?;
            let opts = SummaryOptions {
                pileup: filters.options(true, input.reference.as_deref()),
                min_alt_fraction,
                min_variable_depth: min_depth,
                max_variable_sites: max_sites,
            };
            let summary = ggv_query::summarize(&mut bam, &region, reference.as_deref(), &opts)?;

            match format {
                Format::Json => {
                    serde_json::to_writer_pretty(&mut out, &summary)?;
                    writeln!(out)?;
                }
                // One object, on one line.
                Format::Jsonl => {
                    serde_json::to_writer(&mut out, &summary)?;
                    writeln!(out)?;
                }
                Format::Tsv => {
                    writeln!(out, "region\t{}\nlength\t{}", summary.region, summary.length)?;
                    writeln!(out, "{}", reads_line(&summary.reads))?;
                    writeln!(
                        out,
                        "depth\tmean={} median={} min={} max={}",
                        summary.mean_depth, summary.median_depth, summary.min_depth, summary.max_depth
                    )?;
                    let breadth: Vec<String> =
                        summary.breadth.iter().map(|(k, v)| format!("{k}:{v}")).collect();
                    writeln!(out, "breadth\t{}", breadth.join(" "))?;
                    writeln!(
                        out,
                        "variable_sites\t{} (showing {})",
                        summary.variable_site_count,
                        summary.variable_sites.len()
                    )?;
                    if !summary.variable_sites.is_empty() {
                        writeln!(out, "#pos\tref\tdepth\talt\talt_count\talt_fwd\talt_rev\talt_fraction\tlength")?;
                        for v in &summary.variable_sites {
                            writeln!(
                                out,
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                v.pos,
                                v.ref_base.unwrap_or('.'),
                                v.depth,
                                v.alt,
                                v.alt_count,
                                v.alt_fwd,
                                v.alt_rev,
                                v.alt_fraction,
                                v.length
                            )?;
                        }
                    }
                }
            }
        }
    }

    out.flush()?;
    Ok(())
}

/// `query summary`'s `reads` line: every count, including reads the `--filter` removed.
fn reads_line(r: &ggv_query::ReadCounts) -> String {
    format!(
        "reads\tseen={} used={} duplicate={} secondary={} supplementary={} qc_fail={} unmapped={} low_mapq={} other={} filtered={}",
        r.seen, r.used, r.duplicate, r.secondary, r.supplementary, r.qc_fail, r.unmapped, r.low_mapping_quality, r.other_flags, r.filtered_out
    )
}

/// Image sizes: a width of at least 100 px, and a pixel ratio above 0 (and sane).
fn check_image_size(width: u32, scale: f32) -> Result<()> {
    if width < 100 {
        bail!("--width must be at least 100 pixels (got {width})");
    }
    if !(scale.is_finite() && scale > 0.0 && scale <= 8.0) {
        bail!("--scale must be above 0 and at most 8 (got {scale})");
    }
    Ok(())
}

fn render(args: &RenderArgs) -> Result<()> {
    // Only formats this writes: `-o view.txt` used to get PNG bytes under a text name.
    let ext = args.output.extension().map(|e| e.to_string_lossy().to_ascii_lowercase());
    anyhow::ensure!(
        matches!(ext.as_deref(), Some("png" | "svg")),
        "{}: give a .png or .svg file name",
        args.output.display()
    );
    check_image_size(args.width, args.scale)?;
    let mut bam = open_bam(&args.input.bam, args.input.reference.as_deref())?;
    let mut region = bam.resolve_str(&args.input.locus)?;

    // A single base is shown in context, outlined, the way IGV centres on a position.
    let mut highlight = None;
    if region.len() == 1 {
        highlight = Some((region.start, region.end));
        let contig_len = bam
            .contigs()
            .iter()
            .find(|c| c.name == region.contig)
            .map_or(region.end, |c| c.length);
        region.start = region.start.saturating_sub(args.flank);
        region.end = (region.end + args.flank).min(contig_len);
    }

    let reference = match &args.input.reference {
        Some(path) => Some(ReferenceFasta::open(path)?.fetch(&region)?),
        None => None,
    };

    let info = bam.info();
    let label = info
        .read_groups
        .iter()
        .find_map(|rg| rg.sample.clone())
        .unwrap_or_else(|| {
            args.input.bam.file_name().unwrap_or_default().to_string_lossy().into_owned()
        });

    let view_opts = ggv_layout::ViewOptions {
        filters: args.filters.options(false, args.input.reference.as_deref()),
        bins: args.width as usize,
        max_read_window: args.read_window,
        max_rows: args.max_rows,
        tag: tag_needed(args.color_by, args.group, args.tag)?,
        ..Default::default()
    };
    let mut view = ggv_layout::load_view(&mut bam, &region, reference, label, &view_opts)?;
    let base = highlight.map(|(s, _)| s);
    let sort = match (args.sort, base) {
        (Some(by), Some(pos)) => Some((pos, sort_of(by))),
        (Some(_), None) => bail!("--sort needs a single-base locus"),
        (None, _) => None,
    };
    if args.group != GroupArg::None || sort.is_some() || args.link {
        ggv_layout::repack(&mut view, args.max_rows, group_of(args.group, base, args.tag)?, sort, args.link);
    }

    let render_opts = ggv_render::RenderOptions {
        width: args.width,
        color_by: color_of(args.color_by),
        link: args.link,
        highlight,
        sashimi: args.sashimi,
        min_junction_reads: args.min_junction_reads,
        ..Default::default()
    };
    let scene = ggv_render::build_scene(&view, &render_opts);

    let is_svg = args
        .output
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("svg"));
    let bytes = if is_svg {
        scene.to_svg().into_bytes()
    } else {
        scene.to_png(args.scale).map_err(anyhow::Error::msg)?
    };
    // Report the size actually written: the PNG's own header, or the SVG's.
    let (w, h) = if is_svg || bytes.len() < 24 {
        (scene.width, scene.height)
    } else {
        (u32::from_be_bytes(bytes[16..20].try_into().expect("4 bytes")), u32::from_be_bytes(bytes[20..24].try_into().expect("4 bytes")))
    };
    std::fs::write(&args.output, &bytes)
        .map_err(|e| anyhow::anyhow!("{}: {e}", args.output.display()))?;

    eprintln!("wrote {} ({w}x{h}, {})", args.output.display(), region);
    Ok(())
}

/// Opens an alignment file, warning on stderr when it looks cut short.
/// What `ggv index` did for one file.
enum Indexed {
    Written(PathBuf),
    AlreadyThere(PathBuf),
}

/// Builds `file`'s index unless it has one (`force` builds it regardless).
fn build_missing_index(file: &std::path::Path, force: bool) -> Result<Indexed> {
    let kind = ggv_io::index_kind(file).ok_or_else(|| anyhow::anyhow!("{}: ggv indexes BAM, CRAM, bgzipped VCF (.vcf.gz) and FASTA files, told by their names", file.display()))?;
    if !file.is_file() {
        anyhow::bail!("{}: no such file", file.display());
    }
    // Where an index of this kind may already sit: `x.bam.bai`, a `.csi` in its place, `x.bai`.
    let mut places = vec![file.with_added_extension(kind), file.with_extension(kind)];
    if matches!(kind, "bai" | "tbi") {
        places.push(file.with_added_extension("csi"));
    }
    if !force && let Some(found) = places.into_iter().find(|p| p.is_file()) {
        return Ok(Indexed::AlreadyThere(found));
    }
    Ok(Indexed::Written(ggv_io::build_index(file)?))
}

pub(crate) fn open_bam(path: &std::path::Path, reference: Option<&std::path::Path>) -> Result<BamFile> {
    let bam = BamFile::open_with_reference(path, reference)?;
    if let Some(w) = bam.warning()
        && first_time(&w)
    {
        eprintln!("ggv: warning: {w}");
    }
    Ok(bam)
}

/// Whether this warning is new to the process: batch workers each open every file, and a file
/// may be given twice, but one warning per file is enough.
fn first_time(warning: &str) -> bool {
    static SEEN: std::sync::Mutex<std::collections::BTreeSet<String>> = std::sync::Mutex::new(std::collections::BTreeSet::new());
    SEEN.lock().unwrap_or_else(std::sync::PoisonError::into_inner).insert(warning.to_owned())
}

/// The commands `ggv demo` suggests, as typed in a shell (paths with spaces quoted).
fn demo_commands(d: &ggv_io::demo::Demo) -> Vec<String> {
    // Quoted when anything in it would break a pasted command (spaces, `'`, `&`, parentheses as
    // in "Program Files (x86)"); double quotes work in cmd, PowerShell and sh alike.
    let q = |p: &std::path::Path| {
        let s = p.display().to_string();
        let plain = |c: char| c.is_ascii_alphanumeric() || "_-./:\\".contains(c);
        if s.chars().all(plain) { s } else { format!("\"{s}\"") }
    };
    let vcfs: Vec<String> = d.vcfs.iter().map(|v| format!("-v {}", q(v))).collect();
    let (bam, fa, genes, vcfs) = (q(&d.bam), q(&d.reference), q(&d.genes), vcfs.join(" "));
    vec![
        "ggv-desktop --demo".into(),
        format!("ggv-desktop {bam} -f {fa} {vcfs} -g {genes} {} --group haplotype", ggv_io::demo::STOPS[0].0),
        format!("ggv query pileup {bam} {} -f {fa}", ggv_io::demo::STOPS[1].0.replace(',', "")),
        // `ggv mcp` takes alignments with `-b` only; a bare path was rejected.
        format!("ggv mcp -b {bam} -f {fa} {vcfs} -g {genes}"),
    ]
}

fn open(input: &Input) -> Result<(BamFile, Region, Option<Vec<u8>>)> {
    let bam = open_bam(&input.bam, input.reference.as_deref())?;
    let region = bam.resolve_str(&input.locus)?;
    if region.is_empty() {
        bail!("{} is empty", region);
    }
    let reference = match &input.reference {
        Some(path) => Some(ReferenceFasta::open(path)?.fetch(&region)?),
        None => None,
    };
    Ok((bam, region, reference))
}

const TSV_HEADER: &str = "#contig\tpos\tref\tdepth\tA\tC\tG\tT\tN\tdel\tins\tref_skip\tlow_bq\tfwd\trev\tmean_mapq\tmapq0\tclip_left\tclip_right\talt_fraction\ttop_alt";

fn write_tsv(out: &mut impl Write, s: &PositionStats) -> io::Result<()> {
    let b = &s.bases;
    let fwd = b.a.fwd + b.c.fwd + b.g.fwd + b.t.fwd + b.n.fwd + s.deletions.fwd;
    let rev = b.a.rev + b.c.rev + b.g.rev + b.t.rev + b.n.rev + s.deletions.rev;
    let top = s
        .top_alt
        .as_ref()
        .map_or_else(|| ".".to_owned(), |a| format!("{}:{}", a.allele, a.count.total));
    writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        s.contig,
        s.pos,
        s.ref_base.unwrap_or('.'),
        s.depth,
        b.a.total,
        b.c.total,
        b.g.total,
        b.t.total,
        b.n.total,
        s.deletions.total,
        s.insertions_after.total,
        s.ref_skips,
        s.low_base_quality,
        fwd,
        rev,
        s.mean_mapping_quality,
        s.mapping_quality_zero,
        s.soft_clipped_left,
        s.soft_clipped_right,
        s.alt_fraction,
        top
    )
}

#[cfg(test)]
mod tests {
    #[test]
    fn alt_fractions_must_be_fractions() {
        assert_eq!(super::parse_fraction("0.2"), Ok(0.2));
        assert_eq!(super::parse_fraction("1"), Ok(1.0));
        for bad in ["NaN", "2", "-0.1", "inf", "x"] {
            assert!(super::parse_fraction(bad).is_err(), "{bad}");
        }
    }
}
