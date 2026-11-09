// A windowed app on Windows: opening a file from Explorer shows no console window. Started from a
// terminal, it re-attaches to that terminal's console so --help, --bench and errors still print.
#![cfg_attr(all(windows, not(debug_assertions)), windows_subsystem = "windows")]
//! ggv-desktop: an alignment viewer that stays fast however much is loaded.
//!
//! Reads, coverage, mismatches and indels are uploaded to the GPU once per load as instances
//! in genomic coordinates; panning and zooming only change a uniform. Loading happens on a
//! background thread, a few view widths at a time, so panning rarely touches disk.

mod app;
mod attach;
mod automation;
mod consensus;
mod bookmarks;
mod diff;
mod files;
mod finder;
mod genes;
mod gpu;
mod links;
mod loader;
mod overlays;
mod queue;
mod review;
mod seqsearch;
mod session;
mod store;
mod shell;
mod signals;
mod track;
mod variants;
mod tracks;
mod util;
mod watch;

/// The filter expression language (shared with `ggv`).
use ggv_layout::expr;

use std::path::PathBuf;

/// The system allocator on macOS keeps freed memory mapped and counted against the app: after a
/// few deep spans are decoded and dropped it reported 10 GB with a few hundred MB live. mimalloc
/// returns free pages to the system shortly after they're freed. (Its pages carry VM tag 100,
/// which `footprint` and `vmmap` list as "IOAccelerator": that is heap, not GPU memory.)
#[global_allocator]
static ALLOC: mimalloc::MiMalloc = mimalloc::MiMalloc;

use clap::{Parser, ValueEnum};
use files::{Kind, kind_of};
use eframe::egui;

/// `0.9.0 (commit abc123def456, x86_64-pc-windows-msvc)`: the exact build, for bug reports.
const VERSION: &str = concat!(env!("CARGO_PKG_VERSION"), " (commit ", env!("GGV_COMMIT"), ", ", env!("GGV_TARGET"), ")");

/// GGV desktop viewer.
#[derive(Parser)]
#[command(version = VERSION, about)]
struct Cli {
    /// Indexed BAMs to open, one track each, then optionally a locus
    /// (`ggv-desktop a.bam b.bam chr20:1,895,892`).
    #[arg(value_name = "BAM... [LOCUS]")]
    inputs: Vec<String>,
    /// Locus to show, e.g. `chr20:1,895,892` or `chr20:1,890,000-1,940,000`.
    #[arg(short, long)]
    locus: Option<String>,
    /// Indexed reference FASTA.
    #[arg(short = 'f', long)]
    reference: Option<PathBuf>,
    /// A built-in genome: hg38, hg19, hs1 (T2T), mm39, mm10, sacCer3 (also GRCh38, T2T, mouse...).
    /// Opens its reference and RefSeq genes, downloading them once if needed.
    #[arg(long, value_name = "NAME")]
    genome: Option<String>,
    /// A built-in track for the genome's build: gencode, gencode-basic, mane, refseq, clinvar,
    /// alphamissense, revel (what the Tracks menu lists); repeatable. Downloaded once if needed.
    #[arg(long = "track", value_name = "NAME")]
    tracks: Vec<String>,
    /// Indexed VCF/BCF to show as a variant track; repeatable.
    #[arg(short = 'v', long = "vcf", value_name = "FILE")]
    vcfs: Vec<PathBuf>,
    /// Gene annotation to show (GTF, GFF3 or BED, optionally gzipped); repeatable.
    #[arg(short = 'g', long = "genes", value_name = "FILE")]
    genes: Vec<PathBuf>,
    /// Open a saved session (from Session > Save as…).
    #[arg(long, value_name = "FILE")]
    session: Option<PathBuf>,
    /// Reopen the last session (saved automatically on exit).
    #[arg(long)]
    restore: bool,
    /// Colour reads by.
    #[arg(long, value_enum, default_value_t = ColorArg::None)]
    color: ColorArg,
    /// Group reads by.
    #[arg(long, value_enum, default_value_t = GroupArg::None)]
    group: GroupArg,
    /// Sort reads at the locus (a single base) by.
    #[arg(long, value_enum)]
    sort: Option<SortArg>,
    /// Select a region for statistics, e.g. `chr17:43,090,000-43,100,000`.
    #[arg(long, value_name = "REGION")]
    select: Option<String>,
    /// Outline every alignment of the read with this name (mate, split parts).
    #[arg(long, value_name = "NAME")]
    select_read: Option<String>,
    /// Start reviewing the first VCF (variant triage).
    #[arg(long)]
    review: bool,
    /// Open the demo (real HG002 long reads, 4.3 MB, downloaded once) with a short tour.
    #[arg(long)]
    demo: bool,
    /// Open the evidence panel (exact counts at the outlined base).
    #[arg(long)]
    evidence: bool,
    /// Dark theme.
    #[arg(long)]
    dark: bool,
    /// Show soft-clipped bases.
    #[arg(long)]
    soft_clips: bool,
    /// Start with every transcript shown, not one per gene.
    #[arg(long)]
    all_transcripts: bool,
    /// Save a screenshot of the window once data has loaded, then exit.
    #[arg(long, value_name = "PNG")]
    screenshot: Option<PathBuf>,
    /// Drive panning and zooming, print frame-time statistics, then exit (vsync off).
    #[arg(long)]
    bench: bool,
    /// Zoom all the way out and back in, then print each track's load state (and screenshot if
    /// --screenshot is given). For reproducing zoom-related bugs.
    #[arg(long)]
    zoom_cycle: bool,
    /// Open another pane beside the first at this locus (split view); repeat for more panes.
    #[arg(long, value_name = "LOCUS")]
    split: Vec<String>,
    /// Scroll and zoom split panes together.
    #[arg(long)]
    link_panes: bool,
    /// Step through a list of places (BED, VCF, or one locus per line) and mark each one.
    #[arg(long, value_name = "FILE")]
    queue: Option<PathBuf>,
    /// Filter reads with an expression, e.g. "mapq >= 20 && hp == 1 && !dup".
    #[arg(long, value_name = "EXPR")]
    filter: Option<String>,
    /// Find a sequence in the reference, both strands; append /N to allow N mismatches.
    #[arg(long, value_name = "SEQ[/N]")]
    find_seq: Option<String>,
    /// Find a read by name in the open files.
    #[arg(long, value_name = "NAME")]
    find_read: Option<String>,
    /// Compare the first two alignment tracks (differences in allele fraction).
    #[arg(long)]
    compare: bool,
    /// Outline this base (e.g. 19:11,134,238) without changing the view the locus sets.
    #[arg(long, value_name = "LOCUS")]
    highlight: Option<String>,
    /// Write an evidence report (HTML) for the outlined base once loaded.
    #[arg(long, value_name = "HTML")]
    report: Option<PathBuf>,
    /// Open the command palette with this text (for screenshots).
    #[arg(long, hide = true, value_name = "TEXT")]
    palette: Option<String>,
    /// Show the keyboard shortcuts (for screenshots).
    #[arg(long, hide = true)]
    shortcuts: bool,
    /// Show About GGV and its diagnostics (for screenshots).
    #[arg(long, hide = true)]
    about: bool,
    /// Show reads for views up to this many bases (default 1000000, or what was last chosen in
    /// the app); wider views show coverage only.
    #[arg(long, value_name = "BP")]
    read_window: Option<u64>,
    /// Measure this region (for screenshots).
    #[arg(long, hide = true, value_name = "REGION")]
    measure: Option<String>,
}

#[derive(Clone, Copy, ValueEnum)]
enum ColorArg {
    None,
    Strand,
    Haplotype,
    InsertSize,
    BaseMods,
}

#[derive(Clone, Copy, ValueEnum)]
enum GroupArg {
    None,
    Strand,
    Haplotype,
    ReadGroup,
    /// The base each read shows at the locus (a single base).
    Base,
    /// Haplotype by the alleles at heterozygous sites in view (no HP tags needed).
    Alleles,
}

#[derive(Clone, Copy, ValueEnum)]
enum SortArg {
    Base,
    Strand,
    Mapq,
    Start,
    InsertSize,
}

struct Inputs {
    bams: Vec<PathBuf>,
    vcfs: Vec<PathBuf>,
    genes: Vec<PathBuf>,
    signals: Vec<PathBuf>,
    locus: Option<String>,
}

/// Sorts positional inputs by extension (BAM/CRAM, VCF/BCF, GTF/GFF/BED); the last input is a
/// locus when it is not an existing file and has none of those extensions.
fn split_inputs(mut inputs: Vec<String>, locus: Option<String>) -> Inputs {
    let looks_like_locus = |s: &str| !std::path::Path::new(s).exists() && kind_of(s).is_none();
    let trailing = (locus.is_none() && inputs.last().is_some_and(|s| looks_like_locus(s))).then(|| inputs.pop()).flatten();
    let mut out = Inputs { bams: Vec::new(), vcfs: Vec::new(), genes: Vec::new(), signals: Vec::new(), locus: locus.or(trailing) };
    for s in inputs {
        match kind_of(&s) {
            Some(Kind::Variants) => out.vcfs.push(s.into()),
            Some(Kind::Genes) => out.genes.push(s.into()),
            Some(Kind::Signal) => out.signals.push(s.into()),
            _ => out.bams.push(s.into()),
        }
    }
    out
}

/// Windowed apps start with no console; borrow the parent terminal's if there is one. Output
/// that is redirected (to a file or pipe) is inherited as usual and unaffected.
#[cfg(windows)]
#[allow(unsafe_code)]
fn attach_console() {
    use windows_sys::Win32::System::Console::{ATTACH_PARENT_PROCESS, AttachConsole};
    // SAFETY: AttachConsole takes a process id and has no memory-safety preconditions; it fails
    // harmlessly when there is no parent console or one is already attached.
    unsafe {
        AttachConsole(ATTACH_PARENT_PROCESS);
    }
}

#[cfg(not(windows))]
fn attach_console() {}

/// A `ggv://open?...` link among the inputs, as a session to open.
fn link_session(inputs: &mut Vec<String>) -> Result<Option<session::Session>, String> {
    let Some(i) = inputs.iter().position(|s| links::is_link(s)) else { return Ok(None) };
    let link = links::parse(&inputs.remove(i))?;
    let sorted = files::sort(link.files);
    Ok(Some(session::Session {
        version: session::VERSION,
        bams: sorted.bams,
        vcfs: sorted.vcfs,
        genes: sorted.genes,
        signals: sorted.signals,
        reference: link.reference.or(sorted.reference),
        locus: link.locus,
        highlight: link.highlight,
        ..Default::default()
    }))
}

fn main() -> eframe::Result {
    attach_console();
    let mut cli = Cli::parse();
    let linked_session = match link_session(&mut cli.inputs) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("ggv-desktop: {e}");
            std::process::exit(2);
        }
    };
    // A filter that doesn't parse stops here: starting unfiltered would quietly change every
    // count, screenshot and report.
    if let Some(f) = &cli.filter
        && let Err(e) = expr::intern(f)
    {
        eprintln!("ggv-desktop: --filter {f:?}: {e}");
        std::process::exit(2);
    }
    if let Some(s) = &linked_session
        && let Err(e) = s.validate()
    {
        eprintln!("ggv-desktop: {e}");
        std::process::exit(2);
    }
    app::load_read_window_pref();
    if let Some(bp) = cli.read_window {
        app::set_read_window(bp);
    }
    let mut inputs = split_inputs(cli.inputs, cli.locus);
    inputs.vcfs.extend(cli.vcfs);
    inputs.genes.extend(cli.genes);
    let mut options = eframe::NativeOptions {
        renderer: eframe::Renderer::Wgpu,
        viewport: egui::ViewportBuilder::default().with_inner_size([1400.0, 900.0]).with_title("GGV"),
        ..Default::default()
    };
    if cli.bench {
        // Measure how fast frames can be made, not the display's refresh rate.
        options.wgpu_options.surface.present_mode = eframe::wgpu::PresentMode::AutoNoVsync;
    }
    let result = eframe::run_native(
        "ggv",
        options,
        Box::new(move |cc| {
            let session = match (&cli.session, cli.restore) {
                _ if linked_session.is_some() => linked_session,
                (Some(p), _) => Some(session::Session::load(p)?),
                (None, true) => session::last_session_path().filter(|p| p.exists()).map(|p| session::Session::load(&p)).transpose()?,
                (None, false) => None,
            };
            let split = cli.split.clone();
            let args = app::Args {
                session,
                bams: inputs.bams,
                vcfs: inputs.vcfs,
                genes: inputs.genes,
                signals: inputs.signals,
                all_transcripts: cli.all_transcripts,
                show_evidence: cli.evidence,
                demo: cli.demo,
                review: cli.review,
                select: cli.select,
                select_read: cli.select_read,
                style: tracks::Style {
                    color: match cli.color {
                        ColorArg::None => tracks::ColorBy::None,
                        ColorArg::Strand => tracks::ColorBy::Strand,
                        ColorArg::Haplotype => tracks::ColorBy::Haplotype,
                        ColorArg::InsertSize => tracks::ColorBy::InsertSize,
                        ColorArg::BaseMods => tracks::ColorBy::BaseMods,
                    },
                    group: match cli.group {
                        GroupArg::None => ggv_layout::GroupBy::None,
                        GroupArg::Strand => ggv_layout::GroupBy::Strand,
                        GroupArg::Haplotype => ggv_layout::GroupBy::Haplotype,
                        GroupArg::ReadGroup => ggv_layout::GroupBy::ReadGroup,
                        // Placed at the locus once it is resolved.
                        GroupArg::Base => ggv_layout::GroupBy::BaseAt(u64::MAX),
                        GroupArg::Alleles => ggv_layout::GroupBy::Alleles,
                    },
                    soft_clips: cli.soft_clips,
                    dark: cli.dark,
                    ..Default::default()
                },
                sort: cli.sort.map(|s| match s {
                    SortArg::Base => ggv_layout::SortBy::Base,
                    SortArg::Strand => ggv_layout::SortBy::Strand,
                    SortArg::Mapq => ggv_layout::SortBy::MappingQuality,
                    SortArg::Start => ggv_layout::SortBy::Start,
                    SortArg::InsertSize => ggv_layout::SortBy::InsertSize,
                }),
                reference: cli.reference,
                locus: inputs.locus,
                highlight: cli.highlight,
                genome: cli.genome,
                tracks: cli.tracks,
                automation: automation::Automation::new(cli.screenshot, cli.bench, cli.zoom_cycle),
                bookmarks: None,
                overlays: None,
                queue: cli.queue,
                filter: cli.filter,
                find_seq: cli.find_seq,
                find_read: cli.find_read,
                compare: cli.compare,
                palette: cli.palette,
                shortcuts: cli.shortcuts,
                about: cli.about,
                report: cli.report,
                measure: cli.measure,
            };
            Ok(Box::new(shell::Shell::new(cc, args, split, cli.link_panes)?))
        }),
    );
    if automation::failed() {
        std::process::exit(1);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trailing_locus_is_split_off() {
        let i = split_inputs(vec!["a.bam".into(), "b.bam".into(), "chr1:100".into()], None);
        assert_eq!(i.bams, vec![PathBuf::from("a.bam"), PathBuf::from("b.bam")]);
        assert_eq!(i.locus.as_deref(), Some("chr1:100"));

        let i = split_inputs(vec!["a.bam".into()], Some("chr2".into()));
        assert_eq!(i.bams.len(), 1);
        assert_eq!(i.locus.as_deref(), Some("chr2"));

        let i = split_inputs(vec!["a.bam".into(), "b.BAM".into()], None);
        assert_eq!(i.bams.len(), 2);
        assert_eq!(i.locus, None);
    }

    #[test]
    fn inputs_are_sorted_by_extension() {
        let i = split_inputs(
            vec!["a.bam".into(), "calls.vcf.gz".into(), "x.bcf".into(), "genes.gtf.gz".into(), "t.bed".into(), "BRCA1".into()],
            None,
        );
        assert_eq!(i.bams, vec![PathBuf::from("a.bam")]);
        assert_eq!(i.vcfs, vec![PathBuf::from("calls.vcf.gz"), PathBuf::from("x.bcf")]);
        assert_eq!(i.genes, vec![PathBuf::from("genes.gtf.gz"), PathBuf::from("t.bed")]);
        assert_eq!(i.locus.as_deref(), Some("BRCA1"));
    }
}
