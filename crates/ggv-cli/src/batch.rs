//! `ggv batch`: a picture and an evidence report for every place in a list, with no window.
//!
//! For each place (BED, VCF, or one locus per line): every alignment file drawn into one PNG,
//! an HTML report beside it (consequences, VCF records, allele counts by strand, or a region
//! summary), then `index.html` (a contact sheet linking them all) and `summary.tsv` (one row per
//! place, for spreadsheets and scripts). Places run in parallel, each worker with its own file
//! handles; the annotation is parsed once and shared.

use std::{
    cell::RefCell,
    fmt::Write as _,
    path::{Path, PathBuf},
    sync::{
        Mutex,
        atomic::{AtomicUsize, Ordering},
    },
    time::Instant,
};

use anyhow::{Context, Result};
use ggv_annot::{Annotation, Transcript, coding};
use ggv_io::{BamFile, ReferenceFasta, VariantFile, places::Place};
use ggv_model::Region;
use ggv_query::{PileupOptions, PositionStats, SummaryOptions};
use ggv_render::{
    Prim, Scene,
    report::{Card, Report, Section, contact_sheet, utc_stamp},
};

use crate::transcripts;

pub struct Batch {
    pub list: PathBuf,
    pub bams: Vec<PathBuf>,
    pub reference: Option<PathBuf>,
    pub vcfs: Vec<PathBuf>,
    pub genes: Option<PathBuf>,
    pub out: PathBuf,
    pub width: u32,
    pub scale: f32,
    pub flank: u64,
    pub max_rows: usize,
    pub jobs: usize,
    pub pileup: PileupOptions,
    pub color_by: ggv_render::ColorBy,
    /// Grouping; `BaseAt(u64::MAX)` means at each place's base (single-base places only).
    pub group: ggv_layout::GroupBy,
    /// The tag reads carry for colouring or grouping by it.
    pub tag: Option<[u8; 2]>,
    pub link: bool,
}

/// One place's results (or an unreadable line of the list), for the index and the TSV.
struct Done {
    index: usize,
    label: String,
    locus: String,
    image: String,
    report: String,
    consequence: String,
    per_bam: Vec<String>,
    error: Option<String>,
}

/// Files one worker reads from.
struct Handles {
    bams: Vec<BamFile>,
    labels: Vec<String>,
    fasta: Option<RefCell<ReferenceFasta>>,
    vcfs: Vec<VariantFile>,
}

impl Handles {
    fn open(b: &Batch) -> Result<Self> {
        let mut bams = Vec::new();
        let mut labels = Vec::new();
        for p in &b.bams {
            let bam = crate::open_bam(p, b.reference.as_deref())?;
            let info = bam.info();
            labels.push(info.read_groups.iter().find_map(|rg| rg.sample.clone()).unwrap_or_else(|| file_name(p)));
            bams.push(bam);
        }
        let fasta = b.reference.as_ref().map(ReferenceFasta::open).transpose()?.map(RefCell::new);
        let vcfs = b.vcfs.iter().map(VariantFile::open).collect::<Result<_, _>>()?;
        Ok(Self { bams, labels, fasta, vcfs })
    }
}

fn file_name(p: &Path) -> String {
    p.file_name().map_or_else(|| p.display().to_string(), |n| n.to_string_lossy().into_owned())
}

/// Moves a scene's primitives down by `dy`.
fn shifted(p: &Prim, dy: f32) -> Prim {
    let mut p = p.clone();
    match &mut p {
        Prim::Rect { y, .. } | Prim::StrokeRect { y, .. } | Prim::Text { y, .. } => *y += dy,
        Prim::Line { y1, y2, .. } => {
            *y1 += dy;
            *y2 += dy;
        }
        Prim::Polygon { points, .. } => points.iter_mut().for_each(|pt| pt.1 += dy),
    }
    p
}

/// One image from several scenes, stacked.
pub(crate) fn stack(scenes: Vec<Scene>) -> Scene {
    let width = scenes.iter().map(|s| s.width).max().unwrap_or(1);
    let background = scenes.first().map_or(ggv_render::Color::rgb(255, 255, 255), |s| s.background);
    let mut out = Scene { width, height: 0, background, prims: Vec::new() };
    for s in scenes {
        let dy = out.height as f32;
        out.prims.extend(s.prims.iter().map(|p| shifted(p, dy)));
        out.height += s.height;
    }
    out
}

/// Transcripts over `pos`, trying the annotation's own contig naming.
pub(crate) fn transcripts<'a>(ann: &'a Annotation, contig: &str, pos: u64) -> Vec<&'a Transcript> {
    transcripts_over(ann, contig, pos, pos + 1)
}

/// Like [`transcripts`], for every transcript overlapping `[start, end)`.
pub(crate) fn transcripts_over<'a>(ann: &'a Annotation, contig: &str, start: u64, end: u64) -> Vec<&'a Transcript> {
    let names = std::iter::once(contig.to_owned()).chain(ggv_model::contig_aliases(contig));
    for n in names {
        if ann.contig_index(&n).is_some() {
            let mut v: Vec<&Transcript> = ann.query(&n, start, end.max(start + 1)).collect();
            v.sort_by_key(|t| (t.cds.is_none(), !t.representative));
            return v;
        }
    }
    Vec::new()
}

/// The top non-reference base of a pileup, when it would be a variable site in `ggv query
/// summary` (the same minimum depth and fraction): a single error read is not a change to name.
/// `N` is a read without sequence (`SEQ *`) or an uncalled base, never an allele.
fn credible_alt(s: &PositionStats) -> Option<&'static str> {
    let min = SummaryOptions::default();
    let a = s.top_alt.as_ref()?;
    let enough = s.depth >= min.min_variable_depth && f64::from(a.count.total) >= min.min_alt_fraction * f64::from(s.depth);
    (enough && a.allele.len() == 1 && a.allele != "N").then_some(a.allele)
}

/// Whether `v` is a record of the change at `place`: it starts there, or the place is the
/// bases it changes (a BED row for a deletion, after VCF's padding base). One that only
/// overlaps it (a Mb-scale CNV, a deletion starting upstream) is not.
fn describes(v: &ggv_io::Variant, place: &Region) -> bool {
    v.start == place.start || v.shown() == (place.start, place.end)
}

/// Whether a record is a call of a change: PASS (or unfiltered), with an ALT that changes
/// something (a sequence or a symbolic `<INS>`, not only `*` or `.`), and carried by a sample
/// when samples have genotypes (not a `0/0` or DeepVariant `RefCall`). Only its sequence ALTs
/// name a change (an HGVS one); a symbolic call still takes precedence over the reads.
fn called(v: &ggv_io::Variant) -> bool {
    let pass = v.filters.iter().all(|f| f == "PASS");
    let changes = v.alts.iter().any(|a| a != "*" && a != ".");
    let alt_allele = |gt: &str| gt.split(['/', '|']).any(|x| x.parse::<u32>().is_ok_and(|n| n > 0));
    let typed: Vec<&str> = v.samples.iter().filter(|c| c.class != ggv_io::GenotypeClass::Unknown).map(|c| c.gt.as_str()).collect();
    pass && changes && (typed.is_empty() || typed.iter().any(|gt| alt_allele(gt)))
}

/// Whether an ALT inserts bases after REF: longer sequence, or a symbolic `<INS...>`.
fn inserts(reference: &str, alt: &str) -> bool {
    alt.starts_with("<INS") || (alt.bytes().all(|b| b.is_ascii_alphabetic()) && alt.len() > reference.len())
}

/// A file name part: only characters that need no escaping in a file system or a URL (`#`
/// in a PanSN contig started the link's fragment).
fn file_safe(s: &str) -> String {
    s.chars().map(|c| if c.is_ascii_alphanumeric() || matches!(c, '-' | '.' | '_' | '+') { c } else { '_' }).collect()
}

fn pileup_at(bam: &mut BamFile, region: &Region, base: Option<u8>, opts: &PileupOptions) -> Result<Option<PositionStats>> {
    let mut out = None;
    let seq = base.map(|b| vec![b]);
    ggv_query::pileup(bam, region, seq.as_deref(), opts, |s| {
        out = Some(s);
        Ok(())
    })?;
    Ok(out)
}

fn pileup_section(heading: String, s: &PositionStats) -> Section {
    let mut sec = Section { heading, ..Default::default() };
    let depth = f64::from(s.depth.max(1));
    sec.lines.push(format!("depth {}  ·  reference {}", s.depth, s.ref_base.unwrap_or('.')));
    sec.table.push(["allele", "reads", "%", "forward", "reverse"].map(String::from).to_vec());
    let b = &s.bases;
    for (name, n) in [("A", b.a), ("C", b.c), ("G", b.g), ("T", b.t), ("N", b.n), ("DEL", s.deletions), ("INS after", s.insertions_after)] {
        if n.total > 0 {
            let is_ref = s.ref_base.is_some_and(|r| r.to_string() == name);
            sec.table.push(vec![
                if is_ref { format!("{name} (ref)") } else { name.to_owned() },
                n.total.to_string(),
                format!("{:.1}", 100.0 * f64::from(n.total) / depth),
                n.fwd.to_string(),
                n.rev.to_string(),
            ]);
        }
    }
    sec.lines.push(format!(
        "MAPQ mean {:.1} ({} at 0)  ·  base quality mean {:.1}  ·  soft clips starting {} / ending {}",
        s.mean_mapping_quality, s.mapping_quality_zero, s.mean_base_quality, s.soft_clipped_left, s.soft_clipped_right
    ));
    if let Some(alt) = &s.top_alt {
        sec.lines.push(format!("top alternate {} at {:.1}% ({} forward / {} reverse)", alt.allele, alt.fraction * 100.0, alt.count.fwd, alt.count.rev));
        if alt.count.total >= 4 && (alt.count.fwd == 0 || alt.count.rev == 0) {
            sec.warning = Some("The alternate allele is supported by one strand only — a common artifact.".into());
        }
    }
    sec
}

/// Where `place` is, in the first alignment file that has its contig: the others show "not in
/// this file", whatever order the files were given in. Also the place as that file reads it
/// (a row that is a locus or a BED row on a `:` contig, settled by its contigs).
fn locate<'p>(b: &Batch, h: &Handles, place: &'p Place) -> Result<(Region, Vec<ggv_model::Contig>, &'p Place)> {
    let mut first_error = None;
    for (bam, path) in h.bams.iter().zip(&b.bams) {
        match place.resolve(bam.contigs()) {
            Ok(r) => return Ok((r, bam.contigs().to_vec(), place.settle(bam.contigs()))),
            Err(e) => {
                first_error.get_or_insert(format!("{}: {e}", file_name(path)));
            }
        }
    }
    let e = first_error.unwrap_or_else(|| "no alignment files".into());
    if h.bams.len() > 1 { anyhow::bail!("in none of the alignment files ({e})") } else { anyhow::bail!(e) }
}

fn snapshot(b: &Batch, h: &mut Handles, ann: Option<&Annotation>, index: usize, place: &Place) -> Result<Done> {
    let (target, contigs, place) = locate(b, h, place)?;
    let len_of = |c: &str| contigs.iter().find(|x| x.name == c).map_or(target.end, |x| x.length);
    let single = target.len() == 1;
    let region = if single {
        Region { contig: target.contig.clone(), start: target.start.saturating_sub(b.flank), end: (target.end + b.flank).min(len_of(&target.contig)) }
    } else {
        let pad = (target.len() / 10).max(20);
        Region { contig: target.contig.clone(), start: target.start.saturating_sub(pad), end: (target.end + pad).min(len_of(&target.contig)) }
    };
    let stem = format!("{:03}-{}-{}", index + 1, file_safe(&target.contig), target.start + 1);
    let seq = match &h.fasta {
        Some(f) if region.len() <= 2_000_000 => Some(f.borrow_mut().fetch(&region)?),
        _ => None,
    };
    let base_at = |p: u64| -> Option<u8> {
        if let Some(s) = &seq
            && p >= region.start
            && p < region.end
        {
            return s.get((p - region.start) as usize).copied();
        }
        let f = h.fasta.as_ref()?;
        f.borrow_mut().fetch(&Region { contig: target.contig.clone(), start: p, end: p + 1 }).ok()?.first().copied()
    };

    // The picture: one panel per alignment file.
    let mut scenes = Vec::new();
    for (bam, label) in h.bams.iter_mut().zip(&h.labels) {
        let Ok(local) = bam.resolve_str(&format!("{}:{}-{}", region.contig, region.start + 1, region.end)) else { continue };
        let view_opts = ggv_layout::ViewOptions { filters: b.pileup.clone(), bins: b.width as usize, max_rows: b.max_rows, tag: b.tag, ..Default::default() };
        let mut view = ggv_layout::load_view(bam, &local, seq.clone(), label.clone(), &view_opts)?;
        let group = match b.group {
            ggv_layout::GroupBy::BaseAt(_) if target.len() == 1 => ggv_layout::GroupBy::BaseAt(target.start),
            ggv_layout::GroupBy::BaseAt(_) => ggv_layout::GroupBy::None,
            g => g,
        };
        if group != ggv_layout::GroupBy::None || b.link {
            ggv_layout::repack(&mut view, b.max_rows, group, None, b.link);
        }
        let opts = ggv_render::RenderOptions { width: b.width, color_by: b.color_by, link: b.link, highlight: Some((target.start, target.end)), ..Default::default() };
        scenes.push(ggv_render::build_scene(&view, &opts));
    }
    let png = stack(scenes).to_png(b.scale).map_err(anyhow::Error::msg)?;
    let image = format!("{stem}.png");
    std::fs::write(b.out.join(&image), &png)?;

    let mut sections: Vec<Section> = h
        .bams
        .iter()
        .filter_map(BamFile::warning)
        .map(|w| Section { heading: "Input".into(), warning: Some(w), ..Default::default() })
        .collect();
    let mut per_bam = Vec::new();

    // VCF records here (the list's own row first), and the changes they make. Where a called
    // record is at the place (see `describes`), its alleles are the change, not whatever the
    // reads show; an insertion record there makes the per-file summary count insertions.
    let mut changes: Vec<(u64, String, String)> = Vec::new();
    let mut recorded = false;
    let mut insertion = place.insertion;
    if let Some((r, alts)) = place.alleles.as_deref().and_then(|a| a.split_once('>')) {
        recorded = true;
        for a in alts.split(',') {
            insertion |= inserts(r, a);
            if a.bytes().all(|b| b"ACGTNacgtn".contains(&b)) && !a.is_empty() {
                changes.push((target.start, r.to_owned(), a.to_owned()));
            }
        }
    }
    // The list's own alleles come first in the summary.
    let own = changes.len();
    for (vcf, path) in h.vcfs.iter_mut().zip(&b.vcfs) {
        let records = vcf.query(&target, false)?;
        if records.is_empty() {
            continue;
        }
        // A big SV spanning the place, or a record starting upstream, says nothing about the
        // base the reads show: only one here takes precedence.
        let here: Vec<&ggv_io::Variant> = records.iter().filter(|v| describes(v, &target) && called(v)).collect();
        recorded |= !here.is_empty();
        insertion |= here.iter().any(|v| v.start == target.start && v.alts.iter().any(|a| inserts(&v.reference, a)));
        let mut sec = Section { heading: format!("Variant records · {}", file_name(path)), ..Default::default() };
        sec.table.push(["position", "REF > ALT", "type", "QUAL", "FILTER", "genotypes"].map(String::from).to_vec());
        for v in records.iter().take(30) {
            sec.table.push(vec![
                (v.start + 1).to_string(),
                format!("{} > {}", v.reference, v.alts.join(",")),
                v.kind(),
                v.qual.map_or(".".into(), |q| format!("{q:.1}")),
                if v.filters.is_empty() { ".".into() } else { v.filters.join(",") },
                vcf.samples().iter().zip(&v.samples).map(|(s, c)| format!("{s} {}", c.gt)).collect::<Vec<_>>().join("  "),
            ]);
            // A RefCall, `0/0` or filtered record is listed, but changes nothing.
            if !called(v) {
                continue;
            }
            for a in v.alts.iter().filter(|a| coding::is_sequence(a)) {
                if !changes.iter().any(|c| c.0 == v.start && c.1 == v.reference && c.2 == *a) {
                    changes.push((v.start, v.reference.clone(), a.clone()));
                }
            }
        }
        if records.len() > 30 {
            sec.lines.push(format!("… and {} more", records.len() - 30));
        }
        sections.push(sec);
    }

    // Each alignment file at the base (or across the region).
    for (k, bam) in h.bams.iter_mut().enumerate() {
        let heading = format!("{} · {}", h.labels[k], file_name(&b.bams[k]));
        let Ok(local) = bam.resolve_str(&format!("{}:{}-{}", target.contig, target.start + 1, target.end)) else {
            sections.push(Section { heading, lines: vec!["not in this file".into()], ..Default::default() });
            per_bam.push(String::new());
            continue;
        };
        if single {
            match pileup_at(bam, &local, base_at(target.start), &b.pileup)? {
                Some(s) => {
                    // At an insertion the bases after the anchor are the evidence, not its base.
                    per_bam.push(if insertion {
                        format!("depth {} INS {:.0}%", s.depth, 100.0 * f64::from(s.insertions_after.total) / f64::from(s.depth.max(1)))
                    } else {
                        match &s.top_alt {
                            Some(a) => format!("depth {} {} {:.0}%", s.depth, a.allele, a.fraction * 100.0),
                            None => format!("depth {}", s.depth),
                        }
                    });
                    if !recorded
                        && let (Some(alt), Some(r)) = (credible_alt(&s), s.ref_base)
                        && !changes.iter().any(|c| c.0 == target.start && c.2 == alt)
                    {
                        changes.push((target.start, r.to_string(), alt.to_owned()));
                    }
                    sections.push(pileup_section(heading, &s));
                }
                None => {
                    per_bam.push("depth 0".into());
                    sections.push(Section { heading, lines: vec!["no reads here".into()], ..Default::default() });
                }
            }
        } else {
            let opts = SummaryOptions { pileup: b.pileup.clone(), max_variable_sites: 10, ..Default::default() };
            // A reference contig shorter than the BAM's (another assembly) just runs out, as in
            // `ggv query summary`; slicing past it panicked and lost the whole batch.
            let slice = seq.as_ref().map(|s| {
                let at = |p: u64| ((p - region.start) as usize).min(s.len());
                s[at(target.start)..at(target.end)].to_vec()
            });
            let s = ggv_query::summarize(bam, &local, slice.as_deref(), &opts)?;
            per_bam.push(format!("mean depth {:.1}, {} variable", s.mean_depth, s.variable_site_count));
            let mut sec = Section { heading, ..Default::default() };
            sec.lines.push(format!(
                "reads {} used of {} · depth mean {:.1}, median {}, min {}, max {} · {} variable sites",
                s.reads.used, s.reads.seen, s.mean_depth, s.median_depth, s.min_depth, s.max_depth, s.variable_site_count
            ));
            if !s.variable_sites.is_empty() {
                sec.table.push(["position", "ref", "alt", "%", "forward / reverse", "depth"].map(String::from).to_vec());
                for v in &s.variable_sites {
                    sec.table.push(vec![
                        v.pos.to_string(),
                        v.ref_base.map_or(".".into(), |c| c.to_string()),
                        v.alt.clone(),
                        format!("{:.0}", v.alt_fraction * 100.0),
                        format!("{} / {}", v.alt_fwd, v.alt_rev),
                        v.depth.to_string(),
                    ]);
                }
            }
            sections.push(sec);
        }
    }

    // What the changes do to the genes. The summary is one of them: the list's own alleles,
    // else a change starting at the place (not a record starting upstream that reaches it),
    // then the most severe, then the gene's main transcript; not whichever gene comes first.
    let mut consequence = String::new();
    if let Some(ann) = ann {
        let mut lines = Vec::new();
        let mut best = None;
        for (k, (p, r, a)) in changes.iter().enumerate() {
            let end = *p + r.len() as u64;
            for t in transcripts::pick(transcripts_over(ann, &target.contig, *p, end), *p, end) {
                if let Some(c) = coding::consequence(t, *p, r, a, base_at) {
                    let warn = if c.warnings.is_empty() { String::new() } else { format!("  [warning: {}]", c.warnings.join("; ")) };
                    let line = format!("{} {}  {}{}  {}  ({}){warn}", c.gene, c.transcript, c.hgvs_c, c.hgvs_p.as_deref().map(|p| format!("  {p}")).unwrap_or_default(), c.effect, c.location.feature);
                    let key = (k < own, *p == target.start, transcripts::severity(&c.effect), transcripts::main_rank(t));
                    if best.as_ref().is_none_or(|(b, _)| key > *b) {
                        best = Some((key, format!("{} {}{} {}", c.gene, c.hgvs_c, c.hgvs_p.map(|p| format!(" {p}")).unwrap_or_default(), c.effect)));
                    }
                    lines.push(line);
                }
            }
        }
        if let Some((_, c)) = best {
            consequence = c;
        }
        if lines.is_empty()
            && let Some((t, loc)) = transcripts::pick(transcripts(ann, &target.contig, target.start), target.start, target.start + 1)
                .into_iter()
                .filter_map(|t| Some((t, coding::locate(t, target.start)?)))
                // The first of equals, as in the changes above: `pick` puts coding genes first.
                .rev()
                .max_by_key(|(t, l)| (transcripts::severity(l.kind.name()), transcripts::main_rank(t)))
        {
            lines.push(format!("{} {}  {}  {}", t.gene_name, t.id, loc.feature, loc.kind.name()));
            consequence = format!("{} {}", t.gene_name, loc.kind.name());
        }
        if !lines.is_empty() {
            sections.insert(0, Section { heading: "Consequences".into(), lines, ..Default::default() });
        }
    }

    let locus = if single { format!("{}:{}", target.contig, target.start + 1) } else { format!("{}:{}-{}", target.contig, target.start + 1, target.end) };
    let now = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(0, |d| d.as_secs());
    let mut title = locus.clone();
    if !place.label.is_empty() {
        title.push_str(&format!("  {}", place.label));
    }
    let report = Report {
        title,
        subtitle: format!("{} · {} · place {} of the list {}", utc_stamp(now), b.bams.iter().map(|p| file_name(p)).collect::<Vec<_>>().join(", "), index + 1, file_name(&b.list)),
        png: Some(png),
        sections,
    };
    let report_name = format!("{stem}.html");
    std::fs::write(b.out.join(&report_name), report.html())?;
    Ok(Done { index, label: place.label.clone(), locus, image, report: report_name, consequence, per_bam, error: None })
}

/// What a panic says, for the summary.
fn panic_text(p: &(dyn std::any::Any + Send)) -> &str {
    p.downcast_ref::<&str>().copied().or_else(|| p.downcast_ref::<String>().map(String::as_str)).unwrap_or("unknown cause")
}

/// How a batch went: places written and places that failed.
pub struct Outcome {
    pub written: usize,
    pub failed: usize,
}

pub fn run(b: &Batch) -> Result<Outcome> {
    let t = Instant::now();
    let text = ggv_io::places::read_file(&b.list).with_context(|| b.list.display().to_string())?;
    let entries = ggv_io::places::parse_lines(&text, ggv_io::places::Format::of_path(&b.list));
    // Numbered by entry, so summary rows follow the list and unreadable lines keep their place.
    let places: Vec<(usize, &Place)> = entries.iter().enumerate().filter_map(|(i, e)| Some((i, e.as_ref().ok()?))).collect();
    let unreadable: Vec<(usize, &ggv_io::places::Unreadable)> = entries.iter().enumerate().filter_map(|(i, e)| Some((i, e.as_ref().err()?))).collect();
    if places.is_empty() {
        let first = unreadable.first().map(|(_, u)| format!("; line {}: {}", u.line, u.text)).unwrap_or_default();
        anyhow::bail!("{}: no places found ({} lines not understood{first})", b.list.display(), unreadable.len());
    }
    std::fs::create_dir_all(&b.out).with_context(|| b.out.display().to_string())?;
    let ann = b.genes.as_ref().map(Annotation::open).transpose()?;

    // A line that is not a place fails like a place would, so it is in summary.tsv and the
    // exit status, not just counted.
    let mut failed_lines = Vec::new();
    for (i, u) in &unreadable {
        let error = format!("line {} of the list is not a place (BED, VCF, `chr1:100-200` or `chr1 100`)", u.line);
        let locus = u.text.replace('\t', " ");
        eprintln!("{locus}: {error}");
        failed_lines.push(Done { index: *i, label: String::new(), locus, image: String::new(), report: String::new(), consequence: String::new(), per_bam: Vec::new(), error: Some(error) });
    }

    let next = AtomicUsize::new(0);
    let done: Mutex<Vec<Done>> = Mutex::new(failed_lines);
    // Each file's warning (a truncated BAM), once, though every worker opens every file.
    let warnings: Mutex<Vec<String>> = Mutex::new(Vec::new());
    let jobs = b.jobs.clamp(1, places.len());
    std::thread::scope(|s| -> Result<()> {
        let mut workers = Vec::new();
        for _ in 0..jobs {
            workers.push(s.spawn(|| -> Result<()> {
                let mut h = Handles::open(b)?;
                {
                    let mut w = warnings.lock().unwrap();
                    for x in h.bams.iter().filter_map(BamFile::warning) {
                        if !w.contains(&x) {
                            w.push(x);
                        }
                    }
                }
                loop {
                    let k = next.fetch_add(1, Ordering::Relaxed);
                    let Some(&(i, place)) = places.get(k) else { return Ok(()) };
                    // A bug on one place fails that place, not the run (and its summary).
                    let d = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| snapshot(b, &mut h, ann.as_ref(), i, place)))
                        .unwrap_or_else(|p| {
                            let e = anyhow::anyhow!("internal error: {}", panic_text(p.as_ref()));
                            Handles::open(b).map(|fresh| h = fresh).and(Err(e))
                        });
                    let d = d.unwrap_or_else(|e| Done {
                        index: i,
                        label: place.label.clone(),
                        locus: place.key(),
                        image: String::new(),
                        report: String::new(),
                        consequence: String::new(),
                        per_bam: Vec::new(),
                        error: Some(format!("{e:#}")),
                    });
                    match &d.error {
                        Some(e) => eprintln!("{}: {e}", d.locus),
                        None => eprintln!("{} → {}", d.locus, d.image),
                    }
                    done.lock().unwrap().push(d);
                }
            }));
        }
        for w in workers {
            w.join().map_err(|p| anyhow::anyhow!("internal error: {}", panic_text(p.as_ref())))??;
        }
        Ok(())
    })?;
    let mut done = done.into_inner().unwrap();
    done.sort_by_key(|d| d.index);

    let warnings = warnings.into_inner().unwrap();

    // summary.tsv, after any warnings about the inputs (`##` lines, as in a VCF header).
    let mut tsv: String = warnings.iter().map(|w| format!("## warning: {}\n", w.replace(['\t', '\n'], " "))).collect();
    tsv.push_str("#n\tlocus\tlabel\tconsequence");
    for p in &b.bams {
        let _ = write!(tsv, "\t{}", file_name(p));
    }
    tsv.push_str("\timage\treport\terror\n");
    for d in &done {
        let _ = write!(tsv, "{}\t{}\t{}\t{}", d.index + 1, d.locus, d.label.replace('\t', " "), d.consequence);
        for k in 0..b.bams.len() {
            let _ = write!(tsv, "\t{}", d.per_bam.get(k).map_or("", String::as_str));
        }
        let _ = writeln!(tsv, "\t{}\t{}\t{}", d.image, d.report, d.error.as_deref().unwrap_or(""));
    }
    std::fs::write(b.out.join("summary.tsv"), tsv)?;

    // index.html: a contact sheet.
    let cards: Vec<Card> = done
        .iter()
        .map(|d| Card {
            image: d.image.clone(),
            link: d.report.clone(),
            title: format!("{}. {}", d.index + 1, d.locus),
            lines: [d.label.clone(), d.consequence.clone(), d.per_bam.join(" · "), d.error.clone().unwrap_or_default()]
                .into_iter()
                .filter(|l| !l.is_empty())
                .collect(),
        })
        .collect();
    let errors = done.iter().filter(|d| d.error.is_some()).count();
    let title = format!("ggv batch · {}", file_name(&b.list));
    let mut subtitle = format!("{} places · {} · {:.1} s", done.len(), b.bams.iter().map(|p| file_name(p)).collect::<Vec<_>>().join(", "), t.elapsed().as_secs_f64());
    for w in &warnings {
        subtitle.push_str(&format!(" · warning: {w}"));
    }
    std::fs::write(b.out.join("index.html"), contact_sheet(&title, &subtitle, &cards))?;
    eprintln!("wrote {} snapshots and reports to {} in {:.1} s ({} failed)", done.len() - errors, b.out.display(), t.elapsed().as_secs_f64(), errors);
    Ok(Outcome { written: done.len() - errors, failed: errors })
}

#[cfg(test)]
mod tests {
    use super::*;
    use ggv_query::{AltAllele, StrandCount};

    fn stats(depth: u32, allele: &'static str, n: u32) -> PositionStats {
        let count = StrandCount { total: n, fwd: n, rev: 0 };
        PositionStats { depth, ref_base: Some('A'), top_alt: Some(AltAllele { allele, count, fraction: f64::from(n) / f64::from(depth) }), ..Default::default() }
    }

    #[test]
    fn only_credible_read_alleles_are_changes() {
        assert_eq!(credible_alt(&stats(48, "T", 1)), None); // one error read
        assert_eq!(credible_alt(&stats(5, "T", 5)), None); // too shallow
        assert_eq!(credible_alt(&stats(20, "N", 10)), None); // reads without sequence
        assert_eq!(credible_alt(&stats(20, "DEL", 10)), None);
        assert_eq!(credible_alt(&stats(20, "T", 4)), Some("T"));
    }

    #[test]
    fn insertions_and_file_names() {
        assert!(inserts("C", "CAGATA") && inserts("T", "<INS>") && inserts("T", "<INS:ME>"));
        assert!(!inserts("C", "A") && !inserts("CA", "C") && !inserts("T", "<DEL>") && !inserts("T", "T[chr1:5["));
        assert_eq!(file_safe("HG002#1#chr20%x?y:z"), "HG002_1_chr20_x_y_z");
    }
}
