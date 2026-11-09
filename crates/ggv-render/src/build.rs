use ggv_layout::{AlignedRead, Coverage, ViewData};

use crate::{Anchor, Color, Prim, Scene, text_width};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ColorBy {
    /// IGV's default grey.
    #[default]
    None,
    /// Forward reads red, reverse reads blue.
    Strand,
    /// By `HP` tag: haplotype 1 red, 2 blue, others green; untagged grey.
    Haplotype,
    /// Pairs by orientation and insert size (as IGV): mate elsewhere orange, outward-facing
    /// green, same-strand teal/blue, too long red, too short blue.
    InsertSize,
    /// Each called base by its modification likelihood (MM/ML): red likely modified, blue not.
    BaseMods,
    /// By the value of a tag (load the view with `ViewOptions::tag`): one colour per value.
    Tag,
}

#[derive(Clone, Debug)]
pub struct RenderOptions {
    /// Image width in CSS pixels. Default 1200.
    pub width: u32,
    pub color_by: ColorBy,
    /// A 0-based half-open span to outline across all tracks (e.g. the base an agent asked
    /// about).
    pub highlight: Option<(u64, u64)>,
    /// Coverage track height. Default 50.
    pub coverage_height: f32,
    /// With more rows than this, reads are drawn squished (3 px). Default 40.
    pub squish_above_rows: usize,
    /// Join the parts of split reads that share a row (pack with `link` to put them there).
    pub link: bool,
    /// Draw splice junctions (RNA-seq) as a sashimi band under the coverage: an arc per
    /// intron with its read count, + strand above the line and − below.
    pub sashimi: bool,
    /// Junctions with fewer reads than this are left out of the band. Default 2.
    pub min_junction_reads: u32,
}

impl Default for RenderOptions {
    fn default() -> Self {
        Self {
            width: 1200,
            color_by: ColorBy::None,
            highlight: None,
            coverage_height: 50.0,
            squish_above_rows: 40,
            link: false,
            sashimi: false,
            min_junction_reads: 2,
        }
    }
}

// IGV's palette.
const BG: Color = Color::rgb(255, 255, 255);
const TEXT: Color = Color::rgb(40, 40, 40);
const MUTED: Color = Color::rgb(120, 120, 120);
const RULE: Color = Color::rgb(200, 200, 200);
const READ: Color = Color::rgb(185, 185, 185);
const READ_FWD: Color = Color::rgb(230, 160, 160);
const READ_REV: Color = Color::rgb(160, 160, 230);
const MAPQ0_OUTLINE: Color = Color::rgb(185, 185, 185);
const COVERAGE: Color = Color::rgb(175, 175, 175);
const LINK: Color = Color::rgb(70, 110, 200);
const COVERAGE_MAX: Color = Color::rgb(215, 215, 215);
const DELETION: Color = Color::rgb(30, 30, 30);
const REF_SKIP: Color = Color::rgb(150, 150, 150);
const INSERTION: Color = Color::rgb(138, 94, 161);
const HIGHLIGHT: Color = Color(0, 0, 0, 170);

fn base_color(b: u8) -> Color {
    match b {
        b'A' => Color::rgb(0, 150, 0),
        b'C' => Color::rgb(0, 0, 255),
        b'G' => Color::rgb(209, 113, 5),
        b'T' => Color::rgb(255, 0, 0),
        _ => Color::rgb(100, 100, 100),
    }
}

/// IGV shades mismatches by base quality: faint below 5, solid from 20.
fn quality_alpha(q: u8) -> u8 {
    match q {
        0xff => 255,
        q if q < 5 => 40,
        q if q >= 20 => 255,
        q => (40.0 + (f32::from(q) - 5.0) / 15.0 * 215.0) as u8,
    }
}

const ALT_THRESHOLD: f64 = 0.2;
const MIN_LABELED_DELETION: u32 = 2;
const MIN_INSERTION_ZOOMED_OUT: u32 = 5;

struct Canvas {
    prims: Vec<Prim>,
    start: u64,
    ppb: f32,
    width: f32,
}

impl Canvas {
    fn x(&self, pos0: u64) -> f32 {
        (pos0 as f64 - self.start as f64) as f32 * self.ppb
    }

    fn rect(&mut self, x: f32, y: f32, w: f32, h: f32, fill: Color) {
        let (x0, x1) = (x.max(-10.0), (x + w).min(self.width + 10.0));
        if x1 > x0 {
            self.prims.push(Prim::Rect { x: x0, y, w: x1 - x0, h, fill });
        }
    }

    /// Whether the horizontal span `x0..x1` touches the image (with a small margin).
    fn visible(&self, x0: f32, x1: f32) -> bool {
        x1.max(x0) >= -10.0 && x0.min(x1) <= self.width + 10.0
    }

    fn line(&mut self, x1: f32, y1: f32, x2: f32, y2: f32, color: Color, width: f32) {
        if !self.visible(x1, x2) {
            return;
        }
        let clamp = |x: f32| x.clamp(-10.0, self.width + 10.0);
        self.prims.push(Prim::Line { x1: clamp(x1), y1, x2: clamp(x2), y2, color, width });
    }

    fn polygon(&mut self, points: Vec<(f32, f32)>, fill: Color) {
        let (lo, hi) = points
            .iter()
            .fold((f32::MAX, f32::MIN), |(lo, hi), (x, _)| (lo.min(*x), hi.max(*x)));
        if self.visible(lo, hi) {
            self.prims.push(Prim::Polygon { points, fill });
        }
    }

    fn text(&mut self, x: f32, y: f32, text: impl Into<String>, scale: u8, color: Color, anchor: Anchor) {
        let text = text.into();
        let w = text_width(&text, scale);
        if self.visible(x - w, x + w) {
            self.prims.push(Prim::Text { x, y, text, scale, color, anchor });
        }
    }
}

pub fn build_scene(view: &ViewData, opts: &RenderOptions) -> Scene {
    let region = &view.region;
    let width = opts.width.max(200) as f32;
    let mut c = Canvas {
        prims: Vec::new(),
        start: region.start,
        ppb: width / region.len().max(1) as f32,
        width,
    };

    // Header: the region large, the file name on the right. In a narrow picture the region
    // drops to the small size, and the file name goes if there is still no room for both.
    let title = format!("{}:{}-{}  ({} bp)", region.contig, commas(region.start + 1), commas(region.end), commas(region.len()));
    let label_w = if view.label.is_empty() { 0.0 } else { text_width(&view.label, 1) + 16.0 };
    let big = text_width(&title, 2) + label_w + 8.0 <= width;
    let (scale, title_y) = if big { (2, 5.0) } else { (1, 9.0) };
    c.text(4.0, title_y, title.as_str(), scale, TEXT, Anchor::Start);
    if text_width(&title, scale) + label_w + 8.0 <= width {
        c.text(width - 4.0, 9.0, view.label.clone(), 1, MUTED, Anchor::End);
    }
    let mut y = 26.0;

    // Ruler.
    ruler(&mut c, region.start, region.end, y);
    y += 26.0;
    let tracks_top = y;

    // Reference.
    if let Some(reference) = view.reference.as_deref().filter(|_| c.ppb >= 1.0) {
        let letters = c.ppb >= 7.0;
        for (i, &b) in reference.iter().enumerate() {
            let x = c.x(region.start + i as u64);
            if letters {
                c.text(x + c.ppb / 2.0, y + 1.0, (b as char).to_string(), 1, base_color(b), Anchor::Middle);
            } else {
                c.rect(x, y + 1.0, c.ppb.max(1.0), 7.0, base_color(b));
            }
        }
        y += 12.0;
    }

    // Coverage.
    y = coverage(&mut c, view, y, opts.coverage_height);
    c.line(0.0, y + 0.5, width, y + 0.5, RULE, 1.0);
    y += 4.0;
    if opts.sashimi
        && let Some((reads, _)) = &view.reads
    {
        y = sashimi(&mut c, view, reads, y, opts.min_junction_reads);
        c.line(0.0, y + 0.5, width, y + 0.5, RULE, 1.0);
        y += 4.0;
    }

    // Reads.
    match &view.reads {
        None => {
            let long = "Zoom in to see reads (this region is wider than the read window)";
            let note = if text_width(long, 1) + 8.0 <= width { long } else { "Zoom in to see reads" };
            c.text(4.0, y + 4.0, note, 1, MUTED, Anchor::Start);
            y += 16.0;
        }
        Some((reads, packed)) => {
            let expanded = packed.rows.len() <= opts.squish_above_rows;
            let (row_h, pitch) = if expanded { (10.0, 12.0) } else { (3.0, 4.0) };
            let bounds = if opts.color_by == ColorBy::InsertSize { ggv_layout::classify::insert_bounds(reads) } else { None };
            // Each group is at least as tall as its label (squished groups of a row or two
            // would otherwise run their labels into the next group's).
            let mut row_y = Vec::with_capacity(packed.rows.len() + 1);
            let mut extra = 0.0;
            for r in 0..=packed.rows.len() {
                if let Some(g) = packed.groups.iter().find(|g| g.first_row + g.rows == r && g.rows > 0) {
                    extra += (12.0 - g.rows as f32 * pitch).max(0.0);
                }
                row_y.push(y + r as f32 * pitch + extra);
            }
            for (r, row) in packed.rows.iter().enumerate() {
                let ry = row_y[r];
                for &i in row {
                    draw_read(&mut c, &reads[i], ry, row_h, expanded, opts.color_by, bounds);
                }
                if opts.link {
                    let mid = ry + row_h / 2.0;
                    for w in row.windows(2) {
                        let (a, b) = (&reads[w[0]], &reads[w[1]]);
                        if ggv_layout::same_template(a, b) && b.start > a.end {
                            c.line(c.x(a.end), mid, c.x(b.start), mid, LINK, 1.0);
                        }
                    }
                }
            }
            // Group labels (grouped packings), over each group's first row, as in the app.
            for g in &packed.groups {
                let gy = row_y[g.first_row];
                if g.first_row > 0 {
                    c.line(0.0, gy - pitch / 2.0, width, gy - pitch / 2.0, RULE, 1.0);
                }
                let label = format!("{} ({} read{})", g.label, commas(g.reads as u64), if g.reads == 1 { "" } else { "s" });
                let lw = text_width(&label, 1);
                c.rect(2.0, gy, lw + 6.0, 10.0, BG.with_alpha(220));
                c.text(5.0, gy + 1.5, label, 1, TEXT, Anchor::Start);
            }
            y = row_y[packed.rows.len()] + 4.0;

            let shown = reads.len() - packed.hidden_reads;
            let mut footer = format!("{} reads in {} rows", commas(shown as u64), packed.rows.len());
            if packed.hidden_reads > 0 {
                footer.push_str(&format!(", {} more not shown (row limit)", commas(packed.hidden_reads as u64)));
            }
            let rc = &view.read_counts;
            let excluded: Vec<String> = [
                ("duplicate", rc.duplicate),
                ("secondary", rc.secondary),
                ("supplementary", rc.supplementary),
                ("qc-fail", rc.qc_fail),
                ("low mapq", rc.low_mapping_quality),
                ("other flags", rc.other_flags),
            ]
            .iter()
            .filter(|(_, n)| *n > 0)
            .map(|(k, n)| format!("{} {k}", commas(*n)))
            .collect();
            if !excluded.is_empty() {
                footer.push_str(&format!("; excluded {}", excluded.join(", ")));
            }
            c.text(4.0, y, footer, 1, MUTED, Anchor::Start);
            y += 12.0;
        }
    }

    // Highlight goes on top of everything.
    if let Some((a, b)) = opts.highlight
        && b > region.start
        && a < region.end
    {
        let x0 = c.x(a.max(region.start));
        let x1 = c.x(b.min(region.end)).max(x0 + 1.0);
        c.prims.push(Prim::StrokeRect {
            x: x0 - 0.5,
            y: tracks_top - 2.0,
            w: x1 - x0 + 1.0,
            h: y - tracks_top - 10.0,
            color: HIGHLIGHT,
            width: 1.0,
        });
    }

    Scene { width: width as u32, height: (y + 4.0).ceil() as u32, background: BG, prims: c.prims }
}

fn ruler(c: &mut Canvas, start: u64, end: u64, y: f32) {
    let len = end - start;
    let label_w = text_width(&commas(end), 1) + 16.0;
    let mut step = nice_step(len as f64 / 10.0);
    while (step as f32) * c.ppb < label_w {
        step = nice_step(step as f64 * 1.5);
    }
    let first = (start + 1).div_ceil(step) * step;
    let mut p = first;
    while p <= end {
        let x = c.x(p - 1) + c.ppb / 2.0;
        let label = commas(p);
        let half = text_width(&label, 1) / 2.0;
        if x - half >= 0.0 && x + half <= c.width {
            c.text(x, y, label, 1, TEXT, Anchor::Middle);
        }
        c.line(x, y + 10.0, x, y + 17.0, MUTED, 1.0);
        p += step;
    }
    c.line(0.0, y + 17.5, c.width, y + 17.5, MUTED, 1.0);
}

/// The smallest 1, 2 or 5 × 10ⁿ at least `raw`.
fn nice_step(raw: f64) -> u64 {
    let raw = raw.max(1.0);
    let mag = 10f64.powf(raw.log10().floor());
    let step = [1.0, 2.0, 5.0, 10.0].into_iter().map(|m| m * mag).find(|s| *s >= raw).unwrap_or(10.0 * mag);
    step.round().max(1.0) as u64
}

const JUNCTION_PLUS: Color = Color::rgb(200, 60, 60);
const JUNCTION_MINUS: Color = Color::rgb(60, 90, 200);
const JUNCTION_UNKNOWN: Color = Color::rgb(130, 130, 130);

/// The sashimi band: one arc per intron the reads skip, its height and weight growing with its
/// read count (labelled at the top). An intron whose ends are in view takes its strand from
/// the splice-site motif: + arcs rise above the middle line, − arcs hang below it, and ones
/// whose strand can't be told are grey, above. Returns the band's bottom.
fn sashimi(c: &mut Canvas, view: &ViewData, reads: &[AlignedRead], y: f32, min_reads: u32) -> f32 {
    let mut counts: std::collections::BTreeMap<(u64, u64), u32> = std::collections::BTreeMap::new();
    for r in reads {
        for s in &r.ref_skips {
            *counts.entry((s.pos, s.pos + u64::from(s.len))).or_default() += 1;
        }
    }
    counts.retain(|_, n| *n >= min_reads);
    c.text(4.0, y, format!("Junctions (>= {min_reads} reads)"), 1, MUTED, Anchor::Start);
    if counts.is_empty() {
        c.text(4.0, y + 14.0, "No spliced reads here", 1, MUTED, Anchor::Start);
        return y + 26.0;
    }
    let base = |p: u64| -> Option<u8> {
        let r = view.reference.as_ref()?;
        r.get(p.checked_sub(view.region.start)? as usize).copied()
    };
    let arcs: Vec<(u64, u64, u32, char)> = counts
        .iter()
        .map(|(&(s, e), &n)| {
            let strand = match (base(s), base(s + 1), base(e.saturating_sub(2)), base(e.saturating_sub(1))) {
                (Some(a), Some(b), Some(x), Some(z)) => ggv_model::splice_motif([a, b], [x, z]).1,
                _ => '.',
            };
            (s, e, n, strand)
        })
        .collect();
    // Room above the line only if some arc goes up, below only if some goes down.
    let half = 34.0;
    let up = arcs.iter().any(|a| a.3 != '-');
    let down = arcs.iter().any(|a| a.3 == '-');
    let mid = y + 12.0 + if up { half } else { 4.0 };
    c.line(0.0, mid + 0.5, c.width, mid + 0.5, RULE, 1.0);
    let max = arcs.iter().map(|a| a.2).max().unwrap_or(1).max(1);
    // Labels already placed on each side, to step new ones clear of them.
    let mut placed: [Vec<(f32, f32, f32)>; 2] = [Vec::new(), Vec::new()];
    for &(s, e, n, strand) in &arcs {
        let (color, dir) = match strand {
            '+' => (JUNCTION_PLUS, -1.0),
            '-' => (JUNCTION_MINUS, 1.0),
            _ => (JUNCTION_UNKNOWN, -1.0),
        };
        let weight = (f64::from(n).ln_1p() / f64::from(max).ln_1p()) as f32;
        let h = (8.0 + (half - 10.0) * weight) * dir;
        let (x1, x2) = (c.x(s), c.x(e));
        // A half-sine from one splice site to the other.
        const STEPS: usize = 24;
        let point = |i: usize| {
            let t = i as f32 / STEPS as f32;
            (x1 + (x2 - x1) * t, mid + h * (std::f32::consts::PI * t).sin())
        };
        for i in 0..STEPS {
            let ((ax, ay), (bx, by)) = (point(i), point(i + 1));
            c.line(ax, ay, bx, by, color, 1.0 + 2.0 * weight);
        }
        let label = n.to_string();
        let w = text_width(&label, 1);
        let (cx, cy) = point(STEPS / 2);
        let cx = cx.clamp(w / 2.0 + 2.0, c.width - w / 2.0 - 2.0);
        let side = usize::from(dir > 0.0);
        let mut ly = if dir < 0.0 { cy - 9.0 } else { cy + 2.0 };
        while placed[side].iter().any(|&(a, b, py)| cx + w / 2.0 > a && cx - w / 2.0 < b && (py - ly).abs() < 8.0) {
            ly += 8.0 * dir;
        }
        placed[side].push((cx - w / 2.0, cx + w / 2.0, ly));
        c.text(cx, ly, label, 1, color, Anchor::Middle);
    }
    let lowest = placed[1].iter().map(|p| p.2 + 9.0).fold(mid, f32::max);
    if down { (mid + half).max(lowest) + 4.0 } else { mid + 6.0 }
}

fn coverage(c: &mut Canvas, view: &ViewData, y: f32, height: f32) -> f32 {
    let max = view.coverage.max_depth().max(1);
    let scale = height / max as f32;
    let bottom = y + 10.0 + height;

    match &view.coverage {
        Coverage::Unavailable => {
            c.text(4.0, y + 20.0, "Zoom in to see coverage", 1, MUTED, Anchor::Start);
        }
        Coverage::PerBase(cols) if c.ppb < 1.0 => {
            // Several bases per pixel: one grey bar per pixel at the maximum depth, then the
            // variable columns on top so variants stay visible when zoomed out.
            let px = c.width.ceil() as usize;
            let mut max = vec![0u32; px];
            for (i, col) in cols.iter().enumerate() {
                let p = (c.x(c.start + i as u64) as usize).min(px.saturating_sub(1));
                max[p] = max[p].max(col.depth);
            }
            for (p, d) in max.iter().enumerate().filter(|(_, d)| **d > 0) {
                let h = *d as f32 * scale;
                c.rect(p as f32, bottom - h, 1.0, h, COVERAGE);
            }
            for (i, col) in cols.iter().enumerate().filter(|(_, col)| col.mismatch_fraction() >= ALT_THRESHOLD) {
                let x = c.x(c.start + i as u64).floor();
                stack_alleles(c, col, x, 1.0, bottom, scale);
            }
        }
        Coverage::PerBase(cols) => {
            let w = if c.ppb >= 4.0 { c.ppb - 1.0 } else { c.ppb.max(1.0) };
            for (i, col) in cols.iter().enumerate() {
                if col.depth == 0 {
                    continue;
                }
                let x = c.x(c.start + i as u64);
                let h = col.depth as f32 * scale;
                if col.mismatch_fraction() >= ALT_THRESHOLD {
                    stack_alleles(c, col, x, w, bottom, scale);
                } else {
                    c.rect(x, bottom - h, w, h, COVERAGE);
                }
            }
        }
        Coverage::Binned(bins) | Coverage::Estimated(bins) => {
            let bw = c.width / bins.len().max(1) as f32;
            for (i, b) in bins.iter().enumerate() {
                if b.max == 0 {
                    continue;
                }
                let x = i as f32 * bw;
                c.rect(x, bottom - b.max as f32 * scale, bw.max(1.0), b.max as f32 * scale, COVERAGE_MAX);
                let mh = b.mean as f32 * scale;
                c.rect(x, bottom - mh, bw.max(1.0), mh, COVERAGE);
            }
        }
    }

    // Label last so bars don't cover it.
    c.text(4.0, y, format!("Coverage [0-{}]", commas(u64::from(max))), 1, MUTED, Anchor::Start);
    bottom
}

/// Stacks a column's allele counts in their colours, like IGV's allele-frequency bars; deleted
/// reads go on top in grey.
fn stack_alleles(c: &mut Canvas, col: &ggv_layout::CoverageColumn, x: f32, w: f32, bottom: f32, scale: f32) {
    let mut top = bottom;
    for (k, base) in [b'A', b'C', b'G', b'T', b'N'].into_iter().enumerate() {
        let seg = col.counts[k] as f32 * scale;
        if seg > 0.0 {
            top -= seg;
            c.rect(x, top, w, seg, base_color(base));
        }
    }
    if col.deletions > 0 {
        let seg = col.deletions as f32 * scale;
        c.rect(x, top - seg, w, seg, COVERAGE);
    }
}

fn draw_read(c: &mut Canvas, r: &AlignedRead, y: f32, h: f32, expanded: bool, color_by: ColorBy, bounds: Option<(u32, u32)>) {
    use ggv_layout::classify::{PairClass, pair_class};
    let mid = y + h / 2.0;
    let fill = match color_by {
        ColorBy::None | ColorBy::BaseMods => READ,
        ColorBy::Strand if r.reverse => READ_REV,
        ColorBy::Strand => READ_FWD,
        ColorBy::Haplotype => match r.haplotype {
            Some(1) => Color::rgb(230, 150, 150),
            Some(2) => Color::rgb(150, 150, 230),
            Some(_) => Color::rgb(150, 210, 150),
            None => READ,
        },
        ColorBy::Tag => r.tag.as_deref().map_or(READ, |v| {
            let (red, green, blue) = ggv_layout::classify::tag_rgb(v);
            Color::rgb(red, green, blue)
        }),
        ColorBy::InsertSize => match pair_class(r, bounds) {
            Some(PairClass::OtherContig) => Color::rgb(240, 140, 40),
            Some(PairClass::RightLeft) => Color::rgb(0, 150, 0),
            Some(PairClass::LeftLeft) => Color::rgb(0, 150, 150),
            Some(PairClass::RightRight) => Color::rgb(20, 50, 200),
            Some(PairClass::Large) => Color::rgb(210, 30, 30),
            Some(PairClass::Small) => Color::rgb(30, 30, 210),
            None => READ,
        },
    };
    let mapq0 = r.mapq == 0;

    for s in &r.ref_skips {
        c.line(c.x(s.pos), mid + 0.5, c.x(s.pos + u64::from(s.len)), mid + 0.5, REF_SKIP, 1.0);
    }
    for d in &r.deletions {
        c.line(c.x(d.pos), mid + 0.5, c.x(d.pos + u64::from(d.len)), mid + 0.5, DELETION, 1.0);
    }

    for b in &r.blocks {
        let x0 = c.x(b.ref_start);
        let w = (c.x(b.ref_start + u64::from(b.len)) - x0).max(1.0);
        if !c.visible(x0, x0 + w) {
            continue;
        }
        if mapq0 && expanded {
            c.rect(x0, y, w, h, BG);
            c.prims.push(Prim::StrokeRect { x: x0 + 0.5, y: y + 0.5, w: w - 1.0, h: h - 1.0, color: MAPQ0_OUTLINE, width: 1.0 });
        } else {
            c.rect(x0, y, w, h, if mapq0 { fill.with_alpha(90) } else { fill });
        }
    }

    // Methylation: each called base over the read, by likelihood.
    if color_by == ColorBy::BaseMods {
        for (pos, likelihood) in ggv_layout::classify::mod_likelihoods(r) {
            let x0 = c.x(pos);
            let w = (c.x(pos + 1) - x0).max(1.0);
            if c.visible(x0, x0 + w) {
                let (red, green, blue, alpha) = ggv_layout::classify::mod_rgba(likelihood);
                c.rect(x0, y, w, h, Color::rgb(red, green, blue).with_alpha(alpha));
            }
        }
    }

    // Strand arrowhead.
    if expanded {
        let tip = 4.0f32.min(c.ppb * 3.0).max(2.0);
        let points = if r.reverse {
            let x = c.x(r.start);
            vec![(x, y), (x - tip, mid), (x, y + h)]
        } else {
            let x = c.x(r.end);
            vec![(x, y), (x + tip, mid), (x, y + h)]
        };
        let arrow = if mapq0 { MAPQ0_OUTLINE } else { fill };
        c.polygon(points, arrow);
    }

    // Labels for long deletions, so their size is readable without counting pixels.
    if expanded {
        for d in &r.deletions {
            let (x0, x1) = (c.x(d.pos), c.x(d.pos + u64::from(d.len)));
            let label = d.len.to_string();
            let lw = text_width(&label, 1);
            if d.len >= MIN_LABELED_DELETION && x1 - x0 > lw + 8.0 {
                let cx = (x0 + x1) / 2.0;
                c.rect(cx - lw / 2.0 - 2.0, y + 1.0, lw + 4.0, h - 2.0, BG);
                c.text(cx, y + 1.5, label, 1, DELETION, Anchor::Middle);
            }
        }
    }

    let letters = expanded && c.ppb >= 7.0;
    for m in &r.mismatches {
        let color = base_color(m.base).with_alpha(quality_alpha(m.quality));
        let x = c.x(m.pos);
        if letters {
            // A pale tint of the base colour behind the letter keeps it readable on grey.
            c.rect(x, y, c.ppb, h, BG.with_alpha(color.3 / 2 + 60));
            c.text(x + c.ppb / 2.0, y + 1.5, (m.base as char).to_string(), 1, color, Anchor::Middle);
        } else {
            c.rect(x, y, c.ppb.max(1.0), h, color);
        }
    }

    // Zoomed out, single-base insertions (mostly long-read noise) would hide real signal.
    let min_insertion = if c.ppb < 1.0 { MIN_INSERTION_ZOOMED_OUT } else { 1 };
    for ins in r.insertions.iter().filter(|i| i.len >= min_insertion) {
        let x = c.x(ins.pos);
        c.rect(x - 1.0, y, 2.0, h, INSERTION);
    }
}

pub(crate) fn commas(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    for (i, ch) in s.chars().enumerate() {
        if i > 0 && (s.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(ch);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn commas_and_steps() {
        assert_eq!(commas(0), "0");
        assert_eq!(commas(1234567), "1,234,567");
        assert_eq!(nice_step(7.0), 10);
        assert_eq!(nice_step(13.0), 20);
        assert_eq!(nice_step(0.3), 1);
        assert_eq!(nice_step(4100.0), 5000);
    }

    #[test]
    fn quality_shading() {
        assert_eq!(quality_alpha(2), 40);
        assert_eq!(quality_alpha(30), 255);
        assert!(quality_alpha(12) > 40 && quality_alpha(12) < 255);
    }
}
