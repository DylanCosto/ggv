//! Builds GPU instances for a loaded view. Runs once per data load (on the loader thread) or
//! style change, never per frame: x is in bases, so the shader handles pan and zoom.

use ggv_layout::{AlignedRead, Coverage, CoverageColumn, GroupBy, PackOptions, Packed, SortBy, ViewData};

use crate::gpu::{
    FLAG_COV, FLAG_GAP, FLAG_MARKER, FLAG_SCROLL, Instance, KIND_ARROW_LEFT, KIND_ARROW_RIGHT,
};

/// Vertical layout of one track's section, in points from the section's top. Each track is
/// drawn into its own rect, so these don't depend on how many tracks are stacked above.
pub const HEADER_H: f32 = 16.0;
pub const COVERAGE_TOP: f32 = HEADER_H;
pub const COVERAGE_H: f32 = 56.0;
pub const SEPARATOR_Y: f32 = COVERAGE_TOP + COVERAGE_H + 2.0;
pub const READS_TOP: f32 = SEPARATOR_Y + 4.0;

/// Height of the shared reference row (drawn once, above all tracks).
pub const REF_H: f32 = 10.0;

/// At or above this many points per base, letters replace coloured bars.
pub const LETTERS_PPB: f32 = 7.0;
/// Below one point per base, insertions shorter than this are hidden (long-read noise).
const MIN_INSERTION_ZOOMED_OUT: u32 = 5;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ColorBy {
    /// IGV's grey.
    #[default]
    None,
    /// Forward red, reverse blue.
    Strand,
    /// `HP` tag: 1 red, 2 blue, others green, untagged grey.
    Haplotype,
    /// Insert size and pair orientation (IGV): too large red, too small blue, mate on another
    /// contig orange, unexpected orientation green / teal / dark blue.
    InsertSize,
    /// Base modifications (MM/ML): IGV's two-colour scale, red where a base is likely modified,
    /// blue where it is likely not, stronger the more certain.
    BaseMods,
    /// The value of this tag (`CB`, `PS`...): one colour per value, the same everywhere.
    Tag([u8; 2]),
}

pub use ggv_layout::classify::mod_likelihoods;

/// IGV's two-colour modification scale.
pub fn mod_color(likelihood: u8) -> u32 {
    let (r, g, b, a) = ggv_layout::classify::mod_rgba(likelihood);
    rgba(r, g, b, a)
}

/// Which reads are drawn and counted in coverage. Applied in memory for read-level views (so
/// toggling never touches the disk) and at load time for wide, coverage-only views.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadFilter {
    pub hide_duplicates: bool,
    pub hide_secondary: bool,
    pub hide_supplementary: bool,
    pub hide_qc_fail: bool,
    pub min_mapq: u8,
    /// A filter expression (`mapq >= 20 && hp == 1`), applied where reads are loaded.
    pub expr: Option<crate::expr::ExprId>,
}

impl Default for ReadFilter {
    /// IGV's defaults: duplicates, secondary and QC-fail reads hidden.
    fn default() -> Self {
        Self { hide_duplicates: true, hide_secondary: true, hide_supplementary: false, hide_qc_fail: true, min_mapq: 0, expr: None }
    }
}

impl ReadFilter {
    pub fn exclude_flags(&self) -> u16 {
        0x4 | if self.hide_secondary { 0x100 } else { 0 }
            | if self.hide_qc_fail { 0x200 } else { 0 }
            | if self.hide_duplicates { 0x400 } else { 0 }
            | if self.hide_supplementary { 0x800 } else { 0 }
    }

    /// Whether `r` passes the flag and MAPQ filters and `expr` (this filter's expression,
    /// resolved once by the caller rather than per read).
    pub fn passes(&self, r: &AlignedRead, expr: Option<&crate::expr::Expr>) -> bool {
        r.flags & self.exclude_flags() == 0 && r.mapq >= self.min_mapq && expr.is_none_or(|e| e.matches(r))
    }

    /// The same filter for the pileup engine (evidence, region statistics, binned coverage),
    /// the expression included, so every number agrees with the reads drawn.
    pub fn pileup_options(&self) -> ggv_query::PileupOptions {
        ggv_query::PileupOptions {
            exclude_flags: self.exclude_flags(),
            min_mapping_quality: self.min_mapq,
            read_filter: self.expr.map(|id| crate::expr::predicate(crate::expr::get(id), &crate::expr::text(id))),
            ..Default::default()
        }
    }

    /// [`Self::pileup_options`] for a query run outside a track's loader, with `reference` for
    /// reads that run past the query's own bases (an expression about mismatches).
    pub fn pileup_options_with(&self, reference: Option<&std::path::Path>) -> ggv_query::PileupOptions {
        let mut o = self.pileup_options();
        if o.read_filter.as_ref().is_some_and(|p| p.needs_reference) {
            o.reference = reference.map(|r| ggv_query::ReadReference::from_path(r.to_path_buf()));
        }
        o
    }
}

/// How the coverage track's height is scaled.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum CoverageScale {
    /// To the deepest visible position (IGV's autoscale).
    #[default]
    View,
    /// To the deepest position loaded (steady while panning).
    Loaded,
    /// Log scale, to the deepest visible position.
    Log,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Style {
    pub filter: ReadFilter,
    pub coverage_scale: CoverageScale,
    /// One coverage scale across all alignment tracks, for comparing depth by eye.
    pub shared_scale: bool,
    /// Dark theme (some colours are baked into instances).
    pub dark: bool,
    pub color: ColorBy,
    pub squished: bool,
    /// Splice junctions (RNA-seq) drawn as arcs over the coverage, with their read counts.
    pub junctions: bool,
    pub group: GroupBy,
    /// Sort rows by the read covering this 0-based position.
    pub sort: Option<(u64, SortBy)>,
    pub soft_clips: bool,
    /// Split reads' parts on one row, joined (IGV's "link supplementary alignments").
    pub link_splits: bool,
    /// What the zoom can show (set from the view, not saved): instances that can't show at it
    /// aren't built, and changing band rebuilds from the loaded data.
    pub detail: Detail,
}

/// How much a track's instances carry, by zoom. Loaded data serves every zoom within its span
/// (zooming in doesn't reload), so the instances are rebuilt when the zoom moves to another
/// band rather than built for all of them: below one point per base, a base is a fraction of a
/// pixel and per-base coverage bars pile up hundreds deep, so coverage is drawn in bins of a
/// power of four bases, each a quarter to one pixel wide.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Detail {
    /// Letters on bases (from [`LETTERS_PPB`]): everything.
    #[default]
    Letters,
    /// At least a point per base: per-base coverage, insertions; no letter backings.
    Bases,
    /// Below a point per base: coverage in bins of this many bases (1, 4, 16, ...), and only
    /// insertions long enough to show zoomed out.
    Binned(u32),
}

impl Detail {
    pub fn for_ppb(ppb: f64) -> Self {
        if ppb >= f64::from(LETTERS_PPB) {
            Self::Letters
        } else if ppb >= 1.0 {
            Self::Bases
        } else {
            // The largest power of four with a bin at most a pixel: 4^floor(log4(1 / ppb)).
            let mut bin = 1u32;
            while f64::from(bin) * 4.0 * ppb <= 1.0 && bin < 1 << 20 {
                bin *= 4;
            }
            Self::Binned(bin)
        }
    }

    fn letters(self) -> bool {
        self == Self::Letters
    }

    fn short_insertions(self) -> bool {
        !matches!(self, Self::Binned(_))
    }
}

impl Style {
    /// The tag reads must carry for this style: the grouping's, else the colouring's (one at a
    /// time).
    pub fn tag(&self) -> Option<[u8; 2]> {
        match (self.group, self.color) {
            (GroupBy::Tag(t), _) | (_, ColorBy::Tag(t)) => Some(t),
            _ => None,
        }
    }

    pub fn pack_options(&self) -> PackOptions {
        PackOptions { max_rows: 10_000, min_gap: 2, group: self.group, sort: self.sort, link: self.link_splits }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct RowLayout {
    pub row_h: f32,
    pub pitch: f32,
    pub rows: usize,
}

impl RowLayout {
    pub fn content_height(&self) -> f32 {
        self.rows as f32 * self.pitch
    }
}

pub fn rgba(r: u8, g: u8, b: u8, a: u8) -> u32 {
    u32::from_le_bytes([r, g, b, a])
}

pub fn base_rgb(b: u8) -> [u8; 3] {
    match b {
        b'A' => [0, 150, 0],
        b'C' => [0, 0, 255],
        b'G' => [209, 113, 5],
        b'T' => [255, 0, 0],
        _ => [100, 100, 100],
    }
}

/// [`base_rgb`] for letters drawn straight on the background: brighter on a dark one, where
/// IGV's blue C and dark orange G were nearly invisible.
pub fn base_text_rgb(b: u8, dark: bool) -> [u8; 3] {
    if !dark {
        return base_rgb(b);
    }
    match b {
        b'A' => [90, 210, 90],
        b'C' => [120, 160, 255],
        b'G' => [245, 165, 70],
        b'T' => [255, 100, 100],
        _ => [170, 170, 170],
    }
}

fn base_color(b: u8, a: u8) -> u32 {
    let [r, g, bl] = base_rgb(b);
    rgba(r, g, bl, a)
}

/// IGV shades mismatches by base quality: faint below 5, solid from 20.
pub fn quality_alpha(q: u8) -> u8 {
    match q {
        0xff => 255,
        q if q < 5 => 40,
        q if q >= 20 => 255,
        q => (40.0 + (f32::from(q) - 5.0) / 15.0 * 215.0) as u8,
    }
}

const ALWAYS: (f32, f32) = (0.0, f32::INFINITY);
const ALT_THRESHOLD: f64 = 0.2;

struct Builder {
    out: Vec<Instance>,
    origin: u64,
}

impl Builder {
    #[allow(clippy::too_many_arguments)]
    fn push(&mut self, x0: u64, x1: u64, y0: f32, y1: f32, color: u32, flags: u32, lod: (f32, f32)) {
        self.out.push(Instance {
            x0: (x0 as i64 - self.origin as i64) as f32,
            x1: (x1 as i64 - self.origin as i64) as f32,
            y0,
            y1,
            color,
            flags,
            lod_min: lod.0,
            lod_max: lod.1,
        });
    }
}

pub fn row_layout(packed: &Packed, style: Style) -> RowLayout {
    let rows = packed.rows.len();
    if style.squished {
        RowLayout { row_h: 3.0, pitch: 4.0, rows }
    } else {
        RowLayout { row_h: 10.0, pitch: 12.0, rows }
    }
}

/// The shared reference row: coloured bars from 1 point per base until letters take over.
pub fn build_reference(view: &ViewData) -> Vec<Instance> {
    let mut b = Builder { out: Vec::new(), origin: view.region.start };
    if let Some(reference) = &view.reference {
        for (i, &base) in reference.iter().enumerate() {
            let p = view.region.start + i as u64;
            b.push(p, p + 1, 1.5, REF_H - 1.5, base_color(base, 255), crate::gpu::FLAG_ROUND | crate::gpu::FLAG_GAP, (1.0, LETTERS_PPB));
        }
    }
    b.out
}

/// Packs a view's reads for `style` (sorting and grouping).
/// Indices of the reads that pass the filter.
pub fn filtered(view: &ViewData, style: Style) -> Vec<usize> {
    let expr = style.filter.expr.map(crate::expr::get);
    view.reads.as_ref().map_or_else(Vec::new, |(reads, _)| {
        (0..reads.len()).filter(|&i| style.filter.passes(&reads[i], expr.as_deref())).collect()
    })
}

/// Packs the filtered reads for `style` (sorting and grouping).
pub fn pack(view: &ViewData, style: Style) -> Packed {
    view.reads
        .as_ref()
        .map(|(reads, _)| ggv_layout::pack_indices(reads, filtered(view, style), &style.pack_options()))
        .unwrap_or_default()
}

/// [`pack`], keeping reads in the rows `hints` gives them (the previous span's packing), so
/// panning onto a newly loaded span doesn't move reads already on screen.
pub fn pack_sticky(view: &ViewData, style: Style, hints: &ggv_layout::RowHints) -> Packed {
    view.reads
        .as_ref()
        .map(|(reads, _)| ggv_layout::pack_sticky(reads, filtered(view, style), &style.pack_options(), hints))
        .unwrap_or_default()
}

/// Coverage for `style`: recomputed from the filtered reads when reads are loaded (so filters
/// change coverage too), otherwise the binned depth as loaded.
pub fn coverage_for(view: &ViewData, style: Style) -> Coverage {
    match &view.reads {
        Some((reads, _)) => {
            let keep = filtered(view, style);
            // Every read passes: the coverage counted while loading is the same.
            if keep.len() == reads.len() {
                return view.coverage.clone();
            }
            Coverage::PerBase(ggv_layout::coverage_from_reads(keep.into_iter().map(|i| &reads[i]), &view.region, view.reference.as_deref()))
        }
        None => view.coverage.clone(),
    }
}

/// The height unit for a depth: linear, or log2(1 + depth) for the log scale.
pub fn depth_units(depth: f32, scale: CoverageScale) -> f32 {
    match scale {
        CoverageScale::Log => (1.0 + depth).log2(),
        _ => depth,
    }
}

/// The deepest per-base coverage in each block of [`DepthIndex::BLOCK`] bases, so the deepest
/// in a view takes a few thousand comparisons rather than one per base: over five megabases of
/// per-base coverage, scanning every column took 6-10 ms a track, every frame.
#[derive(Default)]
pub struct DepthIndex {
    blocks: Vec<u32>,
}

impl DepthIndex {
    const BLOCK: usize = 1024;

    pub fn of(coverage: &Coverage) -> Self {
        match coverage {
            Coverage::PerBase(cols) => Self { blocks: cols.chunks(Self::BLOCK).map(|c| c.iter().map(|x| x.depth).max().unwrap_or(0)).collect() },
            _ => Self::default(),
        }
    }

    /// The deepest of `cols[a..b]`: whole blocks from the index, the ends column by column.
    fn max(&self, cols: &[CoverageColumn], a: usize, b: usize) -> u32 {
        if self.blocks.is_empty() || b <= a {
            return cols[a..b.max(a)].iter().map(|c| c.depth).max().unwrap_or(0);
        }
        let (first, last) = (a.div_ceil(Self::BLOCK), b / Self::BLOCK);
        if first >= last {
            return cols[a..b].iter().map(|c| c.depth).max().unwrap_or(0);
        }
        let edges = cols[a..first * Self::BLOCK].iter().chain(&cols[last * Self::BLOCK..b]).map(|c| c.depth).max().unwrap_or(0);
        self.blocks[first..last].iter().copied().max().unwrap_or(0).max(edges)
    }
}

/// The deepest coverage in `[vs, ve)` (or everywhere, for `Loaded`), in height units.
pub fn max_units(coverage: &Coverage, index: &DepthIndex, start: u64, end: u64, vs: f64, ve: f64, scale: CoverageScale) -> f32 {
    let (lo, hi) = match scale {
        CoverageScale::Loaded => (start as f64, end as f64),
        _ => (vs, ve),
    };
    let max = match coverage {
        Coverage::PerBase(cols) => {
            let a = (lo.floor() - start as f64).clamp(0.0, cols.len() as f64) as usize;
            let b = (hi.ceil() - start as f64).clamp(0.0, cols.len() as f64) as usize;
            index.max(cols, a, b)
        }
        Coverage::Binned(bins) | Coverage::Estimated(bins) => {
            let n = bins.len() as f64;
            let len = (end - start).max(1) as f64;
            let a = (((lo - start as f64) / len * n).floor()).clamp(0.0, n) as usize;
            let b = (((hi - start as f64) / len * n).ceil()).clamp(0.0, n) as usize;
            bins[a..b.max(a)].iter().map(|b| b.max).max().unwrap_or(0)
        }
        Coverage::Unavailable => 0,
    };
    depth_units(max as f32, scale).max(1.0)
}

/// One alignment track: coverage and reads, positioned within the track's own section.
pub fn build(view: &ViewData, coverage: &Coverage, packed: &Packed, style: Style) -> Vec<Instance> {
    let mut b = Builder { out: Vec::new(), origin: view.region.start };
    coverage_instances(&mut b, view, coverage, style.coverage_scale, style.detail);

    if let Some((reads, _)) = &view.reads {
        let rl = row_layout(packed, style);
        let sizes = if style.color == ColorBy::InsertSize { InsertSizes::of(reads) } else { None };
        // Zoomed out to several bases a pixel: which columns hold a variant (see `Coarse`).
        let variable: Option<Vec<bool>> = match (style.detail, coverage) {
            (Detail::Binned(n), Coverage::PerBase(cols)) if n >= 4 => Some(cols.iter().map(snv_column).collect()),
            _ => None,
        };
        let coarse = match (style.detail, &variable) {
            (Detail::Binned(n), Some(v)) => Some(Coarse { bases: u64::from(n), variable: v, origin: view.region.start }),
            _ => None,
        };

        for (r, row) in packed.rows.iter().enumerate() {
            let y = READS_TOP + r as f32 * rl.pitch;
            for &i in row {
                read(&mut b, &reads[i], y, rl.row_h, style, sizes.as_ref(), coarse.as_ref());
            }
            // Linked split reads: a thin line from each part to the next.
            if style.link_splits {
                let mid = y + rl.row_h / 2.0;
                for w in row.windows(2) {
                    let (a, c) = (&reads[w[0]], &reads[w[1]]);
                    if ggv_layout::same_template(a, c) && c.start > a.end {
                        b.push(a.end, c.start, mid - 0.5, mid + 0.5, rgba(70, 110, 200, 255), FLAG_SCROLL, ALWAYS);
                    }
                }
            }
        }
    }
    b.out
}

/// The expected insert-size range: the 0.5th to 99.5th percentile of |TLEN| over same-contig
/// pairs in the loaded span (IGV's default). `None` with too few pairs to say.
pub struct InsertSizes {
    pub low: u32,
    pub high: u32,
}

impl InsertSizes {
    pub fn of(reads: &[AlignedRead]) -> Option<Self> {
        ggv_layout::classify::insert_bounds(reads).map(|(low, high)| Self { low, high })
    }
}

/// Pair orientation and insert size, as colours (see `ggv_layout::classify`).
fn pair_color(r: &AlignedRead, sizes: Option<&InsertSizes>) -> Option<u32> {
    use ggv_layout::classify::PairClass::*;
    Some(match ggv_layout::classify::pair_class(r, sizes.map(|s| (s.low, s.high)))? {
        OtherContig => rgba(240, 140, 40, 255),
        RightLeft => rgba(0, 150, 0, 255),
        LeftLeft => rgba(0, 150, 150, 255),
        RightRight => rgba(20, 50, 200, 255),
        Large => rgba(210, 30, 30, 255),
        Small => rgba(30, 30, 210, 255),
    })
}

/// Coverage bars in height units above the baseline (`FLAG_COV`): the shader scales them each
/// frame, so autoscaling to the view never rebuilds anything.
fn coverage_instances(b: &mut Builder, view: &ViewData, coverage: &Coverage, scale: CoverageScale, detail: Detail) {
    let bottom = COVERAGE_TOP + COVERAGE_H;
    let grey = rgba(175, 175, 175, 255);
    let units = |d: f32| depth_units(d, scale);

    match coverage {
        // Zoomed out: a bar per bin at its deepest base, and the variable columns over them as
        // at any zoom, so variants stay visible.
        Coverage::PerBase(cols) if matches!(detail, Detail::Binned(n) if n > 1) => {
            let Detail::Binned(n) = detail else { unreachable!() };
            for (k, bin) in cols.chunks(n as usize).enumerate() {
                let p = view.region.start + k as u64 * u64::from(n);
                let max = bin.iter().map(|c| c.depth).max().unwrap_or(0);
                if max > 0 {
                    b.push(p, p + bin.len() as u64, bottom - units(max as f32), bottom, grey, FLAG_COV, ALWAYS);
                }
                // Several bases a pixel: only columns with a real single-base variant (see
                // `snv_column`), not every column a noisy long-read pileup tips over 20%.
                for (j, col) in bin.iter().enumerate().filter(|(_, c)| if n >= 4 { snv_column(c) } else { c.depth > 0 && c.mismatch_fraction() >= ALT_THRESHOLD }) {
                    stack_alleles(b, col, p + j as u64, bottom, units(col.depth as f32));
                }
            }
        }
        Coverage::PerBase(cols) => {
            for (i, col) in cols.iter().enumerate() {
                if col.depth == 0 {
                    continue;
                }
                let p = view.region.start + i as u64;
                let total = units(col.depth as f32);
                if col.mismatch_fraction() >= ALT_THRESHOLD {
                    stack_alleles(b, col, p, bottom, total);
                } else {
                    b.push(p, p + 1, bottom - total, bottom, grey, FLAG_GAP | FLAG_COV, ALWAYS);
                }
            }
        }
        Coverage::Binned(bins) => {
            let n = bins.len() as u64;
            let len = view.region.len();
            for (i, bin) in bins.iter().enumerate().filter(|(_, b)| b.max > 0) {
                let x0 = view.region.start + i as u64 * len / n;
                let x1 = view.region.start + (i as u64 + 1) * len / n;
                b.push(x0, x1, bottom - units(bin.max as f32), bottom, rgba(215, 215, 215, 255), FLAG_COV, ALWAYS);
                b.push(x0, x1, bottom - units(bin.mean as f32), bottom, grey, FLAG_COV, ALWAYS);
            }
        }
        Coverage::Estimated(bins) => {
            // A lighter, bluer grey: an estimate, not a count.
            let n = bins.len() as u64;
            let len = view.region.len();
            for (i, bin) in bins.iter().enumerate().filter(|(_, b)| b.mean > 0.0) {
                let x0 = view.region.start + i as u64 * len / n;
                let x1 = view.region.start + (i as u64 + 1) * len / n;
                b.push(x0, x1, bottom - units(bin.mean as f32), bottom, rgba(160, 175, 200, 255), FLAG_COV, ALWAYS);
            }
        }
        Coverage::Unavailable => {}
    }
}

/// Allele segments share the column's height in proportion to their counts (so the log scale
/// keeps the proportions readable).
fn stack_alleles(b: &mut Builder, col: &CoverageColumn, p: u64, bottom: f32, total: f32) {
    let per_read = total / col.depth.max(1) as f32;
    let mut top = bottom;
    for (k, base) in [b'A', b'C', b'G', b'T', b'N'].into_iter().enumerate() {
        let seg = col.counts[k] as f32 * per_read;
        if seg > 0.0 {
            b.push(p, p + 1, top - seg, top, base_color(base, 255), FLAG_GAP | FLAG_COV, ALWAYS);
            top -= seg;
        }
    }
    if col.deletions > 0 {
        let seg = col.deletions as f32 * per_read;
        b.push(p, p + 1, top - seg, top, rgba(175, 175, 175, 255), FLAG_GAP | FLAG_COV, ALWAYS);
    }
}

pub fn read_fill(r: &AlignedRead, style: Style, sizes: Option<&InsertSizes>) -> u32 {
    let fill = read_colour(r, style, sizes);
    // MAPQ 0 (the read fits elsewhere equally well): its colour, faded, so strand, haplotype and
    // so on still show but the read reads as unplaced. Alpha, not a mix with white, so it fades
    // into the dark theme's background too. It was flat grey, which hid the colouring.
    if r.mapq == 0 { (fill & 0x00ff_ffff) | (MAPQ0_ALPHA << 24) } else { fill }
}

/// How opaque a MAPQ 0 read is (of 255).
pub const MAPQ0_ALPHA: u32 = 95;

fn read_colour(r: &AlignedRead, style: Style, sizes: Option<&InsertSizes>) -> u32 {
    let grey = rgba(185, 185, 185, 255);
    match style.color {
        ColorBy::None => grey,
        ColorBy::Strand if r.reverse => rgba(160, 160, 230, 255),
        ColorBy::Strand => rgba(230, 160, 160, 255),
        ColorBy::Haplotype => match r.haplotype {
            Some(1) => rgba(230, 150, 150, 255),
            Some(2) => rgba(150, 150, 230, 255),
            Some(_) => rgba(150, 210, 150, 255),
            None => grey,
        },
        ColorBy::InsertSize => pair_color(r, sizes).unwrap_or(grey),
        ColorBy::BaseMods => grey,
        // Tagged reads by value (when the grouping's tag is another, the colouring has none).
        ColorBy::Tag(t) if style.tag() == Some(t) => r.tag.as_deref().map_or(grey, |v| {
            let (red, green, blue) = ggv_layout::classify::tag_rgb(v);
            rgba(red, green, blue, 255)
        }),
        ColorBy::Tag(_) => grey,
    }
}

/// Zoomed out to `bases` per pixel or more: what a read draws is cut to what can show. A long
/// read's aligned blocks (split by every small indel: hundreds for a nanopore read) join into
/// bars wherever the gap is under a pixel; deletions and introns under a pixel aren't drawn;
/// mismatches are drawn only in columns holding a variant (`variable`, from the coverage), not
/// every sequencing error, which at this zoom is a sub-pixel speck. A few megabases of long
/// reads went from ten million instances to a few hundred thousand. Zooming in rebuilds with
/// everything.
struct Coarse<'a> {
    bases: u64,
    variable: &'a [bool],
    origin: u64,
}

fn read(b: &mut Builder, r: &AlignedRead, y: f32, h: f32, style: Style, sizes: Option<&InsertSizes>, coarse: Option<&Coarse>) {
    let s = FLAG_SCROLL;
    let mid = y + h / 2.0;
    let fill = read_fill(r, style, sizes);
    if let Some(c) = coarse {
        return coarse_read(b, r, y, h, style, fill, c);
    }

    // Soft clips: the clipped bases where they would align, paler than the read, with their
    // mismatches against the reference.
    if style.soft_clips {
        // Paler than the read itself, whatever its own opacity (MAPQ 0 reads are faded too).
        let pale = (fill & 0x00ff_ffff) | (((fill >> 24) * 110 / 255) << 24);
        if r.soft_clip_left > 0 {
            b.push(r.start.saturating_sub(u64::from(r.soft_clip_left)), r.start, y, y + h, pale, s, ALWAYS);
        }
        if r.soft_clip_right > 0 {
            b.push(r.end, r.end + u64::from(r.soft_clip_right), y, y + h, pale, s, ALWAYS);
        }
        for m in &r.clip_mismatches {
            let alpha = quality_alpha(m.quality);
            b.push(m.pos, m.pos + 1, y, y + h, base_color(m.base, alpha), s, (0.0, LETTERS_PPB));
            if style.detail.letters() {
                b.push(m.pos, m.pos + 1, y, y + h, rgba(255, 255, 255, alpha / 2 + 60), s, (LETTERS_PPB, f32::INFINITY));
            }
        }
    }

    for skip in &r.ref_skips {
        b.push(skip.pos, skip.pos + u64::from(skip.len), mid, mid + 1.0, rgba(150, 150, 150, 255), s, ALWAYS);
    }
    for d in &r.deletions {
        let line = if style.dark { rgba(230, 230, 230, 255) } else { rgba(30, 30, 30, 255) };
        b.push(d.pos, d.pos + u64::from(d.len), mid, mid + 1.0, line, s, ALWAYS);
    }
    for block in &r.blocks {
        b.push(block.ref_start, block.ref_start + u64::from(block.len), y, y + h, fill, s, ALWAYS);
    }
    if !style.squished {
        let (anchor, kind) = if r.reverse { (r.start, KIND_ARROW_LEFT) } else { (r.end, KIND_ARROW_RIGHT) };
        b.push(anchor, anchor, y, y + h, fill, s | kind, ALWAYS);
    }

    if style.color == ColorBy::BaseMods {
        match style.detail {
            Detail::Binned(n) if n > 1 => mod_bins(b, r, u64::from(n), y, h),
            _ => {
                for (pos, likelihood) in mod_likelihoods(r) {
                    b.push(pos, pos + 1, y, y + h, mod_color(likelihood), s, ALWAYS);
                }
            }
        }
    }

    for m in &r.mismatches {
        let alpha = quality_alpha(m.quality);
        // Bars when zoomed out; when zoomed in, a pale backing for the letter drawn by egui.
        b.push(m.pos, m.pos + 1, y, y + h, base_color(m.base, alpha), s, (0.0, LETTERS_PPB));
        if style.detail.letters() {
            b.push(m.pos, m.pos + 1, y, y + h, rgba(255, 255, 255, alpha / 2 + 60), s, (LETTERS_PPB, f32::INFINITY));
        }
    }

    for ins in r.insertions.iter().filter(|i| style.detail.short_insertions() || i.len >= MIN_INSERTION_ZOOMED_OUT) {
        let lod = if ins.len >= MIN_INSERTION_ZOOMED_OUT { ALWAYS } else { (1.0, f32::INFINITY) };
        b.push(ins.pos, ins.pos, y, y + h, rgba(138, 94, 161, 255), s | FLAG_MARKER, lod);
    }
}

/// Zoomed out, a column shows a variant when one alternate base is in at least this share of the
/// reads: a heterozygous SNV is near half, while a nanopore pileup's systematic errors put one
/// column in 70 over 20% (45,000 columns in 3 Mb of HG002, against about 2,000 real SNVs).
const COARSE_ALT_THRESHOLD: f64 = 0.35;

/// A column holding a single-base variant (see [`COARSE_ALT_THRESHOLD`]). Zoomed out, this,
/// rather than the total mismatch fraction, marks variant columns: deletions and scattered
/// errors in long-read pileups put one column in 25 over 20% non-reference.
fn snv_column(c: &CoverageColumn) -> bool {
    let Some(r) = c.ref_base else { return false };
    let ref_k = [b'A', b'C', b'G', b'T'].iter().position(|&b| b == r.to_ascii_uppercase());
    let top = c.counts[..4].iter().enumerate().filter(|(k, _)| Some(*k) != ref_k).map(|(_, n)| *n).max().unwrap_or(0);
    c.depth > 0 && f64::from(top) >= COARSE_ALT_THRESHOLD * f64::from(c.depth)
}

/// [`read`] zoomed out (see [`Coarse`]).
fn coarse_read(b: &mut Builder, r: &AlignedRead, y: f32, h: f32, style: Style, fill: u32, c: &Coarse) {
    let s = FLAG_SCROLL;
    let mid = y + h / 2.0;
    if style.soft_clips {
        let pale = (fill & 0x00ff_ffff) | (((fill >> 24) * 110 / 255) << 24);
        if r.soft_clip_left > 0 {
            b.push(r.start.saturating_sub(u64::from(r.soft_clip_left)), r.start, y, y + h, pale, s, ALWAYS);
        }
        if r.soft_clip_right > 0 {
            b.push(r.end, r.end + u64::from(r.soft_clip_right), y, y + h, pale, s, ALWAYS);
        }
    }
    for skip in r.ref_skips.iter().filter(|k| u64::from(k.len) >= c.bases) {
        b.push(skip.pos, skip.pos + u64::from(skip.len), mid, mid + 1.0, rgba(150, 150, 150, 255), s, ALWAYS);
    }
    let line = if style.dark { rgba(230, 230, 230, 255) } else { rgba(30, 30, 30, 255) };
    for d in r.deletions.iter().filter(|d| u64::from(d.len) >= c.bases) {
        b.push(d.pos, d.pos + u64::from(d.len), mid, mid + 1.0, line, s, ALWAYS);
    }
    // Blocks joined across gaps under a pixel.
    let mut bar: Option<(u64, u64)> = None;
    for block in &r.blocks {
        let (s0, e0) = (block.ref_start, block.ref_start + u64::from(block.len));
        bar = match bar {
            Some((a, e)) if s0 < e + c.bases => Some((a, e.max(e0))),
            Some((a, e)) => {
                b.push(a, e, y, y + h, fill, s, ALWAYS);
                Some((s0, e0))
            }
            None => Some((s0, e0)),
        };
    }
    if let Some((a, e)) = bar {
        b.push(a, e, y, y + h, fill, s, ALWAYS);
    }
    if !style.squished {
        let (anchor, kind) = if r.reverse { (r.start, KIND_ARROW_LEFT) } else { (r.end, KIND_ARROW_RIGHT) };
        b.push(anchor, anchor, y, y + h, fill, s | kind, ALWAYS);
    }
    if style.color == ColorBy::BaseMods {
        mod_bins(b, r, c.bases, y, h);
    }
    for m in r.mismatches.iter().filter(|m| m.pos.checked_sub(c.origin).and_then(|i| c.variable.get(i as usize)).copied().unwrap_or(false)) {
        b.push(m.pos, m.pos + 1, y, y + h, base_color(m.base, quality_alpha(m.quality)), s, ALWAYS);
    }
    for ins in r.insertions.iter().filter(|i| i.len >= MIN_INSERTION_ZOOMED_OUT) {
        b.push(ins.pos, ins.pos, y, y + h, rgba(138, 94, 161, 255), s | FLAG_MARKER, ALWAYS);
    }
}

/// Zoomed out, a read's modification calls in bins of `n` bases (as coverage is binned), each
/// coloured by its calls' mean likelihood: per base, reads with calls on every A and C made
/// tens of millions of instances for a few megabases.
fn mod_bins(b: &mut Builder, r: &AlignedRead, n: u64, y: f32, h: f32) {
    let mut flush = |bin: u64, sum: u32, calls: u32| {
        let (x0, x1) = ((bin * n).max(r.start), ((bin + 1) * n).min(r.end));
        b.push(x0, x1.max(x0 + 1), y, y + h, mod_color((sum / calls) as u8), FLAG_SCROLL, ALWAYS);
    };
    let mut bin: Option<(u64, u32, u32)> = None;
    for (pos, likelihood) in mod_likelihoods(r) {
        match &mut bin {
            Some((i, sum, calls)) if *i == pos / n => {
                *sum += u32::from(likelihood);
                *calls += 1;
            }
            _ => {
                if let Some((i, sum, calls)) = bin {
                    flush(i, sum, calls);
                }
                bin = Some((pos / n, u32::from(likelihood), 1));
            }
        }
    }
    if let Some((i, sum, calls)) = bin {
        flush(i, sum, calls);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn depth_index_finds_the_same_maximum_as_a_scan() {
        let cols: Vec<CoverageColumn> = (0..10_000u32).map(|i| CoverageColumn { depth: (i * 7919) % 997, ..Default::default() }).collect();
        let cov = Coverage::PerBase(cols.clone());
        let ix = DepthIndex::of(&cov);
        for (a, b) in [(0, 10_000), (5, 6), (1000, 3100), (1023, 1025), (2048, 4096), (9_999, 10_000), (3, 3), (500, 1500)] {
            let scan = cols[a..b].iter().map(|c| c.depth).max().unwrap_or(0);
            assert_eq!(ix.max(&cols, a, b), scan, "{a}..{b}");
        }
    }

    #[test]
    fn zoom_bands() {
        assert_eq!(Detail::for_ppb(8.0), Detail::Letters);
        assert_eq!(Detail::for_ppb(1.0), Detail::Bases);
        // Bins of a power of four, a quarter to one pixel wide.
        assert_eq!(Detail::for_ppb(0.5), Detail::Binned(1));
        assert_eq!(Detail::for_ppb(0.2), Detail::Binned(4));
        // A 310 kb view in 1,400 points: 0.0045 points per base, 64-base bins (0.29 px).
        assert_eq!(Detail::for_ppb(1400.0 / 310_000.0), Detail::Binned(64));
        for ppb in [0.9, 0.3, 0.01, 0.0001] {
            let Detail::Binned(n) = Detail::for_ppb(ppb) else { panic!() };
            let px = f64::from(n) * ppb;
            assert!(px <= 1.0 && (n == 1 || px > 0.25), "{ppb}: bin {n} is {px} px");
        }
    }

    fn paired(start: u64, reverse: bool, mate_pos: i64, mate_reverse: bool, tlen: i32) -> AlignedRead {
        AlignedRead {
            name: String::new(),
            start,
            end: start + 100,
            reverse,
            mapq: 60,
            flags: 0x1,
            blocks: vec![],
            deletions: vec![],
            ref_skips: vec![],
            insertions: vec![],
            mismatches: vec![],
            soft_clip_left: 0,
            soft_clip_right: 0,
            clip_mismatches: vec![],
            haplotype: None,
            read_group: None,
            template_length: tlen,
            mate: Some(ggv_layout::Mate { same_contig: true, contig_id: 0, pos: mate_pos, reverse: mate_reverse }),
            supplementary: vec![],
            mods: vec![],
            low_quality: vec![],
            tag: None,
            bases: None,
        }
    }

    #[test]
    fn pair_orientation_and_insert_size_colours() {
        let sizes = InsertSizes { low: 200, high: 600 };
        let s = Some(&sizes);
        // Normal forward-reverse pair of 400 bp: default grey.
        assert_eq!(pair_color(&paired(1000, false, 1300, true, 400), s), None);
        // Same pair seen from the reverse mate.
        assert_eq!(pair_color(&paired(1300, true, 1000, false, -400), s), None);
        assert_eq!(pair_color(&paired(1000, false, 5000, true, 4100), s), Some(rgba(210, 30, 30, 255)), "too large");
        assert_eq!(pair_color(&paired(1000, false, 1050, true, 150), s), Some(rgba(30, 30, 210, 255)), "too small");
        assert_eq!(pair_color(&paired(1000, true, 1300, false, 400), s), Some(rgba(0, 150, 0, 255)), "RL");
        // Overlapping mates where the reverse one starts first (read-through): not RL.
        assert_eq!(pair_color(&paired(1000, true, 1020, false, -250), s), None, "overlapping read-through");
        assert_eq!(pair_color(&paired(1000, false, 1300, false, 400), s), Some(rgba(0, 150, 150, 255)), "LL");
        let mut other = paired(1000, false, 1300, true, 0);
        other.mate = Some(ggv_layout::Mate { same_contig: false, contig_id: 1, pos: 5, reverse: true });
        assert_eq!(pair_color(&other, s), Some(rgba(240, 140, 40, 255)), "mate elsewhere");
    }

    #[test]
    fn in_memory_filtering_matches_load_time_filtering() {
        // GGV_TEST_BAM + GGV_TEST_REGIONS: loading everything mapped and filtering in memory with
        // the default filter must give the same coverage as filtering at load (which matches
        // samtools), and filters must actually change it.
        let (Ok(bam_path), Ok(regions)) = (std::env::var("GGV_TEST_BAM"), std::env::var("GGV_TEST_REGIONS")) else { return };
        for locus in regions.split(',') {
            let mut bam = ggv_io::BamFile::open(&bam_path).unwrap();
            let region = bam.resolve_str(locus).unwrap();
            let classic = ggv_layout::load_view(&mut bam, &region, None, String::new(), &ggv_layout::ViewOptions::default()).unwrap();
            let opts = ggv_layout::ViewOptions {
                filters: ggv_query::PileupOptions { exclude_flags: 0x4, ..Default::default() },
                ..Default::default()
            };
            let all = ggv_layout::load_view(&mut bam, &region, None, String::new(), &opts).unwrap();
            let style = Style::default();
            assert_eq!(coverage_for(&all, style), classic.coverage, "{locus}: default filter");
            assert_eq!(filtered(&all, style).len(), classic.reads.as_ref().unwrap().0.len(), "{locus}: reads kept");

            let mut keep_dups = style;
            keep_dups.filter.hide_duplicates = false;
            let more = filtered(&all, keep_dups).len();
            let mut strict = style;
            strict.filter.min_mapq = 30;
            let fewer = filtered(&all, strict).len();
            eprintln!("{locus}: {} reads by default, {more} with duplicates, {fewer} at MAPQ >= 30", filtered(&all, style).len());
            assert!(more >= filtered(&all, style).len() && fewer <= filtered(&all, style).len());
        }
    }

    #[test]
    fn expression_filters_reach_the_pileup_engine() {
        // GGV_TEST_BAM + GGV_TEST_REGIONS: with an expression active, the pileup engine (evidence,
        // region statistics, binned coverage) must count exactly the reads the view draws.
        let (Ok(bam_path), Ok(regions)) = (std::env::var("GGV_TEST_BAM"), std::env::var("GGV_TEST_REGIONS")) else { return };
        for text in ["mapq >= 60 && !reverse", "insert > 300 || !paired", "len < 140"] {
            let mut style = Style::default();
            style.filter.expr = Some(crate::expr::intern(text).unwrap());
            for locus in regions.split(',') {
                let mut bam = ggv_io::BamFile::open(&bam_path).unwrap();
                let region = bam.resolve_str(locus).unwrap();
                let opts = ggv_layout::ViewOptions { filters: ggv_query::PileupOptions { exclude_flags: 0x4, ..Default::default() }, ..Default::default() };
                let all = ggv_layout::load_view(&mut bam, &region, None, String::new(), &opts).unwrap();
                let Coverage::PerBase(drawn) = coverage_for(&all, style) else { panic!() };
                let mut counted = vec![0u32; drawn.len()];
                let reads = ggv_query::pileup(&mut bam, &region, None, &style.filter.pileup_options(), |p| {
                    counted[(p.pos - 1 - region.start) as usize] = p.depth;
                    Ok(())
                })
                .unwrap();
                let depth: Vec<u32> = drawn.iter().map(|c| c.depth).collect();
                assert_eq!(depth, counted, "{text} at {locus}");
                assert!(reads.filtered_out > 0, "{text}: the expression should remove some reads");
                eprintln!("{text} at {locus}: {} used, {} removed by the expression", reads.used, reads.filtered_out);
            }
        }
    }

    #[test]
    fn modification_likelihoods_sum_codes_per_base() {
        let mut r = paired(0, false, 0, true, 0);
        let call = |offset, code, prob| ggv_layout::ModCall { offset, code, canonical: b'C', prob };
        r.mods = vec![call(5, b'h', 20), call(5, b'm', 200), call(9, b'm', 10), call(12, b'm', 250), call(12, b'h', 250)];
        assert_eq!(mod_likelihoods(&r).collect::<Vec<_>>(), vec![(5, 220), (9, 10), (12, 255)]);
        // Two-colour scale: strong red, strong blue, faint either side of 50%.
        assert_eq!(mod_color(255), rgba(220, 20, 20, 254));
        assert_eq!(mod_color(0), rgba(20, 20, 220, 255));
        assert_eq!(mod_color(130) >> 24, 4);
        assert_eq!(mod_color(126) >> 24, 4);
    }

    #[test]
    fn modification_calls_bin_zoomed_out() {
        let mut r = paired(10, false, 0, true, 0);
        let call = |pos: u32, prob| ggv_layout::ModCall { offset: pos - 10, code: b'm', canonical: b'C', prob };
        // Bins of 16: [10, 16) gets 12 and 14, [16, 32) gets 20, [96, 110) gets 100.
        r.mods = vec![call(12, 200), call(14, 100), call(20, 0), call(100, 255)];
        let mut b = Builder { out: Vec::new(), origin: 0 };
        mod_bins(&mut b, &r, 16, 0.0, 10.0);
        let spans: Vec<(f32, f32, u32)> = b.out.iter().map(|i| (i.x0, i.x1, i.color)).collect();
        assert_eq!(spans, vec![(10.0, 16.0, mod_color(150)), (16.0, 32.0, mod_color(0)), (96.0, 110.0, mod_color(255))]);
    }

    #[test]
    fn color_packing_is_rgba_little_endian() {
        assert_eq!(rgba(1, 2, 3, 4), 0x0403_0201);
    }

    #[test]
    fn quality_shading_matches_png_renderer() {
        assert_eq!(quality_alpha(2), 40);
        assert_eq!(quality_alpha(30), 255);
        assert_eq!(quality_alpha(0xff), 255);
    }
}
