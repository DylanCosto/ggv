//! A VCF/BCF track: a site row coloured by allele frequency (or by an INFO field: ClinVar's
//! clinical significance) and one genotype row per sample, drawn as GPU instances. Loading and next/previous-variant searches run on the track's own
//! thread.

use std::{
    path::PathBuf,
    sync::{Arc, mpsc},
    thread,
    time::{Duration, Instant},
};

use eframe::egui;
use ggv_io::{GenotypeClass, InfoField, Variant, VariantFile};
use ggv_model::Region;

use crate::{
    gpu::{FLAG_GAP, FLAG_MIN3, FLAG_ROUND, FLAG_SCROLL, Instance, InstanceSet},
    tracks::{HEADER_H, rgba},
};

pub const SITE_TOP: f32 = HEADER_H + 2.0;
pub const SITE_H: f32 = 12.0;
/// Site rows: overlapping records (an SNV inside a deletion, several ALTs at one base) stack,
/// as in IGV's expanded display.
pub const SITE_PITCH: f32 = 14.0;
/// Rows beyond this go into the last one (still overlapping), counted in the header.
const MAX_SITE_ROWS: usize = 12;
pub const GT_PITCH: f32 = 11.0;
pub const GT_H: f32 = 10.0;
/// Sample rows shown before the section scrolls.
const VISIBLE_SAMPLES: usize = 10;
/// Variants load for views up to this wide.
pub const VARIANT_WINDOW: u64 = 10_000_000;
/// Load at least this much around the view.
const MIN_SPAN: f64 = 1_000_000.0;

/// The dark theme's track background, which the dark colours below keep at least 3:1 against.
pub const DARK_BG: [u8; 3] = [24, 24, 24];

/// A genotype's colour (IGV's); on the dark theme het is a lighter blue: IGV's #220CFD was
/// about 2:1 against the background.
pub fn genotype_rgb(c: GenotypeClass, dark: bool) -> Option<[u8; 3]> {
    match c {
        GenotypeClass::HomRef => Some([200, 200, 200]),
        GenotypeClass::Het if dark => Some([100, 130, 255]),
        GenotypeClass::Het => Some([34, 12, 253]),
        GenotypeClass::HomAlt => Some([17, 248, 254]),
        GenotypeClass::NoCall => Some([235, 235, 210]),
        GenotypeClass::Unknown => None,
    }
}

/// The allele-frequency site colours (IGV's): the alternate share and the rest.
fn af_rgb(dark: bool) -> ([u8; 3], [u8; 3]) {
    if dark { ([230, 40, 40], [80, 100, 255]) } else { ([220, 0, 0], [0, 0, 220]) }
}

/// WCAG relative luminance of an sRGB colour.
fn luminance([r, g, b]: [u8; 3]) -> f64 {
    let lin = |c: u8| {
        let c = f64::from(c) / 255.0;
        if c <= 0.039_28 { c / 12.92 } else { ((c + 0.055) / 1.055).powf(2.4) }
    };
    0.2126 * lin(r) + 0.7152 * lin(g) + 0.0722 * lin(b)
}

/// WCAG contrast ratio between two colours (1 to 21).
pub fn contrast(a: [u8; 3], b: [u8; 3]) -> f64 {
    let (la, lb) = (luminance(a), luminance(b));
    (la.max(lb) + 0.05) / (la.min(lb) + 0.05)
}

/// `rgb` mixed toward white just enough to stand at least 3:1 against the dark background.
fn lift_for_dark(rgb: [u8; 3]) -> [u8; 3] {
    (0..=10)
        .map(|step| {
            let t = f64::from(step) / 10.0;
            rgb.map(|c| (f64::from(c) + (255.0 - f64::from(c)) * t).round() as u8)
        })
        .find(|&c| contrast(c, DARK_BG) >= 3.0)
        .unwrap_or([255, 255, 255])
}

/// What colours a site.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SiteColor {
    /// Red for the alternate allele's share (INFO/AF, else from genotypes), blue for the rest.
    AlleleFrequency,
    /// One colour per value of this INFO field: ClinVar's classes for `CLNSIG`, else a palette.
    Info(String),
}

impl SiteColor {
    /// The default for a file: ClinVar's significance when it has one.
    pub fn default_for(fields: &[InfoField]) -> Self {
        if fields.iter().any(|f| f.id == "CLNSIG") { Self::Info("CLNSIG".into()) } else { Self::AlleleFrequency }
    }

    /// As sessions and the MCP tools name it: `af`, or the INFO key.
    pub fn name(&self) -> &str {
        match self {
            Self::AlleleFrequency => "af",
            Self::Info(k) => k,
        }
    }

    pub fn parse(s: &str) -> Self {
        if s.eq_ignore_ascii_case("af") || s.is_empty() { Self::AlleleFrequency } else { Self::Info(s.to_owned()) }
    }
}

/// ClinVar's clinical significance in five classes (UCSC's ClinVar track colours).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Significance {
    Pathogenic,
    Uncertain,
    Benign,
    Conflicting,
    Other,
}

impl Significance {
    pub const ALL: [Self; 5] = [Self::Pathogenic, Self::Uncertain, Self::Benign, Self::Conflicting, Self::Other];

    /// From a CLNSIG value: `Pathogenic`, `Likely_pathogenic`, `Pathogenic/Likely_pathogenic`,
    /// `Benign/Likely_benign`, `Conflicting_classifications_of_pathogenicity`,
    /// `Uncertain_significance|risk_factor` (the part before `|` decides)...
    pub fn of(clnsig: &str) -> Self {
        let main = clnsig.split('|').next().unwrap_or("").to_ascii_lowercase();
        if main.contains("conflicting") {
            Self::Conflicting
        } else if main.starts_with("pathogenic") || main.starts_with("likely_pathogenic") {
            Self::Pathogenic
        } else if main.starts_with("uncertain") {
            Self::Uncertain
        } else if main.starts_with("benign") || main.starts_with("likely_benign") {
            Self::Benign
        } else {
            Self::Other
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Pathogenic => "pathogenic",
            Self::Uncertain => "uncertain",
            Self::Benign => "benign",
            Self::Conflicting => "conflicting",
            Self::Other => "other",
        }
    }

    pub fn rgb(self) -> [u8; 3] {
        self.rgb_in(false)
    }

    /// The class's colour on the light or dark theme. On dark, uncertain (UCSC's navy) and other
    /// (dark grey) were nearly invisible; they get a light violet-blue, still apart from
    /// conflicting's sky blue, and a light grey.
    pub fn rgb_in(self, dark: bool) -> [u8; 3] {
        match self {
            Self::Pathogenic => [210, 30, 30],
            Self::Uncertain if dark => [125, 125, 255],
            Self::Uncertain => [20, 40, 150],
            Self::Benign => [40, 150, 60],
            Self::Conflicting => [100, 170, 230],
            Self::Other if dark => [150, 150, 150],
            Self::Other => [100, 100, 100],
        }
    }
}

/// A site's colour for INFO field `key`'s value (grey when absent), on the light theme.
pub fn info_rgb(key: &str, value: Option<&str>) -> [u8; 3] {
    info_rgb_in(key, value, false)
}

/// A site's colour for INFO field `key`'s value on the light or dark theme.
pub fn info_rgb_in(key: &str, value: Option<&str>, dark: bool) -> [u8; 3] {
    match value {
        None | Some("" | ".") => [190, 190, 190],
        Some(v) if key == "CLNSIG" => Significance::of(v).rgb_in(dark),
        Some(v) => {
            let (r, g, b) = ggv_layout::classify::tag_rgb(v);
            if dark { lift_for_dark([r, g, b]) } else { [r, g, b] }
        }
    }
}

enum Request {
    Load { id: u64, region: Region, color: SiteColor, dark: bool },
    /// Find the nearest variant strictly after (or before) `pos` on the first contig, then on
    /// each following contig in order (in the search direction).
    Seek { id: u64, contigs: Vec<(String, u64)>, pos: u64, forward: bool },
}

enum Reply {
    Loaded { id: u64, region: Region, variants: Arc<Vec<Variant>>, packed: Packed, instances: Arc<Vec<Instance>>, color: SiteColor, sv: Arc<Vec<Option<(&'static str, Option<u64>)>>>, legend: Vec<(String, [u8; 3])>, dark: bool, elapsed: Duration },
    Found { id: u64, hit: Option<(String, u64, u64)> },
    Failed(String),
}

pub struct Loaded {
    pub region: Region,
    pub variants: Arc<Vec<Variant>>,
    /// Each variant's site row, and how many rows there are.
    pub rows: Arc<Vec<u16>>,
    pub site_rows: usize,
    /// Variants that found no free row and overlap others in the last.
    pub overflow: usize,
    pub instances: InstanceSet,
    /// What the instances were coloured by, and for which theme.
    color: SiteColor,
    dark: bool,
    pub elapsed: Duration,
    /// Each variant's structural type and size (`Variant::structural`), worked out once per
    /// load rather than every frame.
    pub sv: Arc<Vec<Option<(&'static str, Option<u64>)>>>,
    /// The colour legend for what's loaded, worked out once per load or recolouring.
    legend: Vec<(String, [u8; 3])>,
}

pub struct VariantTrack {
    pub id: u64,
    pub path: PathBuf,
    pub samples: Vec<String>,
    /// The header's INFO fields with their descriptions.
    pub info_fields: Arc<Vec<InfoField>>,
    pub color: SiteColor,
    /// Colours for the dark theme (set by the pane each frame; a change recolours).
    pub dark: bool,
    /// The header's contig lengths and declared assembly, for checking it matches the rest.
    pub header_contigs: Vec<ggv_model::Contig>,
    pub declared_assembly: Option<String>,
    tx: mpsc::Sender<Request>,
    rx: mpsc::Receiver<Reply>,
    pub current: Option<Loaded>,
    requested: Option<(u64, Region)>,
    seeking: Option<u64>,
    pub scroll_y: f32,
    pub error: Option<String>,
}

impl VariantTrack {
    pub fn open(id: u64, path: PathBuf, ctx: &egui::Context) -> anyhow::Result<Self> {
        let file = VariantFile::open(&path)?;
        let samples = file.samples().to_vec();
        let info_fields = Arc::new(file.info_fields());
        let color = SiteColor::default_for(&info_fields);
        let header_contigs = file.header_contigs();
        let declared_assembly = file.declared_assembly();
        let (tx, worker_rx) = mpsc::channel::<Request>();
        let (worker_tx, rx) = mpsc::channel();
        let ctx = ctx.clone();
        let n = samples.len();
        thread::Builder::new().name("ggv-variants".into()).spawn(move || {
            let mut file = file;
            while let Ok(mut req) = worker_rx.recv() {
                // Loads may be superseded; seeks never are (the user asked for each one).
                if matches!(req, Request::Load { .. }) {
                    while let Ok(newer) = worker_rx.try_recv() {
                        req = newer;
                        if matches!(req, Request::Seek { .. }) {
                            break;
                        }
                    }
                }
                let reply = match req {
                    Request::Load { id, region, color, dark } => {
                        let t = Instant::now();
                        match file.query(&region, false) {
                            Ok(v) => {
                                let packed = pack(&v);
                                let instances = Arc::new(build(&v, &packed, n, region.start, &color, dark));
                                let sv = Arc::new(v.iter().map(|x| x.structural(50)).collect());
                                let legend = legend_of(&color, &v, dark);
                                Reply::Loaded { id, region, variants: Arc::new(v), packed, instances, color, sv, legend, dark, elapsed: t.elapsed() }
                            }
                            Err(e) => Reply::Failed(e.to_string()),
                        }
                    }
                    Request::Seek { id, contigs, pos, forward } => match seek_across(&mut file, &contigs, pos, forward) {
                        Ok(hit) => Reply::Found { id, hit },
                        Err(e) => Reply::Failed(e.to_string()),
                    },
                };
                if worker_tx.send(reply).is_err() {
                    break;
                }
                ctx.request_repaint();
            }
        })?;
        Ok(Self { id, path, samples, info_fields, color, dark: false, header_contigs, declared_assembly, tx, rx, current: None, requested: None, seeking: None, scroll_y: 0.0, error: None })
    }

    pub fn file_name(&self) -> String {
        self.path.file_name().unwrap_or_default().to_string_lossy().into_owned()
    }

    /// A next/previous-site search is under way.
    pub fn is_seeking(&self) -> bool {
        self.seeking.is_some()
    }

    pub fn is_loading(&self) -> bool {
        self.requested.is_some() || self.seeking.is_some()
    }

    /// Site rows in what's loaded (one before anything is).
    pub fn site_rows(&self) -> usize {
        self.current.as_ref().map_or(1, |c| c.site_rows)
    }

    /// Where the genotype rows start, below the site rows.
    pub fn gt_top(&self) -> f32 {
        gt_top(self.site_rows())
    }

    pub fn section_height(&self) -> f32 {
        self.gt_top() + self.samples.len().min(VISIBLE_SAMPLES) as f32 * GT_PITCH + 3.0
    }

    pub fn max_scroll(&self, section_h: f32) -> f32 {
        (self.samples.len() as f32 * GT_PITCH - (section_h - self.gt_top() - 3.0)).max(0.0)
    }

    /// Requests the span around the view if what's loaded (or pending) doesn't cover it.
    pub fn ensure_loaded(&mut self, contig: &str, contig_len: u64, vs: f64, ve: f64, next_id: &mut u64) {
        if ve - vs > VARIANT_WINDOW as f64 {
            return;
        }
        // A span far wider than the view (loaded while zoomed out) is reloaded: instances are
        // f32 bases from the span's start, and millions of bases away a base falls between
        // representable values, so sites drift and jitter against the reads as you zoom.
        let max_span = ((ve - vs) * 8.0).max(MIN_SPAN * 2.0);
        let covers = |r: &Region| {
            r.contig == contig && r.start as f64 <= vs.max(0.0) && r.end as f64 >= ve.min(contig_len as f64) && ((r.end - r.start) as f64) <= max_span
        };
        if self.current.as_ref().is_some_and(|c| covers(&c.region)) || self.requested.as_ref().is_some_and(|(_, r)| covers(r)) {
            return;
        }
        let half = ((ve - vs) * 1.5).max(MIN_SPAN / 2.0);
        let center = (vs + ve) / 2.0;
        let start = (center - half).max(0.0) as u64;
        let end = ((center + half) as u64).min(contig_len).max(start + 1);
        let region = Region { contig: contig.to_owned(), start, end };
        let id = *next_id;
        *next_id += 1;
        let _ = self.tx.send(Request::Load { id, region: region.clone(), color: self.color.clone(), dark: self.dark });
        self.requested = Some((id, region));
    }

    /// `contigs`: the current contig first, then the ones to continue into, in order.
    pub fn seek(&mut self, contigs: Vec<(String, u64)>, pos: u64, forward: bool, next_id: &mut u64) {
        let id = *next_id;
        *next_id += 1;
        let _ = self.tx.send(Request::Seek { id, contigs, pos, forward });
        self.seeking = Some(id);
    }

    /// Takes finished work. Returns a found variant `(contig, start, end)` from a seek, with
    /// `Some(None)` meaning the seek found nothing.
    pub fn receive(&mut self, generation: &mut u64) -> Option<Option<(String, u64, u64)>> {
        let mut found = None;
        for reply in self.rx.try_iter().collect::<Vec<_>>() {
            match reply {
                Reply::Loaded { id, region, variants, packed, instances, color, sv, legend, dark, elapsed } => {
                    if self.requested.as_ref().is_some_and(|(r, _)| *r == id) {
                        self.requested = None;
                    }
                    *generation += 1;
                    self.current = Some(Loaded {
                        region,
                        variants,
                        rows: Arc::new(packed.rows),
                        site_rows: packed.count,
                        overflow: packed.overflow,
                        instances: InstanceSet { generation: *generation, instances },
                        color,
                        dark,
                        elapsed,
                        sv,
                        legend,
                    });
                    self.error = None;
                }
                Reply::Found { id, hit } => {
                    if self.seeking == Some(id) {
                        self.seeking = None;
                        found = Some(hit);
                    }
                }
                Reply::Failed(e) => {
                    self.requested = None;
                    self.seeking = None;
                    self.error = Some(e);
                }
            }
        }
        if self.current.as_ref().is_some_and(|c| c.color != self.color || c.dark != self.dark) {
            self.recolor(generation);
        }
        found
    }

    /// Rebuilds what's loaded for [`VariantTrack::color`] and the theme (a load in flight comes
    /// back in the old colouring and is rebuilt here too).
    fn recolor(&mut self, generation: &mut u64) {
        let n = self.samples.len();
        if let Some(c) = &mut self.current {
            *generation += 1;
            let packed = Packed { rows: c.rows.to_vec(), count: c.site_rows, overflow: c.overflow };
            c.instances = InstanceSet { generation: *generation, instances: Arc::new(build(&c.variants, &packed, n, c.region.start, &self.color, self.dark)) };
            c.color = self.color.clone();
            c.dark = self.dark;
            c.legend = legend_of(&c.color, &c.variants, c.dark);
        }
    }

    /// The legend for the current colouring: each class or value with its colour, the ones seen
    /// in what's loaded (most common first, at most `max`). Empty for allele frequency.
    pub fn legend(&self, max: usize) -> Vec<(String, [u8; 3])> {
        self.current.as_ref().map(|c| c.legend.iter().take(max).cloned().collect()).unwrap_or_default()
    }
}

/// The legend for `variants` coloured by `color` (see [`VariantTrack::legend`]): read once per
/// load, as it reads every record's INFO.
fn legend_of(color: &SiteColor, variants: &[Variant], dark: bool) -> Vec<(String, [u8; 3])> {
    let SiteColor::Info(key) = color else { return Vec::new() };
    if key == "CLNSIG" {
        let mut seen = [false; Significance::ALL.len()];
        for s in variants.iter().filter_map(|v| v.info_value(key)).map(Significance::of) {
            if let Some(i) = Significance::ALL.iter().position(|x| *x == s) {
                seen[i] = true;
            }
        }
        return Significance::ALL.iter().zip(seen).filter(|(_, s)| *s).map(|(s, _)| (s.label().to_owned(), s.rgb_in(dark))).collect();
    }
    let mut counts: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
    for v in variants.iter().filter_map(|v| v.info_value(key)) {
        *counts.entry(v).or_default() += 1;
    }
    let mut counts: Vec<(&str, usize)> = counts.into_iter().collect();
    counts.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(b.0)));
    counts.iter().take(12).map(|(v, _)| (v.replace('_', " "), info_rgb_in(key, Some(v), dark))).collect()
}

/// Seeks on the first contig from `pos`, then from the start (forward) or end (backward) of
/// each following contig.
fn seek_across(file: &mut VariantFile, contigs: &[(String, u64)], pos: u64, forward: bool) -> ggv_io::Result<Option<(String, u64, u64)>> {
    for (k, (contig, len)) in contigs.iter().enumerate() {
        let from = match (k, forward) {
            (0, _) => Some(pos),
            (_, true) => None,
            (_, false) => Some(*len),
        };
        if let Some((s, e)) = seek(file, contig, *len, from, forward)? {
            return Ok(Some((contig.clone(), s, e)));
        }
    }
    Ok(None)
}

/// The nearest variant after/before `pos`, searching outwards in growing windows so a dense
/// file answers from one small query and a sparse one still finds far-off sites.
/// `pos: None` (forward only) means from the contig's start, including a site at position 0.
fn seek(file: &mut VariantFile, contig: &str, contig_len: u64, pos: Option<u64>, forward: bool) -> ggv_io::Result<Option<(u64, u64)>> {
    let after = |v: &ggv_io::Variant| pos.is_none_or(|p| v.start > p);
    let first = pos.map_or(0, |p| p + 1);
    let pos = pos.unwrap_or(0);
    let mut window = 100_000u64;
    loop {
        let (start, end) = if forward {
            (first, (first + window).min(contig_len))
        } else {
            (pos.saturating_sub(window), pos)
        };
        if start >= end {
            return Ok(None);
        }
        let vs = file.query(&Region { contig: contig.to_owned(), start, end }, false)?;
        let hit = if forward {
            vs.iter().filter(|v| after(v)).min_by_key(|v| v.start)
        } else {
            vs.iter().filter(|v| v.start < pos).max_by_key(|v| v.start)
        };
        if let Some(v) = hit {
            return Ok(Some((v.start, v.end)));
        }
        if (forward && end >= contig_len) || (!forward && start == 0) {
            return Ok(None);
        }
        window = window.saturating_mul(4);
    }
}

pub fn gt_top(site_rows: usize) -> f32 {
    SITE_TOP + (site_rows.max(1) - 1) as f32 * SITE_PITCH + SITE_H + 3.0
}

/// The top of site row `row`.
pub fn site_y(row: u16) -> f32 {
    SITE_TOP + f32::from(row) * SITE_PITCH
}

pub struct Packed {
    pub rows: Vec<u16>,
    pub count: usize,
    pub overflow: usize,
}

/// Site rows: the shortest records first, each in the first row where it overlaps nothing, so
/// SNVs keep the top row and the indels and SVs over them stack below (insertions count as one
/// base at their boundary).
pub fn pack(variants: &[Variant]) -> Packed {
    use std::collections::BTreeMap;
    let span = |v: &Variant| {
        let (s, e) = v.shown();
        (s, e.max(s + 1))
    };
    let mut order: Vec<usize> = (0..variants.len()).collect();
    order.sort_by_key(|&i| {
        let (s, e) = span(&variants[i]);
        (e - s, s)
    });
    // Per row, placed spans by start (they never overlap within a row).
    let mut placed: Vec<BTreeMap<u64, u64>> = Vec::new();
    let free = |row: &BTreeMap<u64, u64>, s: u64, e: u64| {
        row.range(..e).next_back().is_none_or(|(_, &end)| end <= s)
    };
    let mut rows = vec![0u16; variants.len()];
    let mut overflow = 0;
    for i in order {
        let (s, e) = span(&variants[i]);
        let r = match placed.iter().position(|row| free(row, s, e)) {
            Some(r) => r,
            None if placed.len() < MAX_SITE_ROWS => {
                placed.push(BTreeMap::new());
                placed.len() - 1
            }
            None => {
                overflow += 1;
                placed.len() - 1
            }
        };
        placed[r].insert(s, e);
        rows[i] = r as u16;
    }
    Packed { rows, count: placed.len().max(1), overflow }
}

fn build(variants: &[Variant], packed: &Packed, samples: usize, origin: u64, color: &SiteColor, dark: bool) -> Vec<Instance> {
    let gt_top = gt_top(packed.count);
    // Rounded, anti-aliased, with a pixel between neighbours once they're wide enough.
    let shape = FLAG_MIN3 | FLAG_ROUND | FLAG_GAP;
    let mut out = Vec::with_capacity(variants.len() * (2 + samples));
    let x = |p: u64| (p as i64 - origin as i64) as f32;
    let mut push = |s: u64, e: u64, y0: f32, y1: f32, color: u32, flags: u32| {
        out.push(Instance { x0: x(s), x1: x(e), y0, y1, color, flags, lod_min: 0.0, lod_max: f32::INFINITY });
    };
    for (v, &row) in variants.iter().zip(&packed.rows) {
        let (start, end) = v.shown();
        let top = site_y(row);
        // Site: red for the alt fraction, blue for the rest (IGV), or the INFO value's colour;
        // faded when filtered.
        let alpha = if v.is_pass() { 255 } else { 90 };
        match color {
            SiteColor::AlleleFrequency => {
                let af = v.af.unwrap_or(0.0).clamp(0.0, 1.0);
                let split = top + SITE_H * (1.0 - af);
                let ([ar, ag, ab], [rr, rg, rb]) = af_rgb(dark);
                push(start, end, top, top + SITE_H, rgba(rr, rg, rb, alpha), shape);
                if af > 0.0 {
                    push(start, end, split, top + SITE_H, rgba(ar, ag, ab, alpha), shape);
                }
            }
            SiteColor::Info(key) => {
                let [r, g, b] = info_rgb_in(key, v.info_value(key), dark);
                push(start, end, top, top + SITE_H, rgba(r, g, b, alpha), shape);
            }
        }
        for (s, call) in v.samples.iter().enumerate() {
            if let Some([r, g, b]) = genotype_rgb(call.class, dark) {
                let y = gt_top + s as f32 * GT_PITCH;
                push(start, end, y, y + GT_H, rgba(r, g, b, 255), shape | FLAG_SCROLL);
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn site(start: u64, reference: &str, alt: &str) -> Variant {
        Variant {
            start,
            end: start + reference.len() as u64,
            id: ".".into(),
            reference: reference.into(),
            alts: vec![alt.into()],
            qual: None,
            filters: vec![],
            info: String::new(),
            af: None,
            samples: vec![],
            ..Default::default()
        }
    }

    #[test]
    fn overlapping_sites_get_their_own_rows() {
        // A 10 bp deletion (bases 101-110), an SNV inside it, an SNV beside it, a second ALT at
        // the first SNV's base, and an insertion at the deletion's edge.
        let vs = [site(100, "AAAAAAAAAAA", "A"), site(105, "C", "T"), site(120, "G", "A"), site(105, "C", "G"), site(110, "A", "AT")];
        let p = pack(&vs);
        // SNVs keep the top row; the second ALT at 105 stacks; the deletion goes under both.
        assert_eq!(p.rows, vec![2, 0, 0, 1, 0]);
        assert_eq!((p.count, p.overflow), (3, 0));
        // Adjacent, not overlapping: one row.
        assert_eq!(pack(&[site(10, "A", "G"), site(11, "C", "T"), site(12, "GA", "G")]).rows, vec![0, 0, 0]);
        // More stacked than the cap: the rest share the last row.
        let same: Vec<Variant> = (0..MAX_SITE_ROWS + 3).map(|_| site(50, "A", "G")).collect();
        let p = pack(&same);
        assert_eq!((p.count, p.overflow), (MAX_SITE_ROWS, 3));
        assert_eq!(pack(&[]).count, 1);
    }

    #[test]
    fn clinvar_significance_classes() {
        for (v, c) in [
            ("Pathogenic", Significance::Pathogenic),
            ("Likely_pathogenic", Significance::Pathogenic),
            ("Pathogenic/Likely_pathogenic|risk_factor", Significance::Pathogenic),
            ("Uncertain_significance", Significance::Uncertain),
            ("Benign/Likely_benign", Significance::Benign),
            ("Likely_benign", Significance::Benign),
            ("Conflicting_classifications_of_pathogenicity", Significance::Conflicting),
            ("Conflicting_interpretations_of_pathogenicity|other", Significance::Conflicting),
            ("not_provided", Significance::Other),
            ("drug_response", Significance::Other),
        ] {
            assert_eq!(Significance::of(v), c, "{v}");
        }
        assert_eq!(info_rgb("CLNSIG", None), [190, 190, 190]);
        assert_eq!(info_rgb("SVTYPE", Some("DEL")), info_rgb("SVTYPE", Some("DEL")), "a value always gets one colour");
        assert_eq!(SiteColor::parse("AF"), SiteColor::AlleleFrequency);
        assert_eq!(SiteColor::parse("CLNSIG").name(), "CLNSIG");
    }

    /// Dark theme: every site and genotype colour stands at least 3:1 against the background
    /// (ClinVar's uncertain was 1.5:1, het 2.1:1), the classes stay apart, and the light theme's
    /// colours are unchanged.
    #[test]
    fn dark_theme_colours_are_readable_and_distinct() {
        assert!(contrast([20, 40, 150], DARK_BG) < 2.0, "the old uncertain colour, for reference");
        let sig: Vec<[u8; 3]> = Significance::ALL.iter().map(|s| s.rgb_in(true)).collect();
        for (s, c) in Significance::ALL.iter().zip(&sig) {
            assert!(contrast(*c, DARK_BG) >= 3.0, "{s:?} {c:?}: {:.2}", contrast(*c, DARK_BG));
        }
        for (i, a) in sig.iter().enumerate() {
            for b in &sig[i + 1..] {
                let d: i32 = a.iter().zip(b).map(|(x, y)| (i32::from(*x) - i32::from(*y)).abs()).sum();
                assert!(d >= 60, "{a:?} and {b:?} look alike");
            }
        }
        let calls = [GenotypeClass::HomRef, GenotypeClass::Het, GenotypeClass::HomAlt, GenotypeClass::NoCall];
        for g in calls {
            let c = genotype_rgb(g, true).unwrap();
            assert!(contrast(c, DARK_BG) >= 3.0, "{g:?} {c:?}");
        }
        let (alt, rest) = af_rgb(true);
        assert!(contrast(alt, DARK_BG) >= 3.0 && contrast(rest, DARK_BG) >= 3.0);
        for v in ["DEL", "INS", "DUP", "INV", "BND", "x", "PASS"] {
            assert!(contrast(info_rgb_in("SVTYPE", Some(v), true), DARK_BG) >= 3.0, "{v}");
        }
        // Light: as before.
        assert_eq!(Significance::Uncertain.rgb(), [20, 40, 150]);
        assert_eq!(Significance::Uncertain.rgb_in(false), [20, 40, 150]);
        assert_eq!(genotype_rgb(GenotypeClass::Het, false), Some([34, 12, 253]));
        assert_eq!(af_rgb(false), ([220, 0, 0], [0, 0, 220]));
        assert_eq!(info_rgb_in("SVTYPE", Some("DEL"), false), info_rgb("SVTYPE", Some("DEL")));
    }

    /// A track switched to the dark theme rebuilds what's loaded in its colours, legend too.
    #[test]
    fn switching_theme_recolours_the_loaded_sites() {
        let data = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio");
        let path = data.join("HG002.revio.deepvariant.phased.vcf.gz");
        if !path.exists() {
            return;
        }
        let ctx = egui::Context::default();
        let mut t = VariantTrack::open(1, path, &ctx).unwrap();
        let mut next = 1;
        let mut generation = 0;
        t.ensure_loaded("chr20", 64_444_167, 31_165_000.0, 31_184_000.0, &mut next);
        let start = Instant::now();
        while t.current.is_none() && start.elapsed() < Duration::from_secs(30) {
            t.receive(&mut generation);
            std::thread::sleep(Duration::from_millis(10));
        }
        let colours = |t: &VariantTrack| t.current.as_ref().unwrap().instances.instances.iter().map(|i| i.color).collect::<std::collections::HashSet<u32>>();
        let light = colours(&t);
        assert!(light.contains(&rgba(34, 12, 253, 255)), "het calls in the light colour");
        t.dark = true;
        t.receive(&mut generation);
        let dark = colours(&t);
        assert!(!dark.contains(&rgba(34, 12, 253, 255)) && dark.contains(&rgba(100, 130, 255, 255)), "het calls recoloured");
        t.dark = false;
        t.receive(&mut generation);
        assert_eq!(colours(&t), light);
    }

    #[test]
    fn info_values_and_flags() {
        let v = Variant {
            start: 0,
            end: 1,
            id: ".".into(),
            reference: "A".into(),
            alts: vec!["G".into()],
            qual: None,
            filters: vec![],
            info: "ALLELEID=1;CLNSIG=Benign;DB;GENEINFO=BRCA1:672".into(),
            af: None,
            samples: vec![],
            ..Default::default()
        };
        assert_eq!(v.info_value("CLNSIG"), Some("Benign"));
        assert_eq!(v.info_value("DB"), Some(""));
        assert_eq!(v.info_value("CLN"), None);
        assert_eq!(v.info_pairs().collect::<Vec<_>>(), vec![("ALLELEID", "1"), ("CLNSIG", "Benign"), ("DB", ""), ("GENEINFO", "BRCA1:672")]);
    }

    #[test]
    fn seek_finds_neighbours() {
        // GGV_TEST_VCF=giab.HG002.vcf.gz (the HG002 chr20 fixture) checks against known sites.
        let Ok(path) = std::env::var("GGV_TEST_VCF") else { return };
        let mut f = VariantFile::open(path).unwrap();
        let len = 64_444_167;
        // 1-based 1,895,892 and 1,896,504 are consecutive GIAB sites.
        assert_eq!(seek(&mut f, "chr20", len, Some(1_895_891), true).unwrap(), Some((1_896_503, 1_896_504)));
        assert_eq!(seek(&mut f, "chr20", len, Some(1_896_503), false).unwrap(), Some((1_895_891, 1_895_892)));
        // Contig naming without `chr` works too, and nothing lies past the end.
        assert_eq!(seek(&mut f, "20", len, Some(1_895_891), true).unwrap(), Some((1_896_503, 1_896_504)));
        assert_eq!(seek(&mut f, "chr20", len, Some(len - 1), true).unwrap(), None);
        // From the contig start, forward finds the first site; backward from it finds nothing.
        let first = seek(&mut f, "chr20", len, None, true).unwrap().unwrap();
        assert_eq!(seek(&mut f, "chr20", len, Some(first.0), false).unwrap(), None);

        // Across contigs: past the end of chr19 (absent from the file) and into chr20; and
        // backward from chr21's end into chr20's last site.
        let contigs = |names: &[&str]| names.iter().map(|n| ((*n).to_owned(), len)).collect::<Vec<_>>();
        let hit = seek_across(&mut f, &contigs(&["chr19", "chr20"]), 1000, true).unwrap().unwrap();
        assert_eq!((hit.0.as_str(), hit.1), ("chr20", first.0));
        let last = seek(&mut f, "chr20", len, Some(len), false).unwrap().unwrap();
        let hit = seek_across(&mut f, &contigs(&["chr21", "chr20"]), 5, false).unwrap().unwrap();
        assert_eq!((hit.0.as_str(), hit.1), ("chr20", last.0));
        assert_eq!(seek_across(&mut f, &contigs(&["chr21", "chr22"]), 0, true).unwrap(), None);
    }
}
