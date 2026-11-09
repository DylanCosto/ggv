//! A gene annotation track: parsed once on a background thread, then drawn like the
//! alignment tracks: exons, UTRs and introns as GPU instances in genomic coordinates; names and
//! strand chevrons per frame for what's visible.

use std::{
    path::PathBuf,
    sync::{Arc, mpsc},
    thread,
    time::{Duration, Instant},
};

use eframe::egui;
use ggv_annot::{Annotation, Transcript};

use crate::{
    gpu::{FLAG_SCROLL, Instance, InstanceSet},
    tracks::{HEADER_H, rgba},
};

/// Rows start here, in points from the section top.
pub const ROWS_TOP: f32 = HEADER_H + 3.0;
/// One row: feature (12 pt) and its name underneath.
pub const PITCH: f32 = 25.0;
pub const FEATURE_H: f32 = 12.0;
/// Rows shown before the section scrolls.
const VISIBLE_ROWS: usize = 4;
/// Room left for a name between neighbours in a row, in points.
const LABEL_ROOM: f64 = 90.0;
/// The most room left for names between neighbours, in bases.
const MAX_GAP: f64 = 100_000.0;
/// Load at least this much around the view.
const MIN_SPAN: f64 = 200_000.0;

pub const CODING: [u8; 3] = [0, 0, 150];
pub const NONCODING: [u8; 3] = [100, 100, 170];

/// Feature colour: navy on light, light blue on dark.
pub fn feature_rgb(coding: bool, dark: bool) -> [u8; 3] {
    match (coding, dark) {
        (true, false) => CODING,
        (false, false) => NONCODING,
        (true, true) => [120, 160, 255],
        (false, true) => [150, 150, 200],
    }
}

enum State {
    Loading(mpsc::Receiver<Result<Annotation, String>>, Instant),
    Ready(Arc<Annotation>, Duration),
    Failed(String),
}

/// Transcripts placed into rows for a span around the view.
pub struct Built {
    pub contig: String,
    pub start: u64,
    pub end: u64,
    all: bool,
    dark: bool,
    ppb: f64,
    /// Transcripts in the span, in start order.
    pub items: Vec<Transcript>,
    /// Indices into `items`, per row, in start order (non-overlapping within a row).
    pub rows: Vec<Vec<usize>>,
    pub instances: InstanceSet,
}

pub struct GeneTrack {
    pub id: u64,
    pub path: PathBuf,
    state: State,
    pub built: Option<Built>,
    pub scroll_y: f32,
}

impl Built {
    pub fn is_all(&self) -> bool {
        self.all
    }
}

impl GeneTrack {
    pub fn open(id: u64, path: PathBuf, ctx: &egui::Context) -> Self {
        let (tx, rx) = mpsc::channel();
        let p = path.clone();
        let ctx = ctx.clone();
        let _ = thread::Builder::new().name("ggv-genes".into()).spawn(move || {
            let _ = tx.send(Annotation::open(&p).map_err(|e| e.to_string()));
            ctx.request_repaint();
        });
        Self { id, path, state: State::Loading(rx, Instant::now()), built: None, scroll_y: 0.0 }
    }

    pub fn file_name(&self) -> String {
        self.path.file_name().unwrap_or_default().to_string_lossy().into_owned()
    }

    pub fn annotation(&self) -> Option<&Arc<Annotation>> {
        match &self.state {
            State::Ready(a, _) => Some(a),
            _ => None,
        }
    }

    pub fn is_loading(&self) -> bool {
        matches!(self.state, State::Loading(..))
    }

    /// A one-line status for the header.
    pub fn status(&self) -> String {
        match &self.state {
            State::Loading(_, t) => format!("loading… {:.1} s", t.elapsed().as_secs_f64()),
            State::Ready(a, d) => format!("{} genes · {} transcripts · parsed in {:.2} s", a.gene_count(), a.len(), d.as_secs_f64()),
            State::Failed(e) => format!("error: {e}"),
        }
    }

    /// Picks up the parsed annotation; returns true when it just became ready.
    pub fn poll(&mut self) -> bool {
        let State::Loading(rx, t) = &self.state else { return false };
        match rx.try_recv() {
            Ok(Ok(a)) => {
                self.state = State::Ready(Arc::new(a), t.elapsed());
                true
            }
            Ok(Err(e)) => {
                self.state = State::Failed(e);
                false
            }
            Err(mpsc::TryRecvError::Empty) => false,
            Err(mpsc::TryRecvError::Disconnected) => {
                self.state = State::Failed("loader stopped".into());
                false
            }
        }
    }

    /// Rebuilds rows and instances when the view leaves the built span, the zoom changes enough
    /// to need different spacing for names, or the transcript mode changes. Cheap: a query and
    /// a few thousand instances, on the UI thread.
    #[allow(clippy::too_many_arguments)]
    pub fn update(&mut self, contig: &str, contig_len: u64, vs: f64, ve: f64, ppb: f64, all: bool, dark: bool, generation: &mut u64) {
        let Some(ann) = self.annotation().cloned() else { return };
        let fresh = self.built.as_ref().is_some_and(|b| {
            b.contig == contig
                && b.all == all
                && b.dark == dark
                && (b.start as f64) <= vs.max(0.0)
                && (b.end as f64) >= ve.min(contig_len as f64)
                && (0.5..=2.0).contains(&(ppb / b.ppb))
        });
        if fresh {
            return;
        }

        let vis = ve - vs;
        let half = (vis * 1.5).max(MIN_SPAN / 2.0);
        let center = (vs + ve) / 2.0;
        let start = (center - half).max(0.0) as u64;
        let end = ((center + half) as u64).min(contig_len).max(start + 1);

        let items: Vec<Transcript> =
            ann.query(contig, start, end).filter(|t| all || t.representative).cloned().collect();
        // Room for names, capped so that zoomed far out genes pack densely (names that would
        // collide are skipped when drawing) instead of spreading over hundreds of rows.
        let gap = (LABEL_ROOM / ppb).clamp(2.0, MAX_GAP) as u64;
        let rows = pack(&items, gap);
        *generation += 1;
        let instances = InstanceSet { generation: *generation, instances: Arc::new(build(&items, &rows, start, dark)) };
        self.built = Some(Built { contig: contig.to_owned(), start, end, all, dark, ppb, items, rows, instances });
    }

    pub fn section_height(&self) -> f32 {
        let rows = self.built.as_ref().map_or(1, |b| b.rows.len()).clamp(1, VISIBLE_ROWS);
        ROWS_TOP + rows as f32 * PITCH + 2.0
    }

    pub fn max_scroll(&self, section_h: f32) -> f32 {
        let rows = self.built.as_ref().map_or(0, |b| b.rows.len());
        (rows as f32 * PITCH - (section_h - ROWS_TOP)).max(0.0)
    }
}

/// First-fit rows, leaving `gap` bases after each transcript for its name.
fn pack(items: &[Transcript], gap: u64) -> Vec<Vec<usize>> {
    let mut rows: Vec<Vec<usize>> = Vec::new();
    let mut ends: Vec<u64> = Vec::new();
    for (i, t) in items.iter().enumerate() {
        match ends.iter().position(|&e| e + gap <= t.start) {
            Some(r) => {
                ends[r] = t.end;
                rows[r].push(i);
            }
            None => {
                ends.push(t.end);
                rows.push(vec![i]);
            }
        }
    }
    rows
}

fn build(items: &[Transcript], rows: &[Vec<usize>], origin: u64, dark: bool) -> Vec<Instance> {
    let mut out = Vec::new();
    let x = |p: u64| (p as i64 - origin as i64) as f32;
    let mut push = |s: u64, e: u64, y0: f32, y1: f32, c: [u8; 3]| {
        out.push(Instance {
            x0: x(s),
            x1: x(e),
            y0,
            y1,
            color: rgba(c[0], c[1], c[2], 255),
            flags: FLAG_SCROLL,
            lod_min: 0.0,
            lod_max: f32::INFINITY,
        });
    };

    for (r, row) in rows.iter().enumerate() {
        let y = ROWS_TOP + r as f32 * PITCH;
        let mid = y + FEATURE_H / 2.0;
        for &i in row {
            let t = &items[i];
            let color = feature_rgb(t.cds.is_some(), dark);
            push(t.start, t.end, mid - 0.5, mid + 0.5, color);
            for &(s, e) in &t.exons {
                match t.cds {
                    // Coding exon: thick inside the CDS, thin (UTR) outside it.
                    Some((cs, ce)) => {
                        let (ts, te) = (s.max(cs), e.min(ce));
                        if ts < te {
                            push(ts, te, y, y + FEATURE_H, color);
                        }
                        if s < ts.min(e) {
                            push(s, ts.min(e), y + 3.0, y + FEATURE_H - 3.0, color);
                        }
                        if te.max(s) < e {
                            push(te.max(s), e, y + 3.0, y + FEATURE_H - 3.0, color);
                        }
                    }
                    None => push(s, e, y + 2.0, y + FEATURE_H - 2.0, color),
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tx(start: u64, end: u64) -> Transcript {
        Transcript {
            contig: 0,
            start,
            end,
            strand: ggv_annot::Strand::Forward,
            id: "t".into(),
            gene_id: "g".into(),
            gene_name: "G".into(),
            biotype: "".into(),
            exons: vec![(start, end)],
            cds: None,
            tags: 0,
            representative: true,
            translation: Default::default(),
            exons_from_span: false,
        }
    }

    #[test]
    fn packing_leaves_room_for_names() {
        let items = vec![tx(0, 100), tx(150, 200), tx(250, 300)];
        // With a 100 bp gap, row 0 (ending at 100) next takes a start ≥ 200: 150 goes to a new
        // row, 250 fits back on row 0.
        assert_eq!(pack(&items, 100), vec![vec![0, 2], vec![1]]);
        assert_eq!(pack(&items, 10), vec![vec![0, 1, 2]]);
    }

    #[test]
    fn coding_exon_splits_into_utr_and_cds() {
        let mut t = tx(0, 1000);
        t.exons = vec![(0, 400), (600, 1000)];
        t.cds = Some((200, 800));
        let inst = build(&[t], &[vec![0]], 0, false);
        // intron line + (UTR 0-200, CDS 200-400) + (CDS 600-800, UTR 800-1000)
        let spans: Vec<(f32, f32, f32)> = inst.iter().map(|i| (i.x0, i.x1, i.y1 - i.y0)).collect();
        assert_eq!(
            spans,
            vec![(0.0, 1000.0, 1.0), (200.0, 400.0, 12.0), (0.0, 200.0, 6.0), (600.0, 800.0, 12.0), (800.0, 1000.0, 6.0)]
        );
    }
}
