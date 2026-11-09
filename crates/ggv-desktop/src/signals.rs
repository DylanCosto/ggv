//! Signal tracks (bigWig, bedGraph, SEG): a bar chart autoscaled to the view, or for SEG one
//! heatmap row per sample (blue loss, red gain).
//!
//! The file opens on a background thread (a bedGraph is read whole); values are then read for
//! the view plus a view's width each side, at about three bins per pixel, and read again only
//! when the view leaves that span or the zoom changes by more than half (a bigWig answers from
//! its zoom levels, so this is a few milliseconds even for a whole chromosome).

use std::{
    path::PathBuf,
    sync::mpsc::{self, Receiver},
};

use eframe::egui;
use ggv_io::signal::{Bin, Signal};

use crate::tracks::HEADER_H;

/// Bar chart height for a single-row track.
pub const CHART_H: f32 = 56.0;
/// Heatmap row height per SEG sample.
pub const SEG_ROW_H: f32 = 10.0;

struct Loaded {
    contig: String,
    start: u64,
    end: u64,
    bp_per_bin: f64,
    rows: Vec<Vec<Option<Bin>>>,
}

pub struct SignalTrack {
    pub path: PathBuf,
    rx: Option<Receiver<Result<Signal, String>>>,
    signal: Option<Signal>,
    pub rows: Vec<String>,
    pub error: Option<String>,
    loaded: Option<Loaded>,
    pub is_seg: bool,
}

impl SignalTrack {
    pub fn open(path: PathBuf, ctx: &egui::Context) -> Self {
        let (tx, rx) = mpsc::channel();
        let (p, ctx) = (path.clone(), ctx.clone());
        let _ = std::thread::Builder::new().name("ggv-signal".into()).spawn(move || {
            let _ = tx.send(Signal::open(&p).map_err(|e| e.to_string()));
            ctx.request_repaint();
        });
        let is_seg = path.to_string_lossy().to_ascii_lowercase().trim_end_matches(".gz").ends_with(".seg");
        Self { path, rx: Some(rx), signal: None, rows: vec![String::new()], error: None, loaded: None, is_seg }
    }

    pub fn file_name(&self) -> String {
        crate::files::display_name(&self.path)
    }

    pub fn is_loading(&self) -> bool {
        self.rx.is_some()
    }

    /// The height its section needs.
    pub fn section_height(&self) -> f32 {
        if self.is_seg {
            HEADER_H + self.rows.len().max(1) as f32 * SEG_ROW_H + 4.0
        } else {
            HEADER_H + CHART_H + 4.0
        }
    }

    /// Takes the file once it has opened; true when it just did.
    pub fn poll(&mut self) -> bool {
        if let Some(rx) = &self.rx
            && let Ok(r) = rx.try_recv()
        {
            self.rx = None;
            match r {
                Ok(s) => {
                    self.rows = s.rows();
                    self.signal = Some(s);
                    return true;
                }
                Err(e) => self.error = Some(e),
            }
        }
        false
    }

    /// The contigs the file has values for (to navigate a window with only signal tracks).
    pub fn contigs(&self) -> Vec<ggv_model::Contig> {
        self.signal.as_ref().map(Signal::contigs).unwrap_or_default()
    }

    /// Takes the opened file, and reads values for the view when needed.
    pub fn update(&mut self, contig: &str, vs: f64, ve: f64, width_px: f32) {
        if let Some(rx) = &self.rx
            && let Ok(r) = rx.try_recv()
        {
            self.rx = None;
            match r {
                Ok(s) => {
                    self.rows = s.rows();
                    self.signal = Some(s);
                }
                Err(e) => self.error = Some(e),
            }
        }
        let Some(sig) = &mut self.signal else { return };
        let (vs, ve) = (vs.max(0.0), ve.max(vs + 1.0));
        let want_bp = (ve - vs) / f64::from(width_px.max(1.0)) / 3.0;
        if let Some(l) = &self.loaded
            && l.contig == contig
            && l.start as f64 <= vs
            && l.end as f64 >= ve
            && (0.5..=2.0).contains(&(l.bp_per_bin / want_bp.max(1e-9)))
        {
            return;
        }
        let w = ve - vs;
        let (start, end) = ((vs - w).max(0.0) as u64, (ve + w) as u64);
        let n = ((end - start) as f64 / want_bp.max(1.0)).clamp(1.0, 20_000.0) as usize;
        let mut rows = Vec::with_capacity(self.rows.len());
        for k in 0..self.rows.len() {
            match sig.bins(k, contig, start, end, n) {
                Ok(b) => rows.push(b),
                Err(e) => {
                    self.error = Some(e.to_string());
                    return;
                }
            }
        }
        self.error = None;
        self.loaded = Some(Loaded { contig: contig.to_owned(), start, end, bp_per_bin: (end - start) as f64 / n as f64, rows });
    }

    /// Per pixel column over `vs..ve` (`n` columns): the bins under it combined (mean of means,
    /// min of mins, max of maxes).
    pub fn columns(&self, row: usize, vs: f64, ve: f64, n: usize) -> Vec<Option<Bin>> {
        let Some(l) = &self.loaded else { return vec![None; n] };
        let Some(bins) = l.rows.get(row) else { return vec![None; n] };
        let per = l.bp_per_bin;
        (0..n)
            .map(|c| {
                let (cs, ce) = (vs + (ve - vs) * c as f64 / n as f64, vs + (ve - vs) * (c + 1) as f64 / n as f64);
                let (a, b) = (((cs - l.start as f64) / per).floor().max(0.0) as usize, (((ce - l.start as f64) / per).ceil().max(0.0) as usize).min(bins.len()));
                let present: Vec<&Bin> = bins.get(a..b.max(a)).unwrap_or(&[]).iter().flatten().collect();
                (!present.is_empty()).then(|| Bin {
                    mean: present.iter().map(|x| x.mean).sum::<f64>() / present.len() as f64,
                    min: present.iter().map(|x| x.min).fold(f64::INFINITY, f64::min),
                    max: present.iter().map(|x| x.max).fold(f64::NEG_INFINITY, f64::max),
                    covered: present.iter().map(|x| x.covered).sum(),
                })
            })
            .collect()
    }

    /// The bin at a position (for the tooltip).
    pub fn at(&self, row: usize, pos: f64) -> Option<Bin> {
        let l = self.loaded.as_ref()?;
        let i = ((pos - l.start as f64) / l.bp_per_bin).floor();
        if i < 0.0 {
            return None;
        }
        *l.rows.get(row)?.get(i as usize)?
    }
}

/// SEG heatmap colour: blue for loss, red for gain, white at 0 (log2 ratios; full at ±1.5).
pub fn seg_color(v: f64) -> egui::Color32 {
    let t = (v / 1.5).clamp(-1.0, 1.0);
    let k = (255.0 * (1.0 - t.abs())) as u8;
    if t >= 0.0 { egui::Color32::from_rgb(255, k, k) } else { egui::Color32::from_rgb(k, k, 255) }
}
