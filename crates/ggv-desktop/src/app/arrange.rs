//! Arranging tracks: drag a track's header to move it, drag the edge between two tracks to
//! give one's height to the other, like a splitter (double-click the edge to go back to the
//! automatic height). With more tracks than fit, the last track's bottom edge resizes it alone.
//! A grip in each header and a notch on each edge show where to take hold.
//!
//! Tracks are known by kind and file, so the arrangement survives tracks opening and closing
//! and is saved with the session. A new track goes next to the last track of its kind (or where
//! its kind goes by default: variants on top, alignments in the middle, genes at the bottom).
//! Room the tracks don't need goes to those with more to show, alignment tracks first (see
//! `arrange`).

use std::{collections::HashMap, path::PathBuf};

use eframe::egui::{self, Pos2, Rect, Stroke, Vec2};
use serde::{Deserialize, Serialize};

use super::{App, Hit, Layout, MIN_SECTION_H, TRACKS_TOP};
use crate::{genes, tracks::READS_TOP, variants};

/// Room under the coverage for a one-line note ("Zoom in to see reads").
const NOTE_H: f32 = 22.0;

/// The track-stack scrollbar's width.
const STACK_BAR_W: f32 = 10.0;

/// The track-stack scrollbar, when the tracks don't all fit.
#[derive(Clone, Copy, Debug)]
pub(super) struct StackBar {
    pub track: Rect,
    pub thumb: Rect,
    pub first: usize,
    pub last_first: usize,
    pub count: usize,
}
use crate::tracks::HEADER_H;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Kind {
    Signals,
    Variants,
    Alignments,
    Genes,
}

pub type Key = (Kind, PathBuf);

/// The order and sizes the person chose.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Arrangement {
    pub order: Vec<Key>,
    pub heights: HashMap<Key, f32>,
}

/// One track as saved in a session: kind, file, and its height if set.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Saved {
    pub kind: Kind,
    pub path: PathBuf,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub height: Option<f32>,
}

/// A drag in progress.
#[derive(Clone, Debug, PartialEq)]
pub enum Arranging {
    /// The edge under `key`: `key` gets taller as the track below (`below`, with its height
    /// when the drag started) gets shorter; with none below, `key` alone changes.
    Resize { key: Key, start_h: f32, start_y: f32, below: Option<(Key, f32)> },
    Move { key: Key, y: f32 },
}

/// Smallest heights a person can drag a track to.
const MIN_ALIGN_H: f32 = 48.0;
const MIN_OTHER_H: f32 = HEADER_H + 12.0;
/// How close to a track's bottom edge the pointer must be to resize it.
const EDGE: f32 = 4.0;

impl Arrangement {
    pub fn saved(&self, present: &[Key]) -> Vec<Saved> {
        present.iter().map(|k| Saved { kind: k.0, path: k.1.clone(), height: self.heights.get(k).copied() }).collect()
    }

    pub fn from_saved(saved: &[Saved]) -> Self {
        Self {
            order: saved.iter().map(|s| (s.kind, s.path.clone())).collect(),
            heights: saved.iter().filter_map(|s| Some(((s.kind, s.path.clone()), s.height?))).collect(),
        }
    }

    /// The tracks in display order: those the person placed, in their order, and the others
    /// next to the last of their kind (or in the default place for the kind).
    pub fn ordered(&self, default: &[Key]) -> Vec<Key> {
        let mut out: Vec<Key> = self.order.iter().filter(|k| default.contains(k)).cloned().collect();
        for k in default {
            if out.contains(k) {
                continue;
            }
            let at = match out.iter().rposition(|x| x.0 == k.0) {
                Some(i) => i + 1,
                None => match k.0 {
                    Kind::Variants => 0,
                    Kind::Alignments | Kind::Signals => out.iter().position(|x| x.0 == Kind::Genes).unwrap_or(out.len()),
                    Kind::Genes => out.len(),
                },
            };
            out.insert(at, k.clone());
        }
        out
    }
}

impl App {
    /// Every open track's key, in the default order.
    pub(super) fn track_keys(&self) -> Vec<Key> {
        let v = self.variants.iter().map(|v| (Kind::Variants, v.path.clone()));
        let a = self.tracks.iter().map(|t| (Kind::Alignments, t.loader.bam_path.clone()));
        let g = self.genes.iter().map(|g| (Kind::Genes, g.path.clone()));
        let s = self.signals.iter().map(|s| (Kind::Signals, s.path.clone()));
        v.chain(a).chain(s).chain(g).collect()
    }

    fn hit_of(&self, key: &Key) -> Option<Hit> {
        match key.0 {
            Kind::Variants => self.variants.iter().position(|v| v.path == key.1).map(Hit::Variants),
            Kind::Alignments => self.tracks.iter().position(|t| t.loader.bam_path == key.1).map(Hit::Track),
            Kind::Genes => self.genes.iter().position(|g| g.path == key.1).map(Hit::Genes),
            Kind::Signals => self.signals.iter().position(|s| s.path == key.1).map(Hit::Signal),
        }
    }

    /// A track's height when the person hasn't set one: what gene and variant tracks need for
    /// their rows; alignment tracks share the rest (`None`).
    fn natural_height(&self, hit: Hit) -> Option<f32> {
        match hit {
            Hit::Variants(i) => Some(self.variants[i].section_height()),
            Hit::Genes(i) => Some(self.genes[i].section_height()),
            Hit::Signal(i) => Some(self.signals[i].section_height()),
            Hit::Track(_) => None,
        }
    }

    /// The least a track can be squeezed to while another is made taller: a height the person
    /// set can shrink to its kind's minimum; others keep what they need (alignment tracks that
    /// share the space, the shared minimum).
    fn floor_of(&self, key: &Key, hit: Hit) -> f32 {
        if self.arrangement.heights.contains_key(key) {
            return min_height(key.0);
        }
        self.natural_height(hit).unwrap_or(MIN_SECTION_H)
    }

    /// The tallest `key` can be so every other track still fits in `avail` (so dragging an
    /// edge, or a session saved on a bigger screen, never pushes tracks out of reach). With more
    /// tracks than fit anyway, the track's minimum.
    pub(super) fn height_cap(&self, key: &Key, avail: f32) -> f32 {
        let others: f32 = self
            .arrangement
            .ordered(&self.track_keys())
            .iter()
            .filter(|k| *k != key)
            .filter_map(|k| Some(self.floor_of(k, self.hit_of(k)?)))
            .sum();
        (avail - others).max(min_height(key.0))
    }

    /// The height of everything a track has to show (all its rows), and whether it shows so
    /// little that it should be no taller than that (an alignment track zoomed out beyond the
    /// read window: coverage and a note). `None` while that isn't known yet (nothing loaded).
    pub(super) fn content_height(&self, hit: Hit) -> Option<(f32, bool)> {
        match hit {
            Hit::Track(i) => {
                let (vs, ve) = self.visible();
                if !self.contigs.is_empty() && ve - vs > super::read_window() as f64 {
                    return Some((READS_TOP + NOTE_H, true));
                }
                let cur = self.tracks[i].current.as_ref()?;
                let (reads, _) = cur.view.reads.as_ref()?;
                // Reads loaded around the view but none on it: coverage and "No reads in view",
                // not a tall empty band sized to rows off screen.
                if !super::reads_in_view(reads, &cur.packed, vs, ve) {
                    return Some((READS_TOP + NOTE_H, true));
                }
                let rl = crate::tracks::row_layout(&cur.packed, cur.style);
                Some((READS_TOP + rl.content_height() + 12.0, false))
            }
            Hit::Genes(i) => {
                let rows = self.genes[i].built.as_ref()?.rows.len().max(1);
                Some((genes::ROWS_TOP + rows as f32 * genes::PITCH + 2.0, false))
            }
            Hit::Variants(i) => {
                let v = &self.variants[i];
                Some((v.gt_top() + v.samples.len() as f32 * variants::GT_PITCH + 3.0, false))
            }
            Hit::Signal(i) => Some((self.signals[i].section_height(), false)),
        }
    }

    /// The least a track is squeezed to so it still reads: its header and one whole row (an
    /// alignment track: its coverage and two rows of reads).
    pub(super) fn readable_height(&self, hit: Hit) -> f32 {
        match hit {
            Hit::Track(_) => READS_TOP + 24.0,
            Hit::Genes(_) => genes::ROWS_TOP + genes::PITCH + 2.0,
            Hit::Variants(i) => {
                let v = &self.variants[i];
                v.gt_top() + v.samples.len().min(1) as f32 * variants::GT_PITCH + 3.0
            }
            Hit::Signal(i) => self.signals[i].section_height().min(HEADER_H + 30.0),
        }
    }

    /// Where every track goes in `rect`.
    ///
    /// Each track starts at the height the person set, or its natural one (alignment tracks
    /// without one: `MIN_SECTION_H`, or less when their reads need less). Then:
    /// - Too little room: alignment tracks give way first (their reads scroll), down to their
    ///   coverage and two rows; then the others, down to their header and one row. Only more
    ///   tracks than fit like that need the track-stack scrollbar.
    /// - Room to spare: it goes to tracks with more to show than fits (alignment tracks first:
    ///   their rows of reads; then gene and sample rows), each up to what it has to show. What
    ///   no track can use stays below the last track (whose bottom edge isn't a handle).
    pub(super) fn arrange(&self, rect: Rect) -> Layout {
        let top = rect.top() + TRACKS_TOP;
        let avail = rect.bottom() - top;
        let keys = self.arrangement.ordered(&self.track_keys());
        let mut items: Vec<(Key, Hit)> = Vec::new();
        let mut fits: Vec<Fit> = Vec::new();
        for k in keys {
            let Some(hit) = self.hit_of(&k) else { continue };
            let content = self.content_height(hit);
            let readable = self.readable_height(hit);
            let mut h = match self.arrangement.heights.get(&k).copied().or_else(|| self.natural_height(hit)) {
                Some(h) => h,
                None => MIN_SECTION_H.min(content.map_or(MIN_SECTION_H, |c| c.0)).max(readable),
            };
            if let Some((c, true)) = content {
                h = h.min(c);
            }
            let want = content.map_or(h, |c| c.0).max(h);
            fits.push(Fit { h, min: readable.min(h), want, align: k.0 == Kind::Alignments });
            items.push((k, hit));
        }
        let total: f32 = fits.iter().map(|f| f.h).sum();
        if total > avail {
            squeeze(&mut fits, total - avail);
        }
        // (After squeezing to whole points, what rounding left.)
        let total: f32 = fits.iter().map(|f| f.h).sum();
        if total < avail {
            grow(&mut fits, avail - total);
        }
        let mut layout = Layout {
            variants: vec![Rect::NOTHING; self.variants.len()],
            sections: vec![Rect::NOTHING; self.tracks.len()],
            genes: vec![Rect::NOTHING; self.genes.len()],
            signals: vec![Rect::NOTHING; self.signals.len()],
            order: Vec::with_capacity(items.len()),
            stack: None,
            avail,
        };
        // Too many to fit (sixteen BAMs, say): scroll by whole tracks, from `stack_first`, with a
        // scrollbar, so every track stays reachable.
        let total: f32 = fits.iter().map(|f| f.h).sum();
        let mut first = 0;
        let mut stack = None;
        if total > avail + 0.5 && items.len() > 1 {
            // The furthest first track that still fills the space below it.
            let mut last_first = items.len() - 1;
            let mut below = fits[last_first].h;
            while last_first > 0 && below + fits[last_first - 1].h <= avail {
                last_first -= 1;
                below += fits[last_first].h;
            }
            first = self.stack_first.min(last_first);
            let above: f32 = fits[..first].iter().map(|f| f.h).sum();
            let track = Rect::from_min_max(Pos2::new(rect.right() - STACK_BAR_W, top), Pos2::new(rect.right(), rect.bottom()));
            let t0 = track.top() + track.height() * above / total;
            let t1 = (t0 + track.height() * avail / total).min(track.bottom());
            stack = Some(StackBar { track, thumb: Rect::from_min_max(Pos2::new(track.left(), t0), Pos2::new(track.right(), t1)), first, last_first, count: items.len() });
            // Scrolled to the end, the last tracks may leave room below them: it goes to the
            // ones shown that have more to show, not to an empty band.
            let shown: f32 = fits[first..].iter().map(|f| f.h).sum();
            if shown < avail {
                grow(&mut fits[first..], avail - shown);
            }
        }
        let rect = if stack.is_some() { Rect::from_min_max(rect.min, Pos2::new(rect.right() - STACK_BAR_W, rect.bottom())) } else { rect };
        let heights: Vec<f32> = fits.iter().map(|f| f.h).collect();
        let mut y = top;
        for (k, (key, hit)) in items.into_iter().enumerate() {
            if k < first {
                continue;
            }
            let r = Rect::from_min_max(Pos2::new(rect.left(), y), Pos2::new(rect.right(), y + heights[k]));
            y = r.bottom();
            match hit {
                Hit::Variants(i) => layout.variants[i] = r,
                Hit::Track(i) => layout.sections[i] = r,
                Hit::Genes(i) => layout.genes[i] = r,
                Hit::Signal(i) => layout.signals[i] = r,
            }
            layout.order.push((key, r));
        }
        layout.stack = stack;
        layout
    }

    /// The track-stack scrollbar: drag or click it, or turn the wheel over it. True when it
    /// took the pointer (so the view doesn't zoom or pan).
    pub(super) fn stack_input(&mut self, ui: &egui::Ui, response: &egui::Response, layout: &Layout) -> bool {
        let Some(bar) = &layout.stack else {
            self.stack_first = 0;
            return false;
        };
        let pointer = response.interact_pointer_pos().or(response.hover_pos());
        let Some(p) = pointer.filter(|p| bar.track.contains(*p)) else { return false };
        let per = bar.track.height() / bar.count as f32;
        if response.is_pointer_button_down_on() || response.clicked() {
            // The track under the pointer goes to the thumb's middle.
            let k = ((p.y - bar.track.top() - bar.thumb.height() / 2.0) / per).round().max(0.0) as usize;
            self.stack_first = k.min(bar.last_first);
        }
        let wheel = ui.input(|i| i.smooth_scroll_delta().y);
        if wheel.abs() > 0.5 {
            self.stack_first = if wheel < 0.0 { (bar.first + 1).min(bar.last_first) } else { bar.first.saturating_sub(1) };
        }
        true
    }

    pub(super) fn draw_stack_bar(&self, painter: &egui::Painter, layout: &Layout) {
        let Some(bar) = &layout.stack else { return };
        let pal = self.pal();
        // A trough that shows on white and on dark, so how far there is to go can be judged.
        painter.rect_filled(bar.track, 0.0, pal.rule);
        painter.vline(bar.track.left() + 0.5, bar.track.y_range(), Stroke::new(1.0, pal.rule_strong));
        painter.rect_filled(bar.thumb.shrink2(egui::vec2(2.0, 1.0)).translate(egui::vec2(0.5, 0.0)), 3.0, pal.muted);
    }

    /// What can be taken hold of: a grip in each header (drag to move the track) and a notch on
    /// each edge between tracks (drag to resize), both stronger under the pointer, and the edge
    /// being dragged drawn as a line.
    pub(super) fn draw_affordances(&self, painter: &egui::Painter, layout: &Layout, hover: Option<Pos2>) {
        let pal = self.pal();
        let accent = egui::Color32::from_rgb(64, 132, 230);
        let resizing = match &self.arranging {
            Some(Arranging::Resize { key, .. }) => Some(key),
            Some(Arranging::Move { .. }) => return,
            None => None,
        };
        let header = hover.and_then(|p| layout.header_at(p)).map(|(k, _)| k);
        for (k, r) in &layout.order {
            let grip = Layout::grip(*r);
            let c = if header == Some(k) { pal.text } else { pal.muted };
            for col in 0..2 {
                for row in 0..3 {
                    let at = Pos2::new(grip.left() + 1.5 + col as f32 * 4.0, grip.top() + 2.5 + row as f32 * 4.0);
                    painter.circle_filled(at, 1.0, c);
                }
            }
        }
        let edge = hover.and_then(|p| layout.edge_at(p)).map(|(k, _)| k);
        for (k, r) in &layout.order[..layout.resize_edges()] {
            let y = r.bottom();
            let active = resizing == Some(k) || (resizing.is_none() && edge == Some(k));
            let notch = Rect::from_center_size(Pos2::new(r.center().x, y), Vec2::new(28.0, if active { 5.0 } else { 3.0 }));
            if active {
                painter.hline(r.x_range(), y, Stroke::new(2.0, accent.gamma_multiply(0.6)));
                painter.rect_filled(notch, 2.0, accent);
            } else {
                painter.rect_filled(notch, 1.5, pal.rule_strong);
            }
        }
    }

    /// Starts, continues or ends a move or resize. Returns whether it used the pointer (so the
    /// view doesn't pan or select as well).
    pub(super) fn arrange_input(&mut self, ui: &egui::Ui, response: &egui::Response, layout: &Layout) -> bool {
        let pointer = response.hover_pos().or(response.interact_pointer_pos());
        let on_edge = |p: Pos2| layout.edge_at(p).cloned();
        let on_header = |p: Pos2| layout.header_at(p).cloned();

        // Cursors say what a drag would do.
        if self.arranging.is_none()
            && let Some(p) = response.hover_pos()
        {
            if on_edge(p).is_some() {
                ui.ctx().set_cursor_icon(egui::CursorIcon::ResizeVertical);
            } else if on_header(p).is_some() {
                ui.ctx().set_cursor_icon(egui::CursorIcon::Grab);
            }
        }

        // Double-click an edge: back to the automatic height.
        if response.double_clicked()
            && let Some((key, _)) = pointer.and_then(on_edge)
        {
            self.arrangement.heights.remove(&key);
            return true;
        }

        // Where the button went down: a drag is only reported once the pointer has moved a few
        // pixels, by when it may have left a thin header or edge.
        if response.drag_started()
            && let Some(p) = ui.input(|i| i.pointer.press_origin()).or(response.interact_pointer_pos())
        {
            if let Some((key, r)) = on_edge(p) {
                // What is drawn becomes what is kept: heights squeezed to fit the window are
                // stored as drawn, so the edge follows the pointer one to one.
                for (k, kr) in &layout.order {
                    if self.arrangement.heights.contains_key(k) {
                        self.arrangement.heights.insert(k.clone(), kr.height().round());
                    }
                }
                let below = layout.order.iter().skip_while(|(k, _)| *k != key).nth(1).filter(|(_, b)| (b.top() - r.bottom()).abs() < 1.0).map(|(k, b)| (k.clone(), b.height()));
                self.arranging = Some(Arranging::Resize { key, start_h: r.height(), start_y: p.y, below });
            } else if let Some((key, _)) = on_header(p) {
                self.arranging = Some(Arranging::Move { key, y: p.y });
            }
        }
        let Some(arranging) = self.arranging.clone() else { return false };
        let p = response.interact_pointer_pos().or(pointer);
        match arranging.clone() {
            Arranging::Resize { key, start_h, start_y, below } => {
                ui.ctx().set_cursor_icon(egui::CursorIcon::ResizeVertical);
                if let Some(p) = p {
                    let dy = p.y - start_y;
                    match below {
                        // A splitter: what one track gains the other loses, each down to its
                        // minimum, and the tracks further away don't move.
                        Some((bkey, bh)) => {
                            let dy = dy.clamp(min_height(key.0) - start_h, bh - min_height(bkey.0)).round();
                            let a_flex = key.0 == Kind::Alignments && !self.arrangement.heights.contains_key(&key);
                            let b_flex = bkey.0 == Kind::Alignments && !self.arrangement.heights.contains_key(&bkey);
                            // An alignment track sharing the free space keeps sharing it when it
                            // alone shares it, so it takes the change exactly; with others
                            // sharing too it would take only its part, and the edge would lag the
                            // pointer: then both sides get their height.
                            let flex = |k: &Key| k.0 == Kind::Alignments && !self.arrangement.heights.contains_key(k);
                            let others = layout.order.iter().filter(|(k, _)| *k != key && *k != bkey && flex(k)).count();
                            let (a_flex, b_flex) = (a_flex && others == 0, b_flex && others == 0);
                            if !a_flex || b_flex {
                                self.arrangement.heights.insert(key.clone(), start_h + dy);
                            }
                            if !b_flex || a_flex {
                                self.arrangement.heights.insert(bkey, bh - dy);
                            }
                        }
                        // The last track: it alone, no taller than leaves the others on screen.
                        None => {
                            let cap = self.height_cap(&key, layout.avail);
                            let h = (start_h + dy).min(cap).max(min_height(key.0)).round();
                            self.arrangement.heights.insert(key, h);
                        }
                    }
                }
            }
            Arranging::Move { key, .. } => {
                ui.ctx().set_cursor_icon(egui::CursorIcon::Grabbing);
                if let Some(p) = p {
                    self.arranging = Some(Arranging::Move { key, y: p.y });
                }
            }
        }
        if response.drag_stopped() {
            if let Arranging::Move { key, y } = arranging {
                let mut order: Vec<Key> = layout.order.iter().map(|(k, _)| k.clone()).collect();
                let target = drop_index(layout, y);
                if let Some(from) = order.iter().position(|k| *k == key) {
                    let k = order.remove(from);
                    order.insert(if target > from { target - 1 } else { target }, k);
                }
                self.arrangement.order = order;
            }
            self.arranging = None;
        }
        true
    }

    /// While moving a track: where it would land.
    pub(super) fn draw_arranging(&self, painter: &egui::Painter, layout: &Layout) {
        let Some(Arranging::Move { key, y }) = &self.arranging else { return };
        let accent = egui::Color32::from_rgb(64, 132, 230);
        if let Some((_, r)) = layout.order.iter().find(|(k, _)| k == key) {
            painter.rect_filled(*r, 0.0, accent.gamma_multiply(0.12));
        }
        let at = drop_index(layout, *y);
        let line_y = match layout.order.get(at) {
            Some((_, r)) => r.top(),
            None => layout.order.last().map_or(0.0, |(_, r)| r.bottom()),
        };
        if let Some((_, r)) = layout.order.first() {
            painter.hline(r.x_range(), line_y, Stroke::new(3.0, accent));
            painter.rect_filled(Rect::from_center_size(Pos2::new(r.left() + 6.0, line_y), Vec2::splat(8.0)), 2.0, accent);
        }
    }
}

impl Layout {
    /// How many tracks (from the top) have a bottom edge that resizes: all but the last, whose
    /// bottom is the end of the stack when all fit (there it would resize nothing visible, and
    /// space below it is never an invisible handle); with more than fit, all.
    pub(super) fn resize_edges(&self) -> usize {
        if self.stack.is_some() { self.order.len() } else { self.order.len().saturating_sub(1) }
    }

    /// The track whose bottom edge `p` is on, if that edge resizes.
    pub(super) fn edge_at(&self, p: Pos2) -> Option<&(Key, Rect)> {
        self.order[..self.resize_edges()].iter().find(|(_, r)| (p.y - r.bottom()).abs() <= EDGE && r.x_range().contains(p.x))
    }

    /// The track whose header (a handle for moving it) `p` is on, its close button aside.
    pub(super) fn header_at(&self, p: Pos2) -> Option<&(Key, Rect)> {
        self.order.iter().find(|(_, r)| (r.top()..r.top() + HEADER_H).contains(&p.y) && r.x_range().contains(p.x) && !Layout::close_button(*r).contains(p))
    }

    /// The grip drawn in a header, just left of the close button.
    pub(super) fn grip(section: Rect) -> Rect {
        Rect::from_min_size(Pos2::new(Layout::close_button(section).left() - 12.0, section.top() + 2.0), Vec2::new(8.0, HEADER_H - 4.0))
    }
}

/// The smallest a person can drag a track of `kind` to.
fn min_height(kind: Kind) -> f32 {
    if kind == Kind::Alignments { MIN_ALIGN_H } else { MIN_OTHER_H }
}

/// One track's height while the stack is arranged: what it has (`h`), the least it reads at
/// (`min`), and what it has to show (`want`).
#[derive(Clone, Copy, Debug, PartialEq)]
struct Fit {
    h: f32,
    min: f32,
    want: f32,
    align: bool,
}

/// Takes `excess` from the tracks, alignment tracks first (their reads scroll), each down to
/// its `min`, in proportion to what it has above that.
fn squeeze(fits: &mut [Fit], mut excess: f32) {
    for align in [true, false] {
        let give: f32 = fits.iter().filter(|f| f.align == align).map(|f| f.h - f.min).sum();
        if excess <= 0.0 || give <= 0.0 {
            continue;
        }
        let take = excess.min(give);
        let keep = 1.0 - take / give;
        for f in fits.iter_mut().filter(|f| f.align == align) {
            f.h = (f.min + (f.h - f.min) * keep).floor().max(f.min);
        }
        excess -= take;
    }
}

/// Shares `spare` among the tracks with more to show than fits, alignment tracks first, each up
/// to its `want`, evenly. What none can use is left over.
fn grow(fits: &mut [Fit], mut spare: f32) {
    for align in [true, false] {
        loop {
            let n = fits.iter().filter(|f| f.align == align && f.want - f.h > 0.5).count();
            if n == 0 || spare < 0.5 {
                break;
            }
            let share = spare / n as f32;
            for f in fits.iter_mut().filter(|f| f.align == align && f.want - f.h > 0.5) {
                let add = share.min(f.want - f.h);
                f.h += add;
                spare -= add;
            }
        }
    }
}

/// Where a track dropped at `y` goes: before the first track whose middle is below `y`.
fn drop_index(layout: &Layout, y: f32) -> usize {
    layout.order.iter().position(|(_, r)| y < r.center().y).unwrap_or(layout.order.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn k(kind: Kind, p: &str) -> Key {
        (kind, PathBuf::from(p))
    }

    #[test]
    fn new_tracks_join_their_kind_and_moved_ones_stay() {
        let default = vec![k(Kind::Variants, "v.vcf"), k(Kind::Alignments, "a.bam"), k(Kind::Alignments, "b.bam"), k(Kind::Genes, "g.gtf")];
        // Nothing arranged: the default order.
        assert_eq!(Arrangement::default().ordered(&default), default);
        // Genes moved to the top; a new BAM joins the other BAMs, a closed file drops out.
        let arranged = Arrangement { order: vec![k(Kind::Genes, "g.gtf"), k(Kind::Variants, "v.vcf"), k(Kind::Alignments, "a.bam"), k(Kind::Alignments, "gone.bam")], ..Default::default() };
        assert_eq!(arranged.ordered(&default), vec![k(Kind::Genes, "g.gtf"), k(Kind::Variants, "v.vcf"), k(Kind::Alignments, "a.bam"), k(Kind::Alignments, "b.bam")]);
        // Saved and restored.
        let mut a = arranged.clone();
        a.heights.insert(k(Kind::Genes, "g.gtf"), 200.0);
        let back = Arrangement::from_saved(&a.saved(&a.ordered(&default)));
        assert_eq!(back.ordered(&default), a.ordered(&default));
        assert_eq!(back.heights.get(&k(Kind::Genes, "g.gtf")), Some(&200.0));
    }
}
