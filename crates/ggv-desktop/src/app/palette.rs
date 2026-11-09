//! The command palette (Ctrl+K): one fuzzy box for loci, genes, contigs, bookmarks, recent
//! files, saved filters and every action. Results are recomputed only when the query changes;
//! gene names are lowercased once per annotation and kept.

use std::{path::PathBuf, sync::Arc};

use eframe::egui::{self, Align2, Key, Vec2};

use super::{App, actions::Action};
use crate::util::commas;

/// How many results to show.
const SHOWN: usize = 12;
/// How many to rank (more than shown, so arrowing down past the first screen works).
const KEPT: usize = 60;

#[derive(Clone, Debug)]
pub enum Item {
    Action(Action),
    Locus(String),
    Gene { name: String, contig: String, start: u64, end: u64 },
    Contig(usize),
    Bookmark(usize),
    Recent(PathBuf),
    Filter(String),
}

#[derive(Default)]
pub struct Palette {
    pub query: String,
    pub(super) results: Vec<(Item, String, String)>,
    /// The query and data (annotations loaded, bookmarks, contigs) the results were ranked for.
    pub(super) computed_for: Option<(String, usize)>,
    pub(super) selected: usize,
    /// Lowercased gene names per annotation, keyed by the annotation's address.
    pub(super) genes: Vec<(usize, Arc<Vec<String>>)>,
    pub(super) focused: bool,
}

/// Fuzzy score of `query` (lowercase) against `text` (lowercase): prefix beats word start beats
/// substring beats scattered letters; `None` if the letters don't all appear in order.
pub fn score(query: &str, text: &str) -> Option<i32> {
    if query.is_empty() {
        return Some(0);
    }
    if text.starts_with(query) {
        return Some(1000 - (text.len() as i32 - query.len() as i32).min(200));
    }
    if let Some(i) = text.find(query) {
        let word_start = i == 0 || !text.as_bytes()[i - 1].is_ascii_alphanumeric();
        return Some(if word_start { 800 } else { 600 } - (i as i32).min(100));
    }
    // Letters in order, scored by how bunched together they are.
    let t = text.as_bytes();
    let mut j = 0;
    let mut gaps = 0;
    let mut last = None;
    for &q in query.as_bytes() {
        while j < t.len() && t[j] != q {
            j += 1;
        }
        if j == t.len() {
            return None;
        }
        if let Some(l) = last {
            gaps += j - l - 1;
        }
        last = Some(j);
        j += 1;
    }
    Some(300 - (gaps as i32).min(250))
}

impl App {
    fn gene_names(&self, p: &mut Palette) {
        for g in &self.genes {
            let Some(a) = g.annotation() else { continue };
            let key = Arc::as_ptr(a) as usize;
            if !p.genes.iter().any(|(k, _)| *k == key) {
                p.genes.push((key, Arc::new(a.genes().iter().map(|s| s.name.to_ascii_lowercase()).collect())));
            }
        }
    }

    fn palette_results(&mut self, p: &mut Palette) {
        let q = p.query.trim().to_ascii_lowercase();
        let mut scored: Vec<(i32, Item, String, String)> = Vec::new();

        // A locus typed out: go there first.
        if !q.is_empty() && !self.contigs.is_empty() && ggv_io::resolve_locus_str(&self.contigs, p.query.trim()).is_ok() && p.query.contains(':') {
            scored.push((5000, Item::Locus(p.query.trim().to_owned()), format!("Go to {}", p.query.trim()), "locus".into()));
        }

        // Genes: prefix and substring matches only (scattered letters over 60k names is noise).
        if q.len() >= 2 {
            self.gene_names(p);
            let mut genes: Vec<(i32, usize, usize)> = Vec::new();
            for (gi, g) in self.genes.iter().enumerate() {
                let Some(a) = g.annotation() else { continue };
                let key = Arc::as_ptr(a) as usize;
                let Some((_, names)) = p.genes.iter().find(|(k, _)| *k == key) else { continue };
                for (i, n) in names.iter().enumerate() {
                    if n.len() >= q.len() && n.contains(&q) {
                        let s = score(&q, n).unwrap_or(0) + if *n == q { 2000 } else { 0 };
                        genes.push((s + 50, gi, i));
                    }
                }
            }
            genes.sort_by_key(|g| std::cmp::Reverse(g.0));
            genes.dedup_by(|a, b| {
                let name = |x: &(i32, usize, usize)| self.genes[x.1].annotation().map(|an| an.genes()[x.2].name.clone());
                name(a) == name(b)
            });
            for (s, gi, i) in genes.into_iter().take(KEPT) {
                let Some(a) = self.genes[gi].annotation() else { continue };
                let span = &a.genes()[i];
                scored.push((
                    s,
                    Item::Gene { name: span.name.clone(), contig: span.contig.clone(), start: span.start, end: span.end },
                    span.name.clone(),
                    format!("gene  {}:{}-{}", span.contig, commas(span.start + 1), commas(span.end)),
                ));
            }
        }

        for (i, c) in self.contigs.iter().enumerate() {
            if let Some(s) = score(&q, &c.name.to_ascii_lowercase()).filter(|_| !q.is_empty()) {
                scored.push((s - 100, Item::Contig(i), c.name.clone(), format!("contig  {} bp", commas(c.length))));
            }
        }
        for (i, b) in self.bookmarks.borrow().iter().enumerate() {
            let text = format!("{} {}", b.note, b.locus()).to_ascii_lowercase();
            if let Some(s) = score(&q, &text) {
                let label = if b.note.is_empty() { b.locus() } else { b.note.clone() };
                scored.push((s + 20, Item::Bookmark(i), label, format!("bookmark  {}", b.locus())));
            }
        }
        for a in Action::all() {
            let name = a.name();
            if let Some(s) = score(&q, &name.to_ascii_lowercase()) {
                let keys = a.keys().map(super::actions::keys_label).unwrap_or_default();
                scored.push((s, Item::Action(a), name, keys));
            }
        }
        for f in self.saved_filters.clone() {
            if let Some(s) = score(&q, &format!("filter {}", f.to_ascii_lowercase())) {
                scored.push((s - 10, Item::Filter(f.clone()), format!("Filter: {f}"), "saved filter".into()));
            }
        }
        for r in &self.recent {
            let name = crate::files::display_name(r);
            if let Some(s) = score(&q, &format!("open {}", name.to_ascii_lowercase())) {
                scored.push((s - 30, Item::Recent(r.clone()), format!("Open {name}"), "recent file".into()));
            }
        }
        // Stable: equal scores keep the order above (genes, contigs, bookmarks, actions…).
        scored.sort_by_key(|x| std::cmp::Reverse(x.0));
        scored.truncate(KEPT);
        p.results = scored.into_iter().map(|(_, i, l, d)| (i, l, d)).collect();
        p.selected = 0;
        p.computed_for = Some((p.query.clone(), self.palette_data_key()));
    }

    /// Changes when something the palette searches finishes loading or is added.
    fn palette_data_key(&self) -> usize {
        self.genes.iter().filter(|g| g.annotation().is_some()).count() * 1_000_003
            + self.bookmarks.borrow().len() * 1_009
            + self.contigs.len()
            + self.saved_filters.len() * 7
    }

    fn pick(&mut self, item: Item, ctx: &egui::Context) {
        match item {
            Item::Action(a) => self.run(a, ctx),
            Item::Locus(l) => self.pending_goto = Some(l),
            Item::Gene { name, contig, start, end } => {
                self.notify(name);
                let aliases = ggv_model::contig_aliases(&contig);
                if let Some(c) = self.contigs.iter().position(|x| x.name == contig || aliases.contains(&x.name)) {
                    self.remember_place();
                    let pad = ((end - start) / 20).max(50);
                    self.contig = c;
                    self.center = (start + end) as f64 / 2.0;
                    self.ppb = f64::from(self.width) / ((end - start + 2 * pad) as f64);
                    self.highlight = None;
                    self.clamp_view();
                } else {
                    self.error = Some(format!("{contig} is not in this genome"));
                }
            }
            Item::Contig(i) => {
                if let Some(c) = self.contigs.get(i) {
                    self.pending_goto = Some(c.name.clone());
                }
            }
            Item::Bookmark(i) => self.go_to_bookmark(i),
            Item::Recent(p) => self.open_paths(vec![p], ctx),
            Item::Filter(f) => {
                if let Some(e) = self.apply_filter_text(&f, None) {
                    self.error = Some(e);
                }
            }
        }
    }

    pub(super) fn palette_window(&mut self, ctx: &egui::Context) {
        let Some(mut p) = self.palette.take() else { return };
        if !self.is_active {
            return; // dropped: only the active pane has one
        }
        if p.computed_for.as_ref().is_none_or(|(q, k)| *q != p.query || *k != self.palette_data_key()) {
            self.palette_results(&mut p);
        }
        let (up, down, enter, esc) = ctx.input_mut(|i| {
            (
                i.consume_key(egui::Modifiers::NONE, Key::ArrowUp),
                i.consume_key(egui::Modifiers::NONE, Key::ArrowDown),
                i.consume_key(egui::Modifiers::NONE, Key::Enter),
                i.consume_key(egui::Modifiers::NONE, Key::Escape),
            )
        });
        let n = p.results.len();
        if n > 0 {
            if down {
                p.selected = (p.selected + 1) % n;
            }
            if up {
                p.selected = (p.selected + n - 1) % n;
            }
        }
        let mut chosen = enter.then_some(p.selected).filter(|&i| i < n);
        let mut close = esc;
        let width = (ctx.content_rect().width() - 32.0).clamp(240.0, 620.0);
        egui::Window::new("palette")
            .id(egui::Id::unique(("palette", self.reference_slot)))
            .title_bar(false)
            .anchor(Align2::CENTER_TOP, Vec2::new(0.0, 60.0))
            .fixed_size(Vec2::new(width, 0.0))
            .resizable(false)
            .show(ctx, |ui| {
                let edit = ui.add(
                    egui::TextEdit::singleline(&mut p.query)
                        .desired_width(f32::INFINITY)
                        .hint_text("Go to a gene, locus, bookmark, or type a command…")
                        .font(egui::TextStyle::Heading),
                );
                if !p.focused {
                    edit.request_focus();
                    p.focused = true;
                }
                ui.separator();
                let first = p.selected.saturating_sub(SHOWN - 1);
                for (i, (_, label, detail)) in p.results.iter().enumerate().skip(first).take(SHOWN) {
                    let selected = i == p.selected;
                    let r = ui
                        .horizontal(|ui| {
                            let resp = ui.selectable_label(selected, label);
                            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                                ui.label(egui::RichText::new(detail).weak().monospace());
                            });
                            resp
                        })
                        .inner;
                    if selected {
                        r.scroll_to_me(None);
                    }
                    if r.clicked() {
                        chosen = Some(i);
                    }
                }
                if p.results.is_empty() {
                    ui.label(egui::RichText::new("No matches").weak());
                } else if n > SHOWN {
                    ui.label(egui::RichText::new(format!("{} of {n} shown · Up/Down to move · Enter to choose · Esc to close", SHOWN.min(n))).weak().small());
                }
            });
        if let Some(i) = chosen {
            let item = p.results[i].0.clone();
            close = true;
            self.pick(item, ctx);
        }
        if !close && self.palette.is_none() {
            self.palette = Some(p);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::score;

    #[test]
    fn ranking_prefers_prefix_then_word_then_scattered() {
        let prefix = score("brca", "brca1").unwrap();
        let word = score("brca", "gene brca1").unwrap();
        let inner = score("rca", "brca1").unwrap();
        let scattered = score("cpl", "command palette").unwrap();
        assert!(prefix > word && word > inner && inner > scattered, "{prefix} {word} {inner} {scattered}");
        assert!(score("zz", "brca1").is_none());
        assert!(score("crm", "copy region as fasta").is_none());
        assert_eq!(score("", "anything"), Some(0));
        // Tighter scatter wins.
        assert!(score("cop", "copy").unwrap() > score("cpy", "copy").unwrap());
        assert!(score("cpy", "copy").unwrap() > score("cpy", "c__p__y").unwrap());
    }

    #[test]
    fn gene_matching_keeps_up_with_typing() {
        // GGV_TEST_GTF (e.g. Ensembl's full GTF): the per-keystroke scan the palette does over
        // every gene name must stay well under a frame.
        let Ok(path) = std::env::var("GGV_TEST_GTF") else { return };
        let a = ggv_annot::Annotation::read(&path).unwrap();
        let names: Vec<String> = a.genes().iter().map(|g| g.name.to_ascii_lowercase()).collect();
        for q in ["br", "brca", "tp53", "znf", "linc"] {
            let t = std::time::Instant::now();
            let mut hits: Vec<(i32, usize)> = names
                .iter()
                .enumerate()
                .filter(|(_, n)| n.contains(q))
                .map(|(i, n)| (score(q, n).unwrap_or(0), i))
                .collect();
            hits.sort_by_key(|h| std::cmp::Reverse(h.0));
            let took = t.elapsed();
            eprintln!("{q:>5}: {} of {} genes in {took:.2?}, first {}", hits.len(), names.len(), a.genes()[hits[0].1].name);
            assert!(took.as_millis() < 16, "{q}: {took:?}");
        }
    }
}
