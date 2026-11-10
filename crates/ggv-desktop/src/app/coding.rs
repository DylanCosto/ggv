//! What a base or a change means for the genes there: codons and amino acids on the gene track
//! when zoomed in, and HGVS consequences (`c.1234A>G p.Lys412Arg missense`) in the evidence
//! panel, tooltips and reports.
//!
//! Reference bases come from the loaded reference row; a codon split across a distant intron
//! falls back to a small FASTA read. Consequences for the evidence base are cached, so the
//! panel costs nothing per frame once computed.

use eframe::egui::{self, Align2, Color32, FontId, Pos2, Rect};
use ggv_annot::{
    Strand, Transcript,
    coding::{self, Consequence},
};

use super::App;

/// Draw codons when each is at least this many points wide.
const CODON_MIN_PT: f64 = 9.0;

impl App {
    /// The reference base at `p` on the view's contig, from the loaded reference row.
    pub(super) fn loaded_base(&self, p: u64) -> Option<u8> {
        let (view, _, _) = self.reference.as_ref()?;
        let r = &view.region;
        if r.contig != self.contigs.get(self.contig)?.name || p < r.start || p >= r.end {
            return None;
        }
        view.reference.as_ref()?.get((p - r.start) as usize).copied()
    }

    /// The reference base at `p` on the view's contig, reading the FASTA when it isn't loaded
    /// (rare: the other part of a codon split across a long intron).
    fn base_anywhere(&self, p: u64) -> Option<u8> {
        self.base_on(&self.contigs.get(self.contig)?.name, p)
    }

    /// The reference base at `p` on `contig`: the loaded row's when the view is there, else read
    /// from the FASTA.
    fn base_on(&self, contig: &str, p: u64) -> Option<u8> {
        if self.contigs.get(self.contig).is_some_and(|c| c.name == contig)
            && let Some(b) = self.loaded_base(p)
        {
            return Some(b);
        }
        let path = self.reference_path.as_ref()?;
        let contig = contig.to_owned();
        let mut cache = self.ref_reader.borrow_mut();
        if cache.as_ref().is_none_or(|(p0, _)| p0 != path) {
            *cache = Some((path.clone(), ggv_io::ReferenceFasta::open(path).ok()?));
        }
        let (_, fasta) = cache.as_mut()?;
        fasta.fetch(&ggv_model::Region { contig, start: p, end: p + 1 }).ok()?.first().copied()
    }

    /// Transcripts overlapping `[start, end)` on `contig` from the gene tracks, coding ones
    /// first, each once: all of them, or one per gene as an agent's `evidence` names them (the
    /// gene's main transcript, or where that one doesn't reach, its best one there, so a coding
    /// exon the main transcript lacks is still a gene). Off the view (a variant's details window
    /// left open) they come from the annotation itself.
    fn transcripts_over(&self, contig: &str, start: u64, end: u64) -> Vec<&Transcript> {
        let end = end.max(start + 1);
        let overlaps = |t: &&Transcript| t.start < end && start < t.end;
        let mut out: Vec<&Transcript> = Vec::new();
        for g in &self.genes {
            let here: Vec<&Transcript> = match (&g.built, g.annotation()) {
                (Some(b), _) if self.all_transcripts && b.contig == contig && b.start <= start && end <= b.end => b.items.iter().filter(overlaps).collect(),
                (_, Some(a)) if self.all_transcripts => a.query(contig, start, end).filter(overlaps).collect(),
                (_, Some(a)) => ggv_annot::one_per_gene(a.query(contig, start, end).filter(overlaps)),
                _ => continue,
            };
            for t in here {
                if !out.iter().any(|x| x.id == t.id) {
                    out.push(t);
                }
            }
        }
        out.sort_by_key(|t| (t.cds.is_none(), !t.representative));
        out
    }

    /// Consequences of `reference` → `alt` (VCF style, at `pos` on `contig`) for every
    /// transcript there.
    pub(super) fn consequences(&self, contig: &str, pos: u64, reference: &str, alt: &str) -> Vec<Consequence> {
        // Every transcript the REF touches, not only the one at its first base.
        self.transcripts_over(contig, pos, pos + reference.len() as u64)
            .into_iter()
            .filter_map(|t| coding::consequence(t, pos, reference, alt, |p| self.base_on(contig, p)))
            .collect()
    }

    /// The changes seen at `pos`: every ALT of the VCF records there, and each track's top
    /// alternate allele (as VCF-style ref/alt), each once.
    fn changes_at(&self, pos: u64) -> Vec<(u64, String, String, String)> {
        let mut out: Vec<(u64, String, String, String)> = Vec::new();
        let mut add = |p: u64, r: String, a: String, source: String| {
            if !out.iter().any(|x| x.0 == p && x.1 == r && x.2 == a) {
                out.push((p, r, a, source));
            }
        };
        for v in &self.variants {
            let Some(cur) = &v.current else { continue };
            for var in cur.variants.iter().filter(|x| x.start == pos) {
                for a in var.alts.iter().filter(|a| !a.starts_with('<') && *a != "*") {
                    add(var.start, var.reference.clone(), a.clone(), v.file_name());
                }
            }
        }
        for t in &self.tracks {
            let Some(s) = t.evidence.as_ref().filter(|e| e.pos == pos).and_then(|e| e.stats.as_ref()) else { continue };
            let Some(alt) = &s.top_alt else { continue };
            let Some(r) = s.ref_base.filter(|b| b.is_ascii_alphabetic()) else { continue };
            if alt.allele == "DEL" {
                if let Some(prev) = pos.checked_sub(1).and_then(|p| self.base_anywhere(p)) {
                    add(pos - 1, format!("{}{r}", prev as char), (prev as char).to_string(), format!("{} reads", t.short_label()));
                }
            } else if alt.allele.len() == 1 {
                add(pos, r.to_string(), alt.allele.to_owned(), format!("{} reads", t.short_label()));
            }
        }
        out
    }

    /// Lines for the evidence panel and reports: one per change and transcript. Cached by base,
    /// changes, loaded genes and reference.
    pub(super) fn consequence_lines(&self, pos: u64) -> Vec<String> {
        let changes = self.changes_at(pos);
        let key = (
            self.contig,
            pos,
            changes.iter().map(|c| format!("{}:{}>{}", c.0, c.1, c.2)).collect::<Vec<_>>().join(";"),
            self.genes.iter().map(|g| g.built.as_ref().map_or(0, |b| b.items.len())).sum::<usize>(),
            // Codons need reference bases: opening a reference changes what's known.
            self.reference_path.clone(),
        );
        if let Some((k, lines)) = &*self.consequence_cache.borrow()
            && *k == key
        {
            return lines.clone();
        }
        let mut lines = Vec::new();
        let contig = self.contigs.get(self.contig).map(|c| c.name.clone()).unwrap_or_default();
        for (p, r, a, source) in &changes {
            for c in self.consequences(&contig, *p, r, a).into_iter().take(4) {
                lines.push(format!(
                    "{} {}  {}{}  {}  ({}, {}){}",
                    c.gene,
                    c.transcript,
                    c.hgvs_c,
                    c.hgvs_p.map(|p| format!("  {p}")).unwrap_or_default(),
                    c.effect,
                    c.location.feature,
                    source,
                    if c.warnings.is_empty() { String::new() } else { format!("  [warning: {}]", c.warnings.join("; ")) }
                ));
            }
        }
        *self.consequence_cache.borrow_mut() = Some((key, lines.clone()));
        lines
    }

    /// The evidence panel's genes at `pos`: transcript (without UCSC's copy suffix, as the
    /// consequences name it), exon or intron, and whether it is coding there.
    pub(super) fn genes_at(&self, pos: u64) -> Vec<String> {
        let Some(contig) = self.contigs.get(self.contig).map(|c| c.name.clone()) else { return Vec::new() };
        // The same transcripts the consequences are worked out for: where a gene's main
        // transcript doesn't reach, its best one here, not nothing.
        self.transcripts_over(&contig, pos, pos + 1)
            .into_iter()
            .map(|t| {
                let exon = t.exons.iter().position(|&(s, e)| (s..e).contains(&pos)).map(|e| {
                    let n = t.exons.len();
                    if t.strand == Strand::Reverse { n - e } else { e + 1 }
                });
                let place = exon.map_or("intron".to_owned(), |e| format!("exon {e}/{}", t.exons.len()));
                let coding = exon.is_some() && t.cds.is_some_and(|(s, e)| (s..e).contains(&pos));
                format!("{} {}  {place}{}", t.gene_name, ggv_annot::without_copy_suffix(&t.id), if coding { ", coding" } else { "" })
            })
            .collect()
    }

    /// For a gene tooltip: the base's `c.` position and, in a coding exon, its codon and amino
    /// acid.
    pub(super) fn base_in_transcript(&self, t: &Transcript, pos: u64) -> Option<String> {
        let loc = coding::locate(t, pos)?;
        let prefix = if t.cds.is_some() { "c." } else { "n." };
        let mut s = format!("{prefix}{}  {}", loc.c, loc.kind.name());
        if let Some(ci) = coding::coding_index(t, pos)
            && let Some(cp) = coding::codon_positions(t, ci / 3)
        {
            let mut codon = cp.map(|p| self.loaded_base(p).unwrap_or(b'N'));
            if t.strand == Strand::Reverse {
                codon = codon.map(coding::complement);
            }
            let aa = coding::amino_acid(t, ci / 3, codon);
            if aa != b'X' {
                s.push_str(&format!("  codon {} {} ({})", ci / 3 + 1, String::from_utf8_lossy(&codon), coding::three_letter(aa)));
            }
        }
        Some(s)
    }

    /// Codons and amino acids over a coding transcript's exons, when zoomed in far enough.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn draw_codons(&self, painter: &egui::Painter, section: Rect, t: &Transcript, y: f32, h: f32, vs: f64, ve: f64) {
        let Some((cs, ce)) = t.cds else { return };
        if 3.0 * self.ppb < CODON_MIN_PT {
            return;
        }
        let dark = self.style.dark;
        // Alternate codons are lightened a little, so their boundaries show.
        let shade = [Color32::TRANSPARENT, Color32::from_rgba_unmultiplied(255, 255, 255, 45)];
        let start_fill = Color32::from_rgb(30, 150, 60);
        let stop_fill = Color32::from_rgb(200, 40, 40);
        let text = if dark { Color32::from_gray(20) } else { Color32::WHITE };
        let font = FontId::monospace((self.ppb as f32 * 2.0).clamp(8.0, 11.0));
        let lo = (vs.floor().max(0.0) as u64).max(cs);
        let hi = (ve.ceil() as u64).min(ce);
        for &(s, e) in &t.exons {
            let (a, b) = (s.max(lo), e.min(hi));
            if a >= b {
                continue;
            }
            let Some(ci0) = coding::coding_index(t, a) else { continue };
            for p in a..b {
                // Coding index along the transcript: up the genome on +, down it on −.
                let ci = if t.strand == Strand::Reverse { ci0 - (p - a) } else { ci0 + (p - a) };
                let k = ci / 3;
                let x0 = self.x_of(section, vs, p as f64);
                let x1 = self.x_of(section, vs, (p + 1) as f64);
                let cell = Rect::from_min_max(Pos2::new(x0, y), Pos2::new(x1, y + h));
                let cp = coding::codon_positions(t, k);
                let aa = cp.map(|cp| {
                    let mut c = cp.map(|q| self.loaded_base(q).unwrap_or(b'N'));
                    if t.strand == Strand::Reverse {
                        c = c.map(coding::complement);
                    }
                    coding::amino_acid(t, k, c)
                });
                let fill = match aa {
                    Some(b'*') => stop_fill,
                    Some(b'M') if k == 0 => start_fill,
                    _ => shade[(k % 2) as usize],
                };
                painter.rect_filled(cell, 0.0, fill);
                // The letter sits on the codon's middle base.
                if ci % 3 == 1
                    && let Some(aa) = aa.filter(|&x| x != b'X')
                {
                    painter.text(cell.center(), Align2::CENTER_CENTER, (aa as char).to_string(), font.clone(), text);
                }
            }
        }
    }
}
