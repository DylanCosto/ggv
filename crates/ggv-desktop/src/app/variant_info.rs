//! A variant's details, as IGV shows them on click: the site, ClinVar's summary when the record
//! is ClinVar's (significance, review stars, condition, gene, molecular consequence), what it
//! does to the genes open here, the genotypes, and every INFO field with the header's
//! description. The window stays until closed; clicking another variant replaces it.

use std::sync::Arc;

use eframe::egui;
use ggv_io::{InfoField, Variant};

use super::App;
use crate::{util::commas, variants::Significance};

pub struct VariantPopup {
    pub file: String,
    pub contig: String,
    pub var: Variant,
    pub fields: Arc<Vec<InfoField>>,
    pub samples: Vec<String>,
    /// Where on screen the variant was clicked, when it was: the window first opens on the
    /// pane's side away from it.
    clicked_x: Option<f32>,
    /// "In the genes open here", and what it was computed from (see `popup_consequences`).
    pub(super) consequences: Option<(ConsequenceKey, Vec<String>)>,
}

/// What a variant's consequences depend on besides the variant: the annotations loaded, one
/// transcript per gene or all, and the reference.
type ConsequenceKey = (usize, bool, Option<std::path::PathBuf>);

/// ClinVar's review status as gold stars (0-4), as ClinVar shows it.
pub fn review_stars(status: &str) -> u8 {
    match status {
        "practice_guideline" => 4,
        "reviewed_by_expert_panel" => 3,
        "criteria_provided,_multiple_submitters,_no_conflicts" => 2,
        s if s.starts_with("criteria_provided") => 1,
        _ => 0,
    }
}

/// ClinVar's text, readable: underscores as spaces, `|` between values.
pub fn readable(v: &str) -> String {
    v.replace('_', " ").replace('|', " | ")
}

/// The ClinVar summary of a record, as (label, value) lines; empty when it isn't ClinVar's.
pub fn clinvar_lines(var: &Variant) -> Vec<(&'static str, String)> {
    let Some(sig) = var.info_value("CLNSIG") else { return Vec::new() };
    let mut out = vec![("Significance", readable(sig))];
    if let Some(s) = var.info_value("CLNREVSTAT") {
        let n = review_stars(s);
        // As text: the default fonts have no star glyphs (they drew as empty boxes).
        out.push(("Review", format!("{n} of 4 stars  ({})", readable(s))));
    }
    if let Some(c) = var.info_value("CLNDN") {
        out.push(("Condition", readable(c)));
    }
    if let Some(g) = var.info_value("GENEINFO") {
        // `BRCA1:672|NBR2:10230`: symbols with their Gene IDs.
        let genes: Vec<&str> = g.split('|').map(|x| x.split(':').next().unwrap_or(x)).collect();
        out.push(("Gene", genes.join(", ")));
    }
    if let Some(mc) = var.info_value("MC") {
        // `SO:0001583|missense_variant,SO:0001819|synonymous_variant`.
        let terms: Vec<String> = mc.split(',').map(|x| x.split('|').nth(1).unwrap_or(x).replace('_', " ")).collect();
        out.push(("Consequence", terms.join(", ")));
    }
    if let Some(h) = var.info_value("CLNHGVS") {
        out.push(("HGVS", h.to_owned()));
    }
    if var.id != "." {
        out.push(("Variation ID", var.id.clone()));
    }
    out
}

impl App {
    /// Shows `var` from variant track `i` until closed or replaced.
    pub(super) fn open_variant_popup(&mut self, i: usize, var: Variant) {
        let Some(t) = self.variants.get(i) else { return };
        let contig = self.contigs.get(self.contig).map(|c| c.name.clone()).unwrap_or_default();
        self.variant_popup =
            Some(VariantPopup { file: t.file_name(), contig, var, fields: t.info_fields.clone(), samples: t.samples.clone(), clicked_x: None, consequences: None });
    }

    fn variant_window_id(&self) -> egui::Id {
        // One per pane (two panes' windows don't collide).
        egui::Id::unique(("variant details", self.reference_slot))
    }

    /// Shows `var` for a click on it at `at`, clear of that spot: the window first opens on the
    /// pane's side away from the click, and one that would come back over it (where it was when
    /// last closed) is placed again. One already open stays where it is.
    pub(super) fn open_variant_popup_at(&mut self, i: usize, var: Variant, at: egui::Pos2) {
        let reopening = self.variant_popup.is_none();
        self.open_variant_popup(i, var);
        let id = self.variant_window_id();
        let Some(p) = &mut self.variant_popup else { return };
        p.clicked_x = Some(at.x);
        if reopening && self.egui_ctx.memory(|m| m.area_rect(id)).is_some_and(|last| last.contains(at)) {
            super::inspector::place_again(&self.egui_ctx, id);
        }
    }

    /// What the window's variant does to the genes open here: at its own contig and position,
    /// wherever the view has gone since (computed again only when genes or the reference
    /// change, not each frame).
    pub(super) fn popup_consequences(&mut self) -> Vec<String> {
        let Some(p) = &self.variant_popup else { return Vec::new() };
        let key = (self.genes.iter().filter(|g| g.annotation().is_some()).count(), self.all_transcripts, self.reference_path.clone());
        if let Some((k, lines)) = &p.consequences
            && *k == key
        {
            return lines.clone();
        }
        let var = &p.var;
        let lines: Vec<String> = var
            .alts
            .iter()
            .filter(|a| !a.starts_with('<') && *a != "*")
            .flat_map(|a| self.consequences(&p.contig, var.start, &var.reference, a).into_iter().take(6))
            .map(|c| format!("{} {}  {}{}  {}", c.gene, c.transcript, c.hgvs_c, c.hgvs_p.map(|p| format!("  {p}")).unwrap_or_default(), c.effect))
            .collect();
        if let Some(p) = &mut self.variant_popup {
            p.consequences = Some((key, lines.clone()));
        }
        lines
    }

    pub(super) fn variant_window(&mut self, ctx: &egui::Context) {
        let consequences = self.popup_consequences();
        let Some(p) = &self.variant_popup else { return };
        let var = &p.var;
        let title = format!("{}:{}  {}", p.contig, commas(var.start + 1), var.kind());
        let text = self.variant_text(p, &consequences);
        let pane = self.pane_rect(ctx);
        // At its pane's top, on the side away from the variant clicked (the right, when it was
        // opened some other way).
        let spot = super::inspector::beside(pane, p.clicked_x.unwrap_or(pane.left()), 300.0, ctx.content_rect());
        let mut open = true;
        let mut copy = false;
        let id = self.variant_window_id();
        super::inspector::placed(ctx, egui::Window::new(title), id, spot, self.pane_settled())
            .open(&mut open)
            .default_width(460.0)
            .show(ctx, |ui| {
                ui.label(egui::RichText::new(&p.file).weak());
                let short = |s: &str| if s.chars().count() > 60 { format!("{}… ({} bp)", s.chars().take(60).collect::<String>(), s.len()) } else { s.to_owned() };
                ui.monospace(format!("{} > {}", short(&var.reference), short(&var.alts.join(","))));
                ui.label(format!(
                    "ID {}   QUAL {}   FILTER {}",
                    var.id,
                    var.qual.map_or_else(|| ".".to_owned(), |q| format!("{q:.1}")),
                    if var.filters.is_empty() { ".".to_owned() } else { var.filters.join(",") }
                ));
                let clinvar = clinvar_lines(var);
                if !clinvar.is_empty() {
                    ui.separator();
                    let sig = var.info_value("CLNSIG").map(Significance::of).unwrap_or(Significance::Other);
                    let [r, g, b] = sig.rgb();
                    ui.label(egui::RichText::new("ClinVar").strong().color(egui::Color32::from_rgb(r, g, b)));
                    egui::Grid::new("clinvar").num_columns(2).striped(true).show(ui, |ui| {
                        for (k, v) in &clinvar {
                            ui.label(egui::RichText::new(*k).weak());
                            ui.add(egui::Label::new(v).wrap());
                            ui.end_row();
                        }
                    });
                }
                if !consequences.is_empty() {
                    ui.separator();
                    ui.label(egui::RichText::new("In the genes open here").strong());
                    for c in &consequences {
                        ui.monospace(c);
                    }
                }
                if !var.samples.is_empty() {
                    ui.separator();
                    ui.label(egui::RichText::new("Genotypes").strong());
                    for (name, call) in p.samples.iter().zip(&var.samples).take(20) {
                        let mut l = format!("{name}: {}", call.gt);
                        if let Some(dp) = call.dp {
                            l.push_str(&format!("  DP {dp}"));
                        }
                        if let Some(gq) = call.gq {
                            l.push_str(&format!("  GQ {gq}"));
                        }
                        if let Some(ad) = call.ad_text() {
                            l.push_str(&format!("  AD {ad}"));
                        }
                        ui.monospace(l);
                    }
                    if var.samples.len() > 20 {
                        ui.label(egui::RichText::new(format!("and {} more samples", var.samples.len() - 20)).weak());
                    }
                }
                // A record with no INFO (`.`) has no section, rather than a heading over nothing.
                let info: Vec<(&str, &str)> = var.info_pairs().collect();
                if !info.is_empty() {
                    ui.separator();
                    ui.label(egui::RichText::new("INFO").strong());
                }
                egui::ScrollArea::vertical().max_height(320.0).show(ui, |ui| {
                    egui::Grid::new("info").num_columns(2).striped(true).show(ui, |ui| {
                        for (k, v) in info {
                            let desc = p.fields.iter().find(|f| f.id == k).map(|f| f.description.as_str()).unwrap_or("");
                            let key = ui.label(egui::RichText::new(k).monospace().strong());
                            if !desc.is_empty() {
                                key.on_hover_text(desc);
                            }
                            ui.vertical(|ui| {
                                ui.add(egui::Label::new(if v.is_empty() { "(flag)".to_owned() } else { v.to_owned() }).wrap());
                                if !desc.is_empty() {
                                    ui.label(egui::RichText::new(desc).weak().small());
                                }
                            });
                            ui.end_row();
                        }
                    });
                });
                ui.separator();
                if ui.button("Copy as text").clicked() {
                    copy = true;
                }
            });
        if copy {
            ctx.copy_text(text);
        }
        if !open {
            self.variant_popup = None;
        }
    }

    /// The window's contents as plain text, for copying.
    fn variant_text(&self, p: &VariantPopup, consequences: &[String]) -> String {
        let var = &p.var;
        let mut out = vec![format!("{}:{} {} {} > {}  ({})", p.contig, var.start + 1, var.id, super::allele(&var.reference), super::alleles(&var.alts.join(",")), p.file)];
        out.extend(clinvar_lines(var).into_iter().map(|(k, v)| format!("{k}: {v}")));
        out.extend(consequences.iter().cloned());
        out.extend(var.info_pairs().map(|(k, v)| {
            let desc = p.fields.iter().find(|f| f.id == k).map(|f| f.description.as_str()).unwrap_or("");
            if desc.is_empty() { format!("{k}={v}") } else { format!("{k}={v}  ({desc})") }
        }));
        out.join("\n")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn clinvar(info: &str) -> Variant {
        Variant {
            start: 43_094_463,
            end: 43_094_464,
            id: "55601".into(),
            reference: "C".into(),
            alts: vec!["T".into()],
            qual: None,
            filters: vec![],
            info: info.into(),
            af: None,
            samples: vec![],
            ..Default::default()
        }
    }

    #[test]
    fn clinvar_summary() {
        let v = clinvar(
            "ALLELEID=70069;CLNDISDB=MedGen:C2676676;CLNDN=Hereditary_breast_ovarian_cancer_syndrome|not_provided;CLNHGVS=NC_000017.11:g.43094464C>T;\
             CLNREVSTAT=reviewed_by_expert_panel;CLNSIG=Pathogenic;GENEINFO=BRCA1:672;MC=SO:0001587|nonsense",
        );
        let lines = clinvar_lines(&v);
        assert_eq!(lines[0], ("Significance", "Pathogenic".to_owned()));
        assert_eq!(lines[1], ("Review", "3 of 4 stars  (reviewed by expert panel)".to_owned()));
        assert_eq!(lines[2], ("Condition", "Hereditary breast ovarian cancer syndrome | not provided".to_owned()));
        assert_eq!(lines[3], ("Gene", "BRCA1".to_owned()));
        assert_eq!(lines[4], ("Consequence", "nonsense".to_owned()));
        assert!(clinvar_lines(&clinvar("DP=30")).is_empty(), "not a ClinVar record");
        assert_eq!(review_stars("criteria_provided,_single_submitter"), 1);
        assert_eq!(review_stars("criteria_provided,_conflicting_classifications"), 1);
        assert_eq!(review_stars("no_assertion_criteria_provided"), 0);
    }
}
