//! Reference and assembly checks (L03): which assembly each open file is, how sure GGV is of
//! that, and where the files disagree, so an hg19 BAM against hg38 annotations or a VCF called
//! on another reference is explained before anyone trusts its mismatch colours or consequence
//! labels.
//!
//! A build is *declared* when a file names it (a VCF's `##reference`, a file name such as
//! `hg38.ncbiRefSeq.gtf.gz`), *inferred* when only chromosome 1's length matches a known build,
//! and *unknown* otherwise. An inferred build is a hint, not a verification: two assemblies can
//! share chr1's length, so GGV also compares every contig the files share and checks VCF REF
//! alleles against the reference. Custom assemblies are fine: unknown is not an error. The
//! identities and contig comparisons are `ggv_io::assembly`, which `ggv mcp` uses too.

use std::{collections::BTreeMap, path::PathBuf};

use eframe::egui::{self, Color32, RichText};
use ggv_io::assembly::{How, Identity, canonical_contig as canonical, conflicts};
use serde_json::{Value, json};

use super::App;

/// REF alleles checked against the reference for one VCF.
#[derive(Clone, Debug, Default)]
pub struct RefCheck {
    /// What was last checked (so each load is counted once).
    seen: std::collections::HashSet<(String, u64)>,
    /// The loaded variants and reference bases last checked (their addresses): the check runs
    /// when either changes, not every frame (every frame over ClinVar cost 100 ms).
    last: Option<(usize, usize)>,
    pub checked: usize,
    pub mismatched: usize,
    /// The first disagreement: locus, REF, reference.
    pub example: Option<(String, String, String)>,
}

/// The checks for the open files.
#[derive(Default)]
pub struct Assembly {
    /// The files (and reference) these were computed for.
    key: Vec<PathBuf>,
    pub identities: Vec<Identity>,
    /// Disagreements between the files, as sentences.
    pub conflicts: Vec<String>,
    pub ref_checks: BTreeMap<PathBuf, RefCheck>,
    /// The reference `ref_checks` were made against.
    checked_against: Option<PathBuf>,
    pub show: bool,
}

impl Assembly {
    /// Forgets the REF checks when the reference changes: they said whether the VCFs match
    /// the old one.
    fn checking_against(&mut self, reference: Option<&std::path::Path>) {
        if self.checked_against.as_deref() != reference {
            self.ref_checks.clear();
            self.checked_against = reference.map(std::path::Path::to_path_buf);
        }
    }
}

use crate::util::commas;

impl App {
    /// Every open file's assembly.
    fn identities(&self) -> Vec<Identity> {
        let lengths = |c: &[ggv_model::Contig]| c.iter().map(|c| (c.name.clone(), c.length)).collect::<Vec<_>>();
        let mut ids = Vec::new();
        let mut push = |kind, path: PathBuf, declared_text: Option<String>, lengths: Vec<(String, u64)>| ids.push(Identity::new(kind, path, declared_text, lengths));
        if let Some(r) = &self.reference_path {
            let l = ggv_io::ReferenceFasta::open(r).map(|f| lengths(f.contigs())).unwrap_or_default();
            push("reference", r.clone(), Some(crate::files::display_name(r)), l);
        }
        for t in &self.tracks {
            push("alignments", t.loader.bam_path.clone(), None, lengths(t.contigs()));
        }
        for v in &self.variants {
            let text = v.declared_assembly.clone().or_else(|| Some(crate::files::display_name(&v.path)));
            push("variants", v.path.clone(), text, lengths(&v.header_contigs));
        }
        for g in &self.genes {
            push("annotation", g.path.clone(), Some(crate::files::display_name(&g.path)), Vec::new());
        }
        for s in &self.signals {
            let bigwig = s.path.extension().is_some_and(|e| e.eq_ignore_ascii_case("bw") || e.eq_ignore_ascii_case("bigwig"));
            let l = if bigwig { lengths(&s.contigs()) } else { Vec::new() };
            push("signal", s.path.clone(), Some(crate::files::display_name(&s.path)), l);
        }
        ids
    }

    /// Recomputes the checks when the open files change, and checks newly loaded VCF REF
    /// alleles against the reference.
    pub(super) fn update_assembly(&mut self) {
        self.assembly.checking_against(self.reference_path.as_deref());
        let mut key: Vec<PathBuf> = self.reference_path.iter().cloned().collect();
        key.extend(self.tracks.iter().map(|t| t.loader.bam_path.clone()));
        key.extend(self.variants.iter().map(|v| v.path.clone()));
        key.extend(self.genes.iter().map(|g| g.path.clone()));
        key.extend(self.signals.iter().map(|s| s.path.clone()));
        if key != self.assembly.key {
            let ids = self.identities();
            self.assembly.conflicts = conflicts(&ids);
            self.assembly.identities = ids;
            self.assembly.ref_checks.retain(|p, _| key.contains(p));
            self.assembly.key = key;
        }
        self.check_refs();
    }

    /// Compares the REF of each loaded variant with the reference bases under it.
    fn check_refs(&mut self) {
        let Some((view, _, _)) = &self.reference else { return };
        let Some(bases) = view.reference.as_ref() else { return };
        let region = &view.region;
        for v in &self.variants {
            let Some(loaded) = &v.current else { continue };
            if canonical(&loaded.region.contig) != canonical(&region.contig) {
                continue;
            }
            let check = self.assembly.ref_checks.entry(v.path.clone()).or_default();
            let now = (std::sync::Arc::as_ptr(&loaded.variants) as usize, std::sync::Arc::as_ptr(view) as usize);
            if check.last == Some(now) {
                continue;
            }
            check.last = Some(now);
            // Variants are in start order: only those that can lie inside the reference bases.
            let from = loaded.variants.partition_point(|x| x.start < region.start);
            let to = loaded.variants.partition_point(|x| x.start < region.end);
            for var in &loaded.variants[from..to] {
                let r = var.reference.as_bytes();
                if r.is_empty() || !r.iter().all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')) {
                    continue;
                }
                let (s, e) = (var.start, var.start + r.len() as u64);
                if s < region.start || e > region.end {
                    continue;
                }
                let site = (region.contig.clone(), var.start);
                if !check.seen.insert(site) {
                    continue;
                }
                let have = &bases[(s - region.start) as usize..(e - region.start) as usize];
                check.checked += 1;
                if !have.eq_ignore_ascii_case(r) {
                    check.mismatched += 1;
                    if check.example.is_none() {
                        check.example = Some((format!("{}:{}", region.contig, commas(s + 1)), var.reference.clone(), String::from_utf8_lossy(have).into_owned()));
                    }
                }
            }
            // Keep the memory small: only the sites needed to avoid recounting nearby loads.
            if check.seen.len() > 200_000 {
                check.seen.clear();
            }
        }
    }

    /// VCFs whose REF alleles disagree with the reference often enough to mean another reference.
    fn ref_conflicts(&self) -> Vec<String> {
        self.assembly
            .ref_checks
            .iter()
            .filter(|(_, c)| c.mismatched >= 3 && c.mismatched * 5 >= c.checked)
            .map(|(p, c)| {
                let eg = c.example.as_ref().map(|(at, r, have)| format!(" (at {at} the VCF says {r}, the reference has {have})")).unwrap_or_default();
                format!(
                    "{} of {} REF alleles checked in {} don't match the reference{eg}: it was probably called on another reference.",
                    c.mismatched,
                    c.checked,
                    crate::files::display_name(p)
                )
            })
            .collect()
    }

    /// Every disagreement now.
    pub(super) fn assembly_conflicts(&self) -> Vec<String> {
        let mut c = self.assembly.conflicts.clone();
        c.extend(self.ref_conflicts());
        c
    }

    /// The open genome's build and how it is known, for the status bar.
    fn assembly_summary(&self) -> (String, String) {
        let ids = &self.assembly.identities;
        let best = ids.iter().find(|i| i.how == How::Declared && i.build.is_some()).or_else(|| ids.iter().find(|i| i.build.is_some()));
        match best {
            Some(i) => {
                let b = i.build.unwrap_or_default();
                let name = ggv_io::genomes::find(b).map_or(b, |g| g.name);
                (name.to_owned(), if i.how == How::Declared { "declared".into() } else { "inferred".into() })
            }
            None if ids.is_empty() => (String::new(), String::new()),
            None => ("custom assembly".into(), "unknown".into()),
        }
    }

    /// The status bar's assembly label, or a warning when the files disagree.
    pub(super) fn assembly_status(&mut self, ui: &mut egui::Ui) {
        let (name, how) = self.assembly_summary();
        if name.is_empty() {
            return;
        }
        let conflicts = self.assembly_conflicts();
        let r = if conflicts.is_empty() {
            ui.add(egui::Button::new(RichText::new(format!("{name} · {how}"))).frame(false))
                .on_hover_text("The open files' assembly. Click for each file's build and how GGV knows it.")
        } else {
            ui.add(egui::Button::new(RichText::new(format!("Assembly mismatch ({})", conflicts.len())).color(Color32::WHITE).strong()).fill(Color32::from_rgb(190, 70, 20)))
                .on_hover_text(conflicts.join("\n"))
        };
        if r.clicked() {
            self.assembly.show = !self.assembly.show;
        }
        ui.separator();
    }

    /// Each file's assembly and the disagreements, explained.
    pub(super) fn assembly_window(&mut self, ctx: &egui::Context) {
        if !self.is_active || !self.assembly.show {
            return;
        }
        let mut open = true;
        let conflicts = self.assembly_conflicts();
        egui::Window::new("Assembly and reference").id(egui::Id::unique(("assembly", self.reference_slot))).open(&mut open).default_width(560.0).show(ctx, |ui| {
            if conflicts.is_empty() {
                ui.label("No disagreements found between the open files.");
                ui.label(
                    RichText::new("GGV compares the contig lengths the files state and checks VCF REF alleles in view against the reference. A build inferred from chr1's length is a hint, not a proof.")
                        .weak()
                        .small(),
                );
            } else {
                ui.colored_label(Color32::from_rgb(190, 70, 20), RichText::new("These files don't agree on the assembly").strong());
                for c in &conflicts {
                    ui.label(format!("• {c}"));
                }
                ui.label(
                    RichText::new(
                        "Until they match, mismatch colours, variant positions, gene models and consequence labels may be wrong. Open files from one assembly (or lift them over), or set the reference the data was aligned to.",
                    )
                    .weak(),
                );
            }
            ui.separator();
            egui::Grid::new("assembly-grid").striped(true).num_columns(4).show(ui, |ui| {
                for h in ["File", "Kind", "Build", "How GGV knows"] {
                    ui.label(RichText::new(h).strong());
                }
                ui.end_row();
                for i in &self.assembly.identities {
                    ui.label(crate::files::display_name(&i.path)).on_hover_text(i.path.display().to_string());
                    ui.label(i.kind);
                    ui.label(i.build.map_or("—".to_owned(), |b| ggv_io::genomes::find(b).map_or(b, |g| g.name).to_owned()));
                    let mut how = i.how.text().to_owned();
                    if let (How::Declared, Some(t)) = (i.how, &i.declared_text) {
                        how = format!("declared ({t})");
                    }
                    if !i.lengths.is_empty() {
                        let n = i.lengths.len();
                        how.push_str(&format!(", {n} contig length{}", if n == 1 { "" } else { "s" }));
                    }
                    ui.label(how);
                    ui.end_row();
                }
            });
            for (p, c) in &self.assembly.ref_checks {
                if c.checked > 0 {
                    ui.label(RichText::new(format!("{}: {} of {} REF alleles in view match the reference", crate::files::display_name(p), c.checked - c.mismatched, c.checked)).small());
                }
            }
        });
        if !open {
            self.assembly.show = false;
        }
    }

    /// The assembly checks for agents.
    pub(super) fn assembly_state(&self) -> Value {
        let (name, how) = self.assembly_summary();
        json!({
            "build": (!name.is_empty()).then_some(name),
            "how": (!how.is_empty()).then_some(how),
            "conflicts": self.assembly_conflicts(),
            "files": self.assembly.identities.iter().map(|i| json!({
                "path": i.path.display().to_string(),
                "kind": i.kind,
                "build": i.build,
                "how": match i.how { How::Declared => "declared", How::Inferred => "inferred", How::Unknown => "unknown" },
            })).collect::<Vec<_>>(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ref_checks_are_redone_for_a_new_reference() {
        let mut a = Assembly::default();
        let (good, bad) = (std::path::Path::new("/d/good.fa"), std::path::Path::new("/d/bad.fa"));
        a.checking_against(Some(good));
        let mismatched = RefCheck { checked: 4, mismatched: 4, ..Default::default() };
        a.ref_checks.insert("/d/v.vcf.gz".into(), RefCheck { checked: 4, ..Default::default() });
        // The same reference: kept.
        a.checking_against(Some(good));
        assert_eq!(a.ref_checks.len(), 1);
        // Swapped either way: checked again from nothing.
        a.checking_against(Some(bad));
        assert!(a.ref_checks.is_empty());
        a.ref_checks.insert("/d/v.vcf.gz".into(), mismatched);
        a.checking_against(Some(good));
        assert!(a.ref_checks.is_empty());
        a.ref_checks.insert("/d/v.vcf.gz".into(), RefCheck::default());
        a.checking_against(None);
        assert!(a.ref_checks.is_empty());
    }
}
