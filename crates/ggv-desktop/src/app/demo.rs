//! The demo (real HG002 long reads, 4.3 MB, downloaded once), the welcome shown when nothing is open,
//! and a short tour: a deletion on one haplotype, Evidence at a base, methylation, variant review,
//! then what your own files need.

use eframe::egui::{self, Align2, Rect, RichText, Vec2};
use ggv_io::demo;

use super::App;
use crate::session::Session;

/// The demo's name among downloads: finishing opens it with the tour.
pub(super) const DEMO_ID: &str = "the demo";

/// The tour's steps: a title, what to look at, and what "Show me" does.
const STEPS: &[(&str, &str)] = &[
    ("A deletion on one haplotype", ""),
    ("Evidence at one base", ""),
    ("Methylation on each haplotype", ""),
    (
        "Review variants",
        "Variant review steps through the calls: mark each real, artifact or unsure (keys 1, 2, 3; N and P step) with an optional note. Decisions are saved as you go, beside the VCF.",
    ),
    (
        "Your own data",
        "Open alignments (BAM or CRAM) with their index (.bai, .csi or .crai; GGV offers to build a missing one), and a reference: a FASTA with its .fai, or a built-in genome from Reference (it downloads once). CRAM needs the reference it was made with. Add variants (VCF or BCF, indexed .tbi or .csi), genes (GTF, GFF3 or BED) and signals (bigWig, bedGraph, SEG). Drag files onto the window, or use the toolbar.",
    ),
];

impl App {
    /// Opens the demo at the tour's first stop, with the tour; downloads it first (once, in the
    /// background, shown in the status bar) if needed.
    pub(super) fn open_demo(&mut self, ctx: &egui::Context) {
        let Some(d) = demo::cached() else {
            self.start_download(DEMO_ID.to_owned(), demo::SIZE, ctx, |progress| demo::fetch(progress).map(|d| d.files()));
            return;
        };
        let s = Session {
            version: crate::session::VERSION,
            bams: vec![d.bam.clone()],
            vcfs: d.vcfs.clone(),
            genes: vec![d.genes.clone()],
            reference: Some(d.reference.clone()),
            locus: Some(demo::STOPS[0].0.to_owned()),
            group: "haplotype".into(),
            ..Default::default()
        };
        self.apply_session(s, false, ctx);
        self.tour = Some(0);
    }

    /// Does what tour step `k` shows.
    fn tour_show(&mut self, k: usize, ctx: &egui::Context) {
        use crate::tracks::ColorBy;
        use ggv_layout::GroupBy;
        match k {
            0 => {
                self.style.group = GroupBy::Haplotype;
                self.style.color = ColorBy::None;
                self.goto(demo::STOPS[0].0);
            }
            1 => {
                self.show_evidence = true;
                self.goto(demo::STOPS[1].0);
            }
            2 => {
                self.style.group = GroupBy::Haplotype;
                self.style.color = ColorBy::BaseMods;
                self.goto(demo::STOPS[2].0);
            }
            3 => {
                self.style.color = ColorBy::None;
                if !self.variants.is_empty() {
                    self.start_review(0);
                }
            }
            _ => self.run(super::actions::Action::AddBam, ctx),
        }
    }

    pub(super) fn tour_window(&mut self, ctx: &egui::Context) {
        let Some(k) = self.tour else { return };
        let k = k.min(STEPS.len() - 1);
        let (title, text) = STEPS[k];
        let text = if text.is_empty() { demo::STOPS[k].1 } else { text };
        let mut open = true;
        let mut go = None;
        let mut show = false;
        egui::Window::new("Tour")
            .id(egui::Id::unique(("tour", self.reference_slot)))
            .anchor(Align2::RIGHT_TOP, Vec2::new(-12.0, 90.0))
            .resizable(false)
            .collapsible(true)
            .open(&mut open)
            .show(ctx, |ui| {
                ui.set_max_width(320.0);
                ui.label(RichText::new(format!("{} of {}  ·  {title}", k + 1, STEPS.len())).strong());
                ui.add_space(4.0);
                ui.label(text);
                ui.add_space(6.0);
                ui.horizontal(|ui| {
                    let label = if k + 1 == STEPS.len() { "Open my files…" } else { "Show me" };
                    if ui.button(label).clicked() {
                        show = true;
                    }
                    if ui.add_enabled(k > 0, egui::Button::new("< Back")).clicked() {
                        go = Some(k - 1);
                    }
                    if ui.add_enabled(k + 1 < STEPS.len(), egui::Button::new("Next >")).clicked() {
                        go = Some(k + 1);
                    }
                });
                ui.label(RichText::new("The data: PacBio Revio HiFi reads of GIAB HG002 on GRCh38 chr20, public.").weak().small());
            });
        if show {
            self.tour_show(k, ctx);
        }
        if let Some(n) = go {
            self.tour = Some(n);
            self.tour_show(n, ctx);
        }
        if !open {
            self.tour = None;
        }
    }

    /// With nothing open: what GGV is, the demo, and what your own files need.
    pub(super) fn welcome(&mut self, ui: &mut egui::Ui, rect: Rect) {
        if !self.tracks.is_empty() || !self.genes.is_empty() || !self.variants.is_empty() || !self.signals.is_empty() || self.pending_genome.is_some() {
            return;
        }
        let ctx = ui.ctx().clone();
        let area = Rect::from_center_size(rect.center(), Vec2::new(rect.width().min(560.0), rect.height().min(360.0)));
        let mut child = ui.new_child(egui::UiBuilder::new().max_rect(area));
        child.vertical_centered(|ui| {
            ui.heading("GGV");
            ui.label("A fast genome viewer for people and agents.");
            ui.add_space(14.0);
            let downloading = self.downloads.iter().any(|d| d.id == DEMO_ID);
            if ui
                .add_enabled(!downloading, egui::Button::new(RichText::new(if downloading { "Downloading the demo…" } else { "Open the demo" }).size(16.0)))
                .on_hover_text("Real PacBio HiFi reads of GIAB HG002 with a short tour (4.3 MB, downloaded once and kept)")
                .clicked()
            {
                self.open_demo(&ctx);
            }
            ui.add_space(10.0);
            ui.horizontal(|ui| {
                let w = 2.0 * 150.0 + ui.spacing().item_spacing.x;
                ui.add_space(((ui.available_width() - w) / 2.0).max(0.0));
                if ui.add_sized([150.0, 22.0], egui::Button::new("Open alignments…")).clicked() {
                    self.run(super::actions::Action::AddBam, &ctx);
                }
                if ui.add_sized([150.0, 22.0], egui::Button::new("Open session…")).clicked() {
                    self.run(super::actions::Action::OpenSession, &ctx);
                }
            });
            ui.add_space(14.0);
            ui.label(RichText::new(STEPS[STEPS.len() - 1].1).weak());
            ui.add_space(6.0);
            ui.label(RichText::new("Ctrl+K: every command  ·  ?: keyboard shortcuts").weak().small());
        });
    }
}
