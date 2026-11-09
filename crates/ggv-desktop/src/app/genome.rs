//! Built-in genomes and tracks in the app: the Reference menu lists the genomes, the Tracks
//! menu the gene annotations, ClinVar and pathogenicity scores for the open genome's build.
//! Choosing one (or `--genome`) opens its files, downloading them once in the background.

use std::{
    path::PathBuf,
    sync::mpsc::{self, Receiver},
};

use eframe::egui;
use ggv_io::catalog::{self, Group, Track};

use super::App;

/// A download in progress.
pub struct Download {
    /// `hg38`, or `hg38/clinvar` for a track.
    pub id: String,
    pub progress: (String, u64, u64),
    rx: Receiver<Msg>,
}

enum Msg {
    Progress(String, u64, u64),
    /// The files to open.
    Done(Result<Vec<PathBuf>, String>),
}

impl App {
    /// Opens a built-in genome's reference and genes; downloads them first if needed.
    pub fn use_genome(&mut self, name: &str, ctx: &egui::Context) {
        let Some(g) = ggv_io::genomes::find(name) else {
            let ids: Vec<&str> = ggv_io::genomes::GENOMES.iter().map(|g| g.id).collect();
            self.error = Some(format!("Unknown genome `{name}`; built in: {}", ids.join(", ")));
            return;
        };
        if let Some(l) = ggv_io::genomes::cached(g) {
            self.open_paths(vec![l.reference, l.genes], ctx);
            return;
        }
        self.start_download(g.id.to_owned(), g.size, ctx, move |progress| {
            ggv_io::genomes::fetch(g, progress).map(|l| vec![l.reference, l.genes])
        });
    }

    /// Opens a built-in track; downloads it first if needed.
    pub fn use_track(&mut self, t: &'static Track, ctx: &egui::Context) {
        if let Some(paths) = catalog::cached(t) {
            self.open_paths(paths, ctx);
            return;
        }
        self.start_download(format!("{}/{}", t.build, t.id), t.size, ctx, move |progress| catalog::fetch(t, progress));
    }

    /// Opens the `--track` names once the genome is known.
    pub(super) fn open_pending_tracks(&mut self, ctx: &egui::Context) {
        if self.pending_tracks.is_empty() || self.contigs.is_empty() {
            return;
        }
        let Some(build) = self.build() else {
            self.error = Some(format!("--track {}: the genome isn't a built-in build", self.pending_tracks.join(", ")));
            self.pending_tracks.clear();
            return;
        };
        for name in std::mem::take(&mut self.pending_tracks) {
            match catalog::find(build, &name) {
                Some(t) => self.use_track(t, ctx),
                None => {
                    let ids: Vec<&str> = catalog::for_build(build).map(|t| t.id).collect();
                    self.error = Some(format!("No built-in track `{name}` for {build}; there are {}", ids.join(", ")));
                }
            }
        }
    }

    pub(super) fn start_download(
        &mut self,
        id: String,
        size: u64,
        ctx: &egui::Context,
        fetch: impl FnOnce(ggv_io::genomes::Progress) -> std::io::Result<Vec<PathBuf>> + Send + 'static,
    ) {
        if self.downloads.iter().any(|d| d.id == id) {
            return;
        }
        let (tx, rx) = mpsc::channel();
        let ctx2 = ctx.clone();
        let spawned = std::thread::Builder::new().name("ggv-download".into()).spawn(move || {
            let mut last = std::time::Instant::now();
            let result = fetch(&mut |file, done, total| {
                if last.elapsed().as_millis() >= 200 {
                    last = std::time::Instant::now();
                    let _ = tx.send(Msg::Progress(file.to_owned(), done, total));
                    ctx2.request_repaint();
                }
            });
            let _ = tx.send(Msg::Done(result.map_err(|e| e.to_string())));
            ctx2.request_repaint();
        });
        match spawned {
            Ok(_) => self.downloads.push(Download { id, progress: (String::new(), 0, size), rx }),
            Err(e) => self.error = Some(format!("Couldn't start the download: {e}")),
        }
    }

    /// Takes download news; opens the files of finished ones.
    pub(super) fn poll_downloads(&mut self, ctx: &egui::Context) {
        let mut finished = Vec::new();
        self.downloads.retain_mut(|job| {
            let mut done = None;
            for m in job.rx.try_iter() {
                match m {
                    Msg::Progress(f, d, t) => job.progress = (f, d, t),
                    Msg::Done(r) => done = Some(r),
                }
            }
            match done {
                Some(r) => {
                    finished.push((job.id.clone(), r));
                    false
                }
                None => true,
            }
        });
        for (id, result) in finished {
            match result {
                // The demo opens as a session, with its tour.
                Ok(_) if id == super::demo::DEMO_ID => {
                    self.open_demo(ctx);
                    self.notify("The demo is ready (kept for next time)");
                }
                Ok(paths) => {
                    self.open_paths(paths, ctx);
                    self.notify(format!("{id} is ready (kept for next time)"));
                }
                Err(e) => self.error = Some(format!("Downloading {id}: {e} (choose it again to resume)")),
            }
        }
        if !self.downloads.is_empty() {
            ctx.request_repaint_after(std::time::Duration::from_millis(250));
        }
    }

    /// The downloads' lines in the status bar.
    pub(super) fn download_status(&self, ui: &mut egui::Ui) {
        for job in &self.downloads {
            let (file, done, total) = &job.progress;
            ui.spinner();
            let text = if *total > 0 {
                format!("Downloading {} {file}: {:.0} / {:.0} MB", job.id, *done as f64 / 1e6, *total as f64 / 1e6)
            } else {
                format!("Downloading {} {file}: {:.0} MB", job.id, *done as f64 / 1e6)
            };
            ui.label(text);
            ui.separator();
        }
    }

    fn downloading(&self, id: &str) -> bool {
        self.downloads.iter().any(|d| d.id == id)
    }

    /// The Reference menu: a FASTA or .2bit from disk, or a built-in genome.
    pub(super) fn reference_menu(&mut self, ui: &mut egui::Ui, ctx: &egui::Context) {
        ui.menu_button("Reference…", |ui| {
            if ui.button("Open FASTA or .2bit…").clicked() {
                ui.close();
                self.ask(rfd::FileDialog::new().add_filter("Reference", &["fa", "fasta", "fna", "gz", "bgz", "2bit"]), super::dialogs::Kind::PickFile, super::dialogs::Then::Reference);
            }
            ui.separator();
            ui.label(egui::RichText::new("Built-in genomes (reference + RefSeq genes)").weak());
            for g in ggv_io::genomes::GENOMES {
                let state = if ggv_io::genomes::cached(g).is_some() {
                    "ready".to_owned()
                } else if self.downloading(g.id) {
                    "downloading…".to_owned()
                } else {
                    size(g.size)
                };
                if ui.button(format!("{}   ({state})", g.name)).clicked() {
                    ui.close();
                    self.use_genome(g.id, ctx);
                }
            }
        });
    }

    /// The open genome's build, from its first chromosome's length.
    pub(super) fn build(&self) -> Option<&'static str> {
        catalog::build_of(self.contigs.iter().map(|c| (c.name.as_str(), c.length)))
    }

    /// The Tracks menu: the catalogue for the open genome's build (or every build, by name,
    /// when it isn't one GGV knows).
    pub(super) fn tracks_menu(&mut self, ui: &mut egui::Ui, ctx: &egui::Context) {
        ui.menu_button("Tracks…", |ui| match self.build() {
            Some(build) => {
                let name = ggv_io::genomes::find(build).map_or(build, |g| g.name);
                ui.label(egui::RichText::new(format!("Built-in tracks for {name}")).weak());
                self.catalog_buttons(ui, build, ctx);
            }
            None => {
                let what = if self.contigs.is_empty() { "No genome open" } else { "This genome isn't a built-in build" };
                ui.label(egui::RichText::new(format!("{what}: pick the build")).weak());
                for g in ggv_io::genomes::GENOMES {
                    ui.menu_button(g.name, |ui| self.catalog_buttons(ui, g.id, ctx));
                }
            }
        });
    }

    fn catalog_buttons(&mut self, ui: &mut egui::Ui, build: &str, ctx: &egui::Context) {
        let mut group = None;
        for t in catalog::for_build(build) {
            if group != Some(t.group) {
                if group.is_some() {
                    ui.separator();
                }
                ui.label(egui::RichText::new(t.group.label()).strong());
                group = Some(t.group);
            }
            let state = if catalog::cached(t).is_some() {
                "ready".to_owned()
            } else if self.downloading(&format!("{}/{}", t.build, t.id)) {
                "downloading…".to_owned()
            } else {
                size(t.size)
            };
            let tip = match t.group {
                Group::Scores => format!("{}. Opens as four signal tracks, A, C, G and T.", t.about),
                _ => t.about.to_owned(),
            };
            if ui.button(format!("{}   ({state})", t.name)).on_hover_text(tip).clicked() {
                ui.close();
                self.use_track(t, ctx);
            }
        }
    }
}

/// A download size for a menu: `125 MB download`, `1.9 GB download`.
fn size(bytes: u64) -> String {
    if bytes >= 1_000_000_000 {
        format!("{:.1} GB download", bytes as f64 / 1e9)
    } else {
        format!("{} MB download", (bytes / 1_000_000).max(1))
    }
}
