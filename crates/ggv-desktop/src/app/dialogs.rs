//! File dialogs off the UI thread. A native dialog run on the UI thread stops the event loop
//! until it closes, and on Windows the window then took seconds to come back, even after
//! Cancel. Here the dialog runs on its own thread; the window keeps drawing, and what was
//! chosen is acted on in the next frame. One dialog at a time.

use std::path::PathBuf;
use std::sync::mpsc::{self, Receiver};

use eframe::egui;

use super::App;

/// What to do with the chosen path(s).
pub(super) enum Then {
    AddTracks,
    AddGenes,
    AddVariants,
    Reference,
    OpenSession,
    SaveSession,
    ExportPng,
    Report,
    OpenQueue,
    /// Write this text to the chosen file.
    SaveText(String),
    /// A session's missing file is the chosen one.
    Locate(PathBuf),
    /// Look for a session's missing files in the chosen folder.
    LookIn,
    /// Keep review and queue decisions in the chosen folder.
    ReviewDir,
}

/// What kind of dialog.
#[derive(Clone, Copy)]
pub(super) enum Kind {
    PickFiles,
    PickFile,
    PickFolder,
    Save,
}

pub(super) struct Pending {
    rx: Receiver<Vec<PathBuf>>,
    then: Then,
}

impl App {
    /// Shows `dialog` on its own thread; `then` happens once something is chosen (nothing does
    /// on Cancel).
    pub(super) fn ask(&mut self, dialog: rfd::FileDialog, kind: Kind, then: Then) {
        if self.dialog.is_some() {
            return;
        }
        let (tx, rx) = mpsc::channel();
        let ctx = self.egui_ctx.clone();
        let spawned = std::thread::Builder::new().name("ggv-dialog".into()).spawn(move || {
            let chosen = match kind {
                Kind::PickFiles => dialog.pick_files().unwrap_or_default(),
                Kind::PickFile => dialog.pick_file().into_iter().collect(),
                Kind::PickFolder => dialog.pick_folder().into_iter().collect(),
                Kind::Save => dialog.save_file().into_iter().collect(),
            };
            let _ = tx.send(chosen);
            ctx.request_repaint();
        });
        match spawned {
            Ok(_) => self.dialog = Some(Pending { rx, then }),
            Err(e) => self.error = Some(format!("couldn't open the file dialog: {e}")),
        }
    }

    /// Acts on a dialog that has closed.
    pub(super) fn poll_dialog(&mut self, ctx: &egui::Context) {
        let Some(p) = &self.dialog else { return };
        let chosen = match p.rx.try_recv() {
            Ok(c) => c,
            Err(mpsc::TryRecvError::Empty) => return,
            Err(mpsc::TryRecvError::Disconnected) => Vec::new(),
        };
        let then = self.dialog.take().expect("checked").then;
        let Some(first) = chosen.first().cloned() else { return };
        match then {
            Then::AddTracks => self.add_tracks(chosen, ctx),
            Then::AddGenes => self.add_genes(chosen, ctx),
            Then::AddVariants => self.add_variants(chosen, ctx),
            Then::Reference => self.set_reference(first, ctx),
            Then::OpenSession => {
                let _ = self.open_session(&first, ctx);
            }
            Then::SaveSession => {
                if let Err(e) = self.session().with_workspace().save(&first) {
                    self.error = Some(format!("saving session: {e:#}"));
                }
            }
            Then::ExportPng => self.capture(ctx, super::tools::Capture::Png(first)),
            Then::Report => self.report_to(first),
            Then::OpenQueue => self.open_queue_file(&first),
            Then::Locate(old) => self.located(&old, &first, ctx),
            Then::LookIn => self.look_in(&first, ctx),
            Then::ReviewDir => self.set_review_dir(Some(first)),
            Then::SaveText(text) => match std::fs::write(&first, text) {
                Ok(()) => self.notify(format!("Saved {}", first.display())),
                Err(e) => self.error = Some(super::path_error(&first, e)),
            },
        }
    }
}
