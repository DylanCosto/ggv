//! Workspace recovery in a pane: files a session names that have moved (Locate them, or a
//! folder to look in), and the pane's place for saving a split window.

use std::path::{Path, PathBuf};

use eframe::egui::{self, RichText};

use super::App;
use crate::session::{PaneView, Session};

/// A session opened with files missing: what it was, and what's still missing. Files opened
/// on their own (the command line, a dialog, a drop, Recent files) that aren't there are listed
/// the same way, in a session of just those files.
pub struct Missing {
    pub session: Session,
    pub paths: Vec<PathBuf>,
    /// Files opened on their own rather than a session: found ones open beside what's open.
    pub opened: bool,
}

/// What a file was opened as, for listing it when it isn't there.
#[derive(Clone, Copy)]
pub(super) enum OpenedAs {
    Alignments,
    Variants,
    Genes,
    Signal,
    Reference,
}

impl App {
    /// This pane's view and outlined base, for saving a split window.
    /// (Called every frame: just these two, not a whole session, which clones every path,
    /// bookmark and overlay.)
    pub fn pane_view(&self) -> PaneView {
        let Some(c) = self.contigs.get(self.contig) else { return PaneView::default() };
        let (vs, ve) = self.visible();
        PaneView {
            locus: Some(format!("{}:{}-{}", c.name, vs.max(0.0) as u64 + 1, ve.min(self.contig_len()) as u64)),
            highlight: self.highlight.and_then(|(ci, pos)| self.contigs.get(ci).map(|c| format!("{}:{}", c.name, pos + 1))),
        }
    }

    /// A session loaded with several panes asks the window to rebuild them.
    pub fn take_pending_panes(&mut self) -> Option<Session> {
        self.pending_panes.take()
    }

    /// Goes to `s`'s place again once this pane has its final width (a place fitted to the
    /// whole window before a split showed half of it).
    pub fn refit(&mut self, s: &Session) {
        if let Some(l) = &s.locus {
            self.pending_goto = Some(l.clone());
            self.pending_highlight = s.highlight.clone();
        }
    }

    /// Opens `s` without its missing files, and lists them to locate.
    pub(super) fn note_missing(&mut self, s: &mut Session) {
        let paths = s.missing();
        if paths.is_empty() {
            self.missing = None;
            return;
        }
        let full = s.clone();
        let gone = |p: &PathBuf| !paths.contains(p);
        s.bams.retain(gone);
        s.vcfs.retain(gone);
        s.genes.retain(gone);
        s.signals.retain(gone);
        if s.reference.as_ref().is_some_and(|r| paths.contains(r)) {
            s.reference = None;
        }
        self.missing = Some(Missing { session: full, paths, opened: false });
    }

    /// `paths` without the files that aren't there, which are listed in Files not found (as a
    /// session's are) rather than opened: a missing VCF or BAM was offered an index build.
    pub(super) fn existing(&mut self, paths: Vec<PathBuf>, kind: OpenedAs) -> Vec<PathBuf> {
        let (here, gone): (Vec<PathBuf>, Vec<PathBuf>) = paths.into_iter().partition(|p| p.exists());
        if gone.is_empty() {
            return here;
        }
        let m = self.missing.get_or_insert_with(|| Missing { session: Session::default(), paths: Vec::new(), opened: true });
        for p in gone {
            if m.paths.contains(&p) {
                continue;
            }
            match kind {
                OpenedAs::Alignments => m.session.bams.push(p.clone()),
                OpenedAs::Variants => m.session.vcfs.push(p.clone()),
                OpenedAs::Genes => m.session.genes.push(p.clone()),
                OpenedAs::Signal => m.session.signals.push(p.clone()),
                OpenedAs::Reference => m.session.reference = Some(p.clone()),
            }
            m.paths.push(p);
        }
        here
    }

    /// `old` is at `new`: so may the other missing files be, in `new`'s folder. Reopens the
    /// session with every file found so far.
    pub(super) fn located(&mut self, old: &Path, new: &Path, ctx: &egui::Context) {
        let Some(m) = &mut self.missing else { return };
        m.session.relocate(old, new);
        if let Some(dir) = new.parent() {
            self.look_in(dir, ctx);
        } else {
            self.reapply_missing(ctx);
        }
    }

    /// Looks for every missing file, by name, in `dir` (and one folder down).
    pub(super) fn look_in(&mut self, dir: &Path, ctx: &egui::Context) {
        let Some(m) = &mut self.missing else { return };
        let subdirs: Vec<PathBuf> = std::fs::read_dir(dir)
            .map(|rd| rd.filter_map(Result::ok).map(|e| e.path()).filter(|p| p.is_dir()).collect())
            .unwrap_or_default();
        for old in m.session.missing() {
            let Some(name) = old.file_name() else { continue };
            if let Some(found) = std::iter::once(dir.to_path_buf()).chain(subdirs.iter().cloned()).map(|d| d.join(name)).find(|p| p.exists()) {
                m.session.relocate(&old, &found);
            }
        }
        self.reapply_missing(ctx);
    }

    fn reapply_missing(&mut self, ctx: &egui::Context) {
        let Some(m) = self.missing.take() else { return };
        let before = m.paths.len();
        if m.opened {
            self.open_found(m.session, ctx);
        } else {
            self.reapply_session(m.session, ctx);
        }
        let left = self.missing.as_ref().map_or(0, |m| m.paths.len());
        if left < before {
            self.notify(format!("Found {} file{}", before - left, if before - left == 1 { "" } else { "s" }));
        }
    }

    /// Opens the files of `s` (files opened on their own) found so far beside what's open; the
    /// others stay listed.
    fn open_found(&mut self, s: Session, ctx: &egui::Context) {
        let here = |v: &[PathBuf]| v.iter().filter(|p| p.exists()).cloned().collect::<Vec<_>>();
        let (bams, vcfs, genes, signals) = (here(&s.bams), here(&s.vcfs), here(&s.genes), here(&s.signals));
        let reference = s.reference.clone().filter(|r| r.exists());
        let mut rest = s;
        let gone = |p: &PathBuf| !p.exists();
        rest.bams.retain(gone);
        rest.vcfs.retain(gone);
        rest.genes.retain(gone);
        rest.signals.retain(gone);
        rest.reference = rest.reference.filter(|r| !r.exists());
        let paths = rest.missing();
        if !paths.is_empty() {
            self.missing = Some(Missing { session: rest, paths, opened: true });
        }
        if let Some(r) = reference {
            self.set_reference(r, ctx);
        }
        self.add_tracks(bams, ctx);
        self.add_variants(vcfs, ctx);
        self.add_genes(genes, ctx);
        self.add_signals(signals, ctx);
    }

    fn reapply_session(&mut self, mut s: Session, ctx: &egui::Context) {
        // Keep the view where it is now.
        let now = self.session();
        s.locus = now.locus.or(s.locus);
        s.highlight = now.highlight.or(s.highlight);
        s.panes.clear();
        self.apply_session(s, false, ctx);
    }

    /// Keeps review and queue decisions in `dir` from now on, moving the open review there.
    pub(super) fn set_review_dir(&mut self, dir: Option<PathBuf>) {
        crate::store::set_output_dir(dir.as_deref());
        if let Some(r) = &mut self.review {
            r.move_to_preferred();
        }
        if let Some(d) = dir {
            self.notify(format!("Review decisions are saved in {}", d.display()));
        }
    }

    /// The window listing files a session couldn't find.
    pub(super) fn missing_window(&mut self, ctx: &egui::Context) {
        if !self.is_active {
            return;
        }
        let Some(m) = &self.missing else { return };
        let paths = m.paths.clone();
        let opened = m.opened;
        let mut locate = None;
        let mut folder = false;
        let mut close = false;
        egui::Window::new("Files not found").id(egui::Id::unique(("missing", self.reference_slot))).collapsible(false).show(ctx, |ui| {
            ui.label(if opened {
                "These files aren't there (moved, renamed or deleted?). Locate one and GGV looks for the others in the same folder."
            } else {
                "The session names files that aren't where it says. Locate one and GGV looks for the others in the same folder."
            });
            ui.add_space(4.0);
            for p in &paths {
                ui.horizontal(|ui| {
                    if ui.button("Locate…").clicked() {
                        locate = Some(p.clone());
                    }
                    ui.label(RichText::new(crate::files::display_name(p)).strong()).on_hover_text(p.display().to_string());
                    ui.label(RichText::new(p.parent().map(|d| d.display().to_string()).unwrap_or_default()).weak().small());
                });
            }
            ui.add_space(4.0);
            ui.horizontal(|ui| {
                if ui.button("Look in a folder…").on_hover_text("Find every missing file by name in a folder (and the folders inside it)").clicked() {
                    folder = true;
                }
                if ui.button("Continue without them").clicked() {
                    close = true;
                }
            });
        });
        if let Some(p) = locate {
            let mut d = rfd::FileDialog::new().set_title(format!("Where is {}?", crate::files::display_name(&p)));
            if let Some(name) = p.file_name() {
                d = d.set_file_name(name.to_string_lossy());
            }
            self.ask(d, super::dialogs::Kind::PickFile, super::dialogs::Then::Locate(p));
        }
        if folder {
            self.ask(rfd::FileDialog::new().set_title("Folder with the missing files"), super::dialogs::Kind::PickFolder, super::dialogs::Then::LookIn);
        }
        if close {
            self.missing = None;
        }
    }
}
