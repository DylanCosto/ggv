//! About GGV and Copy diagnostics: the exact build, the system and graphics, recent errors, and
//! what's open (file names; full paths only when you choose to include them), so a bug report
//! names the build and the failed operation without anyone inspecting the machine. Nothing is
//! sent anywhere: it goes to the clipboard for you to paste.

use std::collections::VecDeque;
use std::sync::OnceLock;

use eframe::egui::{self, RichText};

use super::App;

/// The graphics adapter, recorded at startup.
static GPU: OnceLock<String> = OnceLock::new();

/// Errors kept for diagnostics.
const KEEP_ERRORS: usize = 20;

pub fn set_gpu(info: &eframe::wgpu::AdapterInfo) {
    let _ = GPU.set(format!("{} ({:?}, driver {} {})", info.name, info.backend, info.driver, info.driver_info).replace(" ()", ""));
}

/// `GGV 0.9.0 (commit …, target)`.
pub fn build_line() -> String {
    format!("GGV {} (commit {}, {}, {} build)", env!("CARGO_PKG_VERSION"), env!("GGV_COMMIT"), env!("GGV_TARGET"), env!("GGV_PROFILE"))
}

/// Recent errors, oldest first, each with when it happened.
#[derive(Default)]
pub struct ErrorLog {
    seen: Option<String>,
    pub entries: VecDeque<(u64, String)>,
}

impl ErrorLog {
    /// Records the error showing now, if it's new.
    pub fn note(&mut self, current: Option<&String>) {
        if current == self.seen.as_ref() {
            return;
        }
        self.seen = current.cloned();
        if let Some(e) = current {
            let now = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(0, |d| d.as_secs());
            self.entries.push_back((now, e.clone()));
            while self.entries.len() > KEEP_ERRORS {
                self.entries.pop_front();
            }
        }
    }
}

impl App {
    /// The diagnostics text. Without `paths`, files appear by name only (sample names inside
    /// them aren't included either way).
    pub(super) fn diagnostics(&self, paths: bool) -> String {
        let name = |p: &std::path::Path| if paths { p.display().to_string() } else { crate::files::display_name(p) };
        let mut out = vec![
            build_line(),
            format!("system: {} {}", std::env::consts::OS, std::env::consts::ARCH),
            format!("graphics: {}", GPU.get().map_or("unknown", String::as_str)),
            format!("read window: {}", super::fmt_bp(super::read_window())),
        ];
        out.push(format!("view: {}", self.locus_label().unwrap_or_else(|| "nothing open".into())));
        if let Some(r) = &self.reference_path {
            out.push(format!("reference: {}", name(r)));
        }
        for t in &self.tracks {
            let f = self.track_style(&t.loader.bam_path).filter;
            let mut line = format!("alignments: {}  (exclude flags 0x{:x}, min MAPQ {}", name(&t.loader.bam_path), f.exclude_flags(), f.min_mapq);
            if let Some(e) = f.expr {
                line.push_str(&format!(", filter {}", crate::expr::text(e)));
            }
            line.push(')');
            if let Some(e) = &t.error {
                line.push_str(&format!("  error: {e}"));
            }
            out.push(line);
        }
        for v in &self.variants {
            out.push(format!("variants: {}{}", name(&v.path), v.error.as_ref().map_or(String::new(), |e| format!("  error: {e}"))));
        }
        for g in &self.genes {
            out.push(format!("genes: {}", name(&g.path)));
        }
        for s in &self.signals {
            out.push(format!("signal: {}{}", name(&s.path), s.error.as_ref().map_or(String::new(), |e| format!("  error: {e}"))));
        }
        if self.error_log.entries.is_empty() {
            out.push("recent errors: none".into());
        } else {
            out.push("recent errors (UTC):".into());
            for (at, e) in &self.error_log.entries {
                let e = if paths { e.clone() } else { strip_dirs(e) };
                out.push(format!("  {}  {e}", ggv_render::report::utc_stamp(*at)));
            }
        }
        out.join("\n")
    }

    pub(super) fn about_window(&mut self, ctx: &egui::Context) {
        if !self.show_about {
            return;
        }
        let mut open = true;
        let mut copy = false;
        egui::Window::new("About GGV")
            .id(egui::Id::unique(("about", self.reference_slot)))
            .collapsible(false)
            .resizable(true)
            .default_width(520.0)
            .open(&mut open)
            .show(ctx, |ui| {
                ui.label(RichText::new(build_line()).strong());
                ui.label("A fast genome viewer for people and agents. MIT licensed.");
                ui.add_space(6.0);
                ui.horizontal(|ui| {
                    if ui.button("Copy diagnostics").on_hover_text("For a bug report: the build, system, graphics, what's open and recent errors. Nothing is sent anywhere.").clicked() {
                        copy = true;
                    }
                    ui.checkbox(&mut self.diagnostics_paths, "Include full file paths");
                });
                ui.separator();
                egui::ScrollArea::vertical().max_height(320.0).show(ui, |ui| {
                    ui.add(egui::Label::new(RichText::new(self.diagnostics(self.diagnostics_paths)).monospace().size(11.0)).wrap());
                });
            });
        if copy {
            ctx.copy_text(self.diagnostics(self.diagnostics_paths));
            self.notify("Diagnostics copied: paste them into your bug report");
        }
        if !open {
            self.show_about = false;
        }
    }
}

/// `C:\data\x\a.bam: …` → `a.bam: …`: folders out of error messages unless paths are wanted.
fn strip_dirs(s: &str) -> String {
    s.split(' ')
        .map(|w| {
            let cut = w.rfind(['/', '\\']).filter(|_| w.contains(['/', '\\']) && !w.starts_with("http"));
            cut.map_or_else(|| w.to_owned(), |i| w[i + 1..].to_owned())
        })
        .collect::<Vec<_>>()
        .join(" ")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn folders_come_out_of_errors() {
        assert_eq!(strip_dirs("C:\\data\\x\\a.bam: no index"), "a.bam: no index");
        assert_eq!(strip_dirs("/home/me/s.bam: bad"), "s.bam: bad");
        assert_eq!(strip_dirs("https://x.org/a.bam failed"), "https://x.org/a.bam failed");
    }

    #[test]
    fn errors_are_logged_once_each() {
        let mut log = ErrorLog::default();
        let e = Some("boom".to_owned());
        log.note(e.as_ref());
        log.note(e.as_ref());
        log.note(None);
        log.note(e.as_ref());
        assert_eq!(log.entries.len(), 2);
    }
}
