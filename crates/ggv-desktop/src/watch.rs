//! Watches open files for rewrites (a pipeline re-running) so they can be reloaded in place.
//!
//! A background thread stats each file and its index once a second; the UI is woken only when
//! something changed and then held still for a second (so a file being written is not read
//! half-done). Idle, it costs a few `stat` calls a second and no frames.

use std::{
    path::{Path, PathBuf},
    sync::{Arc, Mutex, mpsc},
    time::{Duration, SystemTime},
};

use eframe::egui;

/// Size and modification time of a file and of the indexes next to it.
pub type Signature = Vec<Option<(u64, SystemTime)>>;

fn stat(p: &Path) -> Option<(u64, SystemTime)> {
    let m = std::fs::metadata(p).ok()?;
    Some((m.len(), m.modified().ok()?))
}

/// The files that change together with `path`: itself and any index beside it.
fn related(path: &Path) -> Vec<PathBuf> {
    let mut out = vec![path.to_owned()];
    for ext in ["bai", "csi", "crai", "tbi", "fai", "gzi"] {
        out.push(path.with_added_extension(ext));
        out.push(path.with_extension(ext));
    }
    out
}

pub fn signature(path: &Path) -> Signature {
    related(path).iter().map(|p| stat(p)).collect()
}

pub struct Watcher {
    paths: Arc<Mutex<Vec<PathBuf>>>,
    pub rx: mpsc::Receiver<PathBuf>,
    pub enabled: Arc<std::sync::atomic::AtomicBool>,
}

impl Watcher {
    pub fn spawn(ctx: &egui::Context) -> Self {
        let paths: Arc<Mutex<Vec<PathBuf>>> = Arc::default();
        let enabled = Arc::new(std::sync::atomic::AtomicBool::new(true));
        let (tx, rx) = mpsc::channel();
        let (p2, e2, ctx) = (paths.clone(), enabled.clone(), ctx.clone());
        let _ = std::thread::Builder::new().name("ggv-watch".into()).spawn(move || {
            // Per path: the signature last accepted, and a newer one waiting to settle.
            let mut seen: Vec<(PathBuf, Signature, Option<Signature>)> = Vec::new();
            loop {
                std::thread::sleep(Duration::from_secs(1));
                if Arc::strong_count(&p2) == 1 {
                    return; // the pane is gone
                }
                if !e2.load(std::sync::atomic::Ordering::Relaxed) {
                    continue;
                }
                let want = p2.lock().unwrap().clone();
                seen.retain(|(p, ..)| want.contains(p));
                for p in want {
                    let now = signature(&p);
                    match seen.iter_mut().find(|(q, ..)| *q == p) {
                        None => seen.push((p, now, None)),
                        Some((_, accepted, pending)) => {
                            if now == *accepted {
                                *pending = None;
                            } else if pending.as_ref() == Some(&now) && now[0].is_some() {
                                // Changed, and unchanged for a second since: reload.
                                *accepted = now;
                                *pending = None;
                                if tx.send(p.clone()).is_err() {
                                    return;
                                }
                                ctx.request_repaint();
                            } else {
                                *pending = Some(now);
                            }
                        }
                    }
                }
            }
        });
        Self { paths, rx, enabled }
    }

    /// The files to watch (cheap to call every frame: only stores the list when it differs).
    pub fn set(&self, files: Vec<PathBuf>) {
        let mut p = self.paths.lock().unwrap();
        if *p != files {
            *p = files;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn signature_sees_rewrites_and_indexes() {
        let dir = std::env::temp_dir().join(format!("ggv-watch-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let f = dir.join("a.bam");
        std::fs::write(&f, b"one").unwrap();
        let before = signature(&f);
        assert!(before[0].is_some());
        std::fs::write(f.with_added_extension("bai"), b"index").unwrap();
        assert_ne!(signature(&f), before, "a new index changes the signature");
        let with_index = signature(&f);
        std::fs::write(&f, b"longer").unwrap();
        assert_ne!(signature(&f), with_index);
        let _ = std::fs::remove_dir_all(dir);
    }
}
