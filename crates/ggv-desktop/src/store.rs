//! Saving work so that a crash, a full disk or a read-only share never costs it.
//!
//! Every save goes to a temporary file in the same folder, is flushed to disk and then renamed
//! over the old file, so the file is always either the old version or the new one, never half
//! of each. The previous version is kept beside it as `<name>.bak`.
//!
//! Review and queue decisions are kept next to the data they are about when that folder can be
//! written, in a folder the person chose (`Variants > Review > Save reviews in…`), or otherwise
//! in GGV's own folder, so a VCF on a read-only share can still be reviewed.

use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
    path::{Path, PathBuf},
};

/// Writes `bytes` to `path` all or nothing, keeping the previous version as `path.bak`.
pub fn write_atomic(path: &Path, bytes: &[u8]) -> std::io::Result<()> {
    write_keeping(path, bytes, true)
}

/// [`write_atomic`]; `rotate`: the previous version becomes the backup (not when it is damaged,
/// so a good backup isn't replaced by it).
fn write_keeping(path: &Path, bytes: &[u8], rotate: bool) -> std::io::Result<()> {
    if let Some(dir) = path.parent().filter(|d| !d.as_os_str().is_empty()) {
        std::fs::create_dir_all(dir)?;
    }
    let tmp = path.with_added_extension(format!("tmp{}", std::process::id()));
    let written = (|| {
        let mut f = std::fs::File::create(&tmp)?;
        f.write_all(bytes)?;
        f.sync_all()
    })();
    if let Err(e) = written {
        let _ = std::fs::remove_file(&tmp);
        return Err(e);
    }
    if rotate && path.exists() {
        let _ = std::fs::copy(path, backup_of(path));
    }
    std::fs::rename(&tmp, path).inspect_err(|_| {
        let _ = std::fs::remove_file(&tmp);
    })
}

/// How many times each decision file was saved in this process, so a pane reviewing the same
/// file as another can tell it changed and take the other's decisions in (`review`, `queue`).
static SAVES: std::sync::Mutex<std::collections::BTreeMap<PathBuf, u64>> = std::sync::Mutex::new(std::collections::BTreeMap::new());

/// The saves of `path` so far.
pub fn saves(path: &Path) -> u64 {
    SAVES.lock().map_or(0, |m| m.get(path).copied().unwrap_or(0))
}

/// Counts a save of `path`; the saves so far.
pub fn saved(path: &Path) -> u64 {
    SAVES.lock().map_or(0, |mut m| {
        let n = m.entry(path.to_path_buf()).or_default();
        *n += 1;
        *n
    })
}

/// `x.json` → `x.json.bak`.
pub fn backup_of(path: &Path) -> PathBuf {
    path.with_added_extension("bak")
}

/// Reads `path` with `parse`; when it exists but can't be read or doesn't parse (a damaged
/// file), its backup. Says which. `parse` should only check that the content is readable: a
/// file it reads but the caller then refuses must not be swapped for an older backup.
pub fn read_recovering<T>(path: &Path, parse: impl Fn(&str) -> anyhow::Result<T>) -> anyhow::Result<(T, bool)> {
    let first = std::fs::read_to_string(path).map_err(anyhow::Error::from).and_then(|t| parse(&t));
    match first {
        Ok(v) => Ok((v, false)),
        Err(e) if !path.exists() => Err(e),
        Err(e) => {
            let bak = backup_of(path);
            match std::fs::read_to_string(&bak).map_err(anyhow::Error::from).and_then(|t| parse(&t)) {
                Ok(v) if bak.exists() => Ok((v, true)),
                _ => Err(e),
            }
        }
    }
}

#[cfg(test)]
thread_local! {
    /// GGV's own folder for the tests on this thread (instead of the person's).
    pub static TEST_CONFIG: std::cell::RefCell<Option<PathBuf>> = const { std::cell::RefCell::new(None) };
}

/// GGV's own folder (where the last session is kept).
pub(crate) fn config_dir() -> Option<PathBuf> {
    #[cfg(test)]
    if let Some(d) = TEST_CONFIG.with_borrow(Clone::clone) {
        return Some(d);
    }
    crate::session::last_session_path().and_then(|p| p.parent().map(Path::to_path_buf))
}

fn prefs_file() -> Option<PathBuf> {
    config_dir().map(|d| d.join("prefs.json"))
}

/// The folder the person chose for review and queue decisions, if any.
pub fn output_dir() -> Option<PathBuf> {
    let text = std::fs::read_to_string(prefs_file()?).ok()?;
    let prefs: serde_json::Value = serde_json::from_str(&text).ok()?;
    prefs.get("review_dir")?.as_str().filter(|s| !s.is_empty()).map(PathBuf::from)
}

/// Where decisions go when neither the data's folder nor a chosen one can be written.
fn fallback_dir() -> Option<PathBuf> {
    config_dir().map(|d| d.join("reviews"))
}

/// A name for `data`'s decisions outside its folder: its file name plus a short hash of its
/// full path, so two `calls.vcf.gz` in different folders don't share decisions.
fn elsewhere(dir: &Path, data: &Path, ext: &str) -> PathBuf {
    let full = std::path::absolute(data).unwrap_or_else(|_| data.to_path_buf());
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for b in full.to_string_lossy().bytes() {
        h = (h ^ u64::from(b)).wrapping_mul(0x100_0000_01b3);
    }
    let name = data.file_name().map_or_else(|| "data".into(), |n| n.to_string_lossy().into_owned());
    dir.join(format!("{name}.{:08x}.{ext}", h as u32))
}

/// Every place `data`'s decisions may be, most preferred first: the chosen folder, next to the
/// data, GGV's own folder.
fn candidates(data: &Path, ext: &str) -> Vec<PathBuf> {
    let mut v = Vec::new();
    if let Some(d) = output_dir() {
        v.push(elsewhere(&d, data, ext));
    }
    v.push(data.with_added_extension(ext));
    if let Some(d) = fallback_dir() {
        v.push(elsewhere(&d, data, ext));
    }
    v
}

/// Where each data file's decisions were last saved in this process, by data file and `ext`:
/// a property of the file, so every pane reviewing it reads and saves the same place, also
/// after a save had to go elsewhere (a read-only folder).
static LOCATED: std::sync::Mutex<BTreeMap<(PathBuf, String), PathBuf>> = std::sync::Mutex::new(BTreeMap::new());

fn located_key(data: &Path, ext: &str) -> (PathBuf, String) {
    (std::path::absolute(data).unwrap_or_else(|_| data.to_path_buf()), ext.to_owned())
}

/// Where `data`'s decisions were last saved in this process, if they were.
pub fn located(data: &Path, ext: &str) -> Option<PathBuf> {
    LOCATED.lock().ok()?.get(&located_key(data, ext)).cloned()
}

fn locate(data: &Path, ext: &str, path: &Path) {
    if let Ok(mut m) = LOCATED.lock() {
        m.insert(located_key(data, ext), path.to_path_buf());
    }
}

/// Where `data`'s decisions (`ext`, e.g. `ggv-review.tsv`) are: where they were last saved in
/// this process, else the first place that has them, or where new ones go.
pub fn sidecar(data: &Path, ext: &str) -> PathBuf {
    if let Some(p) = located(data, ext) {
        return p;
    }
    let all = candidates(data, ext);
    all.iter()
        .find(|p| p.exists())
        .cloned()
        .or_else(|| legacy_candidates(data, ext).into_iter().find(|p| p.exists()))
        .unwrap_or_else(|| all[0].clone())
}

/// Where decisions saved before the rename to GGV may be (`calls.vcf.gz.gv-review.tsv`): read,
/// and kept being saved to, when there are none under the new name.
fn legacy_candidates(data: &Path, ext: &str) -> Vec<PathBuf> {
    ext.strip_prefix("ggv-").map_or_else(Vec::new, |rest| candidates(data, &format!("gv-{rest}")))
}

/// Every place `data`'s decisions may be saved (see [`sidecar`]), then where those from before
/// the rename may be.
pub fn sidecars(data: &Path, ext: &str) -> Vec<PathBuf> {
    let mut v = candidates(data, ext);
    v.extend(legacy_candidates(data, ext));
    v
}

/// Where new decisions for `data` go.
pub fn preferred(data: &Path, ext: &str) -> PathBuf {
    candidates(data, ext).swap_remove(0)
}

/// Sets (or, with `None`, clears) the folder for review and queue decisions.
pub fn set_output_dir(dir: Option<&Path>) {
    let Some(file) = prefs_file() else { return };
    let mut prefs: serde_json::Map<String, serde_json::Value> =
        std::fs::read_to_string(&file).ok().and_then(|t| serde_json::from_str(&t).ok()).unwrap_or_default();
    match dir {
        Some(d) => prefs.insert("review_dir".into(), d.to_string_lossy().into()),
        None => prefs.remove("review_dir"),
    };
    if let Ok(text) = serde_json::to_string_pretty(&prefs) {
        let _ = write_atomic(&file, text.as_bytes());
    }
}

/// Decisions read from a file (`saved`) taken in with the ones in memory, by key: a key only
/// in memory is kept; one in both takes the saved decision (another pane's, the later), except
/// keys decided here and not saved yet (`unsaved`) or, with `memory_wins` (the file was
/// damaged and `saved` is its older backup, or it is a place these decisions are moving to),
/// every key in memory.
pub fn merge_decisions<K: Ord + Clone, D: Clone>(memory: &BTreeMap<K, D>, unsaved: &BTreeSet<K>, mut saved: BTreeMap<K, D>, memory_wins: bool) -> BTreeMap<K, D> {
    for (k, d) in memory {
        if memory_wins || unsaved.contains(k) || !saved.contains_key(k) {
            saved.insert(k.clone(), d.clone());
        }
    }
    saved
}

/// Writes decisions at `path`: `contents(path)` reads what is there first and gives the bytes to
/// write and whether the file there is damaged (then its backup, the last good save, is kept).
pub fn save_merged(path: &Path, contents: &mut impl FnMut(&Path) -> (Vec<u8>, bool)) -> std::io::Result<()> {
    let (bytes, damaged) = contents(path);
    write_keeping(path, &bytes, !damaged)
}

/// Saves decisions for `data` at `*path` (see [`save_merged`]: whatever is at the place written
/// is read and taken in first); when that can't be written, at the next place that can, moving
/// `*path` there and saying so. Where they went is remembered for `data` ([`located`]), so other
/// panes follow.
pub fn save_sidecar(path: &mut PathBuf, data: &Path, ext: &str, mut contents: impl FnMut(&Path) -> (Vec<u8>, bool)) -> Result<Option<String>, String> {
    let first = match save_merged(path, &mut contents) {
        Ok(()) => {
            locate(data, ext, path);
            return Ok(None);
        }
        Err(e) => format!("{}: {e}", path.display()),
    };
    let all = candidates(data, ext);
    let after = all.iter().position(|p| p == path).map_or(0, |i| i + 1);
    for p in &all[after..] {
        if save_merged(p, &mut contents).is_ok() {
            let note = format!("Couldn't save to {} ({first}); decisions are now saved to {}", path.display(), p.display());
            *path = p.clone();
            locate(data, ext, path);
            return Ok(Some(note));
        }
    }
    Err(first)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atomic_writes_keep_the_previous_version() {
        let dir = std::env::temp_dir().join(format!("ggv-store-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let p = dir.join("a.json");
        write_atomic(&p, b"one").unwrap();
        write_atomic(&p, b"two").unwrap();
        assert_eq!(std::fs::read_to_string(&p).unwrap(), "two");
        assert_eq!(std::fs::read_to_string(backup_of(&p)).unwrap(), "one");
        // No temporary files left behind.
        let names: Vec<_> = std::fs::read_dir(&dir).unwrap().map(|e| e.unwrap().file_name()).collect();
        assert_eq!(names.len(), 2, "{names:?}");
        // A damaged file is read from its backup.
        std::fs::write(&p, "{trunc").unwrap();
        let (v, recovered) = read_recovering(&p, |t| Ok(serde_json::from_str::<serde_json::Value>(t)?)).unwrap_or((serde_json::Value::Null, false));
        assert!(v.is_null() && !recovered, "'one' isn't JSON either");
        write_atomic(&p, b"[1]").unwrap();
        write_atomic(&p, b"[1]").unwrap();
        std::fs::write(&p, "[1").unwrap();
        let (v, recovered) = read_recovering(&p, |t| Ok(serde_json::from_str::<serde_json::Value>(t)?)).unwrap();
        assert_eq!((v, recovered), (serde_json::json!([1]), true));
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn names_elsewhere_differ_by_folder() {
        let a = elsewhere(Path::new("/out"), Path::new("/x/calls.vcf.gz"), "ggv-review.tsv");
        let b = elsewhere(Path::new("/out"), Path::new("/y/calls.vcf.gz"), "ggv-review.tsv");
        assert_ne!(a, b);
        assert!(a.file_name().unwrap().to_string_lossy().starts_with("calls.vcf.gz."));
    }
}
