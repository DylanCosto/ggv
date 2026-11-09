//! GGV's per-user folders: settings (`%APPDATA%\ggv`, `$XDG_CONFIG_HOME/ggv` or
//! `~/.config/ggv`) and downloads (`%LOCALAPPDATA%\ggv`, `~/Library/Caches/ggv` or
//! `$XDG_CACHE_HOME/ggv`).
//!
//! Compatibility with GGV's earlier name, gv, kept small on purpose: settings in the old `gv`
//! config folder are copied to the new one the first time it's needed, the old cache folder
//! (downloaded genomes, the demo, tracks, annotation caches) is used as long as there is no new
//! one, and `GV_*` environment variables still work when the `GGV_*` one isn't set.

use std::{
    ffi::OsString,
    path::{Path, PathBuf},
};

/// The folder name, under the per-user config and cache folders.
const NAME: &str = "ggv";
/// The folder name before the rename to GGV.
const OLD_NAME: &str = "gv";

fn config_base() -> Option<PathBuf> {
    std::env::var_os("APPDATA")
        .or_else(|| std::env::var_os("XDG_CONFIG_HOME"))
        .map(PathBuf::from)
        .or_else(|| std::env::var_os("HOME").map(|h| PathBuf::from(h).join(".config")))
}

fn cache_base() -> Option<PathBuf> {
    if cfg!(windows) {
        std::env::var_os("LOCALAPPDATA").map(PathBuf::from)
    } else if cfg!(target_os = "macos") {
        std::env::var_os("HOME").map(|h| PathBuf::from(h).join("Library/Caches"))
    } else {
        std::env::var_os("XDG_CACHE_HOME").map(PathBuf::from).or_else(|| std::env::var_os("HOME").map(|h| PathBuf::from(h).join(".cache")))
    }
}

/// The per-user config folder (settings, the last session, review decisions kept outside the
/// data's folder). Settings from before the rename are copied in first, once.
pub fn config_dir() -> Option<PathBuf> {
    let base = config_base()?;
    let dir = base.join(NAME);
    migrate(&base.join(OLD_NAME), &dir);
    Some(dir)
}

/// The per-user cache folder (downloads and annotation caches): the one from before the rename
/// while there is no new one, so nothing is downloaded again.
pub fn cache_dir() -> Option<PathBuf> {
    cache_base().map(|base| cache_in(&base))
}

/// `base/ggv`, or `base/gv` when only that one exists.
fn cache_in(base: &Path) -> PathBuf {
    let (dir, old) = (base.join(NAME), base.join(OLD_NAME));
    if !dir.exists() && old.is_dir() { old } else { dir }
}

/// The environment variable `GGV_<name>`, else `GV_<name>` (its name before the rename).
pub fn env_var_os(name: &str) -> Option<OsString> {
    std::env::var_os(format!("GGV_{name}")).or_else(|| std::env::var_os(format!("GV_{name}")))
}

/// Copies the config folder `old` to `new` when `new` doesn't exist yet, leaving `old` as it was
/// (an older install may still use it). Endpoint files of running windows (`attach`) aren't
/// copied. The copy is made beside `new` and renamed into place, so a copy cut short is never
/// taken for the real thing; another process doing the same at once is harmless. Returns
/// whether it copied.
fn migrate(old: &Path, new: &Path) -> bool {
    if new.exists() || !old.is_dir() {
        return false;
    }
    let tmp = new.with_extension(format!("copy{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&tmp);
    let copied = copy_dir(old, &tmp, true).and_then(|()| std::fs::rename(&tmp, new));
    if copied.is_err() {
        let _ = std::fs::remove_dir_all(&tmp);
    }
    copied.is_ok()
}

fn copy_dir(from: &Path, to: &Path, top: bool) -> std::io::Result<()> {
    std::fs::create_dir_all(to)?;
    for entry in std::fs::read_dir(from)? {
        let entry = entry?;
        let (src, dest) = (entry.path(), to.join(entry.file_name()));
        let kind = entry.file_type()?;
        if kind.is_dir() {
            if !(top && entry.file_name() == "attach") {
                copy_dir(&src, &dest, false)?;
            }
        } else if kind.is_file() {
            std::fs::copy(&src, &dest)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scratch(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("ggv-folders-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn settings_from_before_the_rename_are_copied_once() {
        let base = scratch("config");
        let (old, new) = (base.join("gv"), base.join("ggv"));
        // Nothing to copy.
        assert!(!migrate(&old, &new));
        assert!(!new.exists());
        std::fs::create_dir_all(old.join("reviews")).unwrap();
        std::fs::create_dir_all(old.join("attach")).unwrap();
        std::fs::write(old.join("last-session.json"), "{\"version\":1}").unwrap();
        std::fs::write(old.join("reviews").join("calls.vcf.gz.0000abcd.gv-review.tsv"), "#contig\n").unwrap();
        std::fs::write(old.join("attach").join("123.json"), "{}").unwrap();
        assert!(migrate(&old, &new));
        assert_eq!(std::fs::read_to_string(new.join("last-session.json")).unwrap(), "{\"version\":1}");
        assert!(new.join("reviews").join("calls.vcf.gz.0000abcd.gv-review.tsv").is_file());
        assert!(!new.join("attach").exists(), "running windows' endpoints stay where they are");
        // The old folder is left as it was, and no partial copy is left behind.
        assert!(old.join("last-session.json").is_file());
        let names: Vec<_> = std::fs::read_dir(&base).unwrap().map(|e| e.unwrap().file_name().into_string().unwrap()).collect();
        assert_eq!(names.len(), 2, "{names:?}");
        // Once the new folder exists, it is GGV's: later changes to the old one aren't copied.
        std::fs::write(new.join("last-session.json"), "{\"version\":2}").unwrap();
        assert!(!migrate(&old, &new));
        assert_eq!(std::fs::read_to_string(new.join("last-session.json")).unwrap(), "{\"version\":2}");
        let _ = std::fs::remove_dir_all(base);
    }

    #[test]
    fn the_old_cache_is_used_until_there_is_a_new_one() {
        let base = scratch("cache");
        assert_eq!(cache_in(&base), base.join("ggv"));
        std::fs::create_dir_all(base.join("gv").join("genomes")).unwrap();
        assert_eq!(cache_in(&base), base.join("gv"));
        std::fs::create_dir_all(base.join("ggv")).unwrap();
        assert_eq!(cache_in(&base), base.join("ggv"));
        let _ = std::fs::remove_dir_all(base);
    }
}
