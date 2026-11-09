// Build identity for GGV's binaries (included as each binary crate's build.rs): the source
// commit and the target, so `--version`, About and bug reports name the exact build.
//
// GGV_COMMIT in the environment wins (packaging from a source archive without .git); otherwise
// `git rev-parse`; otherwise "unknown".

use std::process::Command;

fn main() {
    let commit = std::env::var("GGV_COMMIT").ok().filter(|c| !c.is_empty()).or_else(git_commit).unwrap_or_else(|| "unknown".to_owned());
    let target = std::env::var("TARGET").unwrap_or_default();
    let profile = std::env::var("PROFILE").unwrap_or_default();
    println!("cargo:rustc-env=GGV_COMMIT={commit}");
    println!("cargo:rustc-env=GGV_TARGET={target}");
    println!("cargo:rustc-env=GGV_PROFILE={profile}");
    println!("cargo:rerun-if-env-changed=GGV_COMMIT");
    // A new commit or branch changes these.
    if let Some(git) = git_dir() {
        println!("cargo:rerun-if-changed={}", git.join("HEAD").display());
        println!("cargo:rerun-if-changed={}", git.join("refs").display());
        println!("cargo:rerun-if-changed={}", git.join("packed-refs").display());
    }
}

fn git_commit() -> Option<String> {
    let out = Command::new("git").args(["rev-parse", "--short=12", "HEAD"]).output().ok()?;
    out.status.success().then(|| String::from_utf8_lossy(&out.stdout).trim().to_owned()).filter(|s| !s.is_empty())
}

fn git_dir() -> Option<std::path::PathBuf> {
    let out = Command::new("git").args(["rev-parse", "--absolute-git-dir"]).output().ok()?;
    out.status.success().then(|| std::path::PathBuf::from(String::from_utf8_lossy(&out.stdout).trim()))
}
