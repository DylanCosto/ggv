# Installing GGV

Download the package for your system from the
[Releases](https://github.com/DylanCosto/ggv/releases) page. Each holds `ggv` (command line and
MCP server) and `ggv-desktop` (the viewer), an installer, this README, the LICENSE and THIRD_PARTY_NOTICES.md (the open-source crates
inside). To check
a download, compare `sha256sum <file>` (`Get-FileHash <file>` on Windows) with the release's
`SHA256SUMS`. `ggv --version` prints the exact build:
`ggv 0.9.0 (commit 36ffdb316176, x86_64-pc-windows-msvc)`; give that line in bug reports.

| System | Package | Needs |
|---|---|---|
| Windows 10/11 (x64) | `ggv-<version>-windows-x64.zip` | a GPU driver with DirectX 12 or Vulkan |
| Linux (x86_64) | `ggv-<version>-linux-x86_64.tar.gz` | glibc 2.35+ (Ubuntu 22.04, Debian 12, Fedora 36 or newer), X11 or Wayland, Vulkan or OpenGL drivers |
| macOS 11+ (Apple silicon and Intel) | `ggv-<version>-macos-universal.zip` | nothing else |

**Windows.** Unzip; it makes a `ggv-<version>-windows-x64` folder. In that folder:

```powershell
powershell -ExecutionPolicy Bypass -File install.ps1
```

It copies GGV to `%LOCALAPPDATA%\Programs\ggv` for you alone (no administrator rights), puts `ggv`
on your PATH for new terminals, adds GGV to the Start menu, and adds **Open in GGV** to the
right-click menu of `.bam`, `.cram`, `.bcf`, `.vcf`, `.gtf`, `.gff(3)` and `.bed` files, leaving
your default programs alone (`-NoFileTypes` skips that). `install.ps1 -Uninstall` removes all
of it. Or skip installing: unzip anywhere and run `ggv-desktop.exe`; `register-file-types.ps1`
adds just the menu entries (`-IncludeGz` adds them to `.gz` files, `-DryRun` shows what it
would write, `-Unregister` removes them). `ggv-desktop.exe` opens no console window from
Explorer; started from a terminal it still prints there.

**Linux.** Unpack and install for your user:

```sh
tar xzf ggv-<version>-linux-x86_64.tar.gz
ggv-<version>-linux-x86_64/install.sh
```

Both binaries go to `~/.local/bin`, with a desktop entry and MIME types so file managers offer
GGV for alignment, variant and annotation files and `ggv://` links open it. `install.sh
--uninstall` removes them.

**macOS.** Unzip, then in that folder:

```sh
./install.sh            # GGV.app to ~/Applications, ggv to ~/.local/bin
./install.sh --system   # /Applications and /usr/local/bin, for every user
```

`./install.sh --uninstall` (with `--system` if used) removes them. Or drag `GGV.app` to
Applications yourself. Release builds signed and notarized by Apple open directly; an unsigned
build (ad-hoc signed, as when built from source) needs its quarantine cleared once:
`install.sh` does that, or run `xattr -dr com.apple.quarantine GGV.app ggv`, or right-click the
app and choose Open. Open files from the app, by dragging them onto its window, or with
`open -a ggv --args sample.bam`; Finder's "Open With" isn't wired up yet.

**First run.** With nothing open, GGV shows how to start. **Open the demo** downloads 4.3 MB
of real PacBio HiFi reads of GIAB HG002 (chr20, with calls, genes and the reference) once and
opens it with a five-step tour: a deletion on one haplotype, Evidence at a base, methylation,
variant review, then what your own files need. Also `ggv-desktop --demo`, or `ggv demo` to
download it and print commands to try. **Commands → About GGV** shows the exact build and has
**Copy diagnostics** (build, system, graphics, open files and recent errors, file names only
unless you tick full paths) for bug reports; nothing is sent anywhere. Evidence reports end
with a Provenance section: the build, reference, files and read filters that produced them.

**From source** (Rust 1.91+): `cargo build --release -p ggv-cli -p ggv-desktop`, or package as the
release does: `scripts/windows/package.ps1`, `scripts/linux/package.sh`, `scripts/macos/package.sh`
(`--native` for this Mac's architecture only). From a checkout, `scripts/linux/install.sh` builds
and installs in one step.

**Publishing a release.** Push a tag such as `v0.9.0`: GitHub Actions builds all three packages,
runs each, writes `SHA256SUMS` and attaches everything to that tag's release. To have Apple sign
and notarize the macOS app, add these repository secrets: `MACOS_CERTIFICATE` (a Developer ID
Application certificate as a base64 `.p12`), `MACOS_CERTIFICATE_PASSWORD`, `MACOS_SIGN_IDENTITY`
(e.g. `Developer ID Application: Name (TEAMID)`), `APPLE_ID`, `APPLE_TEAM_ID` and
`APPLE_APP_PASSWORD` (an app-specific password). Without them the macOS app is ad-hoc signed.
