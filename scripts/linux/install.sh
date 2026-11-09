#!/bin/sh
# Installs GGV for the current user: binaries in ~/.local/bin, a desktop entry, and MIME types so
# file managers offer "Open with GGV" for BAM, CRAM, VCF, BCF and annotation files, and ggv:// links
# open GGV.
#
# From the unpacked download (ggv and ggv-desktop next to this script):
#   ./install.sh                 install those binaries
#   ./install.sh --uninstall     remove everything it installed
# From a source checkout (scripts/linux/install.sh):
#   scripts/linux/install.sh             build (cargo, release) and install
#   scripts/linux/install.sh --no-build  install what is already built
#
# What it installs is listed, with each file's SHA-256, in ~/.local/share/ggv/installed. A file
# already at a destination that GGV didn't install is left alone (the script says so and exits
# 1); uninstalling removes only listed files that are unchanged since.
set -eu

here=$(cd "$(dirname "$0")" && pwd)
bin="$HOME/.local/bin"
data="${XDG_DATA_HOME:-$HOME/.local/share}"
apps="$data/applications"
mime="$data/mime"
manifest="$data/ggv/installed"

sha() {
    # From stdin: with a file name, a name holding a backslash changes the output.
    if command -v sha256sum >/dev/null; then sha256sum <"$1"; else shasum -a 256 <"$1"; fi | cut -d' ' -f1
}

# The hash the manifest records for a path (empty if none).
recorded() {
    [ -f "$manifest" ] || return 0
    while IFS= read -r line; do
        if [ "${line#*  }" = "$1" ]; then
            printf '%s\n' "${line%%  *}"
            return 0
        fi
    done <"$manifest"
}

refresh() {
    if command -v update-mime-database >/dev/null; then update-mime-database "$mime" || true; fi
    if command -v update-desktop-database >/dev/null; then update-desktop-database "$apps" || true; fi
}

if [ "${1:-}" = "--uninstall" ]; then
    if [ ! -f "$manifest" ]; then
        echo "nothing to uninstall: no record of an GGV install ($manifest)"
        exit 0
    fi
    kept=0
    while IFS= read -r line; do
        f="${line#*  }"
        [ -e "$f" ] || continue
        if [ "$(sha "$f")" = "${line%%  *}" ]; then
            rm -f "$f"
            echo "removed $f"
        else
            echo "kept $f: it changed after GGV installed it" >&2
            kept=1
        fi
    done <"$manifest"
    rm -f "$manifest"
    rmdir "$(dirname "$manifest")" 2>/dev/null || true
    refresh
    echo "uninstalled"
    exit "$kept"
fi

if [ -x "$here/ggv" ] && [ -x "$here/ggv-desktop" ]; then
    # The release download: binaries ship beside this script.
    target="$here"
else
    root=$(cd "$here/../.." && pwd)
    if [ "${1:-}" != "--no-build" ]; then
        (cd "$root" && cargo build --release -p ggv-cli -p ggv-desktop)
    fi
    target="${CARGO_TARGET_DIR:-$root/target}/release"
fi

case "$bin" in
    *'
'*)
        echo "can't install to $bin: the path has a line break" >&2
        exit 1
        ;;
esac
mkdir -p "$bin" "$apps" "$mime/packages" "$(dirname "$manifest")"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# The desktop entry runs the installed ggv-desktop by its full path, so it works whether or not
# ~/.local/bin is on the desktop session's PATH. The path is quoted and escaped as the Desktop
# Entry spec asks: \\\\ for a backslash, \\" \\` \\$ for those characters, %% for %.
exec_path=$(printf '%s\n' "$bin/ggv-desktop" | sed -e 's/\\/\\\\\\\\/g' -e 's/["`$]/\\\\&/g' -e 's/%/%%/g')
while IFS= read -r line; do
    case "$line" in
        Exec=*) printf 'Exec="%s" %%U\n' "$exec_path" ;;
        *) printf '%s\n' "$line" ;;
    esac
done <"$here/ggv.desktop" >"$work/ggv.desktop"

: >"$work/manifest"
skipped=0
# put SRC DEST MODE: installs SRC at DEST, unless DEST is a file GGV didn't install (one not in
# the manifest, or changed since), which is left alone.
put() {
    want=$(sha "$1")
    if [ -e "$2" ] || [ -L "$2" ]; then
        have=$(sha "$2" 2>/dev/null || echo none)
        was=$(recorded "$2")
        if [ "$have" != "$want" ] && { [ -z "$was" ] || [ "$have" != "$was" ]; }; then
            echo "left $2 alone: GGV didn't install it (move it away and run this again to install there)" >&2
            skipped=1
            return 0
        fi
    fi
    install -m "$3" "$1" "$2"
    printf '%s  %s\n' "$want" "$2" >>"$work/manifest"
}

put "$target/ggv" "$bin/ggv" 755
put "$target/ggv-desktop" "$bin/ggv-desktop" 755
put "$work/ggv.desktop" "$apps/ggv.desktop" 644
put "$here/ggv-mime.xml" "$mime/packages/ggv.xml" 644
mv "$work/manifest" "$manifest"
refresh
if [ "$skipped" = 1 ]; then
    echo "installed the rest to $bin and $data" >&2
    exit 1
fi
echo "installed ggv and ggv-desktop to $bin"
case ":$PATH:" in
    *":$bin:"*) ;;
    *) echo "add $bin to your PATH to run GGV from a terminal" ;;
esac
"$bin/ggv" --version
