#!/bin/sh
# Installs GGV from the unpacked macOS download: GGV.app into ~/Applications (or /Applications
# with --system) and the `ggv` command-line tool into ~/.local/bin (or /usr/local/bin).
#   ./install.sh               install for this user
#   ./install.sh --system      install for every user (asks for your password)
#   ./install.sh --uninstall   remove what it installed (add --system for the system install)
#
# What it installs is listed, with a SHA-256 of each, in ~/.local/share/ggv/installed (or
# /usr/local/share/ggv/installed). Something already at a destination that GGV didn't install is
# left alone (the script says so and exits 1); uninstalling removes only listed items that are
# unchanged since.
set -eu

here=$(cd "$(dirname "$0")" && pwd)
apps="$HOME/Applications"
bin="$HOME/.local/bin"
manifest="$HOME/.local/share/ggv/installed"
sudo=""
uninstall=""
for a in "$@"; do
    case "$a" in
        --system)
            apps="/Applications"
            bin="/usr/local/bin"
            manifest="/usr/local/share/ggv/installed"
            sudo="sudo"
            ;;
        --uninstall) uninstall=1 ;;
    esac
done

# A SHA-256 of a file, or of a folder (each path in it, with its file's hash or link's target).
sha() {
    if [ -d "$1" ] && [ ! -L "$1" ]; then
        (cd "$1" && find . -print | LC_ALL=C sort | while IFS= read -r f; do
            if [ -L "$f" ]; then
                printf 'l %s %s\n' "$f" "$(readlink "$f")"
            elif [ -f "$f" ]; then
                printf 'f %s %s\n' "$f" "$(shasum -a 256 <"$f" | cut -d' ' -f1)"
            else
                printf 'd %s\n' "$f"
            fi
        done) | shasum -a 256 | cut -d' ' -f1
    else
        shasum -a 256 <"$1" | cut -d' ' -f1
    fi
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

if [ -n "$uninstall" ]; then
    if [ ! -f "$manifest" ]; then
        echo "nothing to uninstall: no record of an GGV install ($manifest)"
        exit 0
    fi
    kept=0
    while IFS= read -r line; do
        f="${line#*  }"
        [ -e "$f" ] || [ -L "$f" ] || continue
        if [ "$(sha "$f" 2>/dev/null)" = "${line%%  *}" ]; then
            $sudo rm -rf "$f"
            echo "removed $f"
        else
            echo "kept $f: it changed after GGV installed it" >&2
            kept=1
        fi
    done <"$manifest"
    $sudo rm -f "$manifest"
    $sudo rmdir "$(dirname "$manifest")" 2>/dev/null || true
    exit "$kept"
fi

case "$apps$bin" in
    *'
'*)
        echo "can't install there: the path has a line break" >&2
        exit 1
        ;;
esac

# A copy downloaded through a browser is quarantined; an unnotarized build would then be
# refused. Clearing it here is the same as right-click > Open, for these files only.
xattr -dr com.apple.quarantine "$here/GGV.app" "$here/ggv" 2>/dev/null || true

$sudo mkdir -p "$apps" "$bin" "$(dirname "$manifest")"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
: >"$work/manifest"
skipped=0

# owned SRC DEST: true if DEST is free, or GGV's to replace (listed in the manifest and unchanged,
# or already exactly SRC); otherwise says it is left alone.
owned() {
    [ -e "$2" ] || [ -L "$2" ] || return 0
    have=$(sha "$2" 2>/dev/null || echo none)
    was=$(recorded "$2")
    if [ "$have" = "$(sha "$1")" ] || { [ -n "$was" ] && [ "$have" = "$was" ]; }; then
        return 0
    fi
    echo "left $2 alone: GGV didn't install it (move it away and run this again to install there)" >&2
    skipped=1
    return 1
}

if owned "$here/GGV.app" "$apps/GGV.app"; then
    # Copied beside it first, so a failed copy leaves the old one in place.
    $sudo rm -rf "$apps/.GGV.app.new"
    $sudo cp -R "$here/GGV.app" "$apps/.GGV.app.new"
    $sudo rm -rf "$apps/GGV.app"
    $sudo mv "$apps/.GGV.app.new" "$apps/GGV.app"
    printf '%s  %s\n' "$(sha "$apps/GGV.app")" "$apps/GGV.app" >>"$work/manifest"
fi
if owned "$here/ggv" "$bin/ggv"; then
    $sudo install -m 755 "$here/ggv" "$bin/ggv"
    printf '%s  %s\n' "$(sha "$bin/ggv")" "$bin/ggv" >>"$work/manifest"
fi
$sudo cp "$work/manifest" "$manifest"
$sudo chmod 644 "$manifest"
if [ "$skipped" = 1 ]; then
    exit 1
fi
echo "installed $apps/GGV.app and $bin/ggv"
case ":$PATH:" in
    *":$bin:"*) ;;
    *) echo "add $bin to your PATH to run GGV from a terminal" ;;
esac
"$bin/ggv" --version
