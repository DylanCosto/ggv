#!/bin/sh
# Builds GGV in release mode and packs dist/ggv-<version>-linux-<arch>.tar.gz: ggv, ggv-desktop,
# install.sh (installs them for the current user, with a desktop entry and file types), README
# and LICENSE.
#
#   scripts/linux/package.sh             build and pack
#   scripts/linux/package.sh --no-build  pack what is already built
#
# Build on the oldest Linux you support: the binaries need that glibc or newer (the release
# workflow uses Ubuntu 22.04, glibc 2.35).
set -eu

root=$(cd "$(dirname "$0")/../.." && pwd)
cd "$root"
if [ "${1:-}" != "--no-build" ]; then
    cargo build --release --locked -p ggv-cli -p ggv-desktop
fi
target="${CARGO_TARGET_DIR:-$root/target}/release"
version=$(sed -n 's/^version = "\(.*\)"/\1/p' Cargo.toml | head -1)
arch=$(uname -m)
name="ggv-$version-linux-$arch"
stage="$root/dist/$name"
rm -rf "$stage"
mkdir -p "$stage"
install -m 755 "$target/ggv" "$target/ggv-desktop" "$stage/"
install -m 755 scripts/linux/install.sh "$stage/install.sh"
install -m 644 scripts/linux/ggv.desktop scripts/linux/ggv-mime.xml README.md LICENSE THIRD_PARTY_NOTICES.md "$stage/"
tar -C "$root/dist" -czf "$root/dist/$name.tar.gz" "$name"
echo "packed $root/dist/$name.tar.gz"
