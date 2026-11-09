#!/usr/bin/env bash
# Builds GGV for Apple silicon and Intel, joins them into universal binaries, and packs
# dist/ggv-<version>-macos-universal.zip: GGV.app (the desktop viewer), the `ggv` command-line tool,
# install.sh, README and LICENSE.
#
#   scripts/macos/package.sh            # build both architectures and pack
#   scripts/macos/package.sh --native   # this Mac's architecture only (faster, for testing)
#
# Signing: with MACOS_SIGN_IDENTITY set (a "Developer ID Application: ..." identity in the
# keychain), the app and tool are signed with the hardened runtime; with APPLE_ID, APPLE_TEAM_ID
# and APPLE_APP_PASSWORD (an app-specific password) too, the package is notarized by Apple and
# the ticket stapled, so it opens without any Gatekeeper prompt. Without them the build is
# ad-hoc signed (it runs, but a downloaded copy needs its quarantine cleared once; see README).
set -euo pipefail

root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$root"
target="${CARGO_TARGET_DIR:-$root/target}"
version="$(sed -n 's/^version = "\(.*\)"/\1/p' Cargo.toml | head -1)"

if [[ "${1:-}" == "--native" ]]; then
  cargo build --release --locked -p ggv-cli -p ggv-desktop
  arch="$(uname -m)"
  name="ggv-$version-macos-$arch"
  bins=("$target/release")
else
  targets=(aarch64-apple-darwin x86_64-apple-darwin)
  rustup target add "${targets[@]}"
  for t in "${targets[@]}"; do
    # Runs on macOS 11 and later.
    MACOSX_DEPLOYMENT_TARGET=11.0 cargo build --release --locked -p ggv-cli -p ggv-desktop --target "$t"
  done
  name="ggv-$version-macos-universal"
  bins=("$target/aarch64-apple-darwin/release" "$target/x86_64-apple-darwin/release")
fi

# One binary from one or more architectures.
join() {
  local exe="$1" out="$2" inputs=()
  for b in "${bins[@]}"; do inputs+=("$b/$exe"); done
  if [[ ${#inputs[@]} -eq 1 ]]; then cp "${inputs[0]}" "$out"; else lipo -create -output "$out" "${inputs[@]}"; fi
}

stage="$root/dist/$name"
rm -rf "$stage"
app="$stage/GGV.app/Contents"
mkdir -p "$app/MacOS" "$app/Resources"
join ggv-desktop "$app/MacOS/ggv-desktop"
join ggv "$stage/ggv"
cat > "$app/Info.plist" <<PLIST
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>CFBundleName</key><string>GGV</string>
  <key>CFBundleDisplayName</key><string>GGV</string>
  <key>CFBundleIdentifier</key><string>io.github.dylancosto.ggv</string>
  <key>CFBundleExecutable</key><string>ggv-desktop</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>CFBundleVersion</key><string>$version</string>
  <key>CFBundleShortVersionString</key><string>$version</string>
  <key>LSMinimumSystemVersion</key><string>11.0</string>
  <key>NSHighResolutionCapable</key><true/>
  <key>LSApplicationCategoryType</key><string>public.app-category.developer-tools</string>
</dict>
</plist>
PLIST
cp README.md LICENSE THIRD_PARTY_NOTICES.md "$stage/"
cp scripts/macos/install.sh "$stage/install.sh"
chmod +x "$stage/install.sh"

zip="$root/dist/$name.zip"
pack() {
  rm -f "$zip"
  ditto -c -k --keepParent "$stage" "$zip"
}

if [[ -n "${MACOS_SIGN_IDENTITY:-}" ]]; then
  # Developer ID with the hardened runtime and a secure timestamp: what notarization requires.
  codesign --force --options runtime --timestamp --sign "$MACOS_SIGN_IDENTITY" "$stage/ggv"
  codesign --force --options runtime --timestamp --sign "$MACOS_SIGN_IDENTITY" "$stage/GGV.app"
  codesign --verify --strict --verbose=2 "$stage/GGV.app"
  pack
  if [[ -n "${APPLE_ID:-}" && -n "${APPLE_TEAM_ID:-}" && -n "${APPLE_APP_PASSWORD:-}" ]]; then
    xcrun notarytool submit "$zip" --apple-id "$APPLE_ID" --team-id "$APPLE_TEAM_ID" --password "$APPLE_APP_PASSWORD" --wait
    # The ticket goes on the app (a bare executable can't be stapled; Gatekeeper checks it online).
    xcrun stapler staple "$stage/GGV.app"
    pack
    spctl --assess --type execute --verbose "$stage/GGV.app"
    echo "signed and notarized"
  else
    echo "signed, not notarized (set APPLE_ID, APPLE_TEAM_ID and APPLE_APP_PASSWORD to notarize)"
  fi
else
  # Ad-hoc signature: Apple silicon refuses to run unsigned code, and lipo drops the linker's.
  codesign --force --deep --sign - "$stage/GGV.app" "$stage/ggv"
  pack
  echo "ad-hoc signed (set MACOS_SIGN_IDENTITY to sign for distribution)"
fi

file "$app/MacOS/ggv-desktop" "$stage/ggv"
echo "packed $zip"
