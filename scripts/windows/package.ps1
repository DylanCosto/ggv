<#
.SYNOPSIS
  Builds GGV in release mode and packs ggv.exe, ggv-desktop.exe, install.ps1, the file-type
  script, README and LICENSE into dist\ggv-<version>-windows-x64.zip.

.EXAMPLE
  powershell -ExecutionPolicy Bypass -File scripts\windows\package.ps1
#>
param(
    # Skip the build and pack what is already in the target directory.
    [switch]$NoBuild
)
$ErrorActionPreference = 'Stop'
$root = Resolve-Path (Join-Path $PSScriptRoot '..\..')
Push-Location $root
try {
    if (-not $NoBuild) {
        cargo build --release --locked -p ggv-cli -p ggv-desktop
        if ($LASTEXITCODE -ne 0) { throw "cargo build failed" }
    }
    $target = if ($env:CARGO_TARGET_DIR) { $env:CARGO_TARGET_DIR } else { Join-Path $root 'target' }
    $version = (Select-String -Path (Join-Path $root 'Cargo.toml') -Pattern '^version = "(.+)"' | Select-Object -First 1).Matches[0].Groups[1].Value
    $name = "ggv-$version-windows-x64"
    $dist = Join-Path $root 'dist'
    $stage = Join-Path $dist $name
    if (Test-Path $stage) { Remove-Item -Recurse -Force $stage }
    New-Item -ItemType Directory -Force $stage | Out-Null
    foreach ($exe in 'ggv.exe', 'ggv-desktop.exe') {
        Copy-Item (Join-Path $target "release\$exe") $stage
    }
    Copy-Item (Join-Path $root 'README.md') $stage
    Copy-Item (Join-Path $root 'LICENSE') $stage
    Copy-Item (Join-Path $root 'THIRD_PARTY_NOTICES.md') $stage
    Copy-Item (Join-Path $PSScriptRoot 'register-file-types.ps1') $stage
    Copy-Item (Join-Path $PSScriptRoot 'install.ps1') $stage
    $zip = Join-Path $dist "$name.zip"
    if (Test-Path $zip) { Remove-Item -Force $zip }
    # The folder itself, so unzipping makes one ggv-<version>-windows-x64 folder rather than
    # scattering its files wherever the zip is.
    Compress-Archive -Path $stage -DestinationPath $zip
    Write-Output "packed $zip"
}
finally {
    Pop-Location
}
