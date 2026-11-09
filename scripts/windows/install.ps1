<#
.SYNOPSIS
  Installs GGV from the unzipped download, for the current user (no administrator rights).

.DESCRIPTION
  Copies ggv.exe and ggv-desktop.exe to %LOCALAPPDATA%\Programs\ggv, adds that folder to your
  user PATH (so `ggv` works in new terminals), adds a Start menu shortcut, and adds
  "Open in GGV" to the right-click menu of genomics files (register-file-types.ps1).
  -Uninstall removes all of it. Running it again updates an earlier install.

  The files it copies are listed, with their SHA-256, in installed.txt in that folder. A file
  already there that GGV didn't install is left alone (with a warning); -Uninstall removes only
  listed files that are unchanged since, and then the folder if nothing else is in it.

.EXAMPLE
  powershell -ExecutionPolicy Bypass -File install.ps1
  powershell -ExecutionPolicy Bypass -File install.ps1 -Uninstall
#>
param(
    [switch]$Uninstall,
    # Leave the right-click menu entries alone.
    [switch]$NoFileTypes
)
$ErrorActionPreference = 'Stop'

$dest = Join-Path $env:LOCALAPPDATA 'Programs\ggv'
$files = 'ggv.exe', 'ggv-desktop.exe', 'register-file-types.ps1', 'install.ps1', 'README.md', 'LICENSE', 'THIRD_PARTY_NOTICES.md'
$manifest = Join-Path $dest 'installed.txt'
$desktopExe = Join-Path $dest 'ggv-desktop.exe'
$shortcut = Join-Path ([Environment]::GetFolderPath('Programs')) 'GGV.lnk'
$userPath = [Environment]::GetEnvironmentVariable('Path', 'User')
$parts = @($userPath -split ';' | Where-Object { $_ -and $_.TrimEnd('\') -ne $dest.TrimEnd('\') })

function Get-Sha([string]$path) { (Get-FileHash -Algorithm SHA256 -LiteralPath $path).Hash.ToLowerInvariant() }

# What installed.txt lists: file name -> SHA-256.
function Get-Installed {
    $m = @{}
    if (Test-Path -LiteralPath $manifest -PathType Leaf) {
        foreach ($line in Get-Content -LiteralPath $manifest) {
            $hash, $name = $line -split '  ', 2
            if ($files -contains $name) { $m[$name] = $hash }
        }
    }
    $m
}

# A listed file that is still exactly what was installed.
function Test-Ours([hashtable]$installed, [string]$name) {
    $p = Join-Path $dest $name
    $installed.ContainsKey($name) -and (Test-Path -LiteralPath $p -PathType Leaf) -and ((Get-Sha $p) -eq $installed[$name])
}

# The Start menu shortcut, if it points at the installed ggv-desktop.exe.
function Test-OurShortcut {
    if (-not (Test-Path -LiteralPath $shortcut)) { return $false }
    $target = (New-Object -ComObject WScript.Shell).CreateShortcut($shortcut).TargetPath
    $target -eq $desktopExe
}

if ($Uninstall) {
    $installed = Get-Installed
    if (Test-Ours $installed 'register-file-types.ps1') {
        & (Join-Path $dest 'register-file-types.ps1') -Unregister
    }
    [Environment]::SetEnvironmentVariable('Path', ($parts -join ';'), 'User')
    if (Test-OurShortcut) { Remove-Item -LiteralPath $shortcut }
    foreach ($name in @($installed.Keys)) {
        $p = Join-Path $dest $name
        if (-not (Test-Path -LiteralPath $p)) { continue }
        if (Test-Ours $installed $name) {
            Remove-Item -LiteralPath $p
        } else {
            Write-Warning "kept $p`: it changed after GGV installed it"
        }
    }
    if (Test-Path -LiteralPath $manifest) { Remove-Item -LiteralPath $manifest }
    # Only its own files: the folder goes only if nothing else is in it.
    if ((Test-Path -LiteralPath $dest) -and -not (Get-ChildItem -LiteralPath $dest -Force)) {
        Remove-Item -LiteralPath $dest
    }
    Write-Output "removed GGV from $dest, your PATH and the Start menu"
    return
}

New-Item -ItemType Directory -Force $dest | Out-Null
$installed = Get-Installed
$lines = @()
$skipped = $false
foreach ($f in $files) {
    $src = Join-Path $PSScriptRoot $f
    if (-not (Test-Path -LiteralPath $src -PathType Leaf)) { continue }
    $to = Join-Path $dest $f
    $want = Get-Sha $src
    if ((Test-Path -LiteralPath $to) -and -not (Test-Ours $installed $f)) {
        $same = (Test-Path -LiteralPath $to -PathType Leaf) -and ((Get-Sha $to) -eq $want)
        if (-not $same) {
            Write-Warning "left $to alone: GGV didn't install it (move it away and run this again to install there)"
            $skipped = $true
            continue
        }
    }
    # Running this script from the install folder: it is already in place.
    if ((Resolve-Path -LiteralPath $src).Path -ne $to) { Copy-Item -LiteralPath $src $to -Force }
    $lines += "$want  $f"
}
Set-Content -LiteralPath $manifest -Value $lines -Encoding ascii
$installed = Get-Installed
if (-not ((Test-Ours $installed 'ggv.exe') -and (Test-Ours $installed 'ggv-desktop.exe'))) {
    Write-Warning "GGV was not installed (see above); your PATH and Start menu are unchanged"
    exit 1
}
[Environment]::SetEnvironmentVariable('Path', (($parts + $dest) -join ';'), 'User')

if ((Test-Path -LiteralPath $shortcut) -and -not (Test-OurShortcut)) {
    Write-Warning "left the Start menu shortcut $shortcut alone: it isn't GGV's"
} else {
    $link = (New-Object -ComObject WScript.Shell).CreateShortcut($shortcut)
    $link.TargetPath = $desktopExe
    $link.WorkingDirectory = $dest
    $link.Description = 'GGV genome viewer'
    $link.Save()
}

if (-not $NoFileTypes -and (Test-Ours $installed 'register-file-types.ps1')) {
    & (Join-Path $dest 'register-file-types.ps1') -Exe $desktopExe
}
Write-Output "installed GGV to $dest (new terminals have it on PATH; Start menu: GGV)"
& (Join-Path $dest 'ggv.exe') --version
if ($skipped) {
    Write-Warning "some files were left alone; see above"
    exit 1
}
