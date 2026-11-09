<#
.SYNOPSIS
  Adds "Open in GGV" to the right-click menu of genomics files, and makes ggv:// links open GGV,
  for the current user only.

.DESCRIPTION
  Writes HKCU\Software\Classes\SystemFileAssociations\<ext>\shell\ggv for each extension, and
  HKCU\Software\Classes\ggv for ggv:// links (from Session > Copy link, or a notebook). This
  adds a menu entry; it does not change which program opens a file on double-click, and needs no
  administrator rights. Run it again after moving ggv-desktop.exe; -Unregister removes the entries.

.EXAMPLE
  powershell -ExecutionPolicy Bypass -File register-file-types.ps1
  powershell -ExecutionPolicy Bypass -File register-file-types.ps1 -DryRun
  powershell -ExecutionPolicy Bypass -File register-file-types.ps1 -Unregister
#>
param(
    # ggv-desktop.exe to register; defaults to the one next to this script.
    [string]$Exe = (Join-Path $PSScriptRoot 'ggv-desktop.exe'),
    # Also add the entry to every .gz file (for .vcf.gz, .gtf.gz and .bed.gz).
    [switch]$IncludeGz,
    [switch]$Unregister,
    # Print what would be written or removed, and change nothing.
    [switch]$DryRun
)
$ErrorActionPreference = 'Stop'

$extensions = @('.bam', '.cram', '.bcf', '.vcf', '.gtf', '.gff', '.gff3', '.bed')
if ($IncludeGz) { $extensions += '.gz' }

if (-not $Unregister) {
    if (-not (Test-Path $Exe)) { throw "ggv-desktop.exe not found at $Exe (pass -Exe)" }
    $Exe = (Resolve-Path $Exe).Path
}

foreach ($ext in $extensions) {
    $key = "HKCU:\Software\Classes\SystemFileAssociations\$ext\shell\ggv"
    if ($Unregister) {
        if (Test-Path $key) {
            if ($DryRun) { Write-Output "would remove $key" }
            else { Remove-Item -Recurse -Force $key; Write-Output "removed $key" }
        }
        continue
    }
    $command = "`"$Exe`" `"%1`""
    if ($DryRun) {
        Write-Output "would write $key = 'Open in GGV', command: $command"
        continue
    }
    New-Item -Force -Path "$key\command" | Out-Null
    Set-ItemProperty -Path $key -Name '(Default)' -Value 'Open in GGV'
    Set-ItemProperty -Path $key -Name 'Icon' -Value "`"$Exe`""
    Set-ItemProperty -Path "$key\command" -Name '(Default)' -Value $command
    Write-Output "registered $ext"
}

# ggv:// links.
$proto = 'HKCU:\Software\Classes\ggv'
if ($Unregister) {
    if (Test-Path $proto) {
        if ($DryRun) { Write-Output "would remove $proto" }
        else { Remove-Item -Recurse -Force $proto; Write-Output "removed $proto" }
    }
}
elseif ($DryRun) {
    Write-Output "would write $proto = 'URL:ggv link', command: `"$Exe`" `"%1`""
}
else {
    New-Item -Force -Path "$proto\shell\open\command" | Out-Null
    Set-ItemProperty -Path $proto -Name '(Default)' -Value 'URL:ggv link'
    Set-ItemProperty -Path $proto -Name 'URL Protocol' -Value ''
    Set-ItemProperty -Path "$proto\shell\open\command" -Name '(Default)' -Value "`"$Exe`" `"%1`""
    Write-Output "registered ggv:// links"
}
