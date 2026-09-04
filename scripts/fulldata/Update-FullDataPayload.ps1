<#
.SYNOPSIS
    Refresh the data payload in the private full-data repository.

.DESCRIPTION
    Run this when a curated .mat under data/ has been regenerated (e.g. after a
    DataCuration/ pass) and colleagues need the new version.

    It copies the manifest files from the working repository into the share
    clone, shows what changed, and commits + pushes them to 'main'.

    It does NOT touch code. Code reaches the private repo through the GitHub
    Actions sync workflow and a 'code -> main' pull request - see
    scripts/fulldata/README.md.

    Runs as a DRY RUN unless -Execute.

    NOTE: this file is deliberately pure ASCII. Windows PowerShell 5.1 reads
    .ps1 as ANSI unless there is a BOM, so non-ASCII characters break parsing.

.EXAMPLE
    .\scripts\fulldata\Update-FullDataPayload.ps1
    .\scripts\fulldata\Update-FullDataPayload.ps1 -Execute
#>
[CmdletBinding()]
param(
    [string] $ShareDir   = 'C:\home\git\ATP-share-fulldata',
    [string] $SourceRepo = 'C:\home\git\ATP-depletion-and-heart-failure',
    [switch] $Execute
)

$ErrorActionPreference = 'Stop'
$here = Split-Path -Parent $MyInvocation.MyCommand.Path

function Say  { param($m) Write-Host $m }
function Step { param($m) Write-Host "`n== $m" -ForegroundColor Cyan }
function Run {
    param([string[]] $GitArgs, [string] $WorkDir = $ShareDir)
    $shown = "[$WorkDir] git " + ($GitArgs -join ' ')
    if (-not $Execute) { Say "  DRYRUN  $shown"; return }
    Say "  $shown"
    & git -C $WorkDir @GitArgs
    if ($LASTEXITCODE -ne 0) { throw "git failed ($LASTEXITCODE): $shown" }
}

if (-not $Execute) {
    Write-Host "`n*** DRY RUN - nothing will be written or pushed. Re-run with -Execute. ***" -ForegroundColor Yellow
}

# --- Preconditions --------------------------------------------------------
Step 'Preconditions'
if (-not (Test-Path $ShareDir)) {
    throw "$ShareDir not found. Run Seed-FullDataRepo.ps1 first."
}

# The share clone must be pointed ONLY at the private repository. If a public
# remote ever reappears here, stop: pushing data through it would publish it.
$remotes = & git -C $ShareDir remote -v
$publicHit = $remotes | Where-Object { $_ -match 'ATP-depletion-and-heart-failure\.git' }
if ($publicHit) {
    Write-Host 'The share clone has a remote pointing at the PUBLIC repository:' -ForegroundColor Red
    $publicHit | ForEach-Object { Write-Host "  $_" -ForegroundColor Red }
    throw "Remove it before continuing:  git -C '$ShareDir' remote remove THATNAME"
}
Say '  remotes (private only):'
$remotes | ForEach-Object { Say "    $_" }

Run @('checkout', 'main')
Run @('pull', '--ff-only', 'origin', 'main')

# --- Payload --------------------------------------------------------------
Step 'Payload'
$lines = Get-Content (Join-Path $here 'payload-manifest.txt') |
         ForEach-Object { $_.Trim() } |
         Where-Object { $_ -and -not $_.StartsWith('#') }
$globs    = @($lines | Where-Object { -not $_.StartsWith('!') })
$excludes = @($lines | Where-Object { $_.StartsWith('!') } | ForEach-Object { $_.Substring(1).Trim() })

$files = @()
foreach ($g in $globs) {
    $files += Get-ChildItem -Path (Join-Path $SourceRepo $g) -File -ErrorAction SilentlyContinue
}
foreach ($x in $excludes) {
    $drop = @(Get-ChildItem -Path (Join-Path $SourceRepo $x) -File -ErrorAction SilentlyContinue)
    if ($drop) {
        $dropNames = $drop | ForEach-Object { $_.FullName }
        $files = $files | Where-Object { $dropNames -notcontains $_.FullName }
    }
}
$files = @($files | Sort-Object FullName -Unique)
if (-not $files) { throw "Manifest matched no files under $SourceRepo." }

$totalMB = [math]::Round((($files | Measure-Object Length -Sum).Sum / 1MB), 2)
Say ("  {0} files, {1} MB" -f $files.Count, $totalMB)
if ($excludes) { Say '  excluded by manifest:' ; $excludes | ForEach-Object { Say "    $_" } }
if ($totalMB -gt 50) {
    throw ("Payload is {0} MB - far larger than the expected ~2.5 MB. Check payload-manifest.txt before continuing." -f $totalMB)
}

Step 'Copy'
$changed  = 0
$relPaths = @()
if ($Execute) { New-Item -ItemType Directory -Force -Path (Join-Path $ShareDir 'data') | Out-Null }
foreach ($f in $files) {
    $relPaths += 'data/' + $f.Name
    $dest = Join-Path $ShareDir ('data\' + $f.Name)
    $differs = $true
    if (Test-Path $dest) {
        $a = (Get-FileHash $f.FullName -Algorithm SHA256).Hash
        $b = (Get-FileHash $dest       -Algorithm SHA256).Hash
        $differs = ($a -ne $b)
    }
    if ($differs) {
        $changed++
        Say ("  {0}" -f $f.Name)
        if ($Execute) { Copy-Item $f.FullName $dest -Force }
    }
}
if ($changed -eq 0) {
    Write-Host "`nNothing changed - the payload is already current." -ForegroundColor Green
    return
}
Say ("  {0} file(s) differ" -f $changed)

# --- Commit + push --------------------------------------------------------
Step 'Stage'
# -f is required: .gitignore excludes 'data/' (line 28) and '*.mat' (line 6).
# Explicit per-file paths, not globs, so the '!' exclusions cannot leak back in.
# NB this stages additions and modifications only. If a .mat was DELETED from
# the source repo, remove it from the share clone by hand (git rm) - deliberate,
# so a half-populated data/ can never silently delete a colleague's file.
Run (@('add', '-f') + $relPaths)

Step 'Review'
Run @('status', '--short')
Run @('diff', '--cached', '--stat')

if (-not $Execute) {
    Write-Host "`nDRY RUN complete. Re-run with -Execute to perform the above." -ForegroundColor Yellow
    return
}

Write-Host "`nAbout to COMMIT and PUSH the above to the private repo's main." -ForegroundColor Yellow
$ans = Read-Host 'Type PUSH to continue, anything else to stop'
if ($ans -ne 'PUSH') { Write-Host 'Stopped. Nothing pushed.'; return }

$stamp = Get-Date -Format 'yyyy-MM-dd'
Run @('commit', '-m', "FULL-DATA: refresh curated protocol datasets ($stamp)")
Run @('push', 'origin', 'main')

Write-Host "`nDone. Collaborators get it with a plain 'git pull'." -ForegroundColor Green
