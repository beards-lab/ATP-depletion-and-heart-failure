<#
.SYNOPSIS
    One-time setup of the private full-data repository and its local share clone.

.DESCRIPTION
    Creates a SEPARATE CLONE (not a worktree) whose only remote is the private
    repository, then pushes two branches to it:

        code   pure mirror of the public repo's Further-optim, no data
        main   code + the data payload + the sync workflow  <- what peers clone

    A separate clone rather than a worktree is deliberate: worktrees share
    .git/config, so a worktree cannot have its own 'origin', and GitHub Desktop
    (which only drives 'origin') would push the data commit to the PUBLIC repo.
    Here the public remote is removed outright, so that is impossible.

    Read scripts/fulldata/README.md first. Runs as a DRY RUN unless -Execute.

    NOTE: this file is deliberately pure ASCII. Windows PowerShell 5.1 reads
    .ps1 as ANSI unless there is a BOM, so non-ASCII characters break parsing.

.PARAMETER Execute
    Actually clone, commit and push. Without it, nothing is written or pushed.

.EXAMPLE
    .\scripts\fulldata\Seed-FullDataRepo.ps1
    .\scripts\fulldata\Seed-FullDataRepo.ps1 -Execute
#>
[CmdletBinding()]
param(
    [string] $ShareDir    = 'C:\home\git\ATP-share-fulldata',
    [string] $SourceRepo  = 'C:\home\git\ATP-depletion-and-heart-failure',
    [string] $PublicUrl   = 'https://github.com/beards-lab/ATP-depletion-and-heart-failure.git',
    [string] $PrivateUrl  = 'https://github.com/beards-lab/ATP-depletion-and-heart-failure-full-data.git',
    [string] $Branch      = 'Further-optim',
    [string] $CloneFrom   = '',
    [switch] $SkipPush,
    [switch] $Yes,
    [switch] $Execute
)

# Clone source. Default is the public GitHub URL (anonymous clone, no auth
# needed). Pass -CloneFrom <local repo path> to seed from a local checkout
# instead - use that when the commits you want are not yet pushed to GitHub.
if (-not $CloneFrom) { $CloneFrom = $PublicUrl }

$ErrorActionPreference = 'Stop'
$here = Split-Path -Parent $MyInvocation.MyCommand.Path

function Say  { param($m) Write-Host $m }
function Step { param($m) Write-Host "`n== $m" -ForegroundColor Cyan }
function Run {
    param([string[]] $GitArgs, [string] $WorkDir)
    $shown = 'git ' + ($GitArgs -join ' ')
    if ($WorkDir) { $shown = "[$WorkDir] $shown" }
    if (-not $Execute) { Say "  DRYRUN  $shown"; return }
    Say "  $shown"
    if ($WorkDir) { & git -C $WorkDir @GitArgs } else { & git @GitArgs }
    if ($LASTEXITCODE -ne 0) { throw "git failed ($LASTEXITCODE): $shown" }
}

if (-not $Execute) {
    Write-Host "`n*** DRY RUN - nothing will be written or pushed. Re-run with -Execute. ***" -ForegroundColor Yellow
}

# --- Resolve the payload from the manifest --------------------------------
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
if (-not $files) { throw "Manifest matched no files under $SourceRepo. Is the data present?" }

# Explicit relative paths, not globs, are what gets staged: git then adds
# exactly the set printed here, and the '!' exclusions cannot leak back in.
$relPaths = $files | ForEach-Object { 'data/' + $_.Name }

$totalMB = [math]::Round((($files | Measure-Object Length -Sum).Sum / 1MB), 2)
Say ("  {0} files, {1} MB" -f $files.Count, $totalMB)
Say '  include:'
$globs | ForEach-Object { Say "    $_" }
if ($excludes) { Say '  exclude:' ; $excludes | ForEach-Object { Say "    $_" } }

# A wrong payload cannot be un-published: force-pushes rewrite refs but GitHub
# keeps unreachable objects around. So refuse anything implausibly large here
# rather than after the fact.
if ($totalMB -gt 50) {
    throw ("Payload is {0} MB - far larger than the expected ~2.5 MB. Check payload-manifest.txt before continuing." -f $totalMB)
}

# --- Preconditions --------------------------------------------------------
Step 'Preconditions'
if (Test-Path $ShareDir) {
    throw "$ShareDir already exists. This script is for FIRST-TIME setup; use Update-FullDataPayload.ps1 instead."
}
Say "  share clone target: $ShareDir  (does not exist yet - good)"
Say "  private remote:     $PrivateUrl"

# --- 1. Clone, then sever the public remote -------------------------------
Step "1. Clone $Branch from: $CloneFrom"
# --no-hardlinks: when seeding from a local path, git would otherwise share
# object files between the two repositories. This clone is pushed elsewhere and
# must stand on its own.
Run @('clone', '--no-hardlinks', '--branch', $Branch, $CloneFrom, $ShareDir)

Step '2. Remove the source remote and point origin at the PRIVATE repository'
# This is the safety property of the whole design: after this, the share clone
# has no route to the public repo at all, so GitHub Desktop cannot leak data.
Run @('remote', 'remove', 'origin') $ShareDir
Run @('remote', 'add', 'origin', $PrivateUrl) $ShareDir

# --- 3. Create the code mirror branch -------------------------------------
Step '3. Create the code mirror branch (no data on it)'
Run @('branch', '-f', 'code', $Branch) $ShareDir
if (-not $SkipPush) {
    Run @('push', '-u', 'origin', 'code:code') $ShareDir
} else {
    Say '  (-SkipPush) local branch only; push it from GitHub Desktop'
}

# --- 4. Build main = code + workflow + data -------------------------------
Step '4. Create main and add the sync workflow'
Run @('checkout', '-B', 'main', $Branch) $ShareDir

$wfDir = Join-Path $ShareDir '.github\workflows'
if ($Execute) {
    New-Item -ItemType Directory -Force -Path $wfDir | Out-Null
    Copy-Item (Join-Path $here 'sync-code.yml') (Join-Path $wfDir 'sync-code.yml') -Force
}
Say '  copy sync-code.yml -> .github/workflows/sync-code.yml'

Step '5. Copy the data payload in'
if ($Execute) {
    New-Item -ItemType Directory -Force -Path (Join-Path $ShareDir 'data') | Out-Null
    foreach ($f in $files) {
        Copy-Item $f.FullName (Join-Path $ShareDir ('data\' + $f.Name)) -Force
    }
}
Say ("  {0} files copied into {1}\data" -f $files.Count, $ShareDir)

Step '6. Commit'
# -f is required: .gitignore excludes both 'data/' (line 28) and '*.mat' (line 6).
# NOTE the explicit per-file list is deliberate - do NOT substitute a bare
# 'data/' here, and do not go back to globs (they would re-add the exclusions).
$addArgs = @('add', '-f') + $relPaths + @('.github/workflows/sync-code.yml')
Run $addArgs $ShareDir
Run @('commit', '-m', 'FULL-DATA: curated protocol datasets (not for public release) + sync workflow') $ShareDir

Step '7. Review before pushing'
Run @('show', '--stat', '--oneline', 'HEAD') $ShareDir

if (-not $Execute) {
    Write-Host "`nDRY RUN complete. Re-run with -Execute to perform the above." -ForegroundColor Yellow
    return
}

if ($SkipPush) {
    Write-Host "`nLocal setup complete - NOTHING WAS PUSHED (-SkipPush)." -ForegroundColor Green
    Write-Host "The share clone is ready at: $ShareDir" -ForegroundColor Green
    Write-Host "Its only remote is the PRIVATE repo, so pushing from there is safe:" -ForegroundColor Green
    Write-Host "  $PrivateUrl" -ForegroundColor Green
    Write-Host "`nIn GitHub Desktop: File -> Add local repository -> $ShareDir" -ForegroundColor Green
    Write-Host "then publish/push BOTH branches: 'code' and 'main'." -ForegroundColor Green
    return
}

Write-Host "`nAbout to PUSH the data payload to:" -ForegroundColor Yellow
Write-Host "  $PrivateUrl  (branch main)" -ForegroundColor Yellow
if (-not $Yes) {
    $ans = Read-Host 'Type PUSH to continue, anything else to stop'
    if ($ans -ne 'PUSH') { Write-Host 'Stopped. Nothing pushed.'; return }
}

Run @('push', '-u', 'origin', 'main') $ShareDir

Write-Host "`nDone. Two things left, both in the private repo's web UI:" -ForegroundColor Green
Write-Host "  1. Settings -> General -> Default branch:  set it to 'main'" -ForegroundColor Green
Write-Host '  2. Settings -> Collaborators:              add your peers' -ForegroundColor Green
Write-Host "`nAfter that the routine update is browser-only:" -ForegroundColor Green
Write-Host "  Actions -> 'Sync code from public repo' -> Run workflow" -ForegroundColor Green
Write-Host '  Pull requests -> merge the PR it opens' -ForegroundColor Green
