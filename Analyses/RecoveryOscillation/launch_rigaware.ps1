# launch_rigaware.ps1 -- start the rig-aware campaign on a Windows box.
#
# Four independent MATLAB sessions, one per search profile, each with its own
# tag so they never touch each other's state, snapshot or iteration gallery
# (that isolation is what optimizeFeatures' cfg.tag is for).
#
#   .\launch_rigaware.ps1                     # dry run: pre-flight all four, no search
#   .\launch_rigaware.ps1 -Hours 30 -Go       # the real thing
#   .\launch_rigaware.ps1 -Hours 30 -Go -Profiles a,c
#
# -Campaign picks what to start from (campaignSpec.m). It defaults to
# `floored`, which is the one to run. `rigaware` already ran and is kept only
# so its result reproduces -- conclusions.md sec. 7 is the record of it putting
# 110 % of its gain into force-velocity.
#
# Each session parallelises INSIDE one objective evaluation over the five slack
# segments, so it saturates at ~5-6 workers; four sessions want ~20 cores and
# ~4 x 6 GB. Resume is automatic -- re-running the same command after a crash
# or a reboot continues from params/rigaware_<p>_state.mat.
#
# Logs -> Analyses/RecoveryOscillation/results/rigaware_<p>.log
# Then: matlab -batch "run('Analyses/RecoveryOscillation/RunRigAwareReport.m')"

param(
    [double]   $Hours     = 30,
    [int]      $Workers   = 5,
    [string[]] $Profiles  = @('a','b','c','d'),
    [ValidateSet('floored','rigaware')]
    [string]   $Campaign  = 'floored',
    [switch]   $Go,
    [string]   $Matlab    = "C:\Program Files\MATLAB\R2023a\bin\matlab.exe"
)

$root = Resolve-Path (Join-Path $PSScriptRoot "..\..")
$res  = Join-Path $PSScriptRoot "results"
if (-not (Test-Path $res)) { New-Item -ItemType Directory $res | Out-Null }

if (-not (Test-Path $Matlab)) { throw "MATLAB not found at $Matlab -- pass -Matlab <path>" }

$dry = if ($Go) { "false" } else { "true" }
if (-not $Go) {
    Write-Host "DRY RUN: pre-flighting $($Profiles.Count) profile(s). Add -Go to start the search." -ForegroundColor Yellow
}

Write-Host "campaign $Campaign"

foreach ($p in $Profiles) {
    $log = Join-Path $res "${Campaign}_$p.log"
    $cmd = "cd('$root'); addpath(genpath('.')); CAMPAIGN='$Campaign'; PROFILE='$p'; " +
           "TIME_HRS=$Hours; NWORKERS=$Workers; DRYRUN=$dry; " +
           "run('Analyses/RecoveryOscillation/RunOptimRigAware.m')"

    Write-Host "profile $p -> $log"
    Start-Process -FilePath $Matlab `
        -ArgumentList @('-batch', $cmd) `
        -RedirectStandardOutput $log `
        -RedirectStandardError  "$log.err" `
        -WorkingDirectory $root `
        -NoNewWindow:$false | Out-Null

    # Stagger: four pools starting at the same instant fight over cores, and
    # optimizeFeatures calls rng('shuffle') per round -- simultaneous starts
    # could draw the same first subset.
    Start-Sleep -Seconds 20
}

Write-Host ""
Write-Host "launched. follow one with:  Get-Content -Wait '$res\${Campaign}_a.log'"
Write-Host "check progress:             matlab -batch `"addpath(genpath('.')); campaignStatus('$Campaign')`""
Write-Host "consolidate when done:      matlab -batch `"CAMPAIGN='$Campaign'; run('Analyses/RecoveryOscillation/RunRigAwareReport.m')`""
