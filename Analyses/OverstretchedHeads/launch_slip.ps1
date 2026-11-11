# launch_slip.ps1 -- run the A1-slip campaign arms in parallel on a Windows box.
# The twin of launch_slip.sh; see that file for sizing (one session per arm,
# ~5 workers each: 16 cores = 3 arms) and RunSlipCampaign.m for the arms.
#
#   .\Analyses\OverstretchedHeads\launch_slip.ps1                    # dry run: pre-flight every arm
#   .\Analyses\OverstretchedHeads\launch_slip.ps1 -Hours 20 -Go      # the real thing
#   .\Analyses\OverstretchedHeads\launch_slip.ps1 -Profiles a -Smoke # one tiny round, then clean up
#
# Resume is automatic: re-run the same command to continue from params/slip_<p>_state.mat.
param(
    [double]   $Hours    = 20,
    [int]      $Workers  = 5,
    [string[]] $Profiles = @('a','b','c'),
    [switch]   $Go,
    [switch]   $Smoke,
    [string]   $Matlab   = "C:\Program Files\MATLAB\R2023a\bin\matlab.exe"
)
$root = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$res  = Join-Path $PSScriptRoot "results"
if (-not (Test-Path $Matlab)) { throw "MATLAB not found at $Matlab -- pass -Matlab <path>" }
$mode = if ($Smoke) { "SMOKE=true;" } elseif ($Go) { "" } else { "DRYRUN=true;" }
if (-not $Go -and -not $Smoke) { Write-Host "DRY RUN: pre-flighting only. Add -Go to start the search." -ForegroundColor Yellow }

$procs = @()
foreach ($p in $Profiles) {
    $log = if ($Smoke) { Join-Path $res "slip_${p}_smoke.log" } else { Join-Path $res "slip_$p.log" }
    $cmd = "cd('$root'); addpath(genpath('.')); PROFILE='$p'; TIME_HRS=$Hours; NWORKERS=$Workers; $mode " +
           "run('Analyses/OverstretchedHeads/RunSlipCampaign.m')"
    Write-Host "arm $p -> $log"
    $procs += Start-Process -FilePath $Matlab -ArgumentList @('-batch', "`"$cmd`"") -RedirectStandardOutput $log `
              -RedirectStandardError "$log.err" -NoNewWindow -PassThru
    Start-Sleep -Seconds 30      # stagger the pool starts
}
$procs | Wait-Process
foreach ($p in $Profiles) {
    $log = if ($Smoke) { Join-Path $res "slip_${p}_smoke.log" } else { Join-Path $res "slip_$p.log" }
    Write-Host "--- ${p}:"; Select-String -Path $log -Pattern 'PRE-FLIGHT OK','here vs','IMPROVED','CAMPAIGN DONE','rror' | Select-Object -Last 5 | ForEach-Object { $_.Line }
}
if ($Smoke) { Remove-Item -Recurse -Force (Join-Path $root "params\slip_*_smoke*") -ErrorAction SilentlyContinue; Write-Host "smoke artefacts removed" }
