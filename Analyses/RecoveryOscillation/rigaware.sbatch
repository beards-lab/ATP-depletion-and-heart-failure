#!/bin/bash
# rigaware.sbatch -- the 30-hour campaign as a SLURM batch job (UM Great Lakes).
#
#   sbatch --account=<your_slurm_account> Analyses/RecoveryOscillation/rigaware.sbatch
#
# Submit from the REPO ROOT. Override anything without editing the file:
#   sbatch --account=beards1 --time=36:00:00 \
#          --export=ALL,HOURS=30,PROFILES="a c" Analyses/RecoveryOscillation/rigaware.sbatch
#
# CAMPAIGN defaults to `floored` -- the objective with a one-preparation-width
# deadband on every data term (conclusions.md sec. 9). The older `rigaware`
# campaign already ran and put 110 % of its gain into force-velocity (sec. 7);
# --export=ALL,CAMPAIGN=rigaware reproduces it, but do not relaunch it to make
# progress.
#
# WHY A BATCH JOB AND NOT THE LAUNCHER DIRECTLY. launch_rigaware.sh runs the
# four sessions in the FOREGROUND and waits. In an interactive salloc that ties
# the run to your SSH session -- drop the connection, or hit the allocation's
# wall clock, and 30 hours are gone with nothing written but the per-round
# snapshots. Batch survives logout. (The snapshots ARE the safety net: the
# search checkpoints params/<tag>_state.mat every round and resumes
# automatically, so a killed job loses at most one round.)
#
# SIZING. Four profiles x 5 workers = 20 busy cores; 24 leaves headroom for the
# four client processes. Parallelism is INSIDE one objective evaluation (over
# the five slack segments), so a session saturates at 5-6 workers -- more cores
# per session buys nothing, more sessions buys profiles.
#
# --time is the campaign budget PLUS margin: MATLAB startup, the pre-flight,
# and the final write. HOURS=30 under --time=32:00:00. The search stops itself
# at HOURS; SLURM killing it at --time is the ungraceful path.

#SBATCH --job-name=rigaware
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=48g
#SBATCH --time=32:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=Analyses/RecoveryOscillation/results/rigaware_slurm_%j.out
#SBATCH --error=Analyses/RecoveryOscillation/results/rigaware_slurm_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p Analyses/RecoveryOscillation/results

echo "host      $(hostname)"
echo "job       ${SLURM_JOB_ID:-none} on ${SLURM_CPUS_PER_TASK:-?} cpus"
echo "repo      $PWD"
echo "commit    $(git rev-parse --short HEAD 2>/dev/null || echo 'not a git checkout')"

module load "${MATLAB_MODULE:-matlab}"
export MATLAB="${MATLAB:-matlab}"
echo "matlab    $(command -v "$MATLAB")"

# Pre-flight the DATA before spending the allocation. data/ is gitignored, so a
# fresh clone has the code and the frozen target but none of the recordings,
# and the failure mode is a job that dies on its first evaluation.
echo
echo "=== data pre-flight ==="
"$MATLAB" -nodisplay -nosplash -batch \
  "addpath(genpath('.')); ok = checkCampaignData(); exit(~ok)" \
  || { echo "DATA PRE-FLIGHT FAILED -- copy the required files into data/ and resubmit" >&2; exit 2; }

CAMPAIGN="${CAMPAIGN:-floored}"

echo
echo "=== campaign $CAMPAIGN ==="
GO=1 \
CAMPAIGN="$CAMPAIGN" \
HOURS="${HOURS:-30}" \
WORKERS="${WORKERS:-5}" \
PROFILES="${PROFILES:-a b c d}" \
bash Analyses/RecoveryOscillation/launch_rigaware.sh

echo
echo "=== consolidating ==="
"$MATLAB" -nodisplay -nosplash -batch \
  "addpath(genpath('.')); CAMPAIGN='$CAMPAIGN'; run('Analyses/RecoveryOscillation/RunRigAwareReport.m')"

echo "done. results in Analyses/RecoveryOscillation/results/"
