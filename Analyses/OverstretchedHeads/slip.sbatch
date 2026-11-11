#!/bin/bash
# slip.sbatch -- the A1-slip campaign as one SLURM job (all arms on one node).
#
#   sbatch --account=<acct> Analyses/OverstretchedHeads/slip.sbatch            # submit from the REPO ROOT
#   sbatch --account=<acct> --export=ALL,HOURS=20,PROFILES="a b" Analyses/OverstretchedHeads/slip.sbatch
#
# Three arms x 5 workers + 3 clients = 18 busy processes -> ask for 20 cpus.
# --time = HOURS + ~2 h margin (startup, pre-flight, last round finishing past
# the budget: the budget is checked BETWEEN rounds and a round is 30 evaluations).
#SBATCH --job-name=slip
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=40g
#SBATCH --time=23:00:00
#SBATCH --output=Analyses/OverstretchedHeads/results/slip_slurm_%j.out
#SBATCH --error=Analyses/OverstretchedHeads/results/slip_slurm_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
module load "${MATLAB_MODULE:-matlab}" 2>/dev/null || true
echo "host $(hostname) | job ${SLURM_JOB_ID:-none} | ${SLURM_CPUS_PER_TASK:-?} cpus | commit $(git rev-parse --short HEAD 2>/dev/null || echo '?')"
echo "=== pre-flight ==="
PROFILES="${PROFILES:-a b c}" bash Analyses/OverstretchedHeads/launch_slip.sh
echo "=== campaign ==="
GO=1 HOURS="${HOURS:-20}" WORKERS="${WORKERS:-5}" PROFILES="${PROFILES:-a b c}" bash Analyses/OverstretchedHeads/launch_slip.sh
