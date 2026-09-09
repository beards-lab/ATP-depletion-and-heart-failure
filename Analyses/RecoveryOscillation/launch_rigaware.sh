#!/usr/bin/env bash
# launch_rigaware.sh -- start the rig-aware campaign on a Linux box or cluster
# node. The POSIX twin of launch_rigaware.ps1; see that file's header for the
# reasoning, and RunOptimRigAware.m for what each profile does.
#
#   ./launch_rigaware.sh                  # dry run: pre-flight all four
#   GO=1 HOURS=30 ./launch_rigaware.sh    # the real thing
#   GO=1 HOURS=30 PROFILES="a c" ./launch_rigaware.sh
#
# Under SLURM, ask for one node with >= 20 cores and >= 32 GB and run this as
# the job step -- four MATLAB sessions x ~5 workers each. Do NOT submit the
# four as separate jobs unless they share a filesystem: they write to the same
# params/ directory (different tags, so no collision, but the consolidation
# report expects to find all four).
#
# Logs -> results/rigaware_<p>.log

set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"
res="$here/results"; mkdir -p "$res"

MATLAB="${MATLAB:-matlab}"
HOURS="${HOURS:-30}"
WORKERS="${WORKERS:-5}"
PROFILES="${PROFILES:-a b c d}"
if [ "${GO:-0}" = "1" ]; then DRY=false; else DRY=true
  echo "DRY RUN: pre-flighting only. Set GO=1 to start the search." >&2
fi

command -v "$MATLAB" >/dev/null || { echo "matlab not on PATH; set MATLAB=<path>" >&2; exit 1; }

pids=()
for p in $PROFILES; do
    log="$res/rigaware_$p.log"
    echo "profile $p -> $log"
    "$MATLAB" -nodisplay -nosplash -batch \
      "cd('$root'); addpath(genpath('.')); PROFILE='$p'; TIME_HRS=$HOURS; \
       NWORKERS=$WORKERS; DRYRUN=$DRY; \
       run('Analyses/RecoveryOscillation/RunOptimRigAware.m')" \
      > "$log" 2>&1 &
    pids+=($!)
    sleep 20      # stagger: see launch_rigaware.ps1
done

echo
echo "launched ${#pids[@]} session(s): ${pids[*]}"
echo "follow one with:  tail -f $res/rigaware_a.log"
wait "${pids[@]}"
echo "all sessions finished. consolidate with:"
echo "  $MATLAB -batch \"run('Analyses/RecoveryOscillation/RunRigAwareReport.m')\""
