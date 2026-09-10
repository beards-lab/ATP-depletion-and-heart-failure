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

# On an HPC system MATLAB is almost never on the default PATH -- it lives
# behind the module system, and `module` is a shell function that a
# non-interactive script cannot rely on inheriting. So say what to do rather
# than just what is wrong.
if ! command -v "$MATLAB" >/dev/null; then
    {
        echo "matlab not on PATH."
        if command -v module >/dev/null 2>&1 || [ -n "${MODULESHOME:-}" ]; then
            echo
            echo "This looks like a module-based system. Load it first, in THIS shell:"
            echo "    module avail matlab        # see what is installed"
            echo "    module load matlab         # or e.g. matlab/R2023a"
            echo "    GO=1 HOURS=$HOURS bash $0"
        fi
        echo
        echo "Or point at it directly:  MATLAB=/path/to/matlab GO=1 bash $0"
    } >&2
    exit 1
fi

# Record the release in the log. The campaign was developed on R2023a; the
# selector syntax it depends on (rsK80[1:4]) is plain regex and string work, so
# older releases are expected to be fine -- but if a run behaves oddly this is
# the first line to check, and it costs one banner to have it.
echo "matlab   $(command -v "$MATLAB")"
"$MATLAB" -nodisplay -nosplash -batch "fprintf('release  %s\n', version('-release'))" 2>/dev/null || true

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
