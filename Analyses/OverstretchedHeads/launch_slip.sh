#!/usr/bin/env bash
# launch_slip.sh -- run the A1-slip campaign arms in parallel on one Linux box / node.
#
#   bash Analyses/OverstretchedHeads/launch_slip.sh                   # dry run: pre-flight every arm
#   GO=1 HOURS=20 bash Analyses/OverstretchedHeads/launch_slip.sh     # the real thing
#   GO=1 HOURS=20 PROFILES="a b" bash Analyses/OverstretchedHeads/launch_slip.sh
#   SMOKE=1 PROFILES=a bash Analyses/OverstretchedHeads/launch_slip.sh  # one tiny round, then clean up
#
# SIZING. One MATLAB session per arm, WORKERS (default 5) workers each. The
# parallelism is inside one objective evaluation (the five slack segments), so a
# session saturates at 5-6 workers: 16 cores run 3 arms (15 workers + 3 clients,
# slightly oversubscribed, fine), 20+ cores run them comfortably. More cores buy
# more arms, not faster arms. ~5 GB RAM per session.
#
# Resume is automatic: re-running the same command continues every arm from
# params/slip_<p>_state.mat. Logs -> Analyses/OverstretchedHeads/results/slip_<p>.log
# Afterwards: matlab -batch "SNAPS={'yieldA_opt','yieldB_opt','slip_a_opt','slip_b_opt','slip_c_opt'}; \
#             SKIP_OBJ=false; run('Analyses/OverstretchedHeads/EvaluateYieldCampaigns.m')"
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"
res="$here/results"; mkdir -p "$res"

MATLAB="${MATLAB:-matlab}"
HOURS="${HOURS:-20}"
WORKERS="${WORKERS:-5}"
PROFILES="${PROFILES:-a b c}"
if [ "${SMOKE:-0}" = "1" ]; then MODE="SMOKE=true;"; elif [ "${GO:-0}" = "1" ]; then MODE=""; else
    MODE="DRYRUN=true;"; echo "DRY RUN: pre-flighting only. Set GO=1 to start the search." >&2; fi

if ! command -v "$MATLAB" >/dev/null; then
    echo "matlab not on PATH. On a module system: module load matlab (e.g. matlab/R2023a), or MATLAB=/path/to/matlab" >&2
    exit 1
fi
echo "matlab   $(command -v "$MATLAB")   repo $root   commit $(git -C "$root" rev-parse --short HEAD 2>/dev/null || echo '?')"

pids=()
for p in $PROFILES; do
    log="$res/slip_$p.log"; [ "${SMOKE:-0}" = "1" ] && log="$res/slip_${p}_smoke.log"
    echo "arm $p -> $log"
    "$MATLAB" -nodisplay -nosplash -batch \
      "cd('$root'); addpath(genpath('.')); PROFILE='$p'; TIME_HRS=$HOURS; NWORKERS=$WORKERS; $MODE \
       run('Analyses/OverstretchedHeads/RunSlipCampaign.m')" > "$log" 2>&1 &
    pids+=($!)
    sleep "${STAGGER:-30}"   # stagger the pool starts: simultaneous parpool() calls on one node collide
done
rc=0; for pid in "${pids[@]}"; do wait "$pid" || rc=1; done
for p in $PROFILES; do
    log="$res/slip_$p.log"; [ "${SMOKE:-0}" = "1" ] && log="$res/slip_${p}_smoke.log"
    echo "--- $p:"; grep -E "PRE-FLIGHT OK|seed .* here|IMPROVED|CAMPAIGN DONE|Error|error" "$log" | tail -5 || true
done
if [ "${SMOKE:-0}" = "1" ]; then rm -rf "$root"/params/slip_*_smoke*; echo "smoke artefacts removed"; fi
exit $rc
