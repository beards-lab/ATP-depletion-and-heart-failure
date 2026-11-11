#!/usr/bin/env bash
# build_dashboard.sh -- assemble dashboard/index.html from
#   results/checkin_<TAG>.json          (RunCheckin.m: seed -> best, per feature)
#   results/checkin_history_<TAG>.jsonl (one row per checkpoint)
#   dashboard/decisions.jsonl           (what was changed and why, by hand)
#   the campaign log                    (is it alive, which round, health)
# and copy the latest figures next to it, ready to publish as an Artifact.
#
#   bash Analyses/KtrSegmentSlope/dashboard/build_dashboard.sh /tmp/c6_leg5.log [TAG]
set -u
LOG="${1:-/tmp/c6_leg5.log}"
TAG="${2:-refit6}"
ROOT="/c/home/git/ATP-depletion-and-heart-failure"
DIR="$ROOT/Analyses/KtrSegmentSlope"
OUT="$DIR/dashboard"
RES="$DIR/results"

# ---- campaign health, from the log -----------------------------------------
nproc_ml=$(ps -W 2>/dev/null | grep -c "MATLAB.exe")
round=$(grep -oE "Round [0-9]+/1000" "$LOG" 2>/dev/null | tail -1 | grep -oE "[0-9]+" | head -1)
sent=$(grep -cE "1e\+09|1e\+06" "$LOG" 2>/dev/null)
leaks=$(grep -c "Leaking probs" "$LOG" 2>/dev/null)
wdog=$(grep -c "took longer than MaxRunTime" "$LOG" 2>/dev/null)
if grep -qE "CAMPAIGN DONE|=== Done\." "$LOG" 2>/dev/null; then
  state="done"; label="leg finished"
elif grep -qE "Error using|\{Error|Access violation" "$LOG" 2>/dev/null; then
  state="stopped"; label="crashed"
elif [ "$nproc_ml" -gt 0 ]; then
  state="running"; label="running"
else
  state="stopped"; label="not running"
fi

if [ -f "$OUT/state_override.txt" ]; then state="done"; label=$(cat "$OUT/state_override.txt"); fi

# ---- figures (only the ones that exist) ------------------------------------
figs=""
add_fig () {  # src caption alt
  if [ -f "$RES/$1" ]; then
    cp -f "$RES/$1" "$OUT/$1"
    ts=$(date -r "$RES/$1" +%H:%M)
    [ -n "$figs" ] && figs="$figs,"
    figs="$figs{\"src\":\"$1\",\"caption\":\"$2 (rendered $ts)\",\"alt\":\"$3\"}"
  fi
}
PAIR=$(cat "$OUT/figpair.txt" 2>/dev/null || echo "seed vs best")
add_fig "${TAG}_protocols.png" "${PAIR}. Every protocol against its target. Panel G is the ktr ladder across the slack segments; D-F are the held-out ktr protocol on each day." "refit6 protocols figure"
add_fig "${TAG}_bars.png"      "${PAIR}. Per-feature weighted cost, seed against best, and what moved." "refit6 per-feature cost bars"
add_fig "${TAG}_features.png"  "${PAIR}. Every scored feature: pooled data against seed and best." "refit6 scored features"

# ---- JSON pieces ------------------------------------------------------------
checkin=$(cat "$RES/checkin_${TAG}.json" 2>/dev/null || echo "{}")
history="[$(paste -sd, "$RES/checkin_history_${TAG}.jsonl" 2>/dev/null)]"
decisions="[$(paste -sd, "$OUT/decisions.jsonl" 2>/dev/null)]"
render=$(cat "$RES/render_${TAG}.json" 2>/dev/null || echo "{}")
render=$(cat "$RES/render_${TAG}.json" 2>/dev/null || echo "{}")
status="{\"state\":\"$state\",\"label\":\"$label\",\"round\":${round:-null},\"sentinels\":${sent:-0},\"leaks\":${leaks:-0},\"watchdog\":${wdog:-0},\"procs\":$nproc_ml}"

DATA="$OUT/data.js"
printf 'window.DATA = {"checkin":%s,"history":%s,"log":%s,"status":%s,"figures":[%s],"render":%s};\n' \
       "$checkin" "$history" "$decisions" "$status" "$figs" "$render" > "$DATA"

awk -v datafile="$DATA" '
  /\/\*__DATA__\*\// { while ((getline line < datafile) > 0) print line; next }
  { print }
' "$OUT/template.html" > "$OUT/index.html"

echo "built $OUT/index.html  state=$state round=${round:-?} sentinels=$sent leaks=$leaks watchdog=$wdog"
