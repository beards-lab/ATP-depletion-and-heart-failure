#!/usr/bin/env bash
# watchdog.sh -- wait for MATLAB to become startable again, then launch the
# refit6 12 h run by itself.
#
# WHY. From 2026-09-22 ~21:10 every `matlab -batch` on this machine blocks at
# startup with ~0 CPU and no output. The cause is the online entitlement call:
# each stuck session holds a TCP connection to login.mathworks.com
# (54.75.159.26 / 52.19.196.186 :443) in CloseWait -- the server closed it and
# MATLAB never times out. TCP reachability to those IPs is fine, so it is the
# HTTPS/login session, not the network. Passing -c with the cached .mwlicx
# fails FAST (License Manager Error -2, wrong format for that flag), which
# confirms the core starts and only the online path hangs.
#
# If that is server-side or a token refresh, it may clear on its own; if not it
# needs a reboot or one interactive MATLAB sign-in. Either way this script
# starts the campaign the moment a health check passes, so no time is lost
# waiting for a human.
#
# SAFETY. Every failed health check leaves a ~556 MB session that Windows
# reports as exited but will not release, so the interval is long and the loop
# aborts if the process count gets dangerous.
#
#   nohup bash Analyses/KtrSegmentSlope/watchdog.sh > /tmp/watchdog.log 2>&1 &

ROOT="/c/home/git/ATP-depletion-and-heart-failure"
MLB="/c/Program Files/MATLAB/R2023a/bin/matlab.exe"
SP="C:/Users/fjezek/AppData/Local/Temp/claude/C--home-git-ATP-depletion-and-heart-failure/221c6c43-630e-43cd-845d-701d0090a187/scratchpad"
PROBE="$SP/wd_probe.m"
FLAG="$SP/wd_ok.txt"
INTERVAL=1800          # 30 min between attempts
MAX_MATLAB=12          # abort rather than exhaust memory with dead sessions
MAX_TRIES=26           # ~13 h

cat > "$PROBE" <<'EOF'
f = fopen('C:\Users\fjezek\AppData\Local\Temp\claude\C--home-git-ATP-depletion-and-heart-failure\221c6c43-630e-43cd-845d-701d0090a187\scratchpad\wd_ok.txt','w');
fprintf(f, 'alive %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS'));
fclose(f);
EOF

log() { echo "[$(date +%H:%M:%S)] $*"; }

log "watchdog up; probing every $((INTERVAL/60)) min, max $MAX_TRIES tries"

for ((i=1; i<=MAX_TRIES; i++)); do
    n=$(ps -W 2>/dev/null | grep -c "MATLAB.exe")
    if [ "$n" -ge "$MAX_MATLAB" ]; then
        log "ABORT: $n MATLAB processes already stuck (limit $MAX_MATLAB). Needs a reboot."
        exit 2
    fi

    rm -f "$FLAG"
    log "try $i/$MAX_TRIES (stuck sessions: $n)"
    timeout 150 "$MLB" -batch "run('$PROBE')" >/dev/null 2>&1

    if [ -s "$FLAG" ]; then
        log "MATLAB IS BACK: $(cat "$FLAG")"
        cd "$ROOT" || exit 1

        # if the earlier legs left a state file, resume; otherwise start clean
        if [ -f "$ROOT/params/refit6_state.mat" ]; then
            log "resuming from existing refit6_state.mat"
        else
            log "fresh start from refit6_seed"
        fi

        nohup "$MLB" -batch "TIME_HRS=2; run('Analyses/KtrSegmentSlope/RunRefit6Optim.m')" \
              > /tmp/c6_leg1.log 2>&1 &
        log "leg 1 launched (2 h). Watchdog exiting; the campaign owns it now."
        exit 0
    fi

    log "still blocked; sleeping $((INTERVAL/60)) min"
    sleep "$INTERVAL"
done

log "gave up after $MAX_TRIES tries"
exit 3
