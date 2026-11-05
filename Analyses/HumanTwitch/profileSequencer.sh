#!/bin/bash
# Case-B profile-likelihood sequencer: tB optimum, 4 Ca-error refits, and profile points over r_net, r_SRX, r_kSE.
# Keeps at most MAXJ unfinished jobs; a job is finished when results/<name>_fit.mat exists. Log: results/profile_sequencer.log
cd "$(dirname "$0")/results" || exit 1
ML="/mnt/c/Program Files/MATLAB/R2023a/bin/matlab.exe"; MAXJ=8
jobs=(tB tB_caE1 tB_caE2 tB_caE3 tB_caE4 tB__r_net__1 tB__r_net__2 tB__r_net__3 tB__r_net__4 tB__r_net__5 tB__r_net__6 tB__r_net__7 tB__r_net__8 tB__r_SRX__1 tB__r_SRX__2 tB__r_SRX__3 tB__r_SRX__4 tB__r_SRX__5 tB__r_kSE__1 tB__r_kSE__2 tB__r_kSE__3 tB__r_kSE__4 tB__r_kSE__5)
running=()
for j in "${jobs[@]}"; do
  while true; do n=0; for r in "${running[@]}"; do [ -f "${r}_fit.mat" ] || n=$((n+1)); done; [ $n -lt $MAXJ ] && break; sleep 30; done
  echo "$(date +%T) launch $j" >> profile_sequencer.log
  nohup "$ML" -batch "cd('C:\\home\\git\\ATP-wt-human'); addpath(genpath('C:\\home\\git\\ATP-wt-human')); runTier('$j')" -logfile "C:\\home\\git\\ATP-wt-human\\Analyses\\HumanTwitch\\results\\$j.console.log" > /dev/null 2>&1 &
  running+=("$j"); sleep 2
done
for r in "${running[@]}"; do until [ -f "${r}_fit.mat" ]; do sleep 30; done; done
echo "$(date +%T) all profile jobs done" >> profile_sequencer.log
