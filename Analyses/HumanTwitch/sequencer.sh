#!/bin/bash
# Stage-2 ladder sequencer: launches tiers as their warm-start tier finishes. Logs to results/sequencer.log.
cd "$(dirname "$0")/results" || exit 1
ML="/mnt/c/Program Files/MATLAB/R2023a/bin/matlab.exe"
launch() { echo "$(date +%T) launch $1" >> sequencer.log
  nohup "$ML" -batch "cd('C:\\home\\git\\ATP-wt-human'); addpath(genpath('C:\\home\\git\\ATP-wt-human')); runTier('$1')" -logfile "C:\\home\\git\\ATP-wt-human\\Analyses\\HumanTwitch\\results\\$1.console.log" > /dev/null 2>&1 & }
waitfor() { for t in "$@"; do until [ -f "${t}_fit.mat" ]; do sleep 45; done; echo "$(date +%T) done $t" >> sequencer.log; done; }
waitfor t1_isoT
launch t2a_Pi; launch t2b_SRX; launch t2c_SE
waitfor t2a_Pi t2b_SRX t2c_SE
launch t2d_pas; launch t2e_SL; launch t2_all
waitfor t2_all
echo "$(date +%T) stage 2 complete; s3_free is launched by hand (pool via MCP if the box is free)" >> sequencer.log
