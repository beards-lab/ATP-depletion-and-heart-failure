#!/bin/bash
# T-World Ca-source ladder: tw1_isoT and tw1_isoT_xb now (de Lange transient replaced by T-World's beat),
# tw2_caH (T-World Ca-handling levers freed) once tw1_isoT has finished. Logs to results/sequencer.log.
cd "$(dirname "$0")/results" || exit 1
ML="/mnt/c/Program Files/MATLAB/R2023a/bin/matlab.exe"
launch() { echo "$(date +%T) launch $1" >> sequencer.log
  nohup "$ML" -batch "cd('C:\\home\\git\\ATP-wt-human'); addpath(genpath('C:\\home\\git\\ATP-wt-human')); runTier('$1')" -logfile "C:\\home\\git\\ATP-wt-human\\Analyses\\HumanTwitch\\results\\$1.console.log" > /dev/null 2>&1 & }
waitfor() { for t in "$@"; do until [ -f "${t}_fit.mat" ]; do sleep 45; done; echo "$(date +%T) done $t" >> sequencer.log; done; }
launch tw1_isoT; launch tw1_isoT_xb
waitfor tw1_isoT
launch tw2_caH
waitfor tw1_isoT_xb tw2_caH
echo "$(date +%T) T-World ladder complete" >> sequencer.log
