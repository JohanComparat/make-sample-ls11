#!/bin/bash
# Submit the whole pipeline on CC-IN2P3 (Slurm, partition htc) with dependencies:
#   1 select (array) -> 2 randoms+maps ; 1 -> 3 kcorr (array) ; 1 -> 4 masses (array) ; 2,3,4 -> 5 vlim -> 6 export
# Usage: source cc_in2p3/env_ccin2p3.sh && bash cc_in2p3/submit_all.sh [--dry-run]
set -euo pipefail
: "${LS11_REPO:?source cc_in2p3/env_ccin2p3.sh first}"
NARRAY=${NARRAY:-100}          # array tasks; each processes sweeps i, i+NARRAY, ...
CODE=${CODE:-}                 # stellar-mass code (default: sed.code of the config)
LOGS="$LS11_OUT/logs"; mkdir -p "$LOGS"
SB="sbatch --parsable -p htc -L sps"
[[ "${1:-}" == "--dry-run" ]] && SB="echo sbatch -p htc -L sps"
JOB="$LS11_REPO/cc_in2p3/job_step.sh"
ARR="--array=0-$((NARRAY - 1))"

j1=$($SB -J ls11_select $ARR -t 0-02:00 -c 4 --mem 16G -o "$LOGS/select_%A_%a.log" "$JOB" 01_select.py)
j2=$($SB -J ls11_randoms --dependency=afterok:$j1 -t 1-00:00 -c 8 --mem 96G -o "$LOGS/randoms_%j.log" "$JOB" 02_randoms.py)
j3=$($SB -J ls11_kcorr --dependency=afterok:$j1 $ARR -t 0-04:00 -c 4 --mem 16G -o "$LOGS/kcorr_%A_%a.log" "$JOB" 03_kcorr_absmag.py)
j4=$($SB -J ls11_mstar --dependency=afterok:$j1 $ARR -t 1-00:00 -c 8 --mem 32G -o "$LOGS/mstar_%A_%a.log" "$JOB" 04_stellar_mass.py ${CODE:+--code $CODE})
j5=$($SB -J ls11_vlim --dependency=afterok:$j2:$j3:$j4 -t 0-12:00 -c 4 --mem 128G -o "$LOGS/vlim_%j.log" "$JOB" 05_vlim.py ${CODE:+--code $CODE})
j6=$($SB -J ls11_export --dependency=afterok:$j5 -t 0-02:00 -c 1 --mem 16G -o "$LOGS/export_%j.log" "$JOB" 06_export.py)
echo "submitted: select=$j1 randoms=$j2 kcorr=$j3 mstar=$j4 vlim=$j5 export=$j6"
