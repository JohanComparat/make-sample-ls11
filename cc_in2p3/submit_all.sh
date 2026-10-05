#!/bin/bash
# Submit the whole pipeline on CC-IN2P3 (Slurm, partition htc) with dependencies:
#   1 select (array) -> 2 randoms+maps
#   1 -> 3 kcorr (array)
#   prepare (LePhare libraries, CIGALE filters; once) -> 4 masses, one array per code
#   2, 3, 4 -> 5 vlim -> 6 export
# Usage: source cc_in2p3/env_ccin2p3.sh && bash cc_in2p3/submit_all.sh [--dry-run]
set -euo pipefail
: "${LS11_REPO:?source cc_in2p3/env_ccin2p3.sh first}"
NARRAY=${NARRAY:-100}                       # array tasks; each processes sweeps i, i+NARRAY, ...
CODES=${CODES:-"lephare cigale kcorrect"}   # stellar-mass codes (sed.codes of the config)
LOGS="$LS11_OUT/logs"; mkdir -p "$LOGS"
SB="sbatch --parsable -p htc -L sps"
[[ "${1:-}" == "--dry-run" ]] && SB="echo sbatch -p htc -L sps"
JOB="$LS11_REPO/cc_in2p3/job_step.sh"
ARR="--array=0-$((NARRAY - 1))"

j1=$($SB -J ls11_select $ARR -t 0-02:00 -c 4 --mem 16G -o "$LOGS/select_%A_%a.log" "$JOB" 01_select.py)
j2=$($SB -J ls11_randoms --dependency=afterok:$j1 -t 1-00:00 -c 8 --mem 96G -o "$LOGS/randoms_%j.log" "$JOB" 02_randoms.py)
j3=$($SB -J ls11_kcorr --dependency=afterok:$j1 $ARR -t 0-04:00 -c 4 --mem 16G -o "$LOGS/kcorr_%A_%a.log" "$JOB" 03_kcorr_absmag.py)
jp=$($SB -J ls11_prepare -t 0-04:00 -c 16 --mem 32G -o "$LOGS/prepare_%j.log" "$JOB" 04_stellar_mass.py --prepare --code "${CODES// /,}")
deps="$j2:$j3"
for code in $CODES; do
  case $code in
    lephare) res="-t 1-00:00 -c 16 --mem 48G" ;;    # 2 GB library read per task + OpenMP fit
    cigale)  res="-t 1-00:00 -c 16 --mem 32G" ;;
    *)       res="-t 0-04:00 -c 2 --mem 8G" ;;
  esac
  j4=$($SB -J "ls11_mstar_$code" --dependency=afterok:$j1:$jp $ARR $res -o "$LOGS/mstar_${code}_%A_%a.log" "$JOB" 04_stellar_mass.py --code "$code")
  deps="$deps:$j4"
done
j5=$($SB -J ls11_vlim --dependency=afterok:$deps -t 0-12:00 -c 4 --mem 128G -o "$LOGS/vlim_%j.log" "$JOB" 05_vlim.py)
j6=$($SB -J ls11_export --dependency=afterok:$j5 -t 0-02:00 -c 1 --mem 16G -o "$LOGS/export_%j.log" "$JOB" 06_export.py)
echo "submitted: select=$j1 randoms=$j2 kcorr=$j3 prepare=$jp vlim=$j5 export=$j6 (masses: $deps)"
