#!/bin/bash
# Submit the pipeline on CC-IN2P3 (Slurm, partition htc, licence sps) with dependencies:
#   select (array) -> randoms (+ maps)
#   prepare (LePhare libraries; data fetched beforehand on the login node with 04_stellar_mass.py --fetch)
#   select + prepare -> fit, one array per code (lephare, cigale, kcorrect)
#   [vlim -> export, with STEPS="... vlim export"]
# Usage: source cc_in2p3/env_ccin2p3.sh && bash cc_in2p3/submit_all.sh [--dry-run]
set -euo pipefail
: "${LS11_REPO:?source cc_in2p3/env_ccin2p3.sh first}"
NARRAY=${NARRAY:-100}                         # array tasks; each processes sweeps i, i+NARRAY, ...
CODES=${CODES:-"lephare cigale kcorrect"}     # stellar-mass codes (sed.codes of the config)
STEPS=${STEPS:-"select randoms prepare fit"}  # add "vlim export" once the masses are checked
LOGS="$LS11_DIR/$LS11_REGION/$(basename "$LS11_CONFIG" .yaml)_logs"; mkdir -p "$LOGS"
SB="sbatch --parsable -p htc -L sps"
[[ "${1:-}" == "--dry-run" ]] && SB="echo sbatch -p htc -L sps"
JOB="$LS11_REPO/cc_in2p3/job_step.sh"
ARR="--array=0-$((NARRAY - 1))"
has() { [[ " $STEPS " == *" $1 "* ]]; }
dep() { local d=""; for j in "$@"; do [[ -n "$j" ]] && d="$d:$j"; done; [[ -n "$d" ]] && echo "--dependency=afterok$d"; }

j1=""; j2=""; jp=""; jf=""
has select && j1=$($SB -J ls11_select $ARR -t 0-02:00 -c 4 --mem 16G -o "$LOGS/select_%A_%a.log" "$JOB" 01_select.py)
has randoms && j2=$($SB -J ls11_randoms $(dep $j1) -t 1-00:00 -c 4 --mem 32G -o "$LOGS/randoms_%j.log" "$JOB" 02_randoms.py)
has prepare && jp=$($SB -J ls11_prepare -t 0-06:00 -c 16 --mem 32G -o "$LOGS/prepare_%j.log" "$JOB" 04_stellar_mass.py --prepare --code "${CODES// /,}")
if has fit; then
  for code in $CODES; do
    case $code in
      lephare) res="-t 1-00:00 -c 16 --mem 32G" ;;    # 2 GB library read per task + OpenMP fit
      cigale)  res="-t 1-00:00 -c 16 --mem 32G" ;;
      *)       res="-t 0-06:00 -c 2 --mem 8G" ;;
    esac
    j=$($SB -J "ls11_fit_$code" $(dep $j1 $jp) $ARR $res -o "$LOGS/fit_${code}_%A_%a.log" "$JOB" 04_stellar_mass.py --code "$code")
    jf="$jf $j"
  done
fi
j5=""
has vlim && j5=$($SB -J ls11_vlim $(dep $j2 $jf) -t 0-12:00 -c 4 --mem 128G -o "$LOGS/vlim_%j.log" "$JOB" 05_vlim.py)
has export && $SB -J ls11_export $(dep $j5) -t 0-02:00 -c 1 --mem 16G -o "$LOGS/export_%j.log" "$JOB" 06_export.py >/dev/null
echo "submitted: select=$j1 randoms=$j2 prepare=$jp fit=$jf vlim=$j5 (logs: $LOGS)"
