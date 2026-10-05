#!/bin/bash
# One pipeline step as a Slurm job (or job-array task): job_step.sh <script> [options]
# Job arrays slice the sweeps with SLURM_ARRAY_TASK_ID / SLURM_ARRAY_TASK_COUNT (ls11samples.cli).
set -euo pipefail
source "${LS11_REPO:?source cc_in2p3/env_ccin2p3.sh first}/cc_in2p3/env_ccin2p3.sh"
eval "$LS11_ACTIVATE"
export LS11_NPROC=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
cd "$LS11_REPO"
echo "$(date -Is) $(hostname) step $1 task ${SLURM_ARRAY_TASK_ID:-0}/${SLURM_ARRAY_TASK_COUNT:-1}"
python "scripts/$1" "${@:2}"
echo "$(date -Is) done"
