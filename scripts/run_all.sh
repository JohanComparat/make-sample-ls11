#!/bin/bash
# Run the whole pipeline locally (laptop defaults of the LS11_* variables, see README).
# Steps 1-2 need numpy/astropy/fitsio/healpy; steps 3-5 also need kcorrect (and the SED code
# chosen in the config for step 4).  PY: python to use (default: the dev-full env).
set -euo pipefail
cd "$(dirname "$0")/.."
PY=${PY:-$HOME/software/miniforge3/envs/dev-full/bin/python}
export PYTHONPATH="$PWD${PYTHONPATH:+:$PYTHONPATH}"
"$PY" scripts/01_select.py "$@"
"$PY" scripts/02_randoms.py "$@"
"$PY" scripts/03_kcorr_absmag.py "$@"
"$PY" scripts/04_stellar_mass.py "$@"
"$PY" scripts/05_vlim.py "$@"
"$PY" scripts/06_export.py
