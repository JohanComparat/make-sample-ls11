#!/bin/bash
# Run the whole pipeline locally (laptop defaults of the LS11_* variables, see README).
# Steps 1-2 need numpy/astropy/fitsio/healpy; steps 4-5 also need the SED codes of sed.codes
# (kcorrect, lephare, pcigale).  PY: python to use (default: the dev-full env).
set -euo pipefail
cd "$(dirname "$0")/.."
PY=${PY:-$HOME/software/miniforge3/envs/dev-full/bin/python}
export PYTHONPATH="$PWD${PYTHONPATH:+:$PYTHONPATH}"
"$PY" scripts/01_select.py "$@"
"$PY" scripts/02_randoms.py "$@"
"$PY" scripts/04_stellar_mass.py "$@"          # every code of sed.codes
"$PY" scripts/05_vlim.py "$@"
"$PY" scripts/06_export.py
