# Environment of the LS DR11 sample pipeline on CC-IN2P3 (https://doc.cc.in2p3.fr).
#   source cc_in2p3/env_ccin2p3.sh
# The home directory is small: code, environment, caches and code data live on /sps.
export LS11_SPS=${LS11_SPS:-/sps/lsst/users/$USER}
export LS11_REPO=${LS11_REPO:-$LS11_SPS/software/make-sample-ls11}
export LS11_DIR=${LS11_DIR:-/sps/lsst/datasets/desi/legacysurveys/dr11}   # <region>/sweep/11.0, 11.0-photo-z, randoms
export LS11_REGION=${LS11_REGION:-south}
export LS11_SWEEP_VER=${LS11_SWEEP_VER:-11.0}
# products: <LS11_DIR>/<region>/sweep/11.0-<tag|code>/ (default LS11_SWEEP_OUT) and <LS11_DIR>/<region>/<tag>/
# 2 random files: 5000 deg^-2
export LS11_RANDOMS=${LS11_RANDOMS:-"randoms-${LS11_REGION}-1-[01].fits"}
export LS11_CONFIG=${LS11_CONFIG:-$LS11_REPO/config/bgs_r21_dr10bits.yaml}
export LS11_GAIA_MAPS=${LS11_GAIA_MAPS:-$LS11_SPS/gaia_maps}
export LEPHAREDIR=${LEPHAREDIR:-$LS11_SPS/lephare/data}
export LEPHAREWORK=${LEPHAREWORK:-$LS11_SPS/lephare/work}
export CONDA_PKGS_DIRS=${CONDA_PKGS_DIRS:-$LS11_SPS/.cache/conda_pkgs}
export PIP_CACHE_DIR=${PIP_CACHE_DIR:-$LS11_SPS/.cache/pip}
# python environment built from environment.yml (prefix on /sps)
export LS11_ENV=${LS11_ENV:-$LS11_SPS/envs/ls11}
export LS11_ACTIVATE=${LS11_ACTIVATE:-"source $HOME/miniforge3/bin/activate $LS11_ENV"}
