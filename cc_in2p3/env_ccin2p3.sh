# Environment of the LS DR11 sample pipeline on CC-IN2P3 (https://doc.cc.in2p3.fr).
# Fill in the <...> paths once, then:  source cc_in2p3/env_ccin2p3.sh
export LS11_REPO=${LS11_REPO:-$HOME/software/make-sample-ls11}
export LS11_DIR=${LS11_DIR:-/sps/<group>/<user>/legacysurvey/dr11}       # <region>/sweep/11.0, 11.0-photo-z, randoms
export LS11_REGION=${LS11_REGION:-south}
export LS11_SWEEP_VER=${LS11_SWEEP_VER:-11.0}
export LS11_RANDOMS=${LS11_RANDOMS:-randoms-${LS11_REGION}-1-0.fits}    # glob; several files for denser randoms
export LS11_OUT=${LS11_OUT:-/sps/<group>/<user>/ls11_samples/${LS11_REGION}}
export LS11_GAIA_MAPS=${LS11_GAIA_MAPS:-/sps/<group>/<user>/legacysurvey/dr10/systematics}
export LEPHAREDIR=${LEPHAREDIR:-/sps/<group>/<user>/lephare/data}
export LEPHAREWORK=${LEPHAREWORK:-/sps/<group>/<user>/lephare/work}
# command that activates the python environment built from environment.yml
export LS11_ACTIVATE=${LS11_ACTIVATE:-"source $HOME/miniforge3/bin/activate ls11"}
