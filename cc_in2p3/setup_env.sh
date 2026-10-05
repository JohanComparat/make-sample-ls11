#!/bin/bash
# One-off set-up of the python environment on CC-IN2P3 (login node: needs internet access).
#   source cc_in2p3/env_ccin2p3.sh && bash cc_in2p3/setup_env.sh
# Creates $LS11_ENV from environment.yml, installs CIGALE from a git-lfs clone (its data are in
# git-lfs; numpy >= 2.4 needs np.trapz -> np.trapezoid), then fetches the LePhare data and
# registers the CIGALE filters.
set -euo pipefail
: "${LS11_REPO:?source cc_in2p3/env_ccin2p3.sh first}"
mkdir -p "$CONDA_PKGS_DIRS" "$PIP_CACHE_DIR" "$LS11_SPS/software" "$LS11_SPS/bin"
source "$HOME/miniforge3/bin/activate"
if [[ ! -d "$LS11_ENV" ]]; then
  (cd "$LS11_REPO" && conda env create -y -p "$LS11_ENV" -f environment.yml)
fi
conda activate "$LS11_ENV"
LFS="$LS11_SPS/bin/git-lfs"
if [[ ! -x "$LFS" ]]; then
  tmp=$(mktemp -d)
  wget -q -O "$tmp/lfs.tgz" https://github.com/git-lfs/git-lfs/releases/download/v3.6.1/git-lfs-linux-amd64-v3.6.1.tar.gz
  tar xzf "$tmp/lfs.tgz" -C "$tmp" && cp "$tmp"/git-lfs-3.6.1/git-lfs "$LFS" && rm -rf "$tmp"
fi
CIG="$LS11_SPS/software/cigale"
if [[ ! -d "$CIG" ]]; then
  PATH="$LS11_SPS/bin:$PATH" git -c filter.lfs.smudge="git-lfs smudge -- %f" \
      -c filter.lfs.process="git-lfs filter-process" -c filter.lfs.required=true \
      clone --depth 1 --branch v2025.1 https://gitlab.lam.fr/cigale/cigale.git "$CIG"
  grep -rl "np\.trapz(" --include=*.py "$CIG" | xargs sed -i "s/np\.trapz(/np.trapezoid(/g"
fi
mkdir -p "$LS11_SPS/.cache/tmp"
python -c "import pcigale" 2>/dev/null || TMPDIR="$LS11_SPS/.cache/tmp" pip install "$CIG"   # 3 GB wheel: not in /tmp
cd "$LS11_REPO"
python scripts/04_stellar_mass.py --fetch --code lephare,cigale,kcorrect
python -c "import kcorrect, lephare, pcigale, ls11samples; print('environment ready')"
