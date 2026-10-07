# make-sample-ls11

[![Documentation Status](https://readthedocs.org/projects/make-sample-ls11/badge/?version=latest)](https://make-sample-ls11.readthedocs.io/en/latest/)
[![tests](https://github.com/JohanComparat/make-sample-ls11/actions/workflows/tests.yml/badge.svg)](https://github.com/JohanComparat/make-sample-ls11/actions/workflows/tests.yml)
[![coverage](https://img.shields.io/badge/coverage-92%25-brightgreen)](https://github.com/JohanComparat/make-sample-ls11/actions/workflows/tests.yml)
[![tested with pytest](https://img.shields.io/badge/tested%20with-pytest-0A9EDC?logo=pytest&logoColor=white)](https://make-sample-ls11.readthedocs.io/en/latest/testing.html)
[![python](https://img.shields.io/badge/python-3.11%20%7C%203.12-3776AB?logo=python&logoColor=white)](https://make-sample-ls11.readthedocs.io/en/latest/installation.html)

Volume-limited galaxy samples and matching randoms from the Legacy Surveys DR11
(<https://www.legacysurvey.org/>), for clustering (`sum_stat`) and imaging systematics
(`sys_mapping`). It succeeds the DR10 BGS-like / `LS10_VLIM_*` samples.

Documentation: <https://make-sample-ls11.readthedocs.io/en/latest/>

## Pipeline

Per-sweep products sit next to the sweeps, one small file per sweep and product:
`<sweep root>/<ver>-<name>/<sweep>-<name>.fits`, alongside `11.0/` and `11.0-photo-z/`.
Positions and photometry are never copied; they are read back from the sweep at `SWEEP_ROW`.

| step | script | product |
|---|---|---|
| 1 | `scripts/01_select.py` | `11.0-<tag>/<sweep>-<tag>.fits`: `LS_ID_DR11`, `SWEEP_ROW`, `BEST_Z`, `BEST_Z_ERR`, `Z_SOURCE`, `STAR_FLAG` (22 B/row); HDU `CUTFLOW` |
| 2 | `scripts/02_randoms.py` | `<LS11_OUT>/<tag>/LS11_<tag>_RAND.fits` (RA, DEC, EBV); `footprint/` and `systematics/<nside>/` HEALPix maps |
| 4 | `scripts/04_stellar_mass.py` | `11.0-<code>/<sweep>-<code>.fits` for code = lephare, cigale, kcorrect: `SWEEP_ROW`, `LOGMSTAR`, `LOGMSTAR_ERR`, `MABS_R`, `MABS_R_ERR` (20 B/row), row-aligned with the selection |
| 5 | `scripts/05_vlim.py` | `<LS11_OUT>/<tag>/vlim/`: Mr and M* volume-limited samples + randoms |
| 6 | `scripts/06_export.py` | `sys_mapping` / `sum_stat` format checks + manifest |
| check | `scripts/validate_run.py` | completeness and row alignment of every per-sweep product, missing values, mass / Mr distributions and errors per code, agreement with `sed.primary`, randoms, disk use → `<LS11_OUT>/<tag>/validation.yaml` (exit 1 if a file is missing or misaligned) |

`scripts/run_all.sh` runs the steps in order. `cc_in2p3/submit_all.sh` submits them on CC-IN2P3
as Slurm job arrays with dependencies. Steps 1 and 4 run per sweep, one array task per group of
sweeps.

## Configurations

| config | tag | selection |
|---|---|---|
| `config/default.yaml` | `bgsl` | DESI BGS Bright-like, 13 < r ≤ 19.5 (below) |
| `config/bgs_r21_dr10bits.yaml` | `bgsr21` | every object with 13 < r ≤ 21 in the footprint, with the DR10 FITBITS, and no other cut |

**About `bgsr21`.**
- **No star/galaxy separation.** `STAR_FLAG` records it instead: bit 0 TYPE=PSF, bit 1 Gaia
  G − r_raw ≤ 0.6.
- **Inheritance.** It inherits `default.yaml` (`base:`) and switches the other cuts off with
  `null`.
- **Footprint.** Its MASKBITS are the DR10 ones (0, 1, 13) plus MEDIUM (11) and GALAXY (12), for
  galaxies and randoms. The DR10 FITBITS cut leaves those areas almost empty of galaxies
  (3.3% of the footprint), so they are removed from the randoms as well.

## Selection (`config/default.yaml`, tag `bgsl`)

All cuts are set in the YAML file. Bits are named there, and their numbers are checked against
the `MBIT_n` / `FBIT_n` keywords of each sweep header.

- **Footprint.** This function is shared by galaxies and randoms:
  - NOBS ≥ 1 and GALDEPTH > 0 in g, r and z;
  - MASKBITS NPRIMARY, BRIGHT, MEDIUM, GALAXY and CLUSTER rejected;
  - E(B−V) < 0.1;
  - the DESI north/south split, keeping the south only.
- **Galaxies only.** These are the DESI BGS Bright-like cuts:
  - TYPE ≠ DUP;
  - FITBITS FIT_BACKGROUND, MEDIUM, WALKER, RUNNER, GAIA_POINTSOURCE and ITERATIVE rejected (as
    for DR10);
  - Gaia star rejection (G − r_raw > 0.6, or no Gaia match);
  - colour sanity cuts;
  - FRACMASKED, FRACIN and FRACFLUX quality cuts;
  - the fibre-magnitude cut and the r_fibtot cut;
  - 13 < r ≤ 19.5, dereddened.
- **Redshift.** `BEST_Z` is the spectroscopic redshift when there is one (from the DR11 photo-z
  files), otherwise `Z_PHOT_MEAN_I`, otherwise `Z_PHOT_MEAN`. Every later step uses this
  redshift.
- **ID check.** The photo-z rows are checked object by object against the sweep IDs
  (`LS_ID_DR11`).

## Stellar masses and absolute magnitudes

`ls11samples/sed/` runs five codes behind one interface. All of them work at fixed z = `BEST_Z`,
with identical filter curves (kcorrect's DECam/WISE responses), dereddened grizW1W2 photometry and
error floor:

| code | models |
|---|---|
| kcorrect | NMF templates |
| LePhare | `config/lephare/LS11_zFIX.para`: BC03 Chabrier, Calzetti, the DR10 set-up |
| CIGALE | delayed-τ + BC03 + nebular + modified starburst |
| eazy | `corr_sfhz_13` |
| DSPS | FSPS SSPs and a jax delayed-τ × Z × A_V grid |

`benchmarks/sed_benchmark.py` compares them; the results are in
`benchmarks/results/local_strip/report.md`:
- every code agrees with LePhare within ±0.05 dex in the median (NMAD 0.09–0.10 dex);
- the photo-z error alone moves masses by 0.07–0.11 dex;
- DR11 LePhare differs from DR10 LePhare by +0.03 dex.

**Production.** LePhare, CIGALE and kcorrect run (`sed.codes`); LePhare is `sed.primary`.
- `MABS_R` is the rest-frame DECam r absolute magnitude (AB, H0 = 67.74).
- `LOGMSTAR` is log10 M*/M☉ (Chabrier).
- The errors are 1σ at fixed z; the photo-z term is not included (it can be derived from
  `BEST_Z_ERR`).

How each code gets its values and errors:

| code | LOGMSTAR | LOGMSTAR_ERR | MABS_R ± MABS_R_ERR |
|---|---|---|---|
| LePhare | `MASS_MED` | half the 68% interval | `MAG_ABS`; error of the observed band closest to rest-frame r (LePhare 1.0 `EMAG_ABS` holds m − M) |
| CIGALE | log of the Bayesian mean | Bayesian error (pcigale floor: 5%, i.e. ≥ 0.022 dex) | rest-frame L_ν(r) (floor ≥ 0.054 mag) |
| kcorrect | best fit | standard deviation over 20 Monte Carlo flux realisations | best fit and the same realisations |

**Volume-limited samples (step 5).** These are built in Mr (thresholds −18 … −22.5) and in M*
(9.0 … 11.5).
- Completeness limits come from the data: K_95(z) of a complete low-z SED set for Mr, Pozzetti
  et al. (2010) for M*.
- Each sample has three files, named `LS11_VLIM_ANY_<lo>_<Mr|Mstar>_<hi>_<zmin>_z_<zmax>_N_<N>`:
  - `_DATA.fits`, the galaxies;
  - `_RAND.fits`, with randoms at `vlim.n_rand_factor` per galaxy and shuffled redshifts;
  - `_COLOUR.fits`.
- *sys_mapping*: `--catalog-dir $LS11_OUT/<tag>/vlim --template-dir $LS11_OUT/<tag>/systematics/<nside>`.
- *sum_stat*: `--survey custom` (`BEST_Z`, `LPH_MASS_BEST`, `RAND.Z`).

## CC-IN2P3

- **Code and data locations:**
  - code: `/sps/lsst/users/$USER/software/make-sample-ls11`, a clone of the bare repository
    `/sps/lsst/users/$USER/git/make-sample-ls11.git`;
  - DR11: `/sps/lsst/datasets/desi/legacysurveys/dr11/south`.
- **Set-up, once, on the login node:**
  `source cc_in2p3/env_ccin2p3.sh && bash cc_in2p3/setup_env.sh`. This builds the env, installs
  CIGALE, fetches the LePhare data and registers the filters.
- **Run:** `source cc_in2p3/env_ccin2p3.sh && bash cc_in2p3/submit_all.sh`. This submits select →
  randoms, and prepare → one fit array per code. `STEPS` selects the steps.
- **Restarts:** existing outputs are skipped, so resubmitting picks up the sweeps that have been
  downloaded since.

## Environment

All paths come from environment variables. On the laptop the defaults are enough; on CC-IN2P3,
set them in `cc_in2p3/env_ccin2p3.sh`.

| variable | default |
|---|---|
| `LS11_DIR` | `~/data/legacysurvey/dr11` |
| `LS11_REGION` | `south` |
| `LS11_SWEEP_VER` | `11.0` |
| `LS11_SWEEPS` | `sweep-*.fits` (glob restricting the sweeps) |
| `LS11_RANDOMS` | `randoms-$LS11_REGION-1-0.fits` (glob, relative to `$LS11_DIR/$LS11_REGION/randoms`) |
| `LS11_OUT` | `$LS11_DIR/$LS11_REGION` (per-run products in `<LS11_OUT>/<tag>/`) |
| `LS11_SWEEP_OUT` | `$LS11_DIR/$LS11_REGION/sweep` (per-sweep products in `<ver>-<name>/`) |
| `LS11_CONFIG` | `config/default.yaml` |
| `LS11_GAIA_MAPS` | `~/data/legacysurvey/dr10/systematics` (full-sky Gaia star-density maps) |
| `LS11_NPROC` | `min(8, ncpu)` |
| `LEPHAREDIR`, `LEPHAREWORK` | LePhare data and libraries (default `~/.cache/lephare/{data,work}`) |
| `EAZY_DATA` | eazy-photoz templates and filters (default `~/.cache/eazy-photoz`, cloned when missing) |
| `DSPS_DRN` | DSPS SSP file `ssp_data_fsps_v3.2_lgmet_age.h5` (default `~/.cache/dsps`) |

## Install

The SED-fitting codes have many dependencies and data files, and some are not on PyPI. They are
not installed with the package: follow the recipes of the
[installation page](https://make-sample-ls11.readthedocs.io/en/latest/installation.html) for the
codes you need.

**Core** (steps 1, 2, 5, 6, the tests):

```bash
conda create -n ls11 -c conda-forge python=3.12 numpy scipy astropy fitsio healpy pyyaml \
    matplotlib pytest pytest-cov configobj
conda activate ls11
git clone https://github.com/JohanComparat/make-sample-ls11.git && cd make-sample-ls11
python -m pip install --no-deps -e .
pytest                                  # about 20 s, no data needed
```

**SED codes** (step 4). Each code's recipe is on the installation page:

| code | role | install | data (env variable) |
|---|---|---|---|
| kcorrect 5.1.9 | needed by every SED step (filter curves, Mr) | `pip install kcorrect==5.1.9` | in the wheel |
| LePhare 1.0.0 | production masses | `pip install lephare==1.0.0`, then `04_stellar_mass.py --fetch` and `--prepare` (libraries, about 25 min) | 82 MB + 2.2 GB (`LEPHAREDIR`, `LEPHAREWORK`) |
| CIGALE 2025.1 | production masses | git-lfs clone, `np.trapz` patch, `pip install ./cigale` (builds its database) | about 10 GB |
| eazy 0.8.7 | benchmark | `pip install eazy==0.8.7` and clone eazy-photoz | 200 MB (`EAZY_DATA`) |
| DSPS 0.4.8 | benchmark | `pip install dsps==0.4.8` and download the SSP file | 60 MB (`DSPS_DRN`) |

`environment.yml` holds the core plus kcorrect and LePhare (it is used on CC-IN2P3 by
`cc_in2p3/setup_env.sh`, which also installs CIGALE).

Maintainer setup: on the laptop, use the shared `dev` / `dev-full` envs (`~/software/dev_env`). Do
not create a new environment.

## Documentation

Online at <https://make-sample-ls11.readthedocs.io>, rebuilt by Read the Docs at every push to
`main` (`.readthedocs.yaml`, which installs only `docs/requirements.txt`). Sphinx sources are in
`docs/` (installation, configuration, every step, methods, file formats, API). To build locally:
`pip install -r docs/requirements.txt` then `make -C docs html`, and open `docs/_build/html/index.html`.
