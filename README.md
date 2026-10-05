# make-sample-ls11

Volume-limited galaxy samples and matching randoms from the Legacy Surveys DR11
(<https://www.legacysurvey.org/>), for clustering (`sum_stat`) and imaging systematics
(`sys_mapping`). It succeeds the DR10 BGS-like / `LS10_VLIM_*` samples.

## Pipeline

| step | script | product (under `$LS11_OUT`) |
|---|---|---|
| 1 | `scripts/01_select.py` | `bgsl/BGSl-<sweep>.fits`: every object with dereddened r ≤ 19.5, `SEL_FLAGS`, photo-z, `BEST_Z`; HDU `CUTFLOW` |
| 2 | `scripts/02_randoms.py` | `LS11_BGSl_DATA.fits` (`SEL_FLAGS == 0`), `LS11_BGSl_RAND.fits`, `footprint/` and `systematics/<nside>/` HEALPix maps |
| 3 | `scripts/03_kcorr_absmag.py` | K-corrections and absolute magnitudes |
| 4 | `scripts/04_stellar_mass.py` | stellar masses (code chosen with `benchmarks/`) |
| 5 | `scripts/05_vlim.py` | Mr and M* volume-limited samples + randoms |
| 6 | `scripts/06_export.py` | `sys_mapping` / `sum_stat` inputs + manifest |

Steps 1–2 are implemented; steps 3–6 are in progress. `scripts/run_all.sh` runs the steps in
order. `cc_in2p3/` holds the batch scripts for the full
footprint.

## Selection (`config/default.yaml`)

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
| `LS11_OUT` | `$LS11_DIR/$LS11_REGION/samples` |
| `LS11_CONFIG` | `config/default.yaml` |
| `LS11_GAIA_MAPS` | `~/data/legacysurvey/dr10/systematics` (full-sky Gaia star-density maps) |
| `LS11_NPROC` | `min(8, ncpu)` |

## Install

```bash
pip install -e .          # numpy, scipy, astropy, fitsio, healpy, pyyaml
pip install -e .[sed]     # SED-fitting codes for steps 3-4
pytest
```

Maintainer setup: on the laptop, use the shared `dev` env (`~/software/dev_env`). Do not create
a new environment.
