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

`scripts/run_all.sh` runs the steps in order. `cc_in2p3/submit_all.sh` submits them on CC-IN2P3
as Slurm job arrays with dependencies. Steps 1, 3 and 4 run per sweep (one sweep per process or
array task), and steps 3–5 join products to the galaxies by `LS_ID_DR11`.

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

## K-corrections, absolute magnitudes, stellar masses

**K-corrections (step 3).** kcorrect v5 (Blanton & Roweis 2007) fits the dereddened grizW1W2
fluxes at z = `BEST_Z`. It produces:
- `KCORR_<b>` and `ABSMAG_<b>` for b = G, R, Z, in rest-frame DECam bands;
- `ABSMAG_R01` in SDSS ^{0.1}r;
- `KC_LOGMSTAR`.

The Mr completeness limit uses the 95th percentile of the r-band K-correction of a complete
low-z SED set, evaluated at every z. This avoids biasing the limit with the K-corrections of
flux-limited high-z galaxies.

**Stellar masses (step 4).** `ls11samples/sed/` runs five codes behind one interface, all at
fixed z and with identical filter curves, photometry and error floor:

| code | models |
|---|---|
| kcorrect | NMF templates |
| LePhare | `config/lephare/LS11_zFIX.para`: BC03 Chabrier, Calzetti, the DR10 set-up |
| CIGALE | delayed-τ + BC03 + nebular + modified starburst |
| eazy | `corr_sfhz_13` templates |
| DSPS | FSPS SSPs and a jax-tabulated delayed-τ × Z × A_V grid |

`benchmarks/sed_benchmark.py` compares the five codes on cost, agreement, photo-z sensitivity
and the DR10 LePhare masses. Results on 10k galaxies of the local strip are in
`benchmarks/results/local_strip/report.md`:
- every code agrees with LePhare within ±0.05 dex in the median (NMAD 0.09–0.10 dex);
- the photo-z error alone moves masses by 0.07–0.11 dex;
- DR11 LePhare differs from DR10 LePhare by +0.03 dex (NMAD 0.09).

**Production masses.** LePhare, CIGALE and kcorrect all run (`sed.codes`). Every DATA file carries
`LOGMSTAR_<CODE>[_LO|_HI]` and `LOGSFR_<CODE>` for each code. LePhare (`sed.primary`) defines
`LOGMSTAR` / `LPH_MASS_BEST` and the M* samples. Running `05_vlim.py --code cigale` builds the M*
samples from CIGALE instead, in `vlim_cigale/`. On CC, `04_stellar_mass.py --prepare` builds the
LePhare libraries (2 GB, about 30 min) and registers the CIGALE filters once, before the job arrays.

**Volume-limited samples (step 5).** These are built in Mr (thresholds −18 … −22.5) and in M*
(9.0 … 11.5). The completeness limits come from the data:
- Mr: the K_95 method above;
- M*: Pozzetti et al. (2010).

Each sample has three files, named
`LS11_VLIM_ANY_<lo>_<Mr|Mstar>_<hi>_<zmin>_z_<zmax>_N_<N>`:
- `_DATA.fits`, the galaxies;
- `_RAND.fits`, with 20× more randoms than galaxies and shuffled redshifts;
- `_COLOUR.fits`.

**Downstream (step 6).**
- *sys_mapping*: `--catalog-dir $LS11_OUT/vlim --template-dir $LS11_OUT/systematics/<nside>`.
- *sum_stat*: `--survey custom --data-file … --rand-file …` (`BEST_Z`, `LPH_MASS_BEST`,
  `RAND.Z`).

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
| `LEPHAREDIR`, `LEPHAREWORK` | LePhare data and libraries (default `~/.cache/lephare/{data,work}`) |
| `EAZY_DATA` | eazy-photoz templates and filters (default `~/.cache/eazy-photoz`, cloned when missing) |
| `DSPS_DRN` | DSPS SSP file `ssp_data_fsps_v3.2_lgmet_age.h5` (default `~/.cache/dsps`) |

## Install

```bash
pip install -e .          # numpy, scipy, astropy, fitsio, healpy, pyyaml
pip install -e .[sed]     # SED-fitting codes for steps 3-4
pytest
```

Maintainer setup: on the laptop, use the shared `dev` env (`~/software/dev_env`). Do not create
a new environment.
