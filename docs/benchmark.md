# Stellar-mass benchmark

`benchmarks/sed_benchmark.py` compares the five SED backends of {doc}`stellar_masses`. Every code
fits the same galaxies, with the same photometry (dereddened grizW1W2 and error floor), filter
curves and redshift (`BEST_Z`). The script measures:

- the cost of each code;
- its agreement with a reference code;
- its sensitivity to the photo-z error, from refits of a photo-z subset at `BEST_Z` ± `BEST_Z_ERR`;
- its offset from the DR10 LePhare masses.

```bash
python benchmarks/sed_benchmark.py --codes kcorrect,lephare,cigale,eazy,dsps --n 10000
```

| option | default | meaning |
|---|---|---|
| `--codes` | all five | codes to run |
| `--ref` | `lephare` | reference code of the comparisons |
| `--n` | 10000 | galaxies: half with a spectroscopic redshift, the rest photo-z |
| `--n-zsens` | 1000 | photo-z galaxies refitted at `BEST_Z` ± `BEST_Z_ERR` |
| `--seed` | 2 | sample seed |
| `--out` | `$LS11_OUT/<tag>/benchmark` | output directory: per-code FITS files, `report.md`, figures |
| `--dr10` | DR10 BGS-like LePhare catalogue | for the DR10 comparison |
| `--overwrite` | | refit codes whose output exists |

The report of the run on the local DR11 strip (0 < RA < 5, −15 < Dec < 0) follows. It is
included from `benchmarks/results/local_strip/report.md`.

:::{note}
This run used pcigale's default `additionalerror` (10% of the flux added to every band). In
production it is 0, and the CIGALE masses are higher by about 0.05 dex than the ones below.
:::

```{include} ../benchmarks/results/local_strip/report.md
:heading-offset: 1
:relative-images:
```
