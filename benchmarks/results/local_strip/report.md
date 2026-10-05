# Stellar-mass benchmark, LS DR11 BGS-like (local strip)

10000 galaxies (5000 with spectroscopic redshift), fixed z = BEST_Z, grizW1W2 with the same filter curves and error floor for every code. Reference: **lephare**.

## Cost

Rerun on 2026-10-05 after the CIGALE fix (`additionalerror` = 0: the same error floor as the other
codes). The laptop was shared and loaded (load ≈ 26 on 24 cores, 8 threads per code), and kcorrect
now includes its 20 Monte Carlo realisations. Idle-machine timings from the first run, per 10k
galaxies: LePhare 118 s on 23 threads, CIGALE 58 s on 16, eazy 25 s on 16, DSPS 1.7 s, kcorrect
0.2 s without Monte Carlo. On CC-IN2P3, LePhare took 52 min for 97.5k objects on 16 cores.

| code | init [s] | fit [s] | threads | ms / galaxy / thread | CPU-h for 20M |
|---|---|---|---|---|---|
| kcorrect | 6.8 | 9.1 | 1 | 0.91 | 5.0 |
| lephare | 0.8 | 506.1 | 8 | 404.85 | 2,249 |
| cigale | 2.1 | 102.3 | 8 | 81.81 | 454 |
| eazy | 0.7 | 77.5 | 8 | 62.00 | 344 |
| dsps | 14.9 | 3.3 | 24 | 7.81 | 43 |

## Agreement with the reference

Δ = log M*(code) − log M*(ref); NMAD = 1.4826 median|Δ − median Δ|.

| code | failed [%] | median χ² | median Δ | NMAD | median Δ (spec-z) | median Δ, g−r>0.8 | median Δ, g−r<0.6 | median ΔMr |
|---|---|---|---|---|---|---|---|---|
| kcorrect | 0.00 | 4.05 | +0.043 | 0.087 | +0.036 | +0.016 | +0.122 | +0.004 |
| lephare | 0.00 | 4.16 | +0.000 | 0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| cigale | 0.23 | 4.81 | -0.016 | 0.066 | -0.013 | -0.007 | -0.029 | +0.004 |
| eazy | 0.00 | 8.66 | +0.030 | 0.089 | +0.023 | +0.009 | +0.122 | +0.003 |
| dsps | 0.00 | 5.02 | +0.046 | 0.091 | +0.049 | +0.056 | +0.021 | -0.003 |

Median Δ in redshift bins:

| code | 0.0–0.1 | 0.1–0.2 | 0.2–0.3 | 0.3–0.4 | 0.4–0.5 | 0.5–0.6 |
|---|---|---|---|---|---|---|
| kcorrect | +0.114 | +0.048 | +0.026 | +0.027 | +0.048 | +0.070 |
| lephare | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| cigale | -0.020 | -0.022 | -0.020 | -0.011 | +0.003 | +0.005 |
| eazy | +0.137 | +0.040 | +0.012 | +0.009 | -0.007 | +0.032 |
| dsps | +0.063 | +0.027 | +0.043 | +0.066 | +0.078 | +0.058 |

## Photo-z sensitivity

Refit of photo-z galaxies at BEST_Z ± BEST_Z_ERR: half the difference of the two log M*, i.e. the mass error from the photo-z alone.

| code | N | median ½\|Δ log M*\| | 84th pct |
|---|---|---|---|
| kcorrect | 1000 | 0.110 | 0.216 |
| lephare | 994 | 0.094 | 0.210 |
| cigale | 999 | 0.081 | 0.198 |
| eazy | 993 | 0.086 | 0.197 |
| dsps | 1000 | 0.077 | 0.212 |

## Comparison with the DR10 LePhare masses

4652 galaxies matched within 1″ with |Δz| < 0.005 between DR10 and DR11 BEST_Z.

| code | median log M*(DR11 code) − DR10 LPH_MASS_BEST | NMAD |
|---|---|---|
| kcorrect | +0.069 | 0.159 |
| lephare | +0.030 | 0.090 |
| cigale | +0.049 | 0.124 |
| eazy | +0.059 | 0.160 |
| dsps | +0.108 | 0.118 |

## Figures

![dlogm_vs_logm.png](dlogm_vs_logm.png)
![dlogm_vs_z.png](dlogm_vs_z.png)
![cost.png](cost.png)
