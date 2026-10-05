"""The one redshift used by every later stage (K-corrections, masses, Mr, VLIM samples, randoms).

BEST_Z = Z_SPEC when a valid spectroscopic redshift exists (Z_SOURCE = 1), else the first valid
photo-z in ``redshift.photoz_priority`` (Z_SOURCE = 2, 3, ...); Z_SOURCE = 0: no redshift.
BEST_Z_ERR is 0 for spectroscopic redshifts and the matching Z_PHOT_STD* otherwise.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

Z_SOURCE_NONE, Z_SOURCE_SPEC = 0, 1


def std_column(zcol: str) -> str:
    """Z_PHOT_MEAN_I -> Z_PHOT_STD_I, Z_PHOT_MEDIAN -> Z_PHOT_STD."""
    return zcol.replace("_MEDIAN", "_STD").replace("_MEAN", "_STD")


def pz_columns(cfg: dict) -> list[str]:
    cols = ["Z_SPEC", "SURVEY"]
    for c in cfg["redshift"]["photoz_priority"]:
        cols += [c, std_column(c)]
    return cols


def valid_spec(pz: Mapping, cfg: dict) -> np.ndarray:
    zmin, zmax = cfg["redshift"]["zspec_range"]
    z = np.asarray(pz["Z_SPEC"], np.float64)
    ok = np.isfinite(z) & (z > zmin) & (z < zmax)
    excl = cfg["redshift"].get("zspec_exclude_surveys") or []
    if excl:
        ok &= ~np.isin(np.char.strip(np.asarray(pz["SURVEY"]).astype("U")), excl)
    return ok


def best_z(pz: Mapping, cfg: dict) -> dict[str, np.ndarray]:
    n = len(pz["Z_SPEC"])
    z = np.full(n, np.nan)
    zerr = np.full(n, np.nan)
    src = np.zeros(n, np.int16)
    for code, col in enumerate(cfg["redshift"]["photoz_priority"][::-1]):
        code = len(cfg["redshift"]["photoz_priority"]) + 1 - code
        zp = np.asarray(pz[col], np.float64)
        ok = np.isfinite(zp) & (zp > -1)
        z[ok], zerr[ok], src[ok] = zp[ok], np.asarray(pz[std_column(col)], np.float64)[ok], code
    spec = valid_spec(pz, cfg)
    z[spec], zerr[spec], src[spec] = np.asarray(pz["Z_SPEC"], np.float64)[spec], 0.0, Z_SOURCE_SPEC
    return {"BEST_Z": z.astype(np.float32), "BEST_Z_ERR": zerr.astype(np.float32), "Z_SOURCE": src}


def z_source_names(cfg: dict) -> dict[int, str]:
    names = {Z_SOURCE_NONE: "NONE", Z_SOURCE_SPEC: "Z_SPEC"}
    for i, col in enumerate(cfg["redshift"]["photoz_priority"]):
        names[i + 2] = col
    return names
