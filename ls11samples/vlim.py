"""Step 5: completeness limits and volume-limited samples in absolute magnitude and stellar mass.

Completeness limits are measured from the data (no external calibration):

* Mr_lim(z) = r_lim - DM(z) - K_p(z), with K_p the p-th percentile (reddest galaxies) of the
  r-band K-correction in redshift bins: every galaxy brighter than Mr_lim(z) is above the flux
  limit at z. Measured on the observed galaxies, K_p is biased low where the flux limit removes
  red galaxies (z >~ 0.35 for r < 19.5); pass ``kcorr_r`` evaluated on a complete low-z SED set
  at every z (:mod:`ls11samples.kcorr`) to avoid it;
* logM*_lim(z) (Pozzetti et al. 2010): for the faintest ``faint_fraction`` of galaxies in each
  redshift bin, the mass each would have at the flux limit, logM + 0.4 (r - r_lim); the p-th
  percentile of it.

Both curves are made monotonic (Mr_lim brightens, M*_lim grows with z); a threshold's z_max is
where the curve crosses it, capped at ``z_cap``. The redshift range is rounded inwards to the 2
decimals of the sample name (z_max down, z_min up), so the name gives exactly the range selected:
sys_mapping reads it from the name.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping

import numpy as np

from . import randoms as rnd
from .cosmo import distmod, shell_volume

log = logging.getLogger(__name__)


def _binned_percentile(z, v, zbins, pct, mask=None, min_count=20):
    zc = 0.5 * (zbins[1:] + zbins[:-1])
    out = np.full(zc.size, np.nan)
    idx = np.digitize(z, zbins) - 1
    for i in range(zc.size):
        sel = (idx == i) & np.isfinite(v)
        if mask is not None:
            sel &= mask
        if sel.sum() >= min_count:
            out[i] = np.percentile(v[sel], pct)
    return zc, out


def mr_limit_curve(z, kcorr_r, r_lim: float, cfg: dict, pct: float = 95, zbins=None):
    """(z grid, Mr_lim(z)) from the r-band K-corrections of the galaxies."""
    zbins = np.arange(0.0, 0.62, 0.02) if zbins is None else zbins
    zc, kp = _binned_percentile(np.asarray(z), np.asarray(kcorr_r), zbins, pct)
    ok = np.isfinite(kp)
    zc, kp = zc[ok], kp[ok]
    m = r_lim - distmod(zc, cfg) - kp
    return zc, np.minimum.accumulate(m)          # brighter (smaller) with z


def mr_limit_from_kcurve(zgrid, k_pct, r_lim: float, cfg: dict):
    """(z grid, Mr_lim) from a K-correction percentile curve evaluated on a complete SED set."""
    zgrid = np.asarray(zgrid, np.float64)
    return zgrid, np.minimum.accumulate(r_lim - distmod(zgrid, cfg) - np.asarray(k_pct))


def mstar_limit_curve(z, logm, rmag, r_lim: float, pct: float = 95, faint_fraction: float = 0.2,
                      zbins=None):
    """(z grid, logM*_lim(z)), Pozzetti et al. (2010) method."""
    z, logm, rmag = (np.asarray(a, np.float64) for a in (z, logm, rmag))
    zbins = np.arange(0.0, 0.62, 0.02) if zbins is None else zbins
    mlim = logm + 0.4 * (rmag - r_lim)
    zc = 0.5 * (zbins[1:] + zbins[:-1])
    out = np.full(zc.size, np.nan)
    idx = np.digitize(z, zbins) - 1
    for i in range(zc.size):
        sel = (idx == i) & np.isfinite(mlim)
        if sel.sum() < 20:
            continue
        faint = sel & (rmag >= np.percentile(rmag[sel], 100 * (1 - faint_fraction)))
        out[i] = np.percentile(mlim[faint], pct)
    ok = np.isfinite(out)
    return zc[ok], np.maximum.accumulate(out[ok])  # grows with z


def zmax_brighter(zc, mlim, threshold: float) -> float:
    """Largest z where Mr_lim(z) is fainter than ``threshold`` (Mr_lim decreasing)."""
    ok = mlim >= threshold
    if not ok.any():
        return np.nan
    if ok.all():
        return float(zc[-1])
    i = np.flatnonzero(~ok)[0]
    if i == 0:
        return float(zc[0])
    # linear interpolation between the last complete and the first incomplete bin
    return float(np.interp(threshold, [mlim[i], mlim[i - 1]], [zc[i], zc[i - 1]]))


def zmax_heavier(zc, mlim, threshold: float) -> float:
    """Largest z where logM*_lim(z) is below ``threshold`` (M*_lim increasing)."""
    return zmax_brighter(zc, -np.asarray(mlim), -threshold)


def name_precision(z: float, up: bool = False) -> float:
    """``z`` rounded down (``up``: up) to the 2 decimals of the sample names, so that the selection
    and the name agree; NaN stays NaN."""
    if not np.isfinite(z):
        return np.nan
    x = round(z * 100, 6)                         # 0.29 * 100 = 28.999999999999996
    return float((np.ceil(x) if up else np.floor(x)) / 100)


def sample_name(kind: str, lo: float, hi: float, z0: float, z1: float, n: int, prefix: str = "LS11") -> str:
    """LS11_VLIM_ANY_<lo>_<Mstar|Mr>_<hi>_<z0>_z_<z1>_N_<N:07d> (DR10 naming, LS10 -> LS11)."""
    return f"{prefix}_VLIM_ANY_{lo:.2f}_{kind}_{hi:.2f}_{z0:.2f}_z_{z1:.2f}_N_{n:07d}"


def define_samples(data: Mapping, cfg: dict, r_lim: float, mr_curve=None) -> tuple[list[dict], dict]:
    """Sample definitions (selection rows + metadata) for every threshold with a valid z range.

    ``mr_curve``: (z grid, Mr_lim) computed elsewhere (e.g. from a complete SED set,
    :func:`mr_limit_from_kcurve`); default: binned percentiles of the observed K-corrections."""
    vc = cfg["vlim"]
    pct = vc.get("completeness_percentile", 95)
    z = np.asarray(data["BEST_Z"], np.float64)
    out = []
    curves = {}
    if "MABS_R" in data and ("KCORR_R" in data or mr_curve is not None):
        c = vc["absmag_r"]
        zc, ml = mr_curve if mr_curve is not None else mr_limit_curve(z, data["KCORR_R"], r_lim, cfg, pct)
        curves["Mr"] = (zc, ml)
        for thr in c["thresholds"]:
            z0 = name_precision(c["z_min"], up=True)
            z1 = name_precision(min(zmax_brighter(zc, ml, thr), c["z_cap"]))
            if not np.isfinite(z1) or z1 <= z0:
                continue
            m = np.asarray(data["MABS_R"])
            sel = (m <= thr) & (m > c["bright"]) & (z > z0) & (z <= z1)
            out.append({"kind": "Mr", "lo": c["bright"], "hi": thr, "z0": z0, "z1": z1, "sel": sel})
    if "LOGMSTAR" in data:
        c = vc["mstar"]
        zc, ml = mstar_limit_curve(z, data["LOGMSTAR"], data["MAG_R"], r_lim, pct)
        curves["Mstar"] = (zc, ml)
        for thr in c["thresholds"]:
            z0 = name_precision(c["z_min"], up=True)
            z1 = name_precision(min(zmax_heavier(zc, ml, thr), c["z_cap"]))
            if not np.isfinite(z1) or z1 <= z0:
                continue
            m = np.asarray(data["LOGMSTAR"])
            sel = (m >= thr) & (m < c["max"]) & (z > z0) & (z <= z1)
            out.append({"kind": "Mstar", "lo": thr, "hi": c["max"], "z0": z0, "z1": z1, "sel": sel})
    for s in out:
        s["n"] = int(s["sel"].sum())
        s["name"] = sample_name(s["kind"], s["lo"], s["hi"], s["z0"], s["z1"], s["n"])
    return out, curves


def make_randoms(rand: Mapping, z_data: np.ndarray, n_factor: float, rng: np.random.Generator) -> dict:
    """Randoms of one sample: a subsample of the footprint randoms with shuffled data redshifts."""
    n_want = int(round(n_factor * len(z_data)))
    idx = rnd.subsample(len(rand["RA"]), n_want, rng)
    return {"RA": np.asarray(rand["RA"])[idx], "DEC": np.asarray(rand["DEC"])[idx],
            "EBV": np.asarray(rand["EBV"])[idx].astype(np.float32),
            "Z": rnd.shuffle_z(z_data, len(idx), rng).astype(np.float32)}


def summary_row(s: dict, data: Mapping, area: float, cfg: dict) -> dict:
    sel = s["sel"]
    vol = shell_volume(s["z0"], s["z1"], area, cfg)
    z = np.asarray(data["BEST_Z"])[sel]
    row = {"NAME": s["name"], "KIND": s["kind"], "LO": s["lo"], "HI": s["hi"], "Z_MIN": s["z0"],
           "Z_MAX": s["z1"], "N_GAL": s["n"], "VOLUME": vol, "N_DENS": s["n"] / vol if vol else np.nan,
           "N_DEG2": s["n"] / area, "Z_MEDIAN": float(np.median(z)) if z.size else np.nan}
    for col in ("MABS_R", "LOGMSTAR"):
        if col in data:
            v = np.asarray(data[col])[sel]
            row[f"{col}_MEDIAN"] = float(np.nanmedian(v)) if v.size else np.nan
    return row
