"""Milky-Way-corrected AB magnitudes and fluxes from the sweep nanomaggies."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np

ZP = 22.5


def mag(flux, mw=None, zp: float = ZP) -> np.ndarray:
    """AB magnitude of a nanomaggy flux (dereddened when ``mw`` is the MW transmission); NaN if <= 0."""
    f = np.asarray(flux, np.float64)
    if mw is not None:
        f = f / np.asarray(mw, np.float64)
    out = np.full(f.shape, np.nan)
    pos = f > 0
    out[pos] = zp - 2.5 * np.log10(f[pos])
    return out


def dereddened_flux(t: Mapping, band: str) -> tuple[np.ndarray, np.ndarray]:
    """(flux, ivar) in nanomaggies, corrected for MW extinction."""
    mw = np.asarray(t[f"MW_TRANSMISSION_{band}"], np.float64)
    f = np.asarray(t[f"FLUX_{band}"], np.float64) / mw
    iv = np.asarray(t[f"FLUX_IVAR_{band}"], np.float64) * mw**2
    return f, iv


def flux_arrays(t: Mapping, bands: Sequence[str], err_floor_mag: Mapping[str, float] | None = None
                ) -> tuple[np.ndarray, np.ndarray]:
    """(flux, flux_err) arrays of shape (n, nband), dereddened nanomaggies, with a fractional error
    floor 0.4 ln10 * floor_mag * |flux| added in quadrature. Missing data (ivar <= 0) -> err = inf."""
    fl, er = [], []
    for b in bands:
        f, iv = dereddened_flux(t, b)
        e = np.full(f.shape, np.inf)
        good = iv > 0
        e[good] = 1 / np.sqrt(iv[good])
        if err_floor_mag and b in err_floor_mag:
            floor = 0.4 * np.log(10) * err_floor_mag[b] * np.abs(f)
            e = np.where(good, np.hypot(e, floor), e)
        fl.append(f)
        er.append(e)
    return np.stack(fl, axis=1), np.stack(er, axis=1)


def depth_mag(ivar, nsigma: float = 5.0) -> np.ndarray:
    """``nsigma`` point/galaxy depth in AB mag from a depth inverse variance (NaN where ivar <= 0)."""
    iv = np.asarray(ivar, np.float64)
    out = np.full(iv.shape, np.nan)
    pos = iv > 0
    out[pos] = ZP - 2.5 * np.log10(nsigma / np.sqrt(iv[pos]))
    return out
