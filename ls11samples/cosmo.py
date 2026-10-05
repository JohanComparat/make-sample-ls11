"""Cosmology from the configuration, with interpolated distance moduli for large arrays."""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from astropy.cosmology import FlatLambdaCDM


def cosmology(cfg: dict) -> FlatLambdaCDM:
    c = cfg["cosmology"]
    return _cosmo(float(c["H0"]), float(c["Om0"]))


@lru_cache(maxsize=8)
def _cosmo(H0: float, Om0: float) -> FlatLambdaCDM:
    return FlatLambdaCDM(H0=H0, Om0=Om0)


@lru_cache(maxsize=8)
def _dm_grid(H0: float, Om0: float, zmax: float = 3.0, n: int = 6001):
    z = np.linspace(1e-4, zmax, n)
    return z, _cosmo(H0, Om0).distmod(z).value


def distmod(z, cfg: dict) -> np.ndarray:
    """Distance modulus (mag); NaN for z <= 0 or non-finite z."""
    c = cfg["cosmology"]
    zg, dm = _dm_grid(float(c["H0"]), float(c["Om0"]))
    z = np.asarray(z, np.float64)
    out = np.full(z.shape, np.nan)
    ok = np.isfinite(z) & (z > 0)
    out[ok] = np.interp(z[ok], zg, dm)
    return out


def shell_volume(z0: float, z1: float, area_deg2: float, cfg: dict) -> float:
    """Comoving volume (Mpc^3) between z0 and z1 over ``area_deg2``."""
    cos = cosmology(cfg)
    fsky = area_deg2 / (4 * np.pi * (180 / np.pi) ** 2)
    return float((cos.comoving_volume(z1) - cos.comoving_volume(z0)).value * fsky)
