"""Shared input preparation for the SED backends."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from ..photometry import flux_arrays

BANDS = ("G", "R", "I", "Z", "W1", "W2")
NANOMAGGY_CGS = 3.631e-29      # erg/s/cm^2/Hz
NANOMAGGY_MJY = 3.631e-3       # mJy


def photometry(t: Mapping, cfg: dict, bands=BANDS) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(z, flux, err) with dereddened nanomaggy fluxes (n, nband) and errors including the
    configured floor; err = inf where a band has no measurement."""
    f, e = flux_arrays(t, bands, cfg["photometry"]["err_floor_mag"])
    return np.asarray(t["BEST_Z"], np.float64), f, e


def empty_result(n: int) -> dict[str, np.ndarray]:
    from . import OUTPUT

    return {k: np.full(n, np.nan, np.float32) for k in OUTPUT}


def valid_z(z: np.ndarray, zmax: float = 1.0) -> np.ndarray:
    return np.isfinite(z) & (z > 0.001) & (z < zmax)
