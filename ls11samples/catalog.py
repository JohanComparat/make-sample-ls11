"""Assemble per-sweep catalogues from the small product files and the sweeps themselves.

The selection file (<ver>-<tag>/<sweep>-<tag>.fits) holds LS_ID_DR11, SWEEP_ROW, BEST_Z, BEST_Z_ERR,
Z_SOURCE and STAR_FLAG; the code files (<ver>-<code>/<sweep>-<code>.fits) hold SWEEP_ROW, LOGMSTAR,
LOGMSTAR_ERR, MABS_R and MABS_R_ERR, row-aligned with the selection. Everything else (positions,
fluxes) is read back from the sweep at SWEEP_ROW, so nothing is duplicated on disk.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np

from . import io
from .env import Paths
from .photometry import mag

CODE_COLUMNS = ("LOGMSTAR", "LOGMSTAR_ERR", "MABS_R", "MABS_R_ERR")


def sweep_rows(sweep: str | Path, rows: np.ndarray, columns: Sequence[str]) -> dict:
    """Columns of ``sweep`` at ``rows`` (sorted, as SWEEP_ROW is)."""
    if len(rows) == 0:
        return {c: np.zeros(0) for c in columns}
    return io.read_table(sweep, list(columns), rows=np.asarray(rows))


def add_mags(t: dict, bands: Sequence[str]) -> dict:
    for b in bands:
        t[f"MAG_{b}"] = mag(t[f"FLUX_{b}"], t[f"MW_TRANSMISSION_{b}"]).astype(np.float32)
    return t


def load_sweep(paths: Paths, sweep: str | Path, tag: str, sweep_columns: Sequence[str] = (),
               codes: Sequence[str] = (), mag_bands: Sequence[str] = ()) -> dict:
    """Selection of one sweep + sweep columns at SWEEP_ROW + <COL>_<CODE> for each code."""
    t = io.read_table(paths.product(sweep, tag))
    cols = list(dict.fromkeys([*sweep_columns, *[f"{p}_{b}" for b in mag_bands for p in ("FLUX", "MW_TRANSMISSION")]]))
    if cols:
        t.update(sweep_rows(sweep, t["SWEEP_ROW"], cols))
    add_mags(t, mag_bands)
    for code in codes:
        c = io.read_table(paths.product(sweep, code))
        if not np.array_equal(c["SWEEP_ROW"], t["SWEEP_ROW"]):
            raise ValueError(f"{paths.product(sweep, code)} is not row-aligned with the {tag} selection")
        for k in CODE_COLUMNS:
            t[f"{k}_{code.upper()}"] = c[k]
    return t


def load(paths: Paths, sweeps: Sequence[Path], tag: str, **kw) -> dict:
    """Concatenation of :func:`load_sweep` over sweeps, with SWEEP_INDEX (position in ``sweeps``)."""
    parts = []
    for i, s in enumerate(sweeps):
        t = load_sweep(paths, s, tag, **kw)
        t["SWEEP_INDEX"] = np.full(len(t["SWEEP_ROW"]), i, np.int32)
        parts.append(t)
    return io.concat(parts)
