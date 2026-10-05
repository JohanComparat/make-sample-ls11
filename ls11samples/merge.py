"""Join per-sweep products (K-corrections, stellar masses) to the BGS-like galaxies by LS_ID_DR11."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np

from . import io


def join_by_id(ids: np.ndarray, files: Sequence[Path], columns: Sequence[str] | None = None,
               prefix: str = "") -> dict:
    """Columns of ``files`` (each with LS_ID_DR11) re-ordered to match ``ids``; NaN / 0 where missing."""
    parts = [io.read_table(f, None if columns is None else ["LS_ID_DR11", *columns]) for f in files if Path(f).exists()]
    prod = io.concat(parts)
    if not prod:
        raise FileNotFoundError(f"none of {len(files)} product files exist")
    order = np.argsort(prod["LS_ID_DR11"])
    sid = prod["LS_ID_DR11"][order]
    if np.any(sid[1:] == sid[:-1]):
        raise ValueError("duplicated LS_ID_DR11 in the products")
    pos = np.clip(np.searchsorted(sid, ids), 0, len(sid) - 1)
    found = sid[pos] == ids
    out = {}
    for k, v in prod.items():
        if k == "LS_ID_DR11":
            continue
        v = v[order][pos]
        if v.dtype.kind == "f":
            v = np.where(found.reshape((-1,) + (1,) * (v.ndim - 1)), v, np.nan)
        out[prefix + k] = v
    out[prefix + "FOUND"] = found
    return out
