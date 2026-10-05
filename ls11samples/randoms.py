"""Step 2: randoms of the BGS-like footprint.

The DR11 random files (desitarget, ``DENSITY`` per deg^2 each) are read in chunks; a random is kept
when it lies in one of the processed sweep boxes and passes :func:`selection.footprint_mask`, the
same function applied to the galaxies. Area = N_kept / (sum of the file densities).
"""

from __future__ import annotations

import logging
import time
from collections.abc import Sequence
from pathlib import Path

import numpy as np

from . import io, selection
from .config import config_hash

log = logging.getLogger(__name__)

class BoxSet:
    """Membership test for sweep boxes; O(1) per point when the boxes are on the 5x5 deg sweep grid."""

    def __init__(self, boxes: Sequence[tuple[float, float, float, float]]):
        self.boxes = list(boxes)
        self.grid = all(np.isclose(ra1 - ra0, 5) and np.isclose(d1 - d0, 5)
                        and np.isclose(ra0 % 5, 0) and np.isclose(d0 % 5, 0)
                        for ra0, ra1, d0, d1 in self.boxes)
        if self.grid:
            self.keys = np.array(sorted({self._key(ra0, d0) for ra0, _, d0, _ in self.boxes}), np.int64)

    @staticmethod
    def _key(ra, dec):
        return (np.floor(np.asarray(ra) / 5).astype(np.int64) * 1000
                + np.floor((np.asarray(dec) + 90) / 5).astype(np.int64))

    def contains(self, ra, dec) -> np.ndarray:
        if not self.grid:
            return io.in_boxes(ra, dec, self.boxes)
        k = self._key(ra, dec)
        i = np.clip(np.searchsorted(self.keys, k), 0, len(self.keys) - 1)
        return self.keys[i] == k

    @property
    def area(self) -> float:
        return io.box_area(self.boxes)


def file_density(path: str | Path, cfg: dict) -> float:
    hdr = io.read_header(path, 1)
    dens = float(hdr.get("DENSITY", cfg["randoms"]["density_per_file"]))
    if not np.isclose(dens, cfg["randoms"]["density_per_file"]):
        log.warning("%s: DENSITY=%s differs from the configured %s", path, dens,
                    cfg["randoms"]["density_per_file"])
    return dens


def select_randoms(files: Sequence[str | Path], cfg: dict, boxes: BoxSet | None = None,
                   keep: Sequence[str] = ("RA", "DEC", "EBV"), maps=None) -> tuple[dict, dict]:
    """(randoms passing the footprint inside ``boxes`` (columns ``keep``), header with area and
    counts). ``maps`` (a :class:`ls11samples.maps.MapAccumulator`) is filled chunk by chunk with
    the kept randoms, so map quantities need not be kept."""
    fcols = selection.footprint_columns(cfg)
    extra = list(dict.fromkeys([*keep, *(maps.columns() if maps is not None else [])]))
    chunk = int(cfg["randoms"]["chunk_rows"])
    parts, n_box, n_pass, dens_total = [], 0, 0, 0.0
    for path in files:
        t0 = time.time()
        dens_total += file_density(path, cfg)
        n_file = 0
        for start, t in io.iter_rows(path, fcols, chunk):
            inb = np.ones(len(t["RA"]), bool) if boxes is None else boxes.contains(t["RA"], t["DEC"])
            sub = io.take(t, inb)
            ok = selection.footprint_mask(sub, cfg)
            n_box += int(inb.sum())
            rows = start + np.flatnonzero(inb)[ok]
            if rows.size == 0:
                continue
            r = io.read_table(path, extra, rows=rows)
            if maps is not None:
                maps.add(r)
            parts.append({k: r[k] for k in keep})
            n_file += rows.size
        n_pass += n_file
        log.info("%s: %d kept (%.0fs)", Path(path).name, n_file, time.time() - t0)
    rand = io.concat(parts) or {k: np.zeros(0) for k in keep}
    for k in ("RA", "DEC"):
        if k in rand:
            rand[k] = rand[k].astype(np.float64)
    if "EBV" in rand:
        rand["EBV"] = rand["EBV"].astype(np.float32)
    area = n_pass / dens_total if dens_total else np.nan
    header = {"NRFILES": len(files), "DENSITY": dens_total, "NINBOX": n_box, "NRAND": n_pass,
              "AREA": area, "BOXAREA": boxes.area if boxes is not None else -1.0, "TAG": cfg["tag"],
              "FPHASH": config_hash(cfg, "footprint"), "CFGHASH": config_hash(cfg),
              "MASKREJ": ",".join(selection.maskbits_reject(cfg))}
    for i, f in enumerate(files):
        header[f"RFILE{i}"] = Path(f).name
    return rand, header


def subsample(n_avail: int, n_want: int, rng: np.random.Generator) -> np.ndarray:
    """Indices of a random subsample (all of them, in random order, if fewer than wanted)."""
    if n_want >= n_avail:
        if n_want > n_avail:
            log.warning("asked for %d randoms, only %d available", n_want, n_avail)
        return rng.permutation(n_avail)
    return np.sort(rng.choice(n_avail, n_want, replace=False))


def shuffle_z(z_data: np.ndarray, n: int, rng: np.random.Generator, weights=None) -> np.ndarray:
    """Redshifts for ``n`` randoms drawn from the data redshifts (shuffle method)."""
    p = None if weights is None else np.asarray(weights, float) / np.sum(weights)
    return rng.choice(np.asarray(z_data), size=n, replace=True, p=p)
