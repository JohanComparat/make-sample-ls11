"""HEALPix maps from the footprint randoms: area fraction and sys_mapping templates.

Layout under $LS11_OUT/<tag>::

  footprint/LS11_FRACAREA_NSIDE_<nside:04d>.fits         fraction of each pixel in the footprint
  systematics/<nside:04d>/LS11_<Q>_NSIDE_<nside:04d>.fits  mean of Q over the footprint randoms
  systematics/<nside:04d>/GAIA_nstar_faint_NSIDE_<nside:05d>.fits   copied full-sky Gaia maps

All maps are RING-ordered, equatorial (COORDSYS='C'), one column named after the quantity, UNSEEN
outside the footprint: the conventions of ~/data/legacysurvey/dr10/systematics read by sys_mapping
(``load_templates_from_dir`` reads every ``*.fits`` of a systematics/<nside> directory).
"""

from __future__ import annotations

import logging
import shutil
from collections.abc import Mapping
from pathlib import Path

import healpy as hp
import numpy as np

from .photometry import depth_mag

log = logging.getLogger(__name__)


def pixel_index(ra, dec, nside: int) -> np.ndarray:
    return hp.ang2pix(nside, np.asarray(ra), np.asarray(dec), lonlat=True)  # RING


def mean_map(pix: np.ndarray, values: np.ndarray, nside: int) -> np.ndarray:
    npix = hp.nside2npix(nside)
    good = np.isfinite(values)
    cnt = np.bincount(pix[good], minlength=npix)
    s = np.bincount(pix[good], weights=values[good], minlength=npix)
    out = np.full(npix, hp.UNSEEN)
    out[cnt > 0] = s[cnt > 0] / cnt[cnt > 0]
    return out


def fracarea_map(pix: np.ndarray, nside: int, density: float) -> np.ndarray:
    npix = hp.nside2npix(nside)
    cnt = np.bincount(pix, minlength=npix).astype(float)
    expected = density * hp.nside2pixarea(nside, degrees=True)
    out = np.full(npix, hp.UNSEEN)
    out[cnt > 0] = cnt[cnt > 0] / expected
    return out


def write_map(path: Path, m: np.ndarray, column: str, extra_header: Mapping | None = None) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    hp.write_map(str(path), m, nest=False, coord="C", column_names=[column], dtype=np.float64,
                 overwrite=True, extra_header=list((extra_header or {}).items()))
    return path


def quantity_values(rand: Mapping, q: str, depth_as_mag: bool) -> np.ndarray:
    v = np.asarray(rand[q], np.float64)
    if depth_as_mag and "DEPTH" in q:
        v = depth_mag(v)
    return v


class MapAccumulator:
    """Per-pixel counts and sums of the footprint randoms, filled chunk by chunk (nothing stored)."""

    def __init__(self, cfg: dict):
        mc = cfg["maps"]
        self.cfg = cfg
        self.nsides = list(mc["nsides"])
        self.quantities = list(mc["quantities"])
        self.as_mag = mc.get("depth_as_mag", True)
        self.count = {n: np.zeros(hp.nside2npix(n)) for n in self.nsides}
        self.sums = {n: {q: np.zeros(hp.nside2npix(n)) for q in self.quantities} for n in self.nsides}
        self.nval = {n: {q: np.zeros(hp.nside2npix(n)) for q in self.quantities} for n in self.nsides}

    def columns(self) -> list[str]:
        return ["RA", "DEC", *self.quantities]

    def add(self, r: Mapping) -> None:
        for n in self.nsides:
            pix = pixel_index(r["RA"], r["DEC"], n)
            npix = hp.nside2npix(n)
            self.count[n] += np.bincount(pix, minlength=npix)
            for q in self.quantities:
                if q not in r:
                    continue
                v = quantity_values(r, q, self.as_mag)
                good = np.isfinite(v)
                self.sums[n][q] += np.bincount(pix[good], weights=v[good], minlength=npix)
                self.nval[n][q] += np.bincount(pix[good], minlength=npix)

    def write(self, out: Path, density: float, gaia_dir: Path | None = None) -> list[Path]:
        written = []
        for n in self.nsides:
            expected = density * hp.nside2pixarea(n, degrees=True)
            frac = np.full(self.count[n].size, hp.UNSEEN)
            seen = self.count[n] > 0
            frac[seen] = self.count[n][seen] / expected
            written.append(write_map(out / "footprint" / f"LS11_FRACAREA_NSIDE_{n:04d}.fits", frac, "FRACAREA",
                                     {"DENSITY": density}))
            sdir = out / "systematics" / f"{n:04d}"
            for q in self.quantities:
                m = np.full(self.count[n].size, hp.UNSEEN)
                ok = self.nval[n][q] > 0
                m[ok] = self.sums[n][q][ok] / self.nval[n][q][ok]
                unit = "AB mag, 5 sigma" if (self.as_mag and "DEPTH" in q) else ""
                written.append(write_map(sdir / f"LS11_{q}_NSIDE_{n:04d}.fits", m, q, {"UNITS": unit}))
            written += _copy_gaia(self.cfg, n, sdir, gaia_dir)
        return written


def _copy_gaia(cfg: dict, nside: int, sdir: Path, gaia_dir: Path | None) -> list[Path]:
    out = []
    if gaia_dir is None:
        return out
    for name in cfg["maps"].get("gaia_templates") or []:
        src = Path(gaia_dir) / f"{nside:04d}" / f"GAIA_{name}_NSIDE_{nside:05d}.fits"
        if src.exists():
            dst = sdir / src.name
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, dst)
            out.append(dst)
        else:
            log.warning("no Gaia map %s", src)
    return out
