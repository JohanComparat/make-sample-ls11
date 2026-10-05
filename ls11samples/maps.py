"""HEALPix maps from the footprint randoms: area fraction and sys_mapping templates.

Layout under $LS11_OUT:

  footprint/LS11_FRACAREA_NSIDE_<nside:04d>.fits         fraction of each pixel in the footprint
  systematics/<nside:04d>/LS11_<Q>_NSIDE_<nside:04d>.fits  mean of Q over the footprint randoms
  systematics/<nside:04d>/GAIA_nstar_faint_NSIDE_<nside:05d>.fits   copied full-sky Gaia maps

All maps are RING-ordered, equatorial (COORDSYS='C'), one column named after the quantity, UNSEEN
outside the footprint: the conventions of ~/data/legacysurvey/dr10/systematics read by sys_mapping
(``load_templates_from_dir`` reads every *.fits of a systematics/<nside> directory).
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


def make_maps(rand: Mapping, cfg: dict, out: Path, density: float, gaia_dir: Path | None = None
              ) -> list[Path]:
    mc = cfg["maps"]
    written = []
    for nside in mc["nsides"]:
        pix = pixel_index(rand["RA"], rand["DEC"], nside)
        written.append(write_map(out / "footprint" / f"LS11_FRACAREA_NSIDE_{nside:04d}.fits",
                                 fracarea_map(pix, nside, density), "FRACAREA",
                                 {"DENSITY": density}))
        sdir = out / "systematics" / f"{nside:04d}"
        for q in mc["quantities"]:
            if q not in rand:
                log.warning("randoms have no %s; map skipped", q)
                continue
            m = mean_map(pix, quantity_values(rand, q, mc.get("depth_as_mag", True)), nside)
            unit = "AB mag, 5 sigma" if (mc.get("depth_as_mag", True) and "DEPTH" in q) else ""
            written.append(write_map(sdir / f"LS11_{q}_NSIDE_{nside:04d}.fits", m, q, {"UNITS": unit}))
        if gaia_dir is not None:
            for name in mc.get("gaia_templates") or []:
                src = Path(gaia_dir) / f"{nside:04d}" / f"GAIA_{name}_NSIDE_{nside:05d}.fits"
                if src.exists():
                    dst = sdir / src.name
                    shutil.copyfile(src, dst)
                    written.append(dst)
                else:
                    log.warning("no Gaia map %s", src)
    return written
