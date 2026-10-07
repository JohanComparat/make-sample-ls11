"""A small synthetic Legacy Surveys DR11 tree for the end-to-end tests.

``build(root)`` writes, under ``root``::

    dr11/south/sweep/11.0/sweep-000m005-005p000.fits      full bit list in the primary header
    dr11/south/sweep/11.0/sweep-005m005-010p000.fits      header without MBIT_19, as 2 DR11 sweeps
    dr11/south/sweep/11.0/sweep-010m005-015p000.fits      no photo-z file: must be skipped
    dr11/south/sweep/11.0-photo-z/<sweep>-pz.fits          row-matched photo-z of the first two
    dr11/south/randoms/randoms-south-1-0.fits             DENSITY randoms over the three boxes
    gaia/<nside:04d>/GAIA_nstar_faint_NSIDE_<nside:05d>.fits

Objects are galaxies (85%) with colours that redden with redshift, and stars (TYPE=PSF, Gaia
match). A few percent carry MASKBITS / FITBITS, NOBS = 0, E(B-V) > 0.1 or r > 21, so every cut has
something to remove.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.io import fits

from ls11samples import bits, io

BANDS = ("G", "R", "I", "Z", "W1", "W2")
EXT = {"G": 3.214, "R": 2.165, "I": 1.592, "Z": 1.211, "W1": 0.184, "W2": 0.113}
SWEEPS = ("sweep-000m005-005p000.fits", "sweep-005m005-010p000.fits", "sweep-010m005-015p000.fits")
DENSITY = 2000.0


def _header(drop_bits=()):
    h = fits.Header()
    for b, n in bits.MASKBITS.items():
        if b not in drop_bits:
            h[f"MBIT_{b}"] = n
    for b, n in bits.FITBITS.items():
        h[f"FBIT_{b}"] = n
    h["DRVERSIO"] = 11010
    return h


def _write(path: Path, table: dict, primary: fits.Header | None = None, extname="SWEEP", hdr1=None):
    from astropy.table import Table

    path.parent.mkdir(parents=True, exist_ok=True)
    t = fits.table_to_hdu(Table({k: np.asarray(v) for k, v in table.items()}))
    t.name = extname
    for k, v in (hdr1 or {}).items():
        t.header[k] = v
    fits.HDUList([fits.PrimaryHDU(header=primary), t]).writeto(path, overwrite=True)


def _objects(rng, n, box, brick0):
    ra0, ra1, d0, d1 = box
    star = rng.random(n) < 0.15
    rmag = rng.uniform(14.0, 21.8, n)
    z = np.clip(0.03 + 0.07 * (rmag - 14) + rng.normal(0, 0.04, n), 0.01, 0.95)
    z[star] = 0.0
    col = {"G": 0.6 + 1.2 * z, "R": 0.0, "I": -(0.25 + 0.5 * z), "Z": -(0.4 + 0.8 * z),
           "W1": -(0.2 + 1.6 * z), "W2": -(-0.2 + 1.6 * z)}          # m_b - m_r
    col_star = {"G": 0.5, "R": 0.0, "I": -0.2, "Z": -0.3, "W1": 0.2, "W2": 0.6}
    ebv = rng.uniform(0.01, 0.12, n)
    t = {"RELEASE": np.full(n, 11010, np.int16),
         "BRICKID": (brick0 + rng.integers(0, 50, n)).astype(np.int32),
         "OBJID": np.arange(n, dtype=np.int32),
         "TYPE": np.where(star, "PSF", rng.choice(["REX", "EXP", "DEV", "SER"], n)),
         "RA": rng.uniform(ra0, ra1, n), "DEC": rng.uniform(d0, d1, n), "EBV": ebv.astype(np.float32)}
    for b in BANDS:
        m = rmag + np.where(star, col_star[b], np.asarray(col[b]) if np.ndim(col[b]) else col[b])
        mw = 10 ** (-0.4 * EXT[b] * ebv)
        f = 10 ** (-0.4 * (m - 22.5)) * mw                          # observed, not dereddened
        sig = 0.02 * f + (0.3 if b in "GRIZ" else 3.0)
        t[f"FLUX_{b}"] = (f + rng.normal(0, sig)).astype(np.float32)
        t[f"FLUX_IVAR_{b}"] = (1 / sig**2).astype(np.float32)
        t[f"MW_TRANSMISSION_{b}"] = mw.astype(np.float32)
    for b in "GRIZ":
        t[f"NOBS_{b}"] = rng.choice([0, 1, 2, 3, 4], n, p=[0.02, 0.2, 0.3, 0.3, 0.18]).astype(np.int16)
        t[f"GALDEPTH_{b}"] = rng.uniform(200, 2000, n).astype(np.float32)
        t[f"PSFDEPTH_{b}"] = rng.uniform(300, 3000, n).astype(np.float32)
        t[f"PSFSIZE_{b}"] = rng.uniform(1.0, 1.6, n).astype(np.float32)
    t["PSFDEPTH_W1"] = rng.uniform(1, 5, n).astype(np.float32)
    t["PSFDEPTH_W2"] = rng.uniform(0.5, 2, n).astype(np.float32)
    for q, lo, hi in (("FRACMASKED", 0, 0.5), ("FRACIN", 0.2, 1.0), ("FRACFLUX", 0, 6)):
        for b in "GRZ":
            t[f"{q}_{b}"] = rng.uniform(lo, hi, n).astype(np.float32)
    t["FIBERFLUX_R"] = (0.4 * t["FLUX_R"]).astype(np.float32)
    t["FIBERTOTFLUX_R"] = (0.45 * t["FLUX_R"]).astype(np.float32)
    t["SHAPE_R"] = rng.uniform(0, 3, n).astype(np.float32)
    t["SERSIC"] = rng.uniform(0.5, 4, n).astype(np.float32)
    t["REF_CAT"] = np.where(star, "G3", "")
    t["GAIA_PHOT_G_MEAN_MAG"] = np.where(star, rmag + 0.3, 0.0).astype(np.float32)
    t["MASKBITS"] = rng.choice([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2048, 4096, 1 << 16],
                               n).astype(np.int32)
    t["FITBITS"] = rng.choice([0] * 30 + [2, 512, 4096], n).astype(np.int16)
    return t, z, star


def _pz(rng, t, z, star):
    n = len(z)
    spec = (~star) & (rng.random(n) < 0.15)
    zp = np.where(star, rng.uniform(0, 2.0, n), np.clip(z + rng.normal(0, 0.03, n), 0.005, 1.5))
    zp_i = np.where(rng.random(n) < 0.05, -99.0, zp + rng.normal(0, 0.005, n))
    ids = io.ls_id(t["RELEASE"], t["BRICKID"], t["OBJID"])
    return {"LS_ID_DR11": ids, "RELEASE": t["RELEASE"], "BRICKID": t["BRICKID"], "OBJID": t["OBJID"],
            "Z_SPEC": np.where(spec, z, -99.0).astype(np.float32),
            "SURVEY": np.where(spec, "DESI", ""),
            "Z_PHOT_MEAN": zp.astype(np.float32), "Z_PHOT_STD": np.full(n, 0.03, np.float32),
            "Z_PHOT_MEAN_I": zp_i.astype(np.float32), "Z_PHOT_STD_I": np.full(n, 0.025, np.float32)}


def _randoms(rng, boxes):
    parts = []
    for ra0, ra1, d0, d1 in boxes:
        area = io.box_area([(ra0, ra1, d0, d1)])
        n = int(DENSITY * area)
        sd0, sd1 = np.sin(np.radians(d0)), np.sin(np.radians(d1))
        parts.append({"RA": rng.uniform(ra0, ra1, n), "DEC": np.degrees(np.arcsin(rng.uniform(sd0, sd1, n)))})
    r = io.concat(parts)
    n = len(r["RA"])
    r["EBV"] = rng.uniform(0.01, 0.12, n).astype(np.float32)
    r["MASKBITS"] = rng.choice([0] * 18 + [2, 2048, 4096, 1 << 16], n).astype(np.int32)
    for b in "GRIZ":
        r[f"NOBS_{b}"] = rng.choice([0, 1, 2, 3], n, p=[0.02, 0.3, 0.4, 0.28]).astype(np.int16)
        r[f"GALDEPTH_{b}"] = rng.uniform(200, 2000, n).astype(np.float32)
        r[f"PSFDEPTH_{b}"] = rng.uniform(300, 3000, n).astype(np.float32)
        r[f"PSFSIZE_{b}"] = rng.uniform(1.0, 1.6, n).astype(np.float32)
    r["PSFDEPTH_W1"] = rng.uniform(1, 5, n).astype(np.float32)
    r["PSFDEPTH_W2"] = rng.uniform(0.5, 2, n).astype(np.float32)
    r["PHOTSYS"] = np.full(n, "S")
    return r


def build(root: str | Path, n_per_sweep: int = 3000, seed: int = 7) -> Path:
    """Write the synthetic tree under ``root``; returns the DR11 root (``root/dr11``)."""
    root = Path(root)
    rng = np.random.default_rng(seed)
    dr11 = root / "dr11"
    sw = dr11 / "south" / "sweep"
    boxes = []
    for i, name in enumerate(SWEEPS):
        box = io.sweep_box(name)
        boxes.append(box)
        t, z, star = _objects(rng, n_per_sweep, box, brick0=1000 * (i + 1))
        _write(sw / "11.0" / name, t, primary=_header(drop_bits=(19,) if i == 1 else ()))
        if i < 2:                                               # the third has no photo-z yet
            _write(sw / "11.0-photo-z" / name.replace(".fits", "-pz.fits"), _pz(rng, t, z, star))
    rhdr = fits.Header()
    rhdr["DENSITY"] = DENSITY
    _write(dr11 / "south" / "randoms" / "randoms-south-1-0.fits", _randoms(rng, boxes),
           extname="RANDOMS", hdr1=rhdr)
    import healpy as hp

    for nside in (8, 16):
        m = np.full(hp.nside2npix(nside), 100.0) + rng.uniform(0, 50, hp.nside2npix(nside))
        path = root / "gaia" / f"{nside:04d}" / f"GAIA_nstar_faint_NSIDE_{nside:05d}.fits"
        path.parent.mkdir(parents=True, exist_ok=True)
        hp.write_map(str(path), m, column_names=["nstar_faint"], coord="C", overwrite=True, dtype=np.float64)
    return dr11


CONFIG = """\
# small configuration of the end-to-end tests (inherits the default)
base: {default}
tag: testsel
galaxy:
  r_range: [13.0, 21.0]
  reject_types: []
  flux_ivar_positive: []
  gaia_star_cut: null
  colour: null
  quality: null
  fiber: null
  rfibtot: null
randoms:
  density_per_file: {density}
  chunk_rows: 20000
maps:
  nsides: [8, 16]
  quantities: [GALDEPTH_R, PSFDEPTH_W1, PSFSIZE_R, NOBS_R, EBV]
kcorr:
  nz: 200
  mc: 3
sed:
  codes: [kcorrect]
  primary: kcorrect
vlim:
  n_rand_factor: 2
  mr_code: kcorrect
  exclude_star_flag: 3
  mstar:
    thresholds: [9.5, 10.0, 10.5, 11.0]
  absmag_r:
    thresholds: [-19.0, -20.0, -21.0, -22.0]
"""
