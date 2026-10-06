"""Step 1: one sweep -> selection file <ver>-<tag>/<sweep>-<tag>.fits (+ HDU CUTFLOW).

Only the selected objects are written, with the minimum needed to find them again and to fit them::

  LS_ID_DR11 (i8)   RELEASE<<42 | BRICKID<<22 | OBJID
  SWEEP_ROW  (i4)   row in the sweep (and in the row-matched photo-z file)
  BEST_Z, BEST_Z_ERR (f4), Z_SOURCE (i1)   see ls11samples.redshift
  STAR_FLAG  (u1)   bit 0 TYPE=PSF, bit 1 Gaia star test (information, never a cut)

Positions and photometry stay in the sweep (read back at SWEEP_ROW, ls11samples.catalog).
"""

from __future__ import annotations

import logging
import time
from pathlib import Path

import numpy as np

from . import bits, io, redshift, selection
from .config import config_hash
from .photometry import mag

log = logging.getLogger(__name__)


def _bit_list(names, table) -> str:
    return ",".join(f"{n}({table[n]})" for n in names)


def select_sweep(sweep: str | Path, cfg: dict) -> tuple[dict, dict, dict]:
    """(selection table, cut-flow table, header) of one sweep."""
    t0 = time.time()
    sweep = Path(sweep)
    bits.check_header(io.read_header(sweep, 0),
                      used={"MBIT": selection.maskbits_reject(cfg), "FBIT": cfg["galaxy"].get("fitbits_reject") or []})
    n_all = io.nrows(sweep)
    r_min, r_max = cfg["galaxy"]["r_range"]

    pre = io.read_table(sweep, ["FLUX_R", "MW_TRANSMISSION_R"])
    with np.errstate(invalid="ignore"):
        r = mag(pre["FLUX_R"], pre["MW_TRANSMISSION_R"])
        rows = np.flatnonzero((r > r_min) & (r <= r_max))
    cols = list(dict.fromkeys(["RELEASE", "BRICKID", "OBJID"] + selection.footprint_columns(cfg)
                              + selection.galaxy_columns(cfg)))
    t = io.read_table(sweep, cols, rows=rows)

    cuts = selection.footprint_cuts(t, cfg, fitbits=t["FITBITS"])
    cuts.update(selection.galaxy_cuts(t, cfg))
    flow = selection.cutflow(cuts, n_total=len(rows))
    ok = np.logical_and.reduce(list(cuts.values())) if cuts else np.ones(len(rows), bool)
    sel = np.flatnonzero(ok)

    out = {"LS_ID_DR11": io.ls_id(t["RELEASE"][sel], t["BRICKID"][sel], t["OBJID"][sel]),
           "SWEEP_ROW": rows[sel].astype(np.int32)}
    pz = io.read_pz(sweep, rows[sel], t["RELEASE"][sel], t["BRICKID"][sel], t["OBJID"][sel],
                    columns=redshift.pz_columns(cfg))
    z = redshift.best_z(pz, cfg)
    out["BEST_Z"] = z["BEST_Z"]
    out["BEST_Z_ERR"] = z["BEST_Z_ERR"]
    out["Z_SOURCE"] = z["Z_SOURCE"].astype(np.int8)
    out["STAR_FLAG"] = selection.star_flag(io.take(t, sel))

    header = {
        "SWEEP": sweep.name, "TAG": cfg["tag"], "NSWEEP": n_all, "NMAGCUT": len(rows), "NSEL": len(sel),
        "RMIN": r_min, "RMAX": r_max, "CFGHASH": config_hash(cfg), "FPHASH": config_hash(cfg, "footprint"),
        "CFGFILE": Path(cfg["_path"]).name,
        "MASKREJ": _bit_list(selection.maskbits_reject(cfg), bits.MASKBIT),
        "FITREJ": _bit_list(cfg["galaxy"].get("fitbits_reject") or [], bits.FITBIT),
        "STARFLAG": "bit0 TYPE=PSF; bit1 GAIA G-r_raw<=0.6",
        **{f"ZSRC{k}": v for k, v in redshift.z_source_names(cfg).items()},
    }
    log.info("%s: %d rows, %d with %.1f<r<=%.1f, %d selected (%.0fs)", sweep.name, n_all, len(rows),
             r_min, r_max, len(sel), time.time() - t0)
    return out, flow, header


def check_existing(path: Path, cfg: dict, overwrite: bool) -> bool:
    """True if ``path`` exists and was made with this configuration (then it is kept). A file made
    with another configuration raises unless ``overwrite``."""
    if not path.exists() or overwrite:
        return False
    h = io.read_header(path, 1)
    if h.get("CFGHASH") != config_hash(cfg):
        raise RuntimeError(f"{path} was made with another configuration (CFGHASH {h.get('CFGHASH')}); "
                           "use --overwrite to replace it")
    return True


def run_sweep(sweep: str | Path, cfg: dict, paths, overwrite: bool = False) -> Path:
    out = paths.product(sweep, cfg["tag"])
    if check_existing(out, cfg, overwrite):
        log.info("%s exists, skipped", out.name)
        return out
    t, flow, header = select_sweep(sweep, cfg)
    io.write_table(out, t, header=header, extname="SELECTION", extra=[("CUTFLOW", flow, None)])
    return out
