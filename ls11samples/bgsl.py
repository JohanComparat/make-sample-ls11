"""Step 1: one sweep -> BGS-like parent table (all objects with dereddened r <= r_max) + cut-flow.

The parent keeps every object brighter than the r limit with its SEL_FLAGS, so the cuts can be
changed later without re-reading the sweeps; the BGS-like sample is SEL_FLAGS == 0.
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

BANDS_ALL = ("G", "R", "I", "Z", "W1", "W2")
KEEP_SWEEP = (
    ["RELEASE", "BRICKID", "BRICKNAME", "OBJID", "TYPE", "RA", "DEC", "EBV"]
    + [f"{p}_{b}" for p in ("FLUX", "FLUX_IVAR", "MW_TRANSMISSION") for b in BANDS_ALL]
    + [f"NOBS_{b}" for b in "GRIZ"]
    + [f"{p}_{b}" for p in ("GALDEPTH", "PSFDEPTH", "PSFSIZE") for b in "GRIZ"]
    + ["PSFDEPTH_W1", "PSFDEPTH_W2"]
    + [f"{p}_{b}" for p in ("FRACMASKED", "FRACIN", "FRACFLUX") for b in "GRZ"]
    + ["FIBERFLUX_R", "FIBERTOTFLUX_R", "SHAPE_R", "SERSIC", "REF_CAT", "GAIA_PHOT_G_MEAN_MAG",
       "MASKBITS", "FITBITS"]
)
PZ_DROP = ("LS_ID_DR11", "RELEASE", "BRICKID", "OBJID")


def output_name(sweep: str | Path) -> str:
    return "BGSl-" + Path(sweep).name


def select_sweep(sweep: str | Path, cfg: dict) -> tuple[dict, dict, dict]:
    """(parent table, cut-flow table, header) of one sweep."""
    t0 = time.time()
    sweep = Path(sweep)
    bits.check_header(io.read_header(sweep, 0))
    n_all = io.nrows(sweep)
    r_max = cfg["galaxy"]["r_range"][1]

    # parent: dereddened r <= r_max (then read the other columns for those rows only)
    pre = io.read_table(sweep, ["FLUX_R", "MW_TRANSMISSION_R"])
    with np.errstate(invalid="ignore"):
        rows = np.flatnonzero(mag(pre["FLUX_R"], pre["MW_TRANSMISSION_R"]) <= r_max)
    cols = list(dict.fromkeys(KEEP_SWEEP + selection.footprint_columns(cfg)
                              + selection.galaxy_columns(cfg)))
    t = io.read_table(sweep, cols, rows=rows)

    cuts = selection.footprint_cuts(t, cfg, fitbits=t["FITBITS"])
    cuts.update(selection.galaxy_cuts(t, cfg))
    t["SEL_FLAGS"] = selection.sel_flags(cuts)
    flow = selection.cutflow(cuts, n_total=len(rows))

    t["LS_ID_DR11"] = io.ls_id(t["RELEASE"], t["BRICKID"], t["OBJID"])
    for b in BANDS_ALL:
        t[f"MAG_{b}"] = mag(t[f"FLUX_{b}"], t[f"MW_TRANSMISSION_{b}"]).astype(np.float32)
    t["RFIB"] = mag(t["FIBERFLUX_R"], t["MW_TRANSMISSION_R"]).astype(np.float32)
    t["RFIBTOT"] = mag(t["FIBERTOTFLUX_R"], t["MW_TRANSMISSION_R"]).astype(np.float32)

    pz = io.read_pz(sweep, rows, t["RELEASE"], t["BRICKID"], t["OBJID"])
    for c, v in pz.items():
        if c not in PZ_DROP:
            t[c] = v
    t.update(redshift.best_z(pz, cfg))

    n_sel = int((t["SEL_FLAGS"] == 0).sum())
    header = {
        "SWEEP": sweep.name, "NSWEEP": n_all, "NPARENT": len(rows), "NSEL": n_sel,
        "RMAX": r_max, "CFGHASH": config_hash(cfg), "FPHASH": config_hash(cfg, "footprint"),
        "CFGFILE": Path(cfg["_path"]).name,
        "MASKREJ": ",".join(selection.maskbits_reject(cfg)),
        "FITREJ": ",".join(cfg["galaxy"]["fitbits_reject"]),
        "PHOTSYS": "S", **selection.sel_header(),
        **{f"ZSRC{k}": v for k, v in redshift.z_source_names(cfg).items()},
    }
    log.info("%s: %d rows, %d parent (r<=%.1f), %d selected (%.0fs)", sweep.name, n_all, len(rows),
             r_max, n_sel, time.time() - t0)
    return t, flow, header


def run_sweep(sweep: str | Path, cfg: dict, outdir: str | Path, overwrite: bool = False) -> Path:
    out = Path(outdir) / output_name(sweep)
    if out.exists() and not overwrite:
        log.info("%s exists, skipped", out.name)
        return out
    t, flow, header = select_sweep(sweep, cfg)
    io.write_table(out, t, header=header, extname="PARENT", extra=[("CUTFLOW", flow, None)])
    return out


def read_selected(files, columns=None) -> dict:
    """Concatenate the BGS-like objects (SEL_FLAGS == 0) of per-sweep parent files."""
    parts = []
    for f in files:
        flags = io.read_table(f, ["SEL_FLAGS"])["SEL_FLAGS"]
        rows = np.flatnonzero(flags == 0)
        parts.append(io.read_table(f, columns, rows=rows))
    return io.concat(parts)
