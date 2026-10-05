#!/usr/bin/env python
"""Step 4: stellar masses and rest-frame r absolute magnitudes, with uncertainties, at fixed
z = BEST_Z, for every code of ``sed.codes``:

    <sweep root>/<ver>-<code>/<sweep>-<code>.fits
    SWEEP_ROW (i4), LOGMSTAR, LOGMSTAR_ERR, MABS_R, MABS_R_ERR (f4)

row-aligned with the selection file <ver>-<tag>/<sweep>-<tag>.fits; the photometry is read from the
sweep at SWEEP_ROW. --prepare only builds what the codes need once (LePhare libraries, CIGALE
filters): run it before the job arrays so that parallel tasks never build the same files.
"""

import logging
import time
import zlib
from pathlib import Path

import numpy as np

from ls11samples import catalog, cli, io
from ls11samples.config import config_hash
from ls11samples.cosmo import cosmology
from ls11samples.sed import PRODUCT, get_backend
from ls11samples.sed.common import BANDS

log = logging.getLogger("ls11samples.step4")
PHOT = [f"{p}_{b}" for p in ("FLUX", "FLUX_IVAR", "MW_TRANSMISSION") for b in BANDS]


def up_to_date(out: Path, selhash: str, overwrite: bool) -> bool:
    if not out.exists() or overwrite:
        return False
    if io.read_header(out, 1).get("SELHASH") != selhash:
        raise RuntimeError(f"{out} belongs to another selection; use --overwrite to replace it")
    return True


def run_sweep(sweep: Path, backend, code: str, cfg: dict, paths, overwrite: bool, chunk: int) -> Path:
    tag = cfg["tag"]
    sel_file = paths.product(sweep, tag)
    selhash = io.read_header(sel_file, 1)["CFGHASH"]
    out = paths.product(sweep, code)
    if up_to_date(out, selhash, overwrite):
        return out
    t0 = time.time()
    sel = io.read_table(sel_file, ["SWEEP_ROW", "BEST_Z"])
    t = {**sel, **catalog.sweep_rows(sweep, sel["SWEEP_ROW"], PHOT)}
    if hasattr(backend, "seed"):
        backend.seed = zlib.crc32(sweep.name.encode())             # reproducible Monte Carlo errors
    n = len(t["BEST_Z"])
    parts = [backend.fit(io.take(t, slice(i, i + chunk))) for i in range(0, n, chunk)] if n else [backend.fit(t)]
    res = io.concat(parts) or {k: np.zeros(0) for k in PRODUCT}
    cos = cosmology(cfg)
    header = {"CODE": code, "TAG": tag, "SELHASH": selhash, "CFGHASH": config_hash(cfg),
              "MABSBAND": "rest-frame DECam r (z0=0), AB", "MASSIMF": "Chabrier",
              "ERRDEF": "1 sigma at fixed z (no photo-z term)",
              "COSMO": f"FlatLCDM H0={cos.H0.value} Om0={cos.Om0}"}
    io.write_table(out, {"SWEEP_ROW": sel["SWEEP_ROW"], **{k: np.asarray(res[k], np.float32) for k in PRODUCT}},
                   header=header, extname=code.upper())
    log.info("%s %s: %d objects (%.0fs)", code, out.name, n, time.time() - t0)
    return out


def main():
    p = cli.parser(__doc__)
    cli.add_slicing(p)
    p.add_argument("--code", default=None, help="comma-separated codes (default: sed.codes in the config)")
    p.add_argument("--chunk", type=int, default=20000)
    p.add_argument("--prepare", action="store_true", help="only build the code libraries / filters")
    args = p.parse_args()
    cfg, paths, _ = cli.setup(args)
    codes = args.code.split(",") if args.code else list(cfg["sed"]["codes"])
    sweeps = [s for s in cli.my_part(paths.sweeps(), args.part, args.nparts)
              if paths.product(s, cfg["tag"]).exists()]
    for code in codes:
        backend = get_backend(code)(cfg)          # builds / checks libraries once per process
        if args.prepare:
            log.info("%s ready", code)
            continue
        paths.product_dir(code).mkdir(parents=True, exist_ok=True)
        for s in sweeps:                          # the codes parallelise internally
            run_sweep(s, backend, code, cfg, paths, args.overwrite, args.chunk)


if __name__ == "__main__":
    main()
