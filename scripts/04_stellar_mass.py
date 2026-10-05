#!/usr/bin/env python
"""Step 4: stellar masses at fixed z = BEST_Z ($LS11_OUT/mstar/<code>/MS-<sweep>.fits)."""

import logging
import time
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

import numpy as np

from ls11samples import bgsl, cli, io
from ls11samples.sed import get_backend
from ls11samples.sed.common import BANDS

log = logging.getLogger("ls11samples.step4")
INPUT = ["LS_ID_DR11", "BEST_Z"] + [f"{p}_{b}" for p in ("FLUX", "FLUX_IVAR", "MW_TRANSMISSION") for b in BANDS]


def run_sweep(parent: Path, cfg: dict, code: str, outdir: Path, overwrite: bool, chunk: int) -> Path:
    out = outdir / ("MS-" + parent.name.removeprefix("BGSl-"))
    if out.exists() and not overwrite:
        return out
    t0 = time.time()
    t = bgsl.read_selected([parent], INPUT) or {c: np.zeros(0) for c in INPUT}
    backend = get_backend(code)(cfg)
    n = len(t["BEST_Z"])
    parts = [backend.fit(io.take(t, slice(i, i + chunk))) for i in range(0, max(n, 1), chunk)]
    res = io.concat(parts) if n else backend.fit(t)
    io.write_table(out, {"LS_ID_DR11": t["LS_ID_DR11"], **res}, header={"SEDCODE": code}, extname="MSTAR")
    log.info("%s %s: %d galaxies (%.0fs)", code, out.name, n, time.time() - t0)
    return out


def main():
    p = cli.parser(__doc__)
    cli.add_slicing(p)
    p.add_argument("--code", default=None, help="SED backend (default: sed.code in the config)")
    p.add_argument("--chunk", type=int, default=20000)
    args = p.parse_args()
    cfg, paths, n = cli.setup(args)
    code = args.code or cfg["sed"]["code"]
    files = [paths.bgsl_dir / bgsl.output_name(s) for s in cli.my_part(paths.sweeps(), args.part, args.nparts)]
    files = [f for f in files if f.exists()]
    run = partial(run_sweep, cfg=cfg, code=code, outdir=paths.outdir("mstar", code),
                  overwrite=args.overwrite, chunk=args.chunk)
    # codes that thread internally (LePhare: OpenMP, DSPS: jax) run one sweep at a time
    if n == 1 or code in ("lephare", "dsps", "cigale"):
        list(map(run, files))
    else:
        with ProcessPoolExecutor(n) as ex:
            list(ex.map(run, files))


if __name__ == "__main__":
    main()
