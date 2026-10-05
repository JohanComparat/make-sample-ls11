#!/usr/bin/env python
"""Step 4: stellar masses at fixed z = BEST_Z for every code of ``sed.codes``
($LS11_OUT/mstar/<code>/MS-<sweep>.fits).

--prepare only builds what the codes need once (LePhare libraries, CIGALE filters): run it before
the job arrays so that parallel tasks never build the same files at the same time.
"""

import logging
import time
from functools import partial
from pathlib import Path

import numpy as np

from ls11samples import bgsl, cli, io
from ls11samples.sed import get_backend
from ls11samples.sed.common import BANDS

log = logging.getLogger("ls11samples.step4")
INPUT = ["LS_ID_DR11", "BEST_Z"] + [f"{p}_{b}" for p in ("FLUX", "FLUX_IVAR", "MW_TRANSMISSION") for b in BANDS]


def run_sweep(parent: Path, backend, code: str, outdir: Path, overwrite: bool, chunk: int) -> Path:
    out = outdir / ("MS-" + parent.name.removeprefix("BGSl-"))
    if out.exists() and not overwrite:
        return out
    t0 = time.time()
    t = bgsl.read_selected([parent], INPUT) or {c: np.zeros(0) for c in INPUT}
    n = len(t["BEST_Z"])
    parts = [backend.fit(io.take(t, slice(i, i + chunk))) for i in range(0, n, chunk)] if n else [backend.fit(t)]
    io.write_table(out, {"LS_ID_DR11": t["LS_ID_DR11"], **io.concat(parts)}, header={"SEDCODE": code},
                   extname="MSTAR")
    log.info("%s %s: %d galaxies (%.0fs)", code, out.name, n, time.time() - t0)
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
    files = [paths.bgsl_dir / bgsl.output_name(s) for s in cli.my_part(paths.sweeps(), args.part, args.nparts)]
    files = [f for f in files if f.exists()]
    for code in codes:
        backend = get_backend(code)(cfg)          # builds / checks libraries once per process
        if args.prepare:
            log.info("%s ready", code)
            continue
        run = partial(run_sweep, backend=backend, code=code, outdir=paths.outdir("mstar", code),
                      overwrite=args.overwrite, chunk=args.chunk)
        # the codes parallelise internally (OpenMP, pcigale cores, jax): one sweep at a time
        list(map(run, files))


if __name__ == "__main__":
    main()
