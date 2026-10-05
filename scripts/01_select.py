#!/usr/bin/env python
"""Step 1: BGS-like parent table + cut-flow for every sweep ($LS11_OUT/bgsl/BGSl-<sweep>.fits)."""

from concurrent.futures import ProcessPoolExecutor
from functools import partial

from ls11samples import bgsl, cli


def main():
    p = cli.parser(__doc__)
    cli.add_slicing(p)
    args = p.parse_args()
    cfg, paths, n = cli.setup(args)
    sweeps = cli.my_part(paths.sweeps(), args.part, args.nparts)
    if not sweeps:
        raise SystemExit(f"no sweeps in {paths.sweep_dir} matching {paths.sweep_glob}")
    run = partial(bgsl.run_sweep, cfg=cfg, outdir=paths.bgsl_dir, overwrite=args.overwrite)
    if n == 1:
        list(map(run, sweeps))
    else:
        with ProcessPoolExecutor(n) as ex:
            list(ex.map(run, sweeps))


if __name__ == "__main__":
    main()
