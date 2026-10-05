#!/usr/bin/env python
"""Step 3: kcorrect v5 K-corrections and absolute magnitudes ($LS11_OUT/kcorr/KC-<sweep>.fits)."""

from concurrent.futures import ProcessPoolExecutor
from functools import partial

from ls11samples import bgsl, cli, kcorr


def main():
    p = cli.parser(__doc__)
    cli.add_slicing(p)
    args = p.parse_args()
    cfg, paths, n = cli.setup(args)
    files = [paths.bgsl_dir / bgsl.output_name(s) for s in cli.my_part(paths.sweeps(), args.part, args.nparts)]
    files = [f for f in files if f.exists()]
    run = partial(kcorr.run_sweep, cfg=cfg, outdir=paths.outdir("kcorr"), overwrite=args.overwrite)
    if n == 1:
        list(map(run, files))
    else:
        with ProcessPoolExecutor(n) as ex:
            list(ex.map(run, files))


if __name__ == "__main__":
    main()
