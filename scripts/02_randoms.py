#!/usr/bin/env python
"""Step 2: footprint randoms ($LS11_OUT/<tag>/LS11_<tag>_RAND.fits: RA, DEC, EBV) and the
footprint / systematics HEALPix maps, accumulated while the random files are read. An existing
random file is kept only if it was made from the same footprint, sweeps and random files."""

import logging

from ls11samples import cli, io, randoms
from ls11samples.config import config_hash
from ls11samples.maps import MapAccumulator

log = logging.getLogger("ls11samples.step2")


def main():
    p = cli.parser(__doc__)
    p.add_argument("--no-maps", action="store_true")
    args = p.parse_args()
    cfg, paths, _ = cli.setup(args)
    tag = cfg["tag"]
    sweeps = paths.sweeps()
    missing = [s.name for s in sweeps if not paths.product(s, tag).exists()]
    if missing:
        raise SystemExit(f"step 1 outputs missing for {len(missing)} sweeps, e.g. {missing[:3]}")
    fphash = config_hash(cfg, "footprint")
    for s in sweeps:
        if io.read_header(paths.product(s, tag), 1)["FPHASH"] != fphash:
            raise SystemExit(f"{paths.product(s, tag).name} was made with another footprint; rerun step 1")
    out = paths.rand_file(tag)
    rfiles = paths.random_files()
    if not rfiles:
        raise SystemExit(f"no random files for {paths.randoms_glob}")
    if out.exists() and not args.overwrite:
        stale = randoms.stale_reason(io.read_header(out, 1), fphash, sweeps, rfiles)
        if not stale:
            log.info("%s is up to date (%d sweeps); use --overwrite to redo it", out, len(sweeps))
            return
        log.info("%s is out of date (%s): redoing it", out, stale)
    boxes = randoms.BoxSet([io.sweep_box(s) for s in sweeps])
    acc = None if args.no_maps else MapAccumulator(cfg)
    rand, hdr = randoms.select_randoms(rfiles, cfg, boxes, maps=acc)
    hdr.update(NSWEEPS=len(sweeps), SWPHASH=randoms.sweeps_hash(sweeps))
    io.write_table(out, rand, header=hdr, extname="RANDOMS")
    log.info("randoms: %d in footprint, area %.2f deg2 (sweep boxes %.2f deg2)", len(rand["RA"]), hdr["AREA"], boxes.area)
    if acc is not None:
        written = acc.write(paths.run_dir(tag), float(hdr["DENSITY"]), paths.gaia_maps)
        log.info("%d maps written", len(written))


if __name__ == "__main__":
    main()
