#!/usr/bin/env python
"""Step 2: merge the BGS-like galaxies (LS11_BGSl_DATA.fits), select the footprint randoms
(LS11_BGSl_RAND.fits) and make the footprint / systematics HEALPix maps."""

import logging

import numpy as np

from ls11samples import bgsl, cli, io, maps, randoms
from ls11samples.config import config_hash

log = logging.getLogger("ls11samples.step2")


def main():
    p = cli.parser(__doc__)
    p.add_argument("--no-maps", action="store_true")
    args = p.parse_args()
    cfg, paths, _ = cli.setup(args)
    sweeps = paths.sweeps()
    files = [paths.bgsl_dir / bgsl.output_name(s) for s in sweeps]
    missing = [f.name for f in files if not f.exists()]
    if missing:
        raise SystemExit(f"step 1 outputs missing: {missing[:5]} ... ({len(missing)})")
    fphash = config_hash(cfg, "footprint")
    for f in files:
        h = io.read_header(f, "PARENT")
        if h["FPHASH"] != fphash:
            raise SystemExit(f"{f.name} was made with another footprint configuration; rerun step 1")

    # randoms
    boxes = randoms.BoxSet([io.sweep_box(s) for s in sweeps])
    rfiles = paths.random_files()
    if not rfiles:
        raise SystemExit(f"no random files for {paths.randoms_glob}")
    if paths.rand_file.exists() and not args.overwrite:
        rand = io.read_table(paths.rand_file)
        rhdr = dict(io.read_header(paths.rand_file, 1))
        log.info("read existing %s", paths.rand_file)
    else:
        rand, rhdr = randoms.select_randoms(rfiles, cfg, boxes)
        io.write_table(paths.rand_file, rand, header=rhdr, extname="RANDOMS")
    area = float(rhdr["AREA"])
    log.info("randoms: %d in footprint, area %.2f deg2 (sweep boxes %.2f deg2)", len(rand["RA"]), area, boxes.area)

    # galaxies
    if not paths.data_file.exists() or args.overwrite:
        data = bgsl.read_selected(files)
        hdr = {"NSWEEPS": len(files), "NGAL": len(data["RA"]), "AREA": area,
               "NDENS": len(data["RA"]) / area, "FPHASH": fphash, "CFGHASH": config_hash(cfg),
               "RANDFILE": paths.rand_file.name}
        hdr.update({k: v for k, v in dict(io.read_header(files[0], "PARENT")).items()
                    if k.startswith(("SELBIT", "ZSRC", "MASKREJ", "FITREJ", "RMAX"))})
        io.write_table(paths.data_file, data, header=hdr, extname="BGSL")
        log.info("galaxies: %d, %.1f per deg2", len(data["RA"]), len(data["RA"]) / area)

    if not args.no_maps:
        written = maps.make_maps(rand, cfg, paths.out, float(rhdr["DENSITY"]), paths.gaia_maps)
        log.info("%d maps written", len(written))


if __name__ == "__main__":
    main()
