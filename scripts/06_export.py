#!/usr/bin/env python
"""Step 6: check the volume-limited samples against the sys_mapping / sum_stat input formats and
write $LS11_OUT/<tag>/vlim/manifest.yaml (with the systematics template directories).

sys_mapping:  run_ls10_analysis.py --catalog-dir $LS11_OUT/<tag>/vlim --sample <name>
                                    --template-dir $LS11_OUT/<tag>/systematics/<nside:04d> --nside <nside>
sum_stat:     measure_joint_sumstat.py --survey custom --data-file <name>_DATA.fits --rand-file <name>_RAND.fits
"""

import logging
import sys

import numpy as np
import yaml

from ls11samples import cli, io

log = logging.getLogger("ls11samples.step6")
REQUIRED = {"DATA": {"RA": "f8", "DEC": "f8", "BEST_Z": "f", "LPH_MASS_BEST": "f8", "WEIGHT_COMP": "f"},
            "RAND": {"RA": "f8", "DEC": "f8", "Z": "f"},
            "COLOUR": {"G_MAG": "f", "Z_MAG": "f", "REDSHIFT": "f"}}


def check_sample(vdir, name) -> list[str]:
    errors = []
    lens = {}
    for kind, cols in REQUIRED.items():
        path = vdir / f"{name}_{kind}.fits"
        if not path.exists():
            errors.append(f"{path.name} missing")
            continue
        t = io.read_table(path, list(cols))
        lens[kind] = len(t["RA" if "RA" in t else "REDSHIFT"])
        for c, kindcode in cols.items():
            a = t[c]
            if kindcode == "f8" and a.dtype != np.float64:
                errors.append(f"{path.name}:{c} is {a.dtype}, not float64")
            if c in ("RA", "DEC", "BEST_Z", "Z", "REDSHIFT") and not np.all(np.isfinite(a)):
                errors.append(f"{path.name}:{c} has non-finite values")
    if "DATA" in lens and "COLOUR" in lens and lens["DATA"] != lens["COLOUR"]:
        errors.append(f"{name}: COLOUR not row-aligned with DATA")
    return errors


def main():
    p = cli.parser(__doc__)
    args = p.parse_args()
    cfg, paths, _ = cli.setup(args)
    tag = cfg["tag"]
    vdir = paths.out / tag / "vlim"
    mpath = vdir / "manifest.yaml"
    manifest = yaml.safe_load(mpath.read_text())
    errors = []
    for s in manifest["samples"]:
        errors += check_sample(vdir, s["NAME"])
    manifest["sys_mapping"] = {"catalog_dir": str(vdir),
                               "template_dirs": {int(n): str(paths.out / tag / "systematics" / f"{n:04d}")
                                                 for n in cfg["maps"]["nsides"]},
                               "fracarea": {int(n): str(paths.out / tag / "footprint" / f"LS11_FRACAREA_NSIDE_{n:04d}.fits")
                                            for n in cfg["maps"]["nsides"]}}
    manifest["sum_stat"] = {"survey": "custom", "data_columns": ["RA", "DEC", "BEST_Z", "LPH_MASS_BEST", "MABS_R"],
                            "rand_columns": ["RA", "DEC", "Z"], "colour_file": "<name>_COLOUR.fits"}
    manifest["checks"] = {"n_samples": len(manifest["samples"]), "errors": errors}
    mpath.write_text(yaml.safe_dump(manifest, sort_keys=False))
    for e in errors:
        log.error(e)
    log.info("%d samples checked, %d problems; manifest %s", len(manifest["samples"]), len(errors), mpath)
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
