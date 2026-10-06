#!/usr/bin/env python
"""Validate a run: every sweep has its selection and code files, row-aligned; missing-value
fractions; mass / Mr distributions and errors per code; agreement with the reference code;
randoms; disk use. Writes $LS11_OUT/<tag>/validation.yaml and exits with 1 if files are missing
or misaligned.

    LS11_CONFIG=config/bgs_r21_dr10bits.yaml python scripts/validate_run.py
    (LS11_SWEEP_LIST=<frozen list> restricts the sweeps, as in the jobs)
"""

import sys

import yaml

from ls11samples import cli
from ls11samples.validate import report, validate


def main():
    p = cli.parser(__doc__)
    p.add_argument("--codes", default=None, help="comma-separated codes (default: sed.codes)")
    p.add_argument("--ref", default=None, help="reference code (default: sed.primary)")
    p.add_argument("--frac", type=float, default=0.02, help="galaxy subsample for the distributions")
    args = p.parse_args()
    cfg, paths, _ = cli.setup(args)
    out = validate(paths, cfg, paths.sweeps(), codes=args.codes.split(",") if args.codes else None,
                   frac=args.frac, ref=args.ref)
    print(report(out))
    path = paths.run_dir(cfg["tag"]) / "validation.yaml"
    path.write_text(yaml.safe_dump(out, sort_keys=False))
    print(f"written {path}")
    sys.exit(1 if out["missing"] or out["misaligned"] else 0)


if __name__ == "__main__":
    main()
