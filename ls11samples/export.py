"""Volume-limited sample files in the formats read by sys_mapping and sum_stat::

    <name>_DATA.fits    RA, DEC (f8), EBV, BEST_Z, BEST_Z_ERR, Z_SOURCE, STAR_FLAG, LS_ID_DR11,
                        MAG_G/R/Z, MABS_R and KCORR_R (``vlim.mr_code``), LOGMSTAR (``sed.primary``),
                        LOGMSTAR[_ERR]_<CODE> and MABS_R[_ERR]_<CODE> for every code run,
                        LPH_MASS_BEST (= LOGMSTAR, the column sum_stat reads), WEIGHT_COMP (= 1;
                        sys_mapping uses it when present)
    <name>_RAND.fits    RA, DEC (f8), EBV, Z (shuffled data redshifts)
    <name>_COLOUR.fits  G_MAG, Z_MAG (dereddened DECam AB), REDSHIFT (= BEST_Z), row-aligned with
                        DATA (sum_stat colour classes); header PARENT, NMATCH, NMISS

The name follows the DR10 convention with an LS11 prefix:
``LS11_VLIM_ANY_<lo>_<Mstar|Mr>_<hi>_<zmin>_z_<zmax>_N_<N:07d>``.
"""

from __future__ import annotations

import datetime as _dt
import subprocess
from collections.abc import Mapping
from pathlib import Path

import numpy as np
import yaml

from . import __version__, io

DATA_COLUMNS = ("RA", "DEC", "EBV", "BEST_Z", "BEST_Z_ERR", "Z_SOURCE", "STAR_FLAG", "LS_ID_DR11", "MAG_G",
                "MAG_R", "MAG_Z", "MABS_R", "KCORR_R", "LOGMSTAR")


PER_CODE = ("LOGMSTAR_", "MABS_R_")        # LOGMSTAR[_ERR]_<CODE>, MABS_R[_ERR]_<CODE>: every code run


def data_table(data: Mapping, sel: np.ndarray) -> dict:
    cols = [c for c in DATA_COLUMNS if c in data]
    cols += [c for c in data if c.startswith(PER_CODE) and c not in cols]
    out = {c: np.asarray(data[c])[sel] for c in cols}
    out["RA"] = out["RA"].astype(np.float64)
    out["DEC"] = out["DEC"].astype(np.float64)
    if "LOGMSTAR" in out:
        out["LPH_MASS_BEST"] = out["LOGMSTAR"].astype(np.float64)
    out["WEIGHT_COMP"] = np.ones(int(np.count_nonzero(sel)), np.float32)
    return out


def colour_table(data: Mapping, sel: np.ndarray) -> dict:
    return {"G_MAG": np.asarray(data["MAG_G"])[sel].astype(np.float32),
            "Z_MAG": np.asarray(data["MAG_Z"])[sel].astype(np.float32),
            "REDSHIFT": np.asarray(data["BEST_Z"])[sel].astype(np.float32)}


def write_sample(outdir: Path, sample: dict, data: Mapping, rand: Mapping, header: Mapping) -> list[Path]:
    name = sample["name"]
    sel = sample["sel"]
    hdr = {"SAMPLE": name, "KIND": sample["kind"], "LO": sample["lo"], "HI": sample["hi"],
           "ZMIN": sample["z0"], "ZMAX": sample["z1"], "NGAL": sample["n"], **header}
    d = io.write_table(outdir / f"{name}_DATA.fits", data_table(data, sel), header=hdr, extname="DATA")
    r = io.write_table(outdir / f"{name}_RAND.fits", rand, header={**hdr, "NRAND": len(rand["RA"])},
                       extname="RAND")
    c = io.write_table(outdir / f"{name}_COLOUR.fits", colour_table(data, sel),
                       header={"PARENT": header.get("PARENT", ""),
                               "NMATCH": sample["n"], "NMISS": 0}, extname="COLOUR")
    return [d, r, c]


def git_commit() -> str:
    try:
        return subprocess.run(["git", "-C", str(Path(__file__).resolve().parents[1]), "rev-parse",
                               "--short", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
    except Exception:  # noqa: BLE001 - not a git checkout
        return "unknown"


def write_manifest(path: Path, cfg: dict, summary: list[dict], extra: Mapping | None = None) -> Path:
    from .config import config_hash

    doc = {"created": _dt.datetime.now().isoformat(timespec="seconds"),
           "package_version": __version__, "git_commit": git_commit(),
           "config_file": cfg.get("_path"), "config_hash": config_hash(cfg),
           "footprint_hash": config_hash(cfg, "footprint"), **(extra or {}),
           "samples": [{k: (float(v) if isinstance(v, (np.floating, float)) else
                            int(v) if isinstance(v, (np.integer, int)) else v)
                        for k, v in row.items()} for row in summary]}
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        yaml.safe_dump(doc, f, sort_keys=False)
    return path
