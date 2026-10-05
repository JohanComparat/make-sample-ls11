"""CIGALE (pcigale >= 2025.0) at fixed redshift, run through the ``pcigale`` command line.

pcigale forces OMP_NUM_THREADS=1 when imported, so it is only ever run in subprocesses here (also
to register the shared DECam/WISE filters, named ``ls11.<kcorrect name>``, in its database).

Model grid (``sed.cigale`` in the configuration, defaults below): delayed-tau SFH, BC03 Chabrier,
nebular emission, Calzetti-like modified starburst attenuation; redshifts rounded to
``redshift_decimals`` so the models are computed once per redshift step. The rest-frame DECam r
luminosity (``restframe_parameters``, L_nu at 10 pc) gives MABS_R and its Bayesian error.
pcigale's ``additionalerror`` (10% of the flux added to every band by default) is set to 0: the
input errors already carry the common error floor. pcigale floors the errors of extensive
properties at 5% of the value (``managers/results.py``, not configurable): LOGMSTAR_ERR >= 0.0217
dex and MABS_R_ERR >= 0.054 mag.

Each fit runs in a temporary directory under $TMPDIR (node-local on the cluster), removed after
the results are read, so nothing accumulates on disk.
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Mapping
from pathlib import Path

import numpy as np

from ..env import nproc
from .common import BANDS, NANOMAGGY_MJY, empty_result, photometry, valid_z
from .filters import KCORRECT_NAMES, read_response

log = logging.getLogger(__name__)
FILTER_NAMES = {b: f"ls11.{KCORRECT_NAMES[b]}" for b in BANDS}
PCIGALE = str(Path(sys.executable).parent / "pcigale")     # the env's script, PATH not required

#: pcigale modules and parameters; ``sed.cigale.grid`` overrides them per module.
DEFAULT_GRID = {
    "sfhdelayed": {"tau_main": "250, 500, 1000, 2000, 4000, 8000",
                   "age_main": "500, 1000, 2000, 3000, 5000, 7000, 9000, 11000",
                   "tau_burst": "50", "age_burst": "20", "f_burst": "0.0", "sfr_A": "1.0",
                   "normalise": "True"},
    "bc03": {"imf": "1", "metallicity": "0.008, 0.02", "separation_age": "10"},
    "nebular": {"logU": "-2.0", "emission": "True"},
    "dustatt_modified_starburst": {"E_BV_lines": "0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8",
                                   "E_BV_factor": "0.44", "uv_bump_amplitude": "0.0",
                                   "powerlaw_slope": "0.0"},
    "restframe_parameters": {"beta_calz94": "False", "Dn4000": "False", "IRX": "False", "EW": "",
                             "luminosity_filters": "ls11.decam_r", "colours_filters": ""},
    "redshifting": {"redshift": ""},
}
LNU_R = "param.restframe_Lnu(ls11.decam_r)"
TO_LUMIN = 1e-29 * 4.0 * np.pi * (3.0856775814913673e17) ** 2     # mJy at 10 pc -> W/Hz (as pcigale)


def _run(args, cwd=None):
    env = {k: v for k, v in os.environ.items()}
    r = subprocess.run(args, cwd=cwd, env=env, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"{' '.join(map(str, args))} failed:\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}")
    return r.stdout


def ensure_filters() -> None:
    """Register the ls11.* filters in the pcigale database (once)."""
    code = ("from pcigale.data import SimpleDatabase as D\n"
            "with D('filters') as db: print(' '.join(db.parameters['name']))")
    names = set(_run([sys.executable, "-c", code]).split())
    todo = [b for b in BANDS if FILTER_NAMES[b] not in names]
    if not todo:
        return
    with tempfile.TemporaryDirectory() as tmp:
        files = []
        for b in todo:
            wl, tr = read_response(b)
            f = Path(tmp) / f"{FILTER_NAMES[b]}.dat"
            with open(f, "w") as fh:
                fh.write(f"# {FILTER_NAMES[b]}\n# photon\n# {KCORRECT_NAMES[b]} (kcorrect v5 curve, LS DR11)\n")
                np.savetxt(fh, np.c_[wl, tr], fmt="%.3f %.6e")
            files.append(str(f))
        _run([sys.executable, "-c", f"from pcigale_filters import add_filters; add_filters({files!r})"])
    log.info("registered %s in the pcigale database", [FILTER_NAMES[b] for b in todo])


class CigaleBackend:
    name = "cigale"

    def __init__(self, cfg: dict, workdir: str | Path | None = None, cores: int | None = None,
                 keep_runs: bool = False):
        self.cfg = cfg
        c = (cfg.get("sed") or {}).get("cigale") or {}
        self.grid = {m: {**p, **(c.get("grid", {}).get(m, {}))} for m, p in DEFAULT_GRID.items()}
        self.decimals = int(c.get("redshift_decimals", 2))
        self.cores = cores or nproc()
        self.keep_runs = keep_runs
        if workdir is None:
            workdir = Path(os.environ.get("TMPDIR", tempfile.gettempdir())) / "ls11_cigale"
        self.workroot = Path(workdir)
        self.workroot.mkdir(parents=True, exist_ok=True)
        ensure_filters()

    def _write_input(self, path: Path, z, f, e) -> None:
        cols = ["id", "redshift"]
        for b in BANDS:
            cols += [FILTER_NAMES[b], FILTER_NAMES[b] + "_err"]
        good = np.isfinite(e) & np.isfinite(f)
        data = [np.arange(len(z)), z]
        for i in range(len(BANDS)):
            data += [np.where(good[:, i], f[:, i] * NANOMAGGY_MJY, np.nan),
                     np.where(good[:, i], e[:, i] * NANOMAGGY_MJY, np.nan)]
        np.savetxt(path, np.column_stack(data), header=" ".join(cols), comments="# ",
                   fmt=["%d", "%.5f"] + ["%.6e"] * (2 * len(BANDS)))

    def _write_config(self, run: Path) -> None:
        from configobj import ConfigObj

        _run([PCIGALE, "init"], cwd=run)
        conf = ConfigObj(str(run / "pcigale.ini"), encoding="UTF8")
        conf["data_file"] = "input.txt"
        conf["sed_modules"] = list(self.grid)
        conf["analysis_method"] = "pdf_analysis"
        conf["cores"] = self.cores
        conf.write()
        _run([PCIGALE, "genconf"], cwd=run)
        conf = ConfigObj(str(run / "pcigale.ini"), encoding="UTF8")
        for mod, pars in self.grid.items():
            for k, v in pars.items():
                vals = [x.strip() for x in str(v).split(",")]
                # configobj writes a list unquoted (a grid); a single value as a scalar
                conf["sed_modules_params"][mod][k] = vals if len(vals) > 1 else vals[0]
        # set after genconf (which writes the default 0.1): pcigale would add 10% of the flux in
        # quadrature to every band and to the errors of extensive properties; the input errors
        # already carry the common error floor
        conf["additionalerror"] = 0.0
        ap = conf["analysis_params"]
        ap["variables"] = ["stellar.m_star", "sfh.sfr", LNU_R]
        ap["bands"] = [FILTER_NAMES[b] for b in BANDS]
        ap["save_best_sed"] = False
        ap["redshift_decimals"] = self.decimals
        ap["blocks"] = 1
        conf.write()

    def fit(self, t: Mapping) -> dict[str, np.ndarray]:
        from astropy.table import Table

        n = len(t["BEST_Z"])
        out = empty_result(n)
        z, f, e = photometry(t, self.cfg)
        ok = valid_z(z)
        if not ok.any():
            return out
        run = Path(tempfile.mkdtemp(prefix="run_", dir=self.workroot))
        self._write_input(run / "input.txt", z[ok], f[ok], e[ok])
        self._write_config(run)
        _run([PCIGALE, "run"], cwd=run)
        res = Table.read(run / "out" / "results.fits")
        idx = np.flatnonzero(ok)[np.asarray(res["id"]).astype(int)]
        m = np.asarray(res["bayes.stellar.m_star"], np.float64)
        me = np.asarray(res["bayes.stellar.m_star_err"], np.float64)
        lnu = np.asarray(res[f"bayes.{LNU_R}"], np.float64)
        lnu_err = np.asarray(res[f"bayes.{LNU_R}_err"], np.float64)
        with np.errstate(divide="ignore", invalid="ignore"):
            out["LOGMSTAR"][idx] = np.log10(m)
            out["LOGMSTAR_ERR"][idx] = me / (m * np.log(10))
            out["LOGMSTAR_LO"][idx] = np.log10(np.clip(m - me, 1e-3, None))
            out["LOGMSTAR_HI"][idx] = np.log10(m + me)
            out["LOGSFR"][idx] = np.log10(np.asarray(res["bayes.sfh.sfr"], np.float64))
            # AB absolute magnitude: flux density at 10 pc in mJy vs 3631 Jy
            out["MABS_R"][idx] = -2.5 * np.log10(lnu / TO_LUMIN / 3.631e6)
            out["MABS_R_ERR"][idx] = 2.5 / np.log(10) * lnu_err / lnu
        out["CHI2"][idx] = np.asarray(res["best.chi_square"], np.float64)
        if self.keep_runs:
            self.last_run = run
        else:
            shutil.rmtree(run, ignore_errors=True)
        return out
