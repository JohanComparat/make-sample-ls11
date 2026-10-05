"""eazy-py at fixed redshift (``fit_at_zbest``) with the FSPS-based ``corr_sfhz_13`` templates;
stellar masses and SFRs from ``sps_parameters`` (template mass-to-light ratios), rest-frame DECam r
from the template fit.

Needs the eazy-photoz templates/filters (cloned to $EAZY_DATA, default ~/.cache/eazy-photoz). The
shared DECam/WISE curves are appended to FILTER.RES.latest as filters N+1..N+6.
"""

from __future__ import annotations

import contextlib
import logging
import os
import subprocess
import tempfile
from collections.abc import Mapping
from pathlib import Path

import numpy as np

from ..cosmo import cosmology, distmod
from ..env import nproc
from .common import BANDS, empty_result, photometry, valid_z
from .filters import KCORRECT_NAMES, pivot_wavelength, read_response

log = logging.getLogger(__name__)
REPO = "https://github.com/gbrammer/eazy-photoz.git"


@contextlib.contextmanager
def _cwd(path: Path):
    old = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


class EazyBackend:
    name = "eazy"

    def __init__(self, cfg: dict, workdir: str | Path | None = None, n_proc: int | None = None,
                 templates: str = "templates/sfhz/corr_sfhz_13.param"):
        import eazy.filters as F

        self.cfg = cfg
        self.n_proc = n_proc or nproc()
        self.templates = templates
        data = Path(os.environ.get("EAZY_DATA", "~/.cache/eazy-photoz")).expanduser()
        if not (data / "templates").exists():
            log.info("cloning eazy-photoz into %s", data)
            subprocess.run(["git", "clone", "-q", "--depth", "1", REPO, str(data)], check=True)
        if workdir is None:
            from ..env import get_paths

            workdir = get_paths().out / "sedwork" / "eazy"
        self.work = Path(workdir)
        self.work.mkdir(parents=True, exist_ok=True)
        link = self.work / "templates"
        if not link.exists():
            link.symlink_to(data / "templates")
        ff = F.FilterFile(str(data / "filters" / "FILTER.RES.latest"))
        n0 = ff.NFILT
        for b in BANDS:
            wl, tr = read_response(b)
            ff.filters.append(F.FilterDefinition(
                name=f"ls11_{KCORRECT_NAMES[b]} lambda_c= {pivot_wavelength(b):.4e}", wave=wl,
                throughput=tr, photon_counter=True))
        ff.write(str(self.work / "FILTER.RES.ls11"), verbose=False)
        self.fnum = {b: n0 + 1 + i for i, b in enumerate(BANDS)}
        (self.work / "zphot.translate").write_text("id id\n")

    def fit(self, t: Mapping) -> dict[str, np.ndarray]:
        from astropy.table import Table
        import eazy.photoz as P

        n = len(t["BEST_Z"])
        out = empty_result(n)
        z, f, e = photometry(t, self.cfg)
        ok = valid_z(z)
        if not ok.any():
            return out
        cat = Table()
        cat["id"] = np.arange(ok.sum())
        cat["z_spec"] = z[ok]
        for i, b in enumerate(BANDS):
            good = np.isfinite(e[ok, i]) & np.isfinite(f[ok, i])
            cat[f"F{self.fnum[b]}"] = np.where(good, f[ok, i], -99.0)
            cat[f"E{self.fnum[b]}"] = np.where(good, e[ok, i], -99.0)
        fd, catfile = tempfile.mkstemp(suffix=".fits", prefix="cat_", dir=self.work)
        os.close(fd)
        cat.write(catfile, overwrite=True)
        params = {"CATALOG_FILE": catfile, "MAIN_OUTPUT_FILE": catfile.replace(".fits", ""),
                  "FILTERS_RES": "FILTER.RES.ls11", "TEMPLATES_FILE": self.templates,
                  "TEMP_ERR_FILE": "templates/TEMPLATE_ERROR.eazy_v1.0", "TEMP_ERR_A2": 0.0,
                  "SYS_ERR": 0.0, "Z_MIN": 0.005, "Z_MAX": 1.0, "Z_STEP": 0.005, "Z_STEP_TYPE": 0,
                  "APPLY_PRIOR": "n", "PRIOR_ABZP": 22.5, "MW_EBV": 0.0, "CAT_HAS_EXTCORR": "y",
                  "FIX_ZSPEC": "n", "N_MIN_COLORS": 3}
        with _cwd(self.work):
            pz = P.PhotoZ(param_file=None, translate_file="zphot.translate", zeropoint_file=None,
                          load_prior=False, params=params, n_proc=self.n_proc, cosmology=cosmology(self.cfg))
            pz.fit_at_zbest(zbest=z[ok], prior=False, beta_prior=False, n_proc=self.n_proc)
            sps = pz.sps_parameters(extra_rf_filters=[self.fnum["R"]], n_proc=self.n_proc)
        os.remove(catfile)
        idx = np.flatnonzero(ok)
        with np.errstate(divide="ignore", invalid="ignore"):
            out["LOGMSTAR"][idx] = np.log10(np.asarray(sps["mass"], np.float64))
            out["LOGSFR"][idx] = np.log10(np.asarray(sps["sfr"], np.float64))
            rf = np.asarray(sps[f"rest{self.fnum['R']}"], np.float64)
            # rest-frame r flux (catalogue units, ZP 22.5) seen from the galaxy's distance
            out["ABSMAG_R"][idx] = 22.5 - 2.5 * np.log10(rf) - distmod(z[ok], self.cfg) + 2.5 * np.log10(1 + z[ok])
        out["LOGMSTAR_LO"][:] = out["LOGMSTAR"]
        out["LOGMSTAR_HI"][:] = out["LOGMSTAR"]
        out["CHI2"][idx] = np.asarray(pz.chi2_best, np.float64)
        return out
