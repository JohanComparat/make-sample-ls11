"""LePhare (python, >= 1.0) at fixed redshift, configured by config/lephare/LS11_zFIX.para.

The first use writes the shared DECam/WISE filters into $LEPHAREDIR/filt/ls11, downloads the
auxiliary files the configuration needs (SEDs, extinction law) and builds the libraries (filters,
sedtolib, mag_gal) in $LEPHAREWORK; later uses reuse them.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Mapping
from pathlib import Path

import numpy as np

from .common import BANDS, NANOMAGGY_CGS, empty_result, photometry, valid_z
from .filters import write_lephare_filters

log = logging.getLogger(__name__)
PARA = Path(__file__).resolve().parents[2] / "config" / "lephare" / "LS11_zFIX.para"
PARA_OUT = PARA.with_name("output_ls11.para")


class LephareBackend:
    name = "lephare"

    def __init__(self, cfg: dict, para: str | Path | None = None, rebuild: bool = False, build: bool = True):
        import lephare as lp

        self.lp = lp
        self.cfg = cfg
        self.config = lp.read_config(str(para or PARA))
        lephare_dir = Path(os.environ.get("LEPHAREDIR", lp.LEPHAREDIR))
        self.config["FILTER_REP"] = lp.keyword("FILTER_REP", str(lephare_dir / "filt"))
        self.config["PARA_OUT"] = lp.keyword("PARA_OUT", str(PARA_OUT))
        write_lephare_filters(lephare_dir)
        self._get_data(lephare_dir)
        work = Path(os.environ["LEPHAREWORK"])
        # the .doc file is written last: an empty one means an interrupted build
        doc = work / "lib_mag" / f"{self.config['GAL_LIB_OUT'].value}.doc"
        if build and (rebuild or not doc.exists() or doc.stat().st_size == 0):
            log.info("building the LePhare libraries in %s", work)
            lp.prepare(self.config)

    def _get_data(self, lephare_dir: Path) -> None:
        dr = self.lp.data_retrieval
        files = [f for f in dr.config_to_required_files(self.config) if "filt/ls11/" not in f]
        missing = [f for f in files if not (lephare_dir / f).exists()]
        if missing:
            import tempfile

            log.info("downloading %d LePhare auxiliary files to %s", len(missing), lephare_dir)
            with tempfile.TemporaryDirectory() as tmp:
                registry = str(Path(tmp) / "data_registry.txt")
                dr.download_registry_from_github(outfile=registry)
                retriever = dr.make_retriever(registry_file=registry, data_path=str(lephare_dir))
                dr.download_all_files(retriever, missing, ignore_registry=False)

    def table(self, t: Mapping):
        from astropy.table import Table

        z, f, e = photometry(t, self.cfg)
        tab = Table()
        tab["id"] = np.arange(len(z))
        good = np.isfinite(e) & np.isfinite(f)
        for i, b in enumerate(BANDS):
            tab[f"f{i}"] = np.where(good[:, i], f[:, i] * NANOMAGGY_CGS, -99.0)
            tab[f"e{i}"] = np.where(good[:, i], e[:, i] * NANOMAGGY_CGS, -99.0)
        tab["context"] = (good * (1 << np.arange(len(BANDS)))).sum(axis=1).astype(np.int64)
        tab["zspec"] = z
        tab["string_input"] = np.full(len(z), "x")
        return tab, valid_z(z)

    def fit(self, t: Mapping) -> dict[str, np.ndarray]:
        n = len(t["BEST_Z"])
        out = empty_result(n)
        tab, ok = self.table(t)
        if not ok.any():
            return out
        res, _ = self.lp.process(self.config, tab[ok])
        idx = np.flatnonzero(ok)

        def col(name):
            v = np.asarray(res[name], np.float64)
            return np.where((v > -99) & np.isfinite(v), v, np.nan)

        med, lo, hi, best = col("MASS_MED"), col("MASS_INF"), col("MASS_SUP"), col("MASS_BEST")
        out["LOGMSTAR"][idx] = np.where(np.isfinite(med), med, best)
        out["LOGMSTAR_LO"][idx] = lo
        out["LOGMSTAR_HI"][idx] = hi
        out["LOGMSTAR_ERR"][idx] = 0.5 * (hi - lo)              # half the 68% interval of the PDF
        out["LOGSFR"][idx] = col("SFR_MED")
        ir = BANDS.index("R")
        mabs = np.asarray(res["MAG_ABS()"], np.float64)[:, ir]  # (n, nfilter), filter order = BANDS
        emabs = np.asarray(res["EMAG_ABS()"], np.float64)[:, ir]
        good = (mabs > -90) & (mabs < 0)
        out["MABS_R"][idx] = np.where(good, mabs, np.nan)
        out["MABS_R_ERR"][idx] = np.where(good & (emabs >= 0) & (emabs < 9), emabs, np.nan)
        out["CHI2"][idx] = col("CHI_BEST")
        self.last = res
        return out
