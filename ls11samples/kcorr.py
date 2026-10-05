"""Step 3: K-corrections and absolute magnitudes with kcorrect v5 (Blanton & Roweis 2007).

The fit uses the dereddened grizW1W2 fluxes (with an error floor) at fixed z = BEST_Z. Outputs per
galaxy:

  KCORR_<b>, ABSMAG_<b>   b = G, R, Z (DECam), K-corrected to the band shifted to ``band_shift``
                          (default 0: rest-frame DECam bands); ABSMAG = MAG - DM(z) - KCORR
  KCORR_R01, ABSMAG_R01   SDSS r shifted to z = 0.1 (^{0.1}r, as in Zehavi et al. 2011), from DECam r
  KC_LOGMSTAR             log10 surviving stellar mass of the template fit (Msun, Chabrier)
  KC_CHI2, KC_COEFFS      fit chi^2 and template coefficients

:func:`k_percentile_curve` evaluates the r-band K-correction of a complete low-z set of fitted
SEDs at every z, so the Mr completeness limit does not inherit the bias of the K distribution of the
observed (flux-limited) galaxies at high z.
"""

from __future__ import annotations

import logging
import time
import warnings
from collections.abc import Mapping
from pathlib import Path

import numpy as np

from . import io
from .cosmo import cosmology, distmod
from .photometry import flux_arrays

log = logging.getLogger(__name__)

RESPONSES = {"G": "decam_g", "R": "decam_r", "I": "decam_i", "Z": "decam_z",
             "W1": "wise_w1", "W2": "wise_w2"}
OUT_BANDS = ("G", "R", "Z")


class KcorrectV5:
    """kcorrect v5 fitter for the LS DR11 bands."""

    def __init__(self, cfg: dict):
        import kcorrect.kcorrect as kk

        kc = cfg["kcorr"]
        self.cfg = cfg
        self.bands = list(kc["bands_in"])
        self.band_shift = float(kc.get("band_shift", 0.0))
        out = [RESPONSES[b] for b in OUT_BANDS] + ["sdss_r0"]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.kc = kk.Kcorrect(responses=[RESPONSES[b] for b in self.bands], responses_out=out,
                                  responses_map=[RESPONSES[b] for b in OUT_BANDS] + ["decam_r"],
                                  redshift_range=list(kc.get("z_range", [0.0, 1.0])),
                                  nredshift=int(kc.get("nz", 500)), cosmo=cosmology(cfg))

    def maggies(self, t: Mapping) -> tuple[np.ndarray, np.ndarray]:
        f, e = flux_arrays(t, self.bands, self.cfg["photometry"]["err_floor_mag"])
        maggies = (f * 1e-9).astype(np.float32)
        ivar = np.where(np.isfinite(e) & (e > 0), 1.0 / (e * 1e-9) ** 2, 0.0).astype(np.float32)
        return maggies, ivar

    def fit(self, t: Mapping) -> dict[str, np.ndarray]:
        z = np.asarray(t["BEST_Z"], np.float64)
        n = z.size
        ok = np.isfinite(z) & (z > 0) & (z < self.kc.redshift_range[1])
        out = {"KC_COEFFS": np.zeros((n, len(self.kc.templates.mremain)), np.float32)}
        for b in OUT_BANDS:
            out[f"KCORR_{b}"] = np.full(n, np.nan, np.float32)
            out[f"ABSMAG_{b}"] = np.full(n, np.nan, np.float32)
        out.update({k: np.full(n, np.nan, np.float32)
                    for k in ("KCORR_R01", "ABSMAG_R01", "KC_LOGMSTAR", "KC_CHI2")})
        if not ok.any():
            return out
        maggies, ivar = self.maggies(t)
        zz = z[ok].astype(np.float32)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            coeffs = self.kc.fit_coeffs(redshift=zz, maggies=maggies[ok], ivar=ivar[ok])
            k = self.kc.kcorrect(redshift=zz, coeffs=coeffs, band_shift=self.band_shift)
            k01 = self.kc.kcorrect(redshift=zz, coeffs=coeffs, band_shift=0.1)
            rec = self.kc.reconstruct(redshift=zz, coeffs=coeffs)
        dm = distmod(zz, self.cfg)
        out["KC_COEFFS"][ok] = coeffs
        for i, b in enumerate(OUT_BANDS):
            out[f"KCORR_{b}"][ok] = k[:, i]
            out[f"ABSMAG_{b}"][ok] = np.asarray(t[f"MAG_{b}"])[ok] - dm - k[:, i]
        out["KCORR_R01"][ok] = k01[:, 3]
        out["ABSMAG_R01"][ok] = np.asarray(t["MAG_R"])[ok] - dm - k01[:, 3]
        mremain = coeffs.dot(self.kc.templates.mremain) * 10.0 ** (0.4 * dm)
        with np.errstate(divide="ignore", invalid="ignore"):
            out["KC_LOGMSTAR"][ok] = np.where(mremain > 0, np.log10(mremain), np.nan)
        out["KC_CHI2"][ok] = np.sum((rec - maggies[ok]) ** 2 * ivar[ok], axis=1)
        return out

    def k_curve(self, coeffs: np.ndarray, zgrid: np.ndarray, pct: float = 95, band: int = 1
                ) -> np.ndarray:
        """p-th percentile over the SED set ``coeffs`` of the K-correction (output band ``band``,
        default DECam r) at each z of ``zgrid``."""
        coeffs = np.asarray(coeffs, np.float32)
        out = np.empty(len(zgrid))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for i, zi in enumerate(zgrid):
                k = self.kc.kcorrect(redshift=np.full(len(coeffs), zi, np.float32), coeffs=coeffs,
                                     band_shift=self.band_shift)
                out[i] = np.nanpercentile(k[:, band], pct)
        return out


def complete_sed_set(kc_out: Mapping, z: np.ndarray, r_lim: float, cfg: dict, z_max: float = 0.1,
                     n_max: int = 20000, seed: int = 0) -> np.ndarray:
    """Coefficients of galaxies forming a complete set of SEDs: z <= z_max and Mr brighter than
    the reddest-galaxy limit at z_max, so no SED type is lost to the flux limit."""
    z = np.asarray(z, np.float64)
    kmax = np.nanpercentile(np.asarray(kc_out["KCORR_R"])[z <= z_max], 99)
    mlim = r_lim - distmod(np.array([z_max]), cfg)[0] - kmax
    sel = (z > 0.01) & (z <= z_max) & (np.asarray(kc_out["ABSMAG_R"]) <= mlim)
    sel &= np.asarray(kc_out["KC_COEFFS"]).sum(axis=1) > 0
    idx = np.flatnonzero(sel)
    if idx.size > n_max:
        idx = np.random.default_rng(seed).choice(idx, n_max, replace=False)
    log.info("complete SED set: %d galaxies with z <= %.2f and Mr <= %.2f", idx.size, z_max, mlim)
    return np.asarray(kc_out["KC_COEFFS"])[idx]


KC_INPUT = ["LS_ID_DR11", "BEST_Z", "MAG_G", "MAG_R", "MAG_Z"] + [
    f"{p}_{b}" for p in ("FLUX", "FLUX_IVAR", "MW_TRANSMISSION") for b in RESPONSES]


def output_name(sweep: str | Path) -> str:
    return "KC-" + Path(sweep).name


def run_sweep(parent_file: str | Path, cfg: dict, outdir: str | Path, overwrite: bool = False) -> Path:
    """K-corrections of the BGS-like galaxies (SEL_FLAGS == 0) of one per-sweep parent file."""
    from .bgsl import read_selected

    parent_file = Path(parent_file)
    out = Path(outdir) / output_name(parent_file.name.removeprefix("BGSl-"))
    if out.exists() and not overwrite:
        return out
    t0 = time.time()
    t = read_selected([parent_file], KC_INPUT)
    if not t:
        t = {c: np.zeros(0) for c in KC_INPUT}
    res = KcorrectV5(cfg).fit(t)
    io.write_table(out, {"LS_ID_DR11": t["LS_ID_DR11"], **res},
                   header={"KCCODE": "kcorrect", "BANDSHFT": float(cfg["kcorr"].get("band_shift", 0.0)),
                           "BANDS": ",".join(cfg["kcorr"]["bands_in"])}, extname="KCORR")
    log.info("%s: %d galaxies (%.0fs)", out.name, len(t["LS_ID_DR11"]), time.time() - t0)
    return out
