"""kcorrect v5 (Blanton & Roweis 2007) for the LS DR11 bands at fixed z = BEST_Z.

:meth:`KcorrectV5.fit` returns per galaxy (NaN when not fitted):

  KCORR_<b>, ABSMAG_<b>   b = G, R, Z (DECam), K-corrected to the band shifted to ``band_shift``
                          (default 0: rest-frame DECam bands); ABSMAG = MAG - DM(z) - KCORR
  KCORR_R01, ABSMAG_R01   SDSS r shifted to z = 0.1 (^{0.1}r), from DECam r
  KC_LOGMSTAR             log10 surviving stellar mass of the template fit (Msun, Chabrier)
  ABSMAG_R_ERR, KC_LOGMSTAR_ERR   standard deviations over ``mc`` Monte Carlo realisations of the
                          fluxes (kcorrect ``fit_coeffs(mc=...)``), 0 realisations -> NaN
  KC_CHI2, KC_COEFFS      fit chi^2 and template coefficients

:meth:`KcorrectV5.k_curve` evaluates the r-band K-correction of a set of fitted SEDs at every z
(completeness limit of the Mr samples, :mod:`ls11samples.vlim`).
"""

from __future__ import annotations

import logging
import warnings
from collections.abc import Mapping

import numpy as np

from .cosmo import cosmology, distmod
from .photometry import flux_arrays, mag

log = logging.getLogger(__name__)

RESPONSES = {"G": "decam_g", "R": "decam_r", "I": "decam_i", "Z": "decam_z",
             "W1": "wise_w1", "W2": "wise_w2"}
OUT_BANDS = ("G", "R", "Z")


class KcorrectV5:
    """kcorrect v5 fitter for the LS DR11 bands."""

    def __init__(self, cfg: dict, mc: int | None = None):
        import kcorrect.kcorrect as kk

        kc = cfg["kcorr"]
        self.cfg = cfg
        self.bands = list(kc["bands_in"])
        self.band_shift = float(kc.get("band_shift", 0.0))
        self.mc = int(kc.get("mc", 20) if mc is None else mc)
        self.seed = 0
        out = [RESPONSES[b] for b in OUT_BANDS] + ["sdss_r0"]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.kc = kk.Kcorrect(responses=[RESPONSES[b] for b in self.bands], responses_out=out,
                                  responses_map=[RESPONSES[b] for b in OUT_BANDS] + ["decam_r"],
                                  redshift_range=list(kc.get("z_range", [0.0, 1.0])),
                                  nredshift=int(kc.get("nz", 500)), cosmo=cosmology(cfg))
        self.ntemplates = len(self.kc.templates.mremain)

    def maggies(self, t: Mapping) -> tuple[np.ndarray, np.ndarray]:
        f, e = flux_arrays(t, self.bands, self.cfg["photometry"]["err_floor_mag"])
        maggies = (f * 1e-9).astype(np.float32)
        ivar = np.where(np.isfinite(e) & (e > 0), 1.0 / (e * 1e-9) ** 2, 0.0).astype(np.float32)
        return maggies, ivar

    def fit(self, t: Mapping) -> dict[str, np.ndarray]:
        z = np.asarray(t["BEST_Z"], np.float64)
        n = z.size
        ok = np.isfinite(z) & (z > 0.001) & (z < self.kc.redshift_range[1])
        out = {"KC_COEFFS": np.zeros((n, self.ntemplates), np.float32)}
        for b in OUT_BANDS:
            out[f"KCORR_{b}"] = np.full(n, np.nan, np.float32)
            out[f"ABSMAG_{b}"] = np.full(n, np.nan, np.float32)
        out.update({k: np.full(n, np.nan, np.float32) for k in
                    ("KCORR_R01", "ABSMAG_R01", "KC_LOGMSTAR", "KC_CHI2", "ABSMAG_R_ERR", "KC_LOGMSTAR_ERR")})
        if not ok.any():
            return out
        maggies, ivar = self.maggies(t)
        zz = z[ok].astype(np.float32)
        mags = {b: mag(t[f"FLUX_{b}"], t[f"MW_TRANSMISSION_{b}"])[ok] for b in OUT_BANDS}
        dm = distmod(zz, self.cfg)
        dfac = 10.0 ** (0.4 * dm)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            np.random.seed(self.seed)                      # kcorrect draws its MC fluxes with np.random
            res = self.kc.fit_coeffs(redshift=zz, maggies=maggies[ok], ivar=ivar[ok], mc=self.mc)
            coeffs, coeffs_mc, maggies_mc = res if self.mc > 0 else (res, None, None)
            k = self.kc.kcorrect(redshift=zz, coeffs=coeffs, band_shift=self.band_shift)
            k01 = self.kc.kcorrect(redshift=zz, coeffs=coeffs, band_shift=0.1)
            rec = self.kc.reconstruct(redshift=zz, coeffs=coeffs)
            if self.mc > 0:
                ir = self.bands.index("R")
                mr_mc = np.empty((zz.size, self.mc))
                lm_mc = np.empty((zz.size, self.mc))
                for i in range(self.mc):
                    k_i = self.kc.kcorrect(redshift=zz, coeffs=coeffs_mc[..., i], band_shift=self.band_shift)
                    with np.errstate(divide="ignore", invalid="ignore"):
                        m_i = np.where(maggies_mc[:, ir, i] > 0, -2.5 * np.log10(maggies_mc[:, ir, i]), np.nan)
                        mr_mc[:, i] = m_i - dm - k_i[:, 1]
                        lm_mc[:, i] = np.log10(coeffs_mc[..., i].dot(self.kc.templates.mremain) * dfac)
                out["ABSMAG_R_ERR"][ok] = np.nanstd(mr_mc, axis=1)
                out["KC_LOGMSTAR_ERR"][ok] = np.nanstd(np.where(np.isfinite(lm_mc), lm_mc, np.nan), axis=1)
        out["KC_COEFFS"][ok] = coeffs
        for i, b in enumerate(OUT_BANDS):
            out[f"KCORR_{b}"][ok] = k[:, i]
            out[f"ABSMAG_{b}"][ok] = mags[b] - dm - k[:, i]
        out["KCORR_R01"][ok] = k01[:, 3]
        out["ABSMAG_R01"][ok] = mags["R"] - dm - k01[:, 3]
        mremain = coeffs.dot(self.kc.templates.mremain) * dfac
        with np.errstate(divide="ignore", invalid="ignore"):
            out["KC_LOGMSTAR"][ok] = np.where(mremain > 0, np.log10(mremain), np.nan)
        out["KC_CHI2"][ok] = np.sum((rec - maggies[ok]) ** 2 * ivar[ok], axis=1)
        return out

    def k_curve(self, coeffs: np.ndarray, zgrid: np.ndarray, pct: float = 95, band: int = 1) -> np.ndarray:
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


def complete_set_rows(z: np.ndarray, mabs_r: np.ndarray, mag_r: np.ndarray, r_lim: float, cfg: dict,
                      z_max: float = 0.1, n_max: int = 20000, seed: int = 0) -> np.ndarray:
    """Indices of galaxies forming a complete set of SEDs: z <= z_max and Mr brighter than the
    limit of the reddest galaxies at z_max (99th percentile of the K-correction there)."""
    z, mabs_r, mag_r = (np.asarray(a, np.float64) for a in (z, mabs_r, mag_r))
    low = (z > 0.01) & (z <= z_max) & np.isfinite(mabs_r)
    kcorr = mag_r - distmod(z, cfg) - mabs_r
    kmax = np.nanpercentile(kcorr[low], 99)
    mlim = r_lim - distmod(np.array([z_max]), cfg)[0] - kmax
    idx = np.flatnonzero(low & (mabs_r <= mlim))
    if idx.size > n_max:
        idx = np.sort(np.random.default_rng(seed).choice(idx, n_max, replace=False))
    log.info("complete SED set: %d galaxies with z <= %.2f and Mr <= %.2f", idx.size, z_max, mlim)
    return idx
