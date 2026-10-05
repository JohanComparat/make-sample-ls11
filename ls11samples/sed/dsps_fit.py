"""DSPS (differentiable SPS, jax) grid fit at fixed redshift.

DSPS provides SSP spectra and photometry kernels, not a fitter; this module builds one:

* SSPs: FSPS v3.2 defaults (``ssp_data_fsps_v3.2_lgmet_age.h5``, Kroupa IMF, nebular emission),
  from $DSPS_DRN (default ~/.cache/dsps);
* models: delayed-tau SFH (tau x age since formation) x metallicity x Calzetti (2000) A_V; the
  rest-frame SED of each model per Msun formed is built once, and its observed AB magnitudes in
  the shared DECam/WISE curves are tabulated on a redshift grid with jax;
* fit: for each galaxy (nearest grid redshift) the mass amplitude of every model is solved
  analytically (linear least squares in flux), models older than the Universe at z are excluded,
  and logM* percentiles come from chi^2 weights exp(-chi^2/2) over the grid;
* surviving mass: formed mass times 1 - 0.05 ln(1 + t / 1.4 Myr) (Behroozi et al. 2013, eq. 14)
  integrated over the SFH; masses converted Kroupa -> Chabrier by -0.034 dex (Madau & Dickinson 2014).
"""

from __future__ import annotations

import logging
import os
from collections.abc import Mapping
from pathlib import Path

import numpy as np

from ..cosmo import cosmology
from .common import BANDS, empty_result, photometry, valid_z
from .filters import read_response

log = logging.getLogger(__name__)
SSP_FILE = "ssp_data_fsps_v3.2_lgmet_age.h5"
KROUPA_TO_CHABRIER = -0.034
Z_SUN = 0.0142

DEFAULT_GRID = {"tau_gyr": [0.3, 0.6, 1.0, 2.0, 4.0, 8.0],
                "age_gyr": [0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0],
                "z_over_zsun": [0.2, 0.5, 1.0, 1.6],
                "av": [0.0, 0.1, 0.2, 0.4, 0.6, 0.9, 1.2, 1.6]}


def surviving_fraction(t_gyr):
    """Behroozi et al. (2013) eq. 14: fraction of the formed mass still in stars after t."""
    return 1.0 - 0.05 * np.log1p(np.asarray(t_gyr) * 1e3 / 1.4)


def age_weights(lg_age_gyr: np.ndarray, tau: float, age: float, n: int = 4000):
    """(weights over the SSP ages, surviving mass fraction, SFR now [Msun/yr]) of a delayed-tau
    SFH that started ``age`` Gyr ago, normalised to 1 Msun formed."""
    lookback = np.linspace(0, age, n + 1)[1:] - age / (2 * n)       # bin centres
    t_since = age - lookback
    sfr = t_since / tau * np.exp(-t_since / tau)
    w = sfr / sfr.sum()
    edges = np.concatenate([[-np.inf], 0.5 * (lg_age_gyr[1:] + lg_age_gyr[:-1]), [np.inf]])
    idx = np.clip(np.digitize(np.log10(np.maximum(lookback, 1e-5)), edges) - 1, 0, lg_age_gyr.size - 1)
    weights = np.bincount(idx, weights=w, minlength=lg_age_gyr.size)
    sfr_now = (age / tau) * np.exp(-age / tau) / (sfr.sum() * age / n) / 1e9
    return weights, float(np.sum(w * surviving_fraction(lookback))), float(sfr_now)


class DspsBackend:
    name = "dsps"

    def __init__(self, cfg: dict, ssp_file: str | Path | None = None, grid: dict | None = None,
                 zgrid: np.ndarray | None = None):
        import jax
        import jax.numpy as jnp
        from dsps import calc_obs_mag, calc_rest_mag, load_ssp_templates
        from dsps.dust.att_curves import calzetti00_att_curve

        jax.config.update("jax_enable_x64", True)
        self.cfg = cfg
        g = {**DEFAULT_GRID, **(grid or ((cfg.get("sed") or {}).get("dsps") or {}).get("grid", {}))}
        drn = Path(os.environ.get("DSPS_DRN", "~/.cache/dsps")).expanduser()
        ssp = load_ssp_templates(fn=str(ssp_file or drn / SSP_FILE))
        wave = np.asarray(ssp.ssp_wave, np.float64)
        flux = np.asarray(ssp.ssp_flux, np.float64)                       # (n_met, n_age, n_wave)
        lgmet = np.asarray(ssp.ssp_lgmet)
        lg_age = np.asarray(ssp.ssp_lg_age_gyr)

        seds, pars = [], []
        for zz in g["z_over_zsun"]:
            im = int(np.argmin(np.abs(lgmet - np.log10(zz * Z_SUN))))
            for tau in g["tau_gyr"]:
                for age in g["age_gyr"]:
                    w, fsurv, sfr_now = age_weights(lg_age, tau, age)
                    sed = w @ flux[im]                                     # Lsun/Hz per Msun formed
                    for av in g["av"]:
                        # calzetti00_att_curve returns A(lambda) in mag, not a transmission
                        trans = 10 ** (-0.4 * np.asarray(calzetti00_att_curve(wave / 1e4, av)))
                        seds.append(sed * trans)
                        pars.append((tau, age, zz, av, fsurv, sfr_now))
        self.seds = np.array(seds)
        self.pars = np.array(pars)
        self.model_age = self.pars[:, 1]
        self.fsurv = self.pars[:, 4]
        self.sfr_now = self.pars[:, 5]
        log.info("DSPS grid: %d models", len(self.seds))

        cos = cosmology(cfg)
        h, om0 = cos.H0.value / 100, cos.Om0
        self.zgrid = np.arange(0.005, 1.0001, 0.005) if zgrid is None else np.asarray(zgrid)
        self.t_univ = cos.age(self.zgrid).value
        filters = [read_response(b) for b in BANDS]
        obs = jax.jit(jax.vmap(lambda s, wf, tf, z: calc_obs_mag(wave, s, wf, tf, z, om0, -1.0, 0.0, h),
                               in_axes=(0, None, None, None)))
        mags = np.empty((self.zgrid.size, len(self.seds), len(BANDS)))
        seds_j = jnp.asarray(self.seds)
        for ib, (wf, tf) in enumerate(filters):
            wf, tf = jnp.asarray(wf), jnp.asarray(tf)
            for iz, z in enumerate(self.zgrid):
                mags[iz, :, ib] = np.asarray(obs(seds_j, wf, tf, z))
        self.model_nmgy = 10 ** (-0.4 * (mags - 22.5))                  # nanomaggies per Msun formed
        wr, tr = filters[BANDS.index("R")]
        rest = jax.jit(jax.vmap(lambda s: calc_rest_mag(wave, s, jnp.asarray(wr), jnp.asarray(tr))))
        self.rest_r = np.asarray(rest(seds_j))                            # abs. mag per Msun formed

    def fit(self, t: Mapping, chunk: int = 2000) -> dict[str, np.ndarray]:
        n = len(t["BEST_Z"])
        out = empty_result(n)
        z, f, e = photometry(t, self.cfg)
        ok = valid_z(z, self.zgrid[-1])
        idx_all = np.flatnonzero(ok)
        iz_all = np.clip(np.searchsorted(self.zgrid, z[ok]), 0, self.zgrid.size - 1)
        for s in range(0, idx_all.size, chunk):
            idx, iz = idx_all[s:s + chunk], iz_all[s:s + chunk]
            ff, ee = f[idx], e[idx]
            ivar = np.where(np.isfinite(ee) & np.isfinite(ff), 1 / ee**2, 0.0)
            ff = np.where(ivar > 0, ff, 0.0)
            m = self.model_nmgy[iz]                                       # (k, n_model, n_band)
            num = np.einsum("kb,kmb->km", ff * ivar, m)
            den = np.einsum("kb,kmb->km", ivar, m**2)
            amp = num / den
            chi2 = np.einsum("kb,kmb->km", ivar, (ff[:, None, :] - amp[..., None] * m) ** 2)
            valid = (amp > 0) & np.isfinite(chi2) & (self.model_age[None, :] <= self.t_univ[iz][:, None])
            chi2 = np.where(valid, chi2, np.inf)
            best = np.argmin(chi2, axis=1)
            chi2min = chi2[np.arange(len(idx)), best]
            with np.errstate(divide="ignore", invalid="ignore"):
                logm = np.log10(amp * self.fsurv[None, :]) + KROUPA_TO_CHABRIER
                w = np.where(valid, np.exp(-0.5 * (chi2 - chi2min[:, None])), 0.0)
            pct = _weighted_percentiles(np.where(valid, logm, 0.0), w, (16, 50, 84))
            out["LOGMSTAR_LO"][idx], out["LOGMSTAR"][idx], out["LOGMSTAR_HI"][idx] = pct
            out["LOGMSTAR_ERR"][idx] = 0.5 * (pct[2] - pct[0])
            out["CHI2"][idx] = chi2min
            ab = amp[np.arange(len(idx)), best]
            with np.errstate(divide="ignore", invalid="ignore"):
                out["MABS_R"][idx] = self.rest_r[best] - 2.5 * np.log10(ab)
                out["LOGSFR"][idx] = np.log10(ab * self.sfr_now[best]) + KROUPA_TO_CHABRIER
        return out


def _weighted_percentiles(values: np.ndarray, weights: np.ndarray, pcts) -> list[np.ndarray]:
    """Row-wise weighted percentiles of ``values`` (k, m) with ``weights`` (k, m)."""
    order = np.argsort(values, axis=1)
    v = np.take_along_axis(values, order, axis=1)
    w = np.take_along_axis(weights, order, axis=1)
    cw = np.cumsum(w, axis=1)
    tot = cw[:, -1:]
    res = []
    for p in pcts:
        with np.errstate(invalid="ignore", divide="ignore"):
            i = np.argmax(cw >= p / 100 * tot, axis=1)
        r = v[np.arange(len(v)), i]
        res.append(np.where(tot[:, 0] > 0, r, np.nan))
    return res
