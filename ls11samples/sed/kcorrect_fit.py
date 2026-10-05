"""kcorrect v5 as an SED backend (NMF templates; mass of the best-fit template mix; errors from
Monte Carlo realisations of the fluxes)."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from ..kcorr import KcorrectV5
from .common import empty_result


class KcorrectBackend:
    name = "kcorrect"

    def __init__(self, cfg: dict):
        self.kc = KcorrectV5(cfg)

    @property
    def seed(self):
        return self.kc.seed

    @seed.setter
    def seed(self, value):
        self.kc.seed = int(value)

    def fit(self, t: Mapping) -> dict[str, np.ndarray]:
        n = len(t["BEST_Z"])
        out = empty_result(n)
        if n == 0:
            return out
        r = self.kc.fit(t)
        out["LOGMSTAR"] = r["KC_LOGMSTAR"]
        out["LOGMSTAR_ERR"] = r["KC_LOGMSTAR_ERR"]
        out["LOGMSTAR_LO"] = r["KC_LOGMSTAR"] - r["KC_LOGMSTAR_ERR"]
        out["LOGMSTAR_HI"] = r["KC_LOGMSTAR"] + r["KC_LOGMSTAR_ERR"]
        out["MABS_R"] = r["ABSMAG_R"]
        out["MABS_R_ERR"] = r["ABSMAG_R_ERR"]
        out["CHI2"] = r["KC_CHI2"]
        return out
