"""kcorrect v5 as an SED backend (NMF templates; mass of the best-fit template mix)."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from ..kcorr import KcorrectV5
from .common import empty_result


class KcorrectBackend:
    name = "kcorrect"

    def __init__(self, cfg: dict):
        self.kc = KcorrectV5(cfg)

    def fit(self, t: Mapping) -> dict[str, np.ndarray]:
        n = len(t["BEST_Z"])
        out = empty_result(n)
        if n == 0:
            return out
        tt = dict(t)
        if "MAG_R" not in tt:
            from ..photometry import mag

            for b in ("G", "R", "Z"):
                tt[f"MAG_{b}"] = mag(t[f"FLUX_{b}"], t[f"MW_TRANSMISSION_{b}"])
        r = self.kc.fit(tt)
        out["LOGMSTAR"] = r["KC_LOGMSTAR"]
        out["LOGMSTAR_LO"] = r["KC_LOGMSTAR"]
        out["LOGMSTAR_HI"] = r["KC_LOGMSTAR"]
        out["ABSMAG_R"] = r["ABSMAG_R"]
        out["CHI2"] = r["KC_CHI2"]
        return out
