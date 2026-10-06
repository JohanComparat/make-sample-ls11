"""BGS-like selection of DR11 sweep objects.

Two kinds of cuts, both configured in ``config/default.yaml``:

* footprint cuts (:func:`footprint_cuts`): NOBS, galaxy depth, MASKBITS, E(B-V), south-only. They
  only use columns the randoms also have and the SAME function is applied to the randoms, so
  galaxies and randoms share one footprint;
* galaxy cuts (:func:`galaxy_cuts`): TYPE, FITBITS, flux ivar, Gaia star rejection, colours,
  DESI BGS quality and fibre-magnitude cuts and the r-band range.

The selected sample passes every cut; :func:`cutflow` counts the objects cut by cut, in the order
of :data:`CUTS`.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from . import bits
from .photometry import mag

#: Cuts applied to galaxies and randoms alike.
FOOTPRINT_CUTS = ("nobs", "galdepth", "maskbits", "ebv", "south")
#: Object-level cuts, galaxies only.
GALAXY_CUTS = ("rmag", "type", "fitbits", "flux_ivar", "gaia", "colour", "quality", "fiber", "rfibtot")
#: All cuts, in cut-flow order.
CUTS = FOOTPRINT_CUTS + GALAXY_CUTS

# ICRS -> Galactic rotation (Hipparcos / astropy), rows are the Galactic x, y, z axes.
_ICRS_TO_GAL = np.array([[-0.0548755604162154, -0.8734370902348850, -0.4838350155487132],
                         [+0.4941094278755837, -0.4448296299600112, +0.7469822444972189],
                         [-0.8676661490190047, -0.1980763734312015, +0.4559837761750669]])
NORTH_DEC_MIN = 32.375          # DESI: north catalogue photometry for Dec >= 32.375 in the NGC


def galactic_b(ra, dec) -> np.ndarray:
    ra, dec = np.radians(ra), np.radians(dec)
    v = np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])
    z = _ICRS_TO_GAL[2] @ v
    return np.degrees(np.arcsin(np.clip(z, -1, 1)))


def is_south(ra, dec) -> np.ndarray:
    """DESI north/south split: north = (Dec >= 32.375) & (b > 0)."""
    dec = np.asarray(dec)
    north = dec >= NORTH_DEC_MIN
    if north.any():
        north[north] = galactic_b(np.asarray(ra)[north], dec[north]) > 0
    return ~north


def footprint_columns(cfg: dict) -> list[str]:
    fp = cfg["footprint"]
    cols = ["RA", "DEC", "EBV", "MASKBITS"]
    cols += [f"NOBS_{b}" for b in fp["nobs_min"]]
    cols += [f"GALDEPTH_{b}" for b in fp["galdepth_positive"]]
    return list(dict.fromkeys(cols))


def maskbits_reject(cfg: dict) -> list[str]:
    fp = cfg["footprint"]
    return list(fp["maskbits_reject"]) + list(fp.get("maskbits_reject_extra") or [])


def footprint_cuts(t: Mapping, cfg: dict, fitbits=None) -> dict[str, np.ndarray]:
    """{cut name: pass} for the footprint cuts. ``fitbits`` (galaxies only) exempts SGA centrals
    (FITBITS LARGEGALAXY) from the MASKBITS GALAXY bit when ``galaxy.keep_sga`` is true."""
    fp = cfg["footprint"]
    n = len(t["RA"])
    cuts = {}
    ok = np.ones(n, bool)
    for b, nmin in fp["nobs_min"].items():
        ok &= np.asarray(t[f"NOBS_{b}"]) >= nmin
    cuts["nobs"] = ok
    ok = np.ones(n, bool)
    for b in fp["galdepth_positive"]:
        d = np.asarray(t[f"GALDEPTH_{b}"])
        ok &= np.isfinite(d) & (d > 0)
    cuts["galdepth"] = ok
    names = maskbits_reject(cfg)
    mb = np.asarray(t["MASKBITS"]).astype(np.int64)
    if fitbits is not None and cfg["galaxy"].get("keep_sga") and "GALAXY" in names:
        others = bits.mask_value([x for x in names if x != "GALAXY"], bits.MASKBIT)
        sga = bits.any_set(fitbits, bits.mask_value(["LARGEGALAXY"], bits.FITBIT))
        cuts["maskbits"] = ((mb & others) == 0) & (((mb & (1 << bits.MASKBIT["GALAXY"])) == 0) | sga)
    else:
        cuts["maskbits"] = (mb & bits.mask_value(names, bits.MASKBIT)) == 0
    ebv_max = fp.get("ebv_max")
    cuts["ebv"] = np.ones(n, bool) if ebv_max is None else np.asarray(t["EBV"]) < ebv_max
    cuts["south"] = is_south(t["RA"], t["DEC"]) if fp.get("south_only") else np.ones(n, bool)
    return cuts


def footprint_mask(t: Mapping, cfg: dict) -> np.ndarray:
    """Pass/fail of all footprint cuts (randoms)."""
    return np.logical_and.reduce(list(footprint_cuts(t, cfg).values()))


def galaxy_columns(cfg: dict) -> list[str]:
    g = cfg["galaxy"]
    cols = ["TYPE", "FITBITS", "FLUX_G", "FLUX_R", "FLUX_Z", "MW_TRANSMISSION_G", "MW_TRANSMISSION_R",
            "MW_TRANSMISSION_Z", "FIBERFLUX_R", "FIBERTOTFLUX_R", "GAIA_PHOT_G_MEAN_MAG"]
    cols += [f"FLUX_IVAR_{b}" for b in g.get("flux_ivar_positive") or []]
    for b in (g.get("quality") or {}).get("bands", []):
        cols += [f"FRACMASKED_{b}", f"FRACIN_{b}", f"FRACFLUX_{b}"]
    return list(dict.fromkeys(cols))


def galaxy_cuts(t: Mapping, cfg: dict) -> dict[str, np.ndarray]:
    """{cut name: pass} for the object-level cuts (dereddened magnitudes). A cut whose configuration
    is null or empty is switched off (all pass), e.g. in config/bgs_r21_dr10bits.yaml."""
    g = cfg["galaxy"]
    n = len(t["RA"])
    allpass = np.ones(n, bool)
    gm = mag(t["FLUX_G"], t["MW_TRANSMISSION_G"])
    rm = mag(t["FLUX_R"], t["MW_TRANSMISSION_R"])
    zm = mag(t["FLUX_Z"], t["MW_TRANSMISSION_Z"])
    cuts = {}
    r0, r1 = g["r_range"]
    cuts["rmag"] = (rm > r0) & (rm <= r1)
    rt = g.get("reject_types") or []
    cuts["type"] = ~np.isin(np.char.strip(np.asarray(t["TYPE"]).astype("U")), rt) if rt else allpass
    fb = g.get("fitbits_reject") or []
    cuts["fitbits"] = ~bits.any_set(t["FITBITS"], bits.mask_value(fb, bits.FITBIT)) if fb else allpass
    ok = np.ones(n, bool)
    for b in g.get("flux_ivar_positive") or []:
        ok &= np.asarray(t[f"FLUX_IVAR_{b}"]) > 0
    cuts["flux_ivar"] = ok
    with np.errstate(invalid="ignore"):
        if g.get("gaia_star_cut") is not None:
            gaia_g = np.asarray(t["GAIA_PHOT_G_MEAN_MAG"], np.float64)
            cuts["gaia"] = (gaia_g == 0) | (gaia_g - mag(t["FLUX_R"]) > g["gaia_star_cut"])
        else:
            cuts["gaia"] = allpass
        if g.get("colour"):
            (gr0, gr1), (rz0, rz1) = g["colour"]["gr"], g["colour"]["rz"]
            cuts["colour"] = (gm - rm > gr0) & (gm - rm < gr1) & (rm - zm > rz0) & (rm - zm < rz1)
        else:
            cuts["colour"] = allpass
    q = g.get("quality")
    ok = np.ones(n, bool)
    if q:
        for b in q["bands"]:
            ok &= np.asarray(t[f"FRACMASKED_{b}"]) < q["fracmasked_max"]
            ok &= np.asarray(t[f"FRACIN_{b}"]) > q["fracin_min"]
            ok &= np.asarray(t[f"FRACFLUX_{b}"]) < q["fracflux_max"]
    cuts["quality"] = ok
    with np.errstate(invalid="ignore"):
        fb_cfg = g.get("fiber")
        if fb_cfg:
            rfib = mag(t["FIBERFLUX_R"], t["MW_TRANSMISSION_R"])
            lim = np.where(rm < fb_cfg["r_pivot"], fb_cfg["rfib_max"] + (rm - fb_cfg["r_pivot"]), fb_cfg["rfib_max"])
            cuts["fiber"] = rfib < lim
        else:
            cuts["fiber"] = allpass
        rt_cfg = g.get("rfibtot")
        if rt_cfg:
            rfibtot = mag(t["FIBERTOTFLUX_R"], t["MW_TRANSMISSION_R"])
            cuts["rfibtot"] = ~((rm > rt_cfg["r_min"]) & (rfibtot < rt_cfg["rfibtot_min"]))
        else:
            cuts["rfibtot"] = allpass
    return cuts


STAR_PSF, STAR_GAIA = 1, 2


def star_flag(t: Mapping, gaia_cut: float = 0.6) -> np.ndarray:
    """Bit 0: TYPE == PSF; bit 1: Gaia match with G - r_raw <= gaia_cut (r not dereddened, the DESI /
    DR10 star test). Information only, never a cut."""
    flag = np.zeros(len(t["TYPE"]), np.uint8)
    flag |= np.where(np.char.strip(np.asarray(t["TYPE"]).astype("U")) == "PSF", STAR_PSF, 0).astype(np.uint8)
    g = np.asarray(t["GAIA_PHOT_G_MEAN_MAG"], np.float64)
    with np.errstate(invalid="ignore"):
        star = (g > 0) & (g - mag(t["FLUX_R"]) <= gaia_cut)
    flag |= np.where(star, STAR_GAIA, 0).astype(np.uint8)
    return flag


def cutflow(cuts: Mapping[str, np.ndarray], n_total: int | None = None) -> dict[str, np.ndarray]:
    """Cumulative counts after each cut (in :data:`CUTS` order) and counts failing each cut alone."""
    names, n_cum, n_fail = [], [], []
    if n_total is not None:
        names.append("all"), n_cum.append(n_total), n_fail.append(0)
    ok = None
    for name in CUTS:
        if name not in cuts:
            continue
        ok = cuts[name].copy() if ok is None else ok & cuts[name]
        names.append(name)
        n_cum.append(int(ok.sum()))
        n_fail.append(int((~cuts[name]).sum()))
    return {"CUT": np.array(names), "N_PASS_CUMUL": np.array(n_cum, np.int64),
            "N_FAIL_ALONE": np.array(n_fail, np.int64)}
