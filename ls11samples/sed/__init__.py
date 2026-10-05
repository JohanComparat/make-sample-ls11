"""Stellar-mass / K-correction backends at fixed redshift, with one interface.

Each backend is a class built from the configuration with ``fit(table) -> dict`` returning, per
galaxy (NaN when the fit failed):

  LOGMSTAR, LOGMSTAR_LO, LOGMSTAR_HI   log10 stellar mass (Msun, Chabrier IMF) and 16/84% bounds
                                       (LO = HI = LOGMSTAR for best-fit-only codes)
  LOGSFR                               log10 SFR (Msun/yr), NaN when the code gives none
  ABSMAG_R                             rest-frame DECam r absolute magnitude (AB)
  CHI2                                 best-fit chi^2

Input table columns: BEST_Z, FLUX_<b>, FLUX_IVAR_<b>, MW_TRANSMISSION_<b> for b in G R I Z W1 W2.
"""

BACKENDS = {
    "kcorrect": "ls11samples.sed.kcorrect_fit:KcorrectBackend",
    "lephare": "ls11samples.sed.lephare_fit:LephareBackend",
    "cigale": "ls11samples.sed.cigale_fit:CigaleBackend",
    "eazy": "ls11samples.sed.eazy_fit:EazyBackend",
    "dsps": "ls11samples.sed.dsps_fit:DspsBackend",
}
OUTPUT = ("LOGMSTAR", "LOGMSTAR_LO", "LOGMSTAR_HI", "LOGSFR", "ABSMAG_R", "CHI2")


def get_backend(name: str):
    import importlib

    mod, cls = BACKENDS[name].split(":")
    return getattr(importlib.import_module(mod), cls)
