"""Stellar-mass / K-correction backends at fixed redshift, with one interface.

Each backend is a class built from the configuration with ``fit(table) -> dict`` returning, per
galaxy (NaN when the fit failed or the code gives no such quantity):

  LOGMSTAR, LOGMSTAR_ERR   log10 stellar mass (Msun, Chabrier IMF) and its 1-sigma uncertainty
  MABS_R, MABS_R_ERR       rest-frame DECam r absolute magnitude (AB, H0 of the config) and error
  LOGMSTAR_LO, LOGMSTAR_HI 16/84% bounds, LOGSFR (Msun/yr), CHI2: benchmark only

PRODUCT lists the columns written to the per-sweep <ver>-<code> files (with SWEEP_ROW). The
uncertainties are at fixed redshift (no photo-z term).

Input table columns: BEST_Z, FLUX_<b>, FLUX_IVAR_<b>, MW_TRANSMISSION_<b> for b in G R I Z W1 W2.
"""

BACKENDS = {
    "kcorrect": "ls11samples.sed.kcorrect_fit:KcorrectBackend",
    "lephare": "ls11samples.sed.lephare_fit:LephareBackend",
    "cigale": "ls11samples.sed.cigale_fit:CigaleBackend",
    "eazy": "ls11samples.sed.eazy_fit:EazyBackend",
    "dsps": "ls11samples.sed.dsps_fit:DspsBackend",
}
PRODUCT = ("LOGMSTAR", "LOGMSTAR_ERR", "MABS_R", "MABS_R_ERR")
OUTPUT = PRODUCT + ("LOGMSTAR_LO", "LOGMSTAR_HI", "LOGSFR", "CHI2")


def get_backend(name: str):
    import importlib

    mod, cls = BACKENDS[name].split(":")
    return getattr(importlib.import_module(mod), cls)
