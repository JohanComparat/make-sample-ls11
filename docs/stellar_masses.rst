Stellar masses and absolute magnitudes
======================================

:mod:`ls11samples.sed` runs five codes behind one interface. Every code fits the same input:

* the same galaxies, at fixed redshift :math:`z=` ``BEST_Z``;
* the same dereddened grizW1W2 photometry, with the same error floor;
* the same filter curves.

Three codes run in production (``sed.codes``): LePhare, CIGALE and kcorrect. LePhare is
``sed.primary``. eazy and DSPS are used in the benchmark only.

Common input
------------

:func:`ls11samples.sed.common.photometry` returns the redshift, fluxes and errors of a table that
has ``BEST_Z`` and, for :math:`b \in` {G, R, I, Z, W1, W2}, ``FLUX_b``, ``FLUX_IVAR_b`` and
``MW_TRANSMISSION_b``:

* fluxes are dereddened, :math:`F/T_\mathrm{MW}`, and inverse variances become
  :math:`\mathrm{ivar}\,T_\mathrm{MW}^2`;
* the error is :math:`\sigma = \sqrt{1/\mathrm{ivar} + \sigma_\mathrm{floor}^2}` with
  :math:`\sigma_\mathrm{floor} = 0.4\ln 10\;\delta m\;|F|` (``photometry.err_floor_mag``: 0.02 mag
  in griz, 0.05 mag in W1, W2);
* a band without a measurement (:math:`\mathrm{ivar} \le 0`) has an infinite error and is left out
  of the fit;
* a galaxy is fitted only when :math:`0.001 < z < 1`.

**Filters.** The filter curves are the kcorrect v5 response curves
(:mod:`ls11samples.sed.filters`): DECam total throughput including the atmosphere, and the WISE
relative response, both per photon. Each code gets them in its own format:

* LePhare: ``$LEPHAREDIR/filt/ls11/*.pb``;
* CIGALE: ``ls11.*`` entries in the pcigale database;
* eazy: entries appended to ``FILTER.RES``;
* DSPS: arrays.

Interface
---------

Each backend is a class built from the configuration, ``get_backend(name)(cfg)``. Its method
``fit(table) -> dict`` returns, per galaxy, the columns below; a value is NaN when the fit failed
or the code gives no such quantity.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - column
     - meaning
   * - ``LOGMSTAR``, ``LOGMSTAR_ERR``
     - :math:`\log_{10} M_\star/M_\odot` (Chabrier IMF) and its 1σ uncertainty
   * - ``MABS_R``, ``MABS_R_ERR``
     - rest-frame DECam r absolute magnitude (AB, :math:`H_0` of the configuration) and its error
   * - ``LOGMSTAR_LO``, ``LOGMSTAR_HI``
     - 16% and 84% bounds (benchmark only)
   * - ``LOGSFR``
     - :math:`\log_{10}` SFR [:math:`M_\odot\,\mathrm{yr}^{-1}`] (benchmark only)
   * - ``CHI2``
     - fit :math:`\chi^2` (benchmark only)

Step 4 writes the first four columns, :data:`ls11samples.sed.PRODUCT`, with ``SWEEP_ROW`` to
``<ver>-<code>/<sweep>-<code>.fits``. These files are row-aligned with the selection file. All
uncertainties are at fixed redshift; the photo-z term is not included.

Production codes
----------------

.. list-table::
   :header-rows: 1
   :widths: 12 26 26 36

   * - code
     - ``LOGMSTAR``
     - ``LOGMSTAR_ERR``
     - ``MABS_R`` ± ``MABS_R_ERR``
   * - LePhare
     - ``MASS_MED`` (``MASS_BEST`` if missing)
     - half the 68% interval, (``MASS_SUP`` − ``MASS_INF``)/2
     - ``MAG_ABS``, ``EMAG_ABS`` in r
   * - CIGALE
     - log of the Bayesian mean
     - Bayesian error, :math:`\sigma_M / (M \ln 10)`
     - rest-frame :math:`L_\nu(r)` and its Bayesian error
   * - kcorrect
     - best fit
     - standard deviation over 20 Monte Carlo flux realisations
     - best fit, and the same realisations

kcorrect
^^^^^^^^

:class:`ls11samples.kcorr.KcorrectV5` wraps kcorrect v5 (Blanton & Roweis 2007), which fits NMF
templates. The backend is :class:`ls11samples.sed.kcorrect_fit.KcorrectBackend`.

* **Mass.** The stellar mass is the surviving mass of the best-fit template mix,
  :math:`\sum_i c_i\,m_{\mathrm{remain},i} \times 10^{0.4\,\mathrm{DM}(z)}`.
* **Absolute magnitude.** :math:`M_r = r - \mathrm{DM}(z) - K_r(z)`, with no band shift.
* **Errors.** kcorrect refits ``kcorr.mc`` = 20 Monte Carlo realisations of the fluxes
  (``fit_coeffs(mc=...)``); the errors are the standard deviations of :math:`\log M_\star` and
  :math:`M_r` over them.

:meth:`~ls11samples.kcorr.KcorrectV5.fit` also returns K-corrections and absolute magnitudes in g,
r and z, and :math:`{}^{0.1}r` (SDSS r shifted to :math:`z = 0.1`). Its
:meth:`~ls11samples.kcorr.KcorrectV5.k_curve` evaluates the r-band K-correction of a set of
fitted SEDs on a redshift grid; step 5 uses it for the :math:`M_r` completeness limit
(:doc:`volume_limited`).

LePhare
^^^^^^^

:class:`ls11samples.sed.lephare_fit.LephareBackend` runs LePhare ≥ 1.0 with
``config/lephare/LS11_zFIX.para``, which reconstructs the DR10 run:

* the BC03 Chabrier composite library;
* Calzetti attenuation with :math:`E(B-V) \le 0.7`;
* physical emission lines;
* the redshift fixed to ``BEST_Z``.

On first use, the backend writes the filters, downloads the auxiliary files the configuration
needs and builds the libraries in ``$LEPHAREWORK``. It rebuilds them if the ``.doc`` file of the
library is missing or empty, which is the sign of an interrupted build.

CIGALE
^^^^^^

:class:`ls11samples.sed.cigale_fit.CigaleBackend` runs ``pcigale`` ≥ 2025.0 through its command
line. pcigale sets ``OMP_NUM_THREADS=1`` when imported, so it is never imported in the main
process.

* **Model grid** (:data:`ls11samples.sed.cigale_fit.DEFAULT_GRID`, overridable in
  ``sed.cigale.grid``): delayed-τ SFH, BC03 Chabrier, nebular emission, and a modified-starburst
  (Calzetti-like) attenuation.
* **Redshift steps.** Redshifts are rounded to ``redshift_decimals`` (default 2), so the models are
  computed once per redshift step.
* **Absolute magnitude.** ``MABS_R`` comes from the Bayesian rest-frame :math:`L_\nu` in DECam r:
  :math:`-2.5\log_{10}(F_{10\,\mathrm{pc}} / 3631\,\mathrm{Jy})`.
* **Errors.** pcigale's ``additionalerror`` is set to 0. Its default would add 10% of the flux in
  quadrature to every band, while the input errors already carry the common error floor.
  pcigale floors the errors of extensive properties at 5% of the value, which cannot be
  configured. As a result, ``LOGMSTAR_ERR`` ≥ 0.0217 dex and ``MABS_R_ERR`` ≥ 0.054 mag.
* **Working directory.** Each fit runs in a temporary directory under ``$TMPDIR``, which is
  node-local on the cluster, and the directory is removed afterwards.

Benchmark-only codes
--------------------

**eazy** (:class:`ls11samples.sed.eazy_fit.EazyBackend`):

* ``fit_at_zbest`` with the FSPS-based ``corr_sfhz_13`` templates;
* masses and SFRs from the template mass-to-light ratios (``sps_parameters``);
* :math:`M_r` from the rest-frame r flux of the fit;
* no uncertainties.

**DSPS** (:class:`ls11samples.sed.dsps_fit.DspsBackend`). DSPS provides SSP spectra and photometry
kernels but no fitter, so the backend builds one:

* **SSPs:** the FSPS v3.2 SSPs (Kroupa IMF, nebular emission).
* **Models:** a delayed-τ × metallicity × Calzetti (2000) :math:`A_V` grid.
* **Photometry:** observed AB magnitudes tabulated with jax on a redshift grid.
* **Fit:** the mass amplitude of each model is solved analytically. Models older than the Universe
  at :math:`z` are excluded. :math:`\log M_\star` percentiles are weighted by
  :math:`e^{-\chi^2/2}`.
* **Surviving mass:** the formed mass times :math:`1 - 0.05\ln(1 + t/1.4\,\mathrm{Myr})`
  (Behroozi et al. 2013).
* **IMF:** converted from Kroupa to Chabrier by −0.034 dex.

Benchmark
---------

``benchmarks/sed_benchmark.py`` fits 10 000 galaxies of the local DR11 strip with the five codes,
half of them with spectroscopic redshifts (:doc:`benchmark`). In summary:

* every code agrees with LePhare within ±0.05 dex in the median, with an NMAD of 0.09–0.10 dex;
* the photo-z error alone moves the masses by 0.07–0.11 dex;
* DR11 LePhare differs from the DR10 LePhare masses by +0.03 dex.

The CIGALE results of the benchmark predate the change of ``additionalerror`` to 0. That change
raises the CIGALE masses by about 0.05 dex.

Cost, in CPU hours for 20 million galaxies: kcorrect ≈ 0.1, DSPS ≈ 22, eazy ≈ 225,
CIGALE ≈ 520, LePhare ≈ 1500.
