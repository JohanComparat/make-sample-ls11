Configuration
=============

Every selection choice is set in one YAML file, read by :func:`ls11samples.config.load_config`.
The file is the one given by ``--config`` or, failing that, by ``$LS11_CONFIG``; the default is
``config/default.yaml``.

Inheritance
-----------

A configuration may start with ``base: <file>``, a path relative to its own directory. The base
file is read first, then this file's keys are deep-merged on top
(:func:`ls11samples.config.deep_merge`):

* mappings are merged key by key;
* everything else (lists, scalars) replaces the base value;
* ``null`` therefore switches a cut off.

``tag`` names the selection. Its per-sweep files go to ``<ver>-<tag>/`` and its per-run products to
``$LS11_OUT/<tag>/``. Two configurations with different tags can share the same sweeps.

.. _hashes:

Configuration hashes
--------------------

:func:`ls11samples.config.config_hash` is the first 12 hexadecimal digits of the SHA-1 of the
configuration, dumped as YAML with sorted keys (``SWPHASH`` hashes sweep names instead). It is computed either over the whole configuration
or over one section. The hashes are written to the output headers and checked before files are
mixed:

.. list-table::
   :header-rows: 1
   :widths: 15 35 50

   * - keyword
     - hash of
     - checked by
   * - ``CFGHASH``
     - the whole configuration
     - step 1: an existing selection file made with another configuration raises, unless
       ``--overwrite``
   * - ``FPHASH``
     - the ``footprint`` section
     - step 2: every selection file must have the footprint of the randoms
   * - ``SELHASH``
     - ``CFGHASH`` of the selection file a mass file was fitted from
     - step 4: a mass file made for another selection raises, unless ``--overwrite``
   * - ``SWPHASH``
     - the sorted names of the sweeps a random file covers
     - step 2: a random file made from other sweeps is redone; step 5: it raises

Sections
--------

``cosmology``
   ``H0`` and ``Om0`` of the flat ΛCDM cosmology (:mod:`ls11samples.cosmo`). It is used for the
   distance moduli, volumes, kcorrect and DSPS; the LePhare ``.para`` file repeats it.

``photometry``
   ``zeropoint`` (22.5, nanomaggies to AB), ``bands`` and ``err_floor_mag``. The error floor is
   added in quadrature to the flux error of every band before an SED fit:
   :math:`\sigma_\mathrm{floor} = 0.4\ln 10\;\delta m\;|F|`, with :math:`\delta m` = 0.02 mag in
   griz and 0.05 mag in W1, W2.

``footprint``
   These cuts apply to galaxies **and** randoms (:ref:`footprint-cuts`):

   * ``nobs_min`` (per band);
   * ``galdepth_positive`` (bands);
   * ``maskbits_reject`` and ``maskbits_reject_extra`` (bit names);
   * ``ebv_max``;
   * ``south_only``.

``galaxy``
   These object-level cuts apply to galaxies only (:ref:`galaxy-cuts`):

   * ``r_range``: dereddened r, :math:`(r_\mathrm{min}, r_\mathrm{max}]`;
   * ``reject_types``;
   * ``fitbits_reject`` and ``keep_sga``;
   * ``flux_ivar_positive``;
   * ``gaia_star_cut``;
   * ``colour``, ``quality``, ``fiber`` and ``rfibtot``.

   Setting a cut to ``null`` (or to an empty list) switches it off.

``redshift``
   How ``BEST_Z`` is chosen (:ref:`best-z`):

   * ``zspec_range``: the open interval in which a spectroscopic redshift is valid;
   * ``zspec_exclude_surveys``: the ``SURVEY`` values whose ``Z_SPEC`` is not used;
   * ``photoz_priority``: the photo-z columns, tried in order.

``kcorr``
   kcorrect v5 settings (:class:`ls11samples.kcorr.KcorrectV5`):

   * ``bands_in``: the fitted bands;
   * ``band_shift``: 0 gives rest-frame DECam bands;
   * ``z_range`` and ``nz``: the redshift grid;
   * ``mc``: Monte Carlo realisations for the errors;
   * ``complete_set_zmax``: maximum redshift of the complete SED set behind the :math:`M_r`
     completeness curve.

``sed``
   ``codes``: the codes step 4 runs. ``primary``: the code giving ``LOGMSTAR`` and the
   :math:`M_\star` volume limits. Two optional subsections override the model grids:

   * ``sed.cigale``: ``grid`` (per pcigale module) and ``redshift_decimals``;
   * ``sed.dsps``: ``grid``.

   See :doc:`stellar_masses`.

``randoms``
   ``density_per_file`` (deg\ :sup:`-2`) is checked against the ``DENSITY`` keyword of each random
   file. ``chunk_rows`` sets how many rows are read at a time.

``maps``
   * ``nsides`` of the HEALPix maps;
   * ``quantities``: the random columns averaged into ``sys_mapping`` templates;
   * ``depth_as_mag``: write depths as 5σ AB magnitudes;
   * ``gaia_templates``: Gaia maps copied from ``$LS11_GAIA_MAPS``.

``vlim``
   * ``n_rand_factor``: randoms per galaxy;
   * ``mr_code``: the code giving ``MABS_R``;
   * ``exclude_star_flag``: a ``STAR_FLAG`` mask of objects left out;
   * ``completeness_percentile``;
   * ``mstar`` and ``absmag_r``: the thresholds, the faint or heavy bound, the minimum redshift
     ``z_min`` and the redshift cap ``z_cap``.

   See :doc:`volume_limited`.

Provided configurations
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 32 10 58

   * - file
     - tag
     - selection
   * - ``config/default.yaml``
     - ``bgsl``
     - DESI BGS Bright-like, :math:`13 < r \le 19.5`
   * - ``config/bgs_r21_dr10bits.yaml``
     - ``bgsr21``
     - every object with :math:`13 < r \le 21` in the footprint, with the DR10 FITBITS cut and no
       other galaxy cut

``bgsr21`` inherits ``default.yaml`` and switches the other galaxy cuts off with ``null``:

* **No star/galaxy separation.** ``STAR_FLAG`` records it instead, and step 5 leaves out both PSF
  and Gaia stars (``exclude_star_flag: 3``).
* **Footprint.** Its MASKBITS are the DR10 ones (NPRIMARY, BRIGHT, CLUSTER) plus MEDIUM and GALAXY,
  for galaxies and randoms alike. The DR10 FITBITS cut leaves the MEDIUM and GALAXY areas almost
  empty of galaxies (3.3% of the footprint). Masking them removes the same area from the randoms,
  where it would otherwise show up as a density deficit.
* **Randoms.** It uses 5 randoms per galaxy, as for DR10.

config/default.yaml
^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../config/default.yaml
   :language: yaml

config/bgs_r21_dr10bits.yaml
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../config/bgs_r21_dr10bits.yaml
   :language: yaml

LePhare configuration
---------------------

``config/lephare/LS11_zFIX.para`` reconstructs the DR10 LePhare run:

* the BC03 Chabrier composite library;
* Calzetti attenuation with :math:`E(B-V) \le 0.7`;
* physical emission lines;
* the redshift fixed to ``BEST_Z``.

``config/lephare/output_ls11.para`` lists the output columns. See :doc:`stellar_masses`.
