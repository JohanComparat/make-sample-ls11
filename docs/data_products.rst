Data products
=============

All products are FITS binary tables. Paths are given relative to ``$LS11_SWEEP_OUT`` (per-sweep
products) or ``$LS11_OUT/<tag>`` (per-run products); see :doc:`environment`. Types are FITS/numpy
codes: ``i8`` is a 64-bit integer, ``f4`` a 32-bit float, and so on.

Selection (step 1)
------------------

``<ver>-<tag>/<sweep>-<tag>.fits``: one file per sweep, holding the selected objects only.

HDU 1, ``SELECTION``:

.. list-table::
   :header-rows: 1
   :widths: 20 8 72

   * - column
     - type
     - content
   * - ``LS_ID_DR11``
     - i8
     - ``RELEASE << 42 | BRICKID << 22 | OBJID``
   * - ``SWEEP_ROW``
     - i4
     - row in the sweep and in the row-matched photo-z file
   * - ``BEST_Z``
     - f4
     - redshift used by every later step (:ref:`best-z`)
   * - ``BEST_Z_ERR``
     - f4
     - 0 for a spectroscopic redshift, else the photo-z standard deviation
   * - ``Z_SOURCE``
     - i1
     - 0 none, 1 ``Z_SPEC``, 2… photo-z columns in ``redshift.photoz_priority`` order
   * - ``STAR_FLAG``
     - u1
     - bit 0 ``TYPE`` = PSF; bit 1 Gaia star test (information, never a cut)

The header keywords are:

``SWEEP``, ``TAG``
   the sweep file name and the selection tag;
``NSWEEP``, ``NMAGCUT``, ``NSEL``
   the numbers of rows in the sweep, inside ``r_range`` and selected;
``RMIN``, ``RMAX``
   the r-band range;
``CFGHASH``, ``FPHASH``, ``CFGFILE``
   the configuration hashes and file name;
``MASKREJ``, ``FITREJ``
   the rejected bits, as ``NAME(number)``;
``STARFLAG``
   the meaning of ``STAR_FLAG``;
``ZSRC<n>``
   the meaning of each ``Z_SOURCE`` code.

HDU 2, ``CUTFLOW``: one row per cut, with ``CUT``, ``N_PASS_CUMUL`` and ``N_FAIL_ALONE``
(see :doc:`selection`).

Stellar masses (step 4)
-----------------------

``<ver>-<code>/<sweep>-<code>.fits``: one file per sweep and code, row-aligned with the selection
file. The HDU is named after the code (``LEPHARE``, ``CIGALE``, ``KCORRECT``).

.. list-table::
   :header-rows: 1
   :widths: 20 8 72

   * - column
     - type
     - content
   * - ``SWEEP_ROW``
     - i4
     - as in the selection file (checked when joined)
   * - ``LOGMSTAR``
     - f4
     - :math:`\log_{10} M_\star/M_\odot`, Chabrier IMF
   * - ``LOGMSTAR_ERR``
     - f4
     - 1σ, at fixed redshift
   * - ``MABS_R``
     - f4
     - rest-frame DECam r absolute magnitude, AB
   * - ``MABS_R_ERR``
     - f4
     - 1σ, at fixed redshift

The header keywords are:

``CODE``, ``TAG``
   the SED code and the selection tag;
``SELHASH``
   the ``CFGHASH`` of the selection file;
``CFGHASH``
   the hash of the configuration;
``MABSBAND``, ``MASSIMF``, ``ERRDEF``, ``COSMO``
   the definitions of the absolute magnitude, IMF, errors and cosmology.

Randoms (step 2)
----------------

``LS11_<tag>_RAND.fits``, HDU ``RANDOMS``: the footprint randoms, with columns ``RA``, ``DEC``
(f8) and ``EBV`` (f4).

The header keywords are:

``NRFILES``, ``RFILE<i>``
   the number and names of the random files read;
``DENSITY``
   the summed density of those files, in deg\ :sup:`-2`;
``NSWEEPS``, ``SWPHASH``
   the number of processed sweeps whose boxes the randoms cover, and a hash of their names;
``NINBOX``, ``NRAND``
   the numbers of randoms inside the sweep boxes and inside the footprint;
``AREA``
   the footprint area, in deg\ :sup:`2`;
``BOXAREA``
   the area of the sweep boxes, in deg\ :sup:`2`;
``TAG``, ``FPHASH``, ``CFGHASH``, ``MASKREJ``
   the selection tag, configuration hashes and rejected MASKBITS.

HEALPix maps (step 2)
---------------------

All maps follow the conventions of the DR10 systematics maps read by ``sys_mapping``:

* RING ordering and equatorial coordinates (``COORDSYS='C'``);
* one column, named after the quantity;
* ``UNSEEN`` outside the footprint.

.. list-table::
   :header-rows: 1
   :widths: 55 45

   * - file
     - content
   * - ``footprint/LS11_FRACAREA_NSIDE_<nside:04d>.fits``
     - fraction of each pixel inside the footprint (column ``FRACAREA``)
   * - ``systematics/<nside:04d>/LS11_<Q>_NSIDE_<nside:04d>.fits``
     - mean of ``Q`` over the footprint randoms of each pixel; depths as 5σ AB magnitudes
   * - ``systematics/<nside:04d>/GAIA_<name>_NSIDE_<nside:05d>.fits``
     - Gaia star-density map copied from ``$LS11_GAIA_MAPS``

``Q`` runs over ``maps.quantities``: ``GALDEPTH_{G,R,I,Z}``, ``PSFDEPTH_W1``,
``PSFSIZE_{G,R,I,Z}``, ``NOBS_{G,R,I,Z}`` and ``EBV``. ``nside`` runs over ``maps.nsides``: 32,
64, 128 and 256. ``sys_mapping`` reads every ``*.fits`` file of a ``systematics/<nside>``
directory as a template.

Volume-limited samples (step 5)
-------------------------------

The files sit in ``vlim/`` (or ``vlim_<code>/``). Each sample ``<name>``
(:ref:`naming <sample-names>`) has three files.

``<name>_DATA.fits``, HDU ``DATA``:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - column
     - content
   * - ``RA``, ``DEC``
     - f8
   * - ``EBV``
     - from the sweep
   * - ``BEST_Z``, ``BEST_Z_ERR``, ``Z_SOURCE``, ``STAR_FLAG``, ``LS_ID_DR11``
     - from the selection file
   * - ``MAG_G``, ``MAG_R``, ``MAG_Z``
     - dereddened DECam AB magnitudes
   * - ``MABS_R``, ``KCORR_R``
     - from ``vlim.mr_code``; :math:`K = r - \mathrm{DM} - M_r`
   * - ``LOGMSTAR``
     - from ``sed.primary``
   * - ``LOGMSTAR_<CODE>``, ``LOGMSTAR_ERR_<CODE>``, ``MABS_R_<CODE>``, ``MABS_R_ERR_<CODE>``
     - one set per code run
   * - ``LPH_MASS_BEST``
     - f8, equal to ``LOGMSTAR``; the column ``sum_stat`` reads
   * - ``WEIGHT_COMP``
     - 1; ``sys_mapping`` uses it when present

``<name>_RAND.fits``, HDU ``RAND``: ``RA``, ``DEC`` (f8), ``EBV`` and ``Z``. ``Z`` holds the
shuffled data redshifts.

``<name>_COLOUR.fits``, HDU ``COLOUR``: ``G_MAG``, ``Z_MAG`` (dereddened DECam AB) and
``REDSHIFT`` (equal to ``BEST_Z``). It is row-aligned with ``DATA`` and gives the ``sum_stat``
colour classes. Its header has ``PARENT``, ``NMATCH`` and ``NMISS``.

The headers of ``DATA`` and ``RAND`` carry:

``SAMPLE``, ``KIND``, ``LO``, ``HI``, ``ZMIN``, ``ZMAX``, ``NGAL``
   the sample definition, as used in the selection and in the name;
``AREA``, ``RLIM``, ``TAG``
   the footprint area, the r-band limit and the selection tag;
``SEDCODE``, ``SEDCODES``, ``MRCODE``
   the code giving ``LOGMSTAR``, all codes run, and the code giving ``MABS_R``;
``COSMO``, ``CFGHASH``, ``FPHASH``, ``PARENT``
   the cosmology, configuration hashes and parent selection.

``RAND`` also has ``NRAND``.

``vlim_boundaries.fits`` has three HDUs:

* ``SAMPLES``: one row per sample, with ``NAME``, ``KIND``, ``LO``, ``HI``, ``Z_MIN``, ``Z_MAX``,
  ``N_GAL``, ``VOLUME`` (comoving, Mpc\ :sup:`3`), ``N_DENS`` (Mpc\ :sup:`-3`), ``N_DEG2``,
  ``Z_MEDIAN``, ``MABS_R_MEDIAN`` and ``LOGMSTAR_MEDIAN``;
* ``LIMIT_MR`` and ``LIMIT_MSTAR``: the completeness curves, with columns ``Z`` and ``LIMIT``.

Manifest (steps 5-6)
--------------------

``vlim/manifest.yaml`` records:

* the creation time, package version and git commit;
* the configuration file and its hashes;
* the area, tag, SED codes and random file;
* one entry per sample, as in ``SAMPLES``.

Step 6 adds three sections:

``sys_mapping``
   the catalogue directory, the template directory per nside and the area-fraction map per nside;
``sum_stat``
   the survey type and the data, random and colour columns;
``checks``
   the number of samples checked and the list of problems found.
