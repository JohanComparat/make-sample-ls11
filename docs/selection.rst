Selection
=========

The selection (:mod:`ls11samples.selection`, applied per sweep by :mod:`ls11samples.bgsl`) has
two kinds of cuts:

* **footprint cuts** use only columns the randoms also have. They are applied by the same
  function to galaxies and randoms, so both share one footprint;
* **galaxy cuts** are object-level cuts that the randoms cannot carry.

Magnitudes are dereddened, :math:`m = 22.5 - 2.5\log_{10}(F/T_\mathrm{MW})`, except where stated.
The values below are those of ``config/default.yaml`` (tag ``bgsl``).

.. _footprint-cuts:

Footprint cuts
--------------

:func:`ls11samples.selection.footprint_cuts` returns one pass/fail array per cut.

.. list-table::
   :header-rows: 1
   :widths: 14 50 36

   * - cut
     - condition to pass
     - configuration
   * - ``nobs``
     - ``NOBS_b`` ≥ 1 in g, r, z
     - ``footprint.nobs_min``
   * - ``galdepth``
     - ``GALDEPTH_b`` finite and > 0 in g, r, z
     - ``footprint.galdepth_positive``
   * - ``maskbits``
     - none of NPRIMARY, BRIGHT, MEDIUM, GALAXY, CLUSTER set in ``MASKBITS``
     - ``footprint.maskbits_reject`` (+ ``maskbits_reject_extra``)
   * - ``ebv``
     - ``EBV`` < 0.1
     - ``footprint.ebv_max``
   * - ``south``
     - not in the DESI north: north means Dec ≥ 32.375 and Galactic latitude :math:`b > 0`
     - ``footprint.south_only``

Notes on these cuts:

* **MEDIUM.** It is masked because, in DR11, 97% of the objects with FITBITS FIT_BACKGROUND (which
  the galaxy cuts reject) lie in MEDIUM-star masks. Masking MEDIUM removes that area from the
  randoms too. Without it, the galaxies around medium-bright stars would be missing while the
  randoms would not.
* **SGA centrals.** When ``galaxy.keep_sga`` is true, an object with FITBITS LARGEGALAXY is exempt
  from the MASKBITS GALAXY bit. The randoms have no FITBITS, so the GALAXY mask still applies to
  them. The default is false.
* **North/south split.** It replaces the ad hoc DR10 cut (100 < RA < 180, Dec < 32). It has no
  effect in the south Galactic cap. :func:`ls11samples.selection.galactic_b` computes :math:`b`
  with the ICRS-to-Galactic rotation matrix.
* **Stricter mask.** ``maskbits_reject_extra`` adds bits to the mask. The stricter mask used by
  ``rema`` is given as a comment in ``config/default.yaml``.

.. _galaxy-cuts:

Galaxy cuts
-----------

:func:`ls11samples.selection.galaxy_cuts` returns one pass/fail array per cut. A cut whose
configuration is ``null`` or empty passes everything.

.. list-table::
   :header-rows: 1
   :widths: 14 56 30

   * - cut
     - condition to pass
     - configuration
   * - ``rmag``
     - :math:`13 < r \le 19.5`
     - ``galaxy.r_range``
   * - ``type``
     - ``TYPE`` ≠ DUP
     - ``galaxy.reject_types``
   * - ``fitbits``
     - none of FIT_BACKGROUND, MEDIUM, WALKER, RUNNER, GAIA_POINTSOURCE, ITERATIVE set in
       ``FITBITS`` (the DR10 bits 1, 6, 10, 11, 12, 13)
     - ``galaxy.fitbits_reject``
   * - ``flux_ivar``
     - ``FLUX_IVAR_b`` > 0 in g, r, z
     - ``galaxy.flux_ivar_positive``
   * - ``gaia``
     - no Gaia match (``GAIA_PHOT_G_MEAN_MAG`` = 0), or :math:`G - r_\mathrm{raw} > 0.6`, with
       :math:`r_\mathrm{raw}` **not** dereddened
     - ``galaxy.gaia_star_cut``
   * - ``colour``
     - :math:`-1 < g-r < 4` and :math:`-1 < r-z < 4` (DESI main BGS)
     - ``galaxy.colour``
   * - ``quality``
     - in g, r, z: ``FRACMASKED`` < 0.4, ``FRACIN`` > 0.3, ``FRACFLUX`` < 5 (DESI SV BGS)
     - ``galaxy.quality``
   * - ``fiber``
     - :math:`r_\mathrm{fib} < 22.9 + (r - 17.8)` if :math:`r < 17.8`, else
       :math:`r_\mathrm{fib} < 22.9`
     - ``galaxy.fiber``
   * - ``rfibtot``
     - not (:math:`r > 12` and :math:`r_\mathrm{fibtot} < 15`), which removes bright-star
       fragments
     - ``galaxy.rfibtot``

The selected sample is the set of objects passing every footprint and galaxy cut. The cut-flow
(HDU ``CUTFLOW`` of the selection file, :func:`ls11samples.selection.cutflow`) counts the objects
cut by cut, in the order of :data:`ls11samples.selection.CUTS`, the one of the two tables above.
Only the selected objects are written to the selection file.

For speed, step 1 first reads ``FLUX_R`` and ``MW_TRANSMISSION_R`` of every row and keeps the
``r_range`` rows. Only for those rows does it read the other columns.

Bit definitions
---------------

:data:`ls11samples.bits.MASKBITS` and :data:`ls11samples.bits.FITBITS` hold the DR11 bit
definitions. The configuration names bits and never gives their numbers.
:func:`ls11samples.bits.check_header` compares the ``MBIT_n`` / ``FBIT_n`` keywords of every sweep
header with these tables and raises on any difference. The randoms carry the same MASKBITS; their
header has no bit list, so the sweep check covers them. The header of the selection file records
the rejected bits with their numbers (``MASKREJ``, ``FITREJ``).

Cut flow
--------

Each selection file has a ``CUTFLOW`` HDU (:func:`ls11samples.selection.cutflow`) with one row per
cut, in the order of :data:`~ls11samples.selection.CUTS`. Its first row, ``all``, counts the rows
inside ``r_range``. The columns are:

``N_PASS_CUMUL``
   the number of objects passing this cut and all the cuts before it;
``N_FAIL_ALONE``
   the number of objects failing this cut, whatever the others.

Summing the ``CUTFLOW`` tables of all sweeps gives the cut flow of the whole footprint.

Star flag
---------

``STAR_FLAG`` (:func:`ls11samples.selection.star_flag`) is information only and never a cut:

* bit 0 (value 1): ``TYPE`` = PSF;
* bit 1 (value 2): a Gaia match with :math:`G - r_\mathrm{raw} \le 0.6`, the DESI / DR10 star
  test.

``bgsl`` already rejects Gaia stars through the ``gaia`` cut, so only bit 0 can be set there. In
``bgsr21`` both bits can be set. Step 5 leaves out the objects with
``STAR_FLAG & vlim.exclude_star_flag``:

* ``bgsl`` sets 2 (Gaia stars);
* ``bgsr21`` sets 3 (Gaia stars and PSF objects).

.. _best-z:

Redshift
--------

The DR11 photo-z files (``<ver>-photo-z/<sweep>-pz.fits``) are row-matched to the sweeps.
:func:`ls11samples.io.read_pz` checks them:

* both files have the same number of rows;
* for every selected object, the photo-z ``RELEASE``, ``BRICKID`` and ``OBJID``, and the
  ``LS_ID_DR11`` column itself, match the sweep's
  :math:`\mathtt{LS\_ID\_DR11} = \mathtt{RELEASE} \ll 42 \;|\; \mathtt{BRICKID} \ll 22 \;|\; \mathtt{OBJID}`.

Any mismatch raises. :func:`ls11samples.redshift.best_z` then sets one redshift per object.
Every later step uses it.

.. list-table::
   :header-rows: 1
   :widths: 12 38 50

   * - ``Z_SOURCE``
     - ``BEST_Z``
     - condition
   * - 1
     - ``Z_SPEC``
     - finite, :math:`0.001 < z < 3`, ``SURVEY`` not COSMOS2015
   * - 2
     - ``Z_PHOT_MEAN_I``
     - no valid ``Z_SPEC``; photo-z finite and > −1
   * - 3
     - ``Z_PHOT_MEAN``
     - neither of the above
   * - 0
     - NaN
     - no valid redshift

The photo-z order is ``redshift.photoz_priority``: code 2 is its first entry, code 3 the second,
and so on. ``BEST_Z_ERR`` is 0 for a spectroscopic redshift and the matching ``Z_PHOT_STD*``
column otherwise. The COSMOS2015 ``Z_SPEC`` values are photometric reference redshifts, which is
why they are excluded. The header keywords ``ZSRC<n>`` of the selection file record the meaning of
each code.
