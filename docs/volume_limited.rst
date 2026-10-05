Volume-limited samples
======================

Step 5 builds volume-limited samples in absolute magnitude :math:`M_r` and in stellar mass
:math:`M_\star` (:mod:`ls11samples.vlim`). Both completeness limits are measured from the data;
there is no external calibration. :math:`r_\mathrm{lim}` is the faint end of ``galaxy.r_range``
(19.5 for ``bgsl``, 21 for ``bgsr21``). DM is the distance modulus.

Input catalogue
---------------

:func:`ls11samples.catalog.load` joins three sources, sweep by sweep:

* the selection files;
* the code files of every code run, as ``LOGMSTAR_<CODE>``, ``LOGMSTAR_ERR_<CODE>``,
  ``MABS_R_<CODE>`` and ``MABS_R_ERR_<CODE>``;
* ``RA``, ``DEC``, ``EBV`` and the g, r, z magnitudes, read from the sweep at ``SWEEP_ROW``.

It checks that every code file is row-aligned with the selection. Then:

* objects with ``STAR_FLAG & vlim.exclude_star_flag`` are left out;
* ``LOGMSTAR`` is taken from ``sed.primary`` and ``MABS_R`` from ``vlim.mr_code`` (kcorrect);
* ``KCORR_R`` is defined as :math:`r - \mathrm{DM}(z) - M_r`.

Absolute-magnitude limit
------------------------

Every galaxy brighter than :math:`M_{r,\mathrm{lim}}(z)` is above the flux limit at :math:`z`:

.. math::

   M_{r,\mathrm{lim}}(z) = r_\mathrm{lim} - \mathrm{DM}(z) - K_p(z),

where :math:`K_p(z)` is the :math:`p`-th percentile (``vlim.completeness_percentile``, 95) of the
r-band K-correction: the K-correction of the reddest galaxies.

Measured on the observed galaxies, :math:`K_p` would be biased low wherever the flux limit removes
red galaxies (:math:`z \gtrsim 0.35` for :math:`r < 19.5`). It is therefore measured on a
complete set of SEDs at low redshift instead:

1. **Complete set** (:func:`ls11samples.kcorr.complete_set_rows`). It takes the galaxies with
   :math:`0.01 < z \le z_c` (``kcorr.complete_set_zmax``, 0.1) and brighter than the limit of the
   reddest galaxies at :math:`z_c`:
   :math:`M_r \le r_\mathrm{lim} - \mathrm{DM}(z_c) - K_{99}`. Here :math:`K_{99}` is the 99th
   percentile of ``KCORR_R`` in that redshift range. At most 20 000 galaxies are kept.
2. **Refit.** These galaxies are refitted with kcorrect at their own redshift, without Monte Carlo.
3. **K-correction curve.** The r-band K-correction of every fitted SED is evaluated at each
   :math:`z` of the grid 0.01, 0.02, …, 0.60 (:meth:`ls11samples.kcorr.KcorrectV5.k_curve`), and
   :math:`K_p(z)` is its :math:`p`-th percentile.

If kcorrect cannot be imported, step 5 falls back on the binned percentiles of the observed
``KCORR_R`` (:func:`ls11samples.vlim.mr_limit_curve`) and logs a warning.

Stellar-mass limit
------------------

The stellar-mass limit follows Pozzetti et al. (2010) (:func:`ls11samples.vlim.mstar_limit_curve`).
Each galaxy's mass is rescaled to the flux limit:

.. math::

   \log M_{\mathrm{lim},i} = \log M_{\star,i} + 0.4\,(r_i - r_\mathrm{lim}).

In each redshift bin of width 0.02 (with at least 20 galaxies), :math:`\log M_{\star,\mathrm{lim}}(z)`
is the :math:`p`-th percentile of :math:`\log M_{\mathrm{lim},i}` over the faintest 20% of the
galaxies.

Sample boundaries
-----------------

Both curves are made monotonic: :math:`M_{r,\mathrm{lim}}` can only brighten with :math:`z`, and
:math:`M_{\star,\mathrm{lim}}` can only grow. For each threshold, the maximum redshift is set as
follows:

1. :math:`z_\mathrm{max}` is where the curve crosses the threshold, interpolated linearly between
   the last complete bin and the first incomplete one
   (:func:`~ls11samples.vlim.zmax_brighter`, :func:`~ls11samples.vlim.zmax_heavier`);
2. it is capped at ``z_cap`` (0.35);
3. it is rounded down to 2 decimals, the precision of the sample name; ``z_min`` is rounded up
   to 2 decimals (:func:`~ls11samples.vlim.name_precision`).

The redshift range in a sample name is therefore exactly the range selected. ``sys_mapping``
reads it from the name. A threshold whose :math:`z_\mathrm{max}` is undefined, or not above
:math:`z_\mathrm{min}`, gives no sample.

.. list-table::
   :header-rows: 1
   :widths: 14 42 44

   * - kind
     - selection
     - thresholds (``config/default.yaml``)
   * - ``Mr``
     - :math:`-25 < M_r \le M_\mathrm{thr}`, :math:`0.02 < z \le z_\mathrm{max}`
     - :math:`M_\mathrm{thr}` = −18.0, −18.5, …, −22.5 (``vlim.absmag_r``)
   * - ``Mstar``
     - :math:`\log M_\mathrm{thr} \le \log M_\star < 12`, :math:`0.05 < z \le z_\mathrm{max}`
     - :math:`\log M_\mathrm{thr}` = 9.0, 9.5, 10.0, 10.25, …, 11.5 (``vlim.mstar``)

The samples are threshold samples: each holds every galaxy brighter, or more massive, than its
threshold, over its own redshift range. The samples therefore overlap.

Randoms
-------

The randoms of a sample (:func:`ls11samples.vlim.make_randoms`) are built in two steps:

1. **Positions.** A random subsample of the footprint randoms of step 2, of size
   ``vlim.n_rand_factor`` × :math:`N_\mathrm{gal}` (20 for ``bgsl``, 5 for ``bgsr21``). If fewer
   randoms are available, all of them are used and a warning is logged.
2. **Redshifts.** Redshifts drawn with replacement from the sample's own ``BEST_Z`` (the shuffle
   method).

The draws use ``--seed`` (default 1), so the samples are reproducible.

.. _sample-names:

Names and outputs
-----------------

Each sample is named, as in DR10 but with an LS11 prefix
(:func:`ls11samples.vlim.sample_name`):

.. code-block:: text

   LS11_VLIM_ANY_<lo>_<Mr|Mstar>_<hi>_<zmin>_z_<zmax>_N_<N:07d>

   LS11_VLIM_ANY_-25.00_Mr_-20.00_0.02_z_0.17_N_0123456     # -25 < Mr <= -20
   LS11_VLIM_ANY_10.50_Mstar_12.00_0.05_z_0.22_N_0123456    # 10.5 <= log M* < 12

The bounds and redshifts in the name are written with 2 decimals. They are the values used in
the selection, which are also in the ``ZMIN`` / ``ZMAX`` header keywords and in
``vlim_boundaries.fits``.

The outputs go to ``$LS11_OUT/<tag>/vlim/`` (:doc:`data_products`):

* ``<name>_DATA.fits``, ``<name>_RAND.fits`` and ``<name>_COLOUR.fits`` for each sample;
* ``vlim_boundaries.fits``: one row per sample with its bounds, redshift range, number of galaxies,
  comoving volume, number density and medians; plus the two completeness curves;
* ``vlim_planes.png``: the :math:`M_r`–:math:`z` and :math:`M_\star`–:math:`z` planes, with the
  completeness limits and the sample boxes;
* ``manifest.yaml``: the version, git commit, configuration and hashes, area and sample list;
  completed by step 6.

``--code <other>`` defines the :math:`M_\star` samples with another code's masses and writes them
to ``vlim_<code>/``, for comparison.
