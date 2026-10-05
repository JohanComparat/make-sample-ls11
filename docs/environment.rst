Paths and environment variables
===============================

All paths come from environment variables, read by :func:`ls11samples.env.get_paths`. On the
laptop the defaults are enough. On CC-IN2P3, ``cc_in2p3/env_ccin2p3.sh`` sets them
(:doc:`cc_in2p3`).

.. list-table::
   :header-rows: 1
   :widths: 22 38 40

   * - variable
     - meaning
     - default
   * - ``LS11_DIR``
     - Legacy Surveys DR11 root
     - ``~/data/legacysurvey/dr11``
   * - ``LS11_REGION``
     - ``south`` or ``north``
     - ``south``
   * - ``LS11_SWEEP_VER``
     - sweep version directory
     - ``11.0``
   * - ``LS11_SWEEPS``
     - glob restricting the sweeps
     - ``sweep-*.fits``
   * - ``LS11_SWEEP_LIST``
     - file with one sweep file name per line; only those sweeps are
       processed
     - unset (no restriction)
   * - ``LS11_RANDOMS``
     - glob of the random files, relative to ``$LS11_DIR/$LS11_REGION/randoms`` or absolute
     - ``randoms-$LS11_REGION-1-0.fits``
   * - ``LS11_OUT``
     - root of the per-run products, in ``<LS11_OUT>/<tag>/``
     - ``$LS11_DIR/$LS11_REGION``
   * - ``LS11_SWEEP_OUT``
     - root of the per-sweep products, in ``<LS11_SWEEP_OUT>/<ver>-<name>/``
     - ``$LS11_DIR/$LS11_REGION/sweep``
   * - ``LS11_CONFIG``
     - selection configuration (:doc:`configuration`); ``--config`` overrides it
     - ``config/default.yaml``
   * - ``LS11_GAIA_MAPS``
     - full-sky Gaia star-density maps, ``<nside:04d>/GAIA_*.fits``
     - ``~/data/legacysurvey/dr10/systematics``
   * - ``LS11_NPROC``
     - worker processes; ``--nproc`` overrides it
     - ``min(8, ncpu)``
   * - ``LEPHAREDIR``, ``LEPHAREWORK``
     - LePhare data and libraries
     - ``~/.cache/lephare/{data,work}``
   * - ``EAZY_DATA``
     - eazy-photoz templates and filters (cloned when missing)
     - ``~/.cache/eazy-photoz``
   * - ``DSPS_DRN``
     - directory of the DSPS SSP file
     - ``~/.cache/dsps``
   * - ``TMPDIR``
     - scratch space for the CIGALE runs
     - system default

.. _sweep-list:

Which sweeps are processed
--------------------------

:meth:`ls11samples.env.Paths.sweeps` returns the sweeps of
``$LS11_DIR/$LS11_REGION/sweep/$LS11_SWEEP_VER`` that match ``LS11_SWEEPS`` **and** have their
photo-z file in ``<ver>-photo-z/``. A sweep whose photo-z file is still downloading is left out,
and picked up by the next run.

When ``LS11_SWEEP_LIST`` names a file, only the sweeps listed in it are considered. On CC-IN2P3,
``submit_all.sh`` writes this list at submission, so every array task slices the same list even
while new sweeps keep arriving (:doc:`cc_in2p3`).

Restricting ``LS11_SWEEPS`` is the way to run a small test, for example:

.. code-block:: bash

   export LS11_SWEEPS='sweep-000m0*.fits'        # the strip 0 < RA < 5, -15 < Dec < 0
   export LS11_OUT=$HOME/scratch/ls11_test
   export LS11_SWEEP_OUT=$HOME/scratch/ls11_test/sweep

The randoms of step 2 are restricted to the boxes of the processed sweeps. The area and the random
catalogue therefore always match the processed galaxies.

File names
----------

:class:`ls11samples.env.Paths` builds every file name:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - method
     - path
   * - ``sweep_dir``
     - ``$LS11_DIR/$LS11_REGION/sweep/<ver>``
   * - ``pz_dir``
     - ``$LS11_DIR/$LS11_REGION/sweep/<ver>-photo-z``
   * - ``product_dir(name)``
     - ``$LS11_SWEEP_OUT/<ver>-<name>``
   * - ``product(sweep, name)``
     - ``$LS11_SWEEP_OUT/<ver>-<name>/<sweep stem>-<name>.fits``
   * - ``run_dir(tag, *parts)``
     - ``$LS11_OUT/<tag>/<parts>`` (created when missing)
   * - ``rand_file(tag)``
     - ``$LS11_OUT/<tag>/LS11_<tag>_RAND.fits``
