Testing
=======

.. code-block:: bash

   pip install -e .[test]       # pytest, pytest-cov, kcorrect, matplotlib, configobj
   pytest                       # about 20 s
   pytest --cov                 # with coverage of ls11samples and scripts (fails below 90%)

The tests are fast and need no real data: ``tests/fake_dr11.py`` writes a small synthetic DR11
tree (two sweeps of 3000 objects with real bit headers, one sweep whose header lacks ``MBIT_19``,
one without photo-z, row-matched photo-z files, a random file and Gaia maps), on which the whole
pipeline runs in a few seconds. LePhare and CIGALE are exercised through stand-ins of the
``lephare`` module and of the ``pcigale`` command line, so their libraries and databases are not
needed. The tests that read a real sweep take the first local sweep
(:meth:`ls11samples.env.Paths.sweeps`) and are skipped when there is none; the DSPS backend test
runs when ``dsps`` and its SSP file are installed; the eazy backend test only with
``LS11_SLOW_TESTS=1``.

GitHub Actions (``.github/workflows/tests.yml``) runs ``pytest --cov`` with Python 3.11 and 3.12 at
every push to ``main`` and on pull requests.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - file
     - what is tested
   * - ``tests/test_bits_io.py``
     - bit-definition check against a header (match and mismatch); bit masks from names;
       ``LS_ID_DR11`` and sweep-name parsing; FITS write/read round trip
   * - ``tests/test_cli.py``
     - the job-array slicing covers every item exactly once
   * - ``tests/test_layout.py``
     - configuration deep merge; per-sweep and per-run file names; joining code files to the
       selection, and the row-alignment check
   * - ``tests/test_selection.py``
     - Galactic latitude against astropy; the north/south split; the same footprint for galaxies
       and randoms; the SGA exemption; ``BEST_Z`` priority; the fast sweep-box test against the
       loop; a random file is redone when the footprint, sweeps or random files change; a monotonic
       cut flow on a real sweep; uniform area fraction; ``bgsr21`` applies only
       the bit and magnitude cuts; ``STAR_FLAG``
   * - ``tests/test_vlim.py``
     - on a flux-limited mock with a known K-correction distribution: the :math:`M_r` limit matches
       the analytic one; the :math:`M_\star` limit is monotonic; samples are defined and written; the
       redshift range in a sample name is the one selected
   * - ``tests/test_pipeline.py``
     - the scripts end to end on the synthetic tree, in-process: selection (the sweep without
       photo-z skipped, existing files kept, another configuration refused), randoms and maps (kept
       when up to date, redone for another sweep list), kcorrect fits row-aligned with the
       selection, volume-limited samples, export checks and run validation; then the failures:
       a sweep without selection, another footprint, randoms of other sweeps, a code file of another
       selection, broken samples, and step 5 without kcorrect or matplotlib
   * - ``tests/test_sed_wrappers.py``
     - the shared photometry and filter curves; LePhare (data download, library build or reuse,
       input table, output mapping, ``MABS_R_ERR`` from the band nearest to rest-frame r) and CIGALE
       (filter registration, generated ``pcigale.ini`` with ``additionalerror = 0``, flux units,
       ``MABS_R`` from the rest-frame luminosity, run directory removed) with stand-ins; DSPS helpers
       and backend; eazy (slow, opt-in)
   * - ``tests/test_edges.py``
     - error and edge branches: sweep names, blank-padded strings, crashed writers, photo-z checks,
       headers without bits, existing selections, off-grid sweep boxes, random densities, missing
       map quantities and Gaia maps, completeness-curve edges, plots without masses, kcorrect without
       usable redshifts, validation of missing or empty sweeps, the subprocess wrapper
