Testing
=======

.. code-block:: bash

   pytest

The tests run in a few seconds and need the core dependencies only. The tests that read a real
sweep take the first local sweep (:meth:`ls11samples.env.Paths.sweeps`) and are skipped when there
is none.

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
