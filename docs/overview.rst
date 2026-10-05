Overview
========

Pipeline
--------

The pipeline has five scripts. Steps 1 and 4 run per sweep and can be split over a job array.
Steps 2, 5 and 6 run once over the whole footprint.

.. list-table::
   :header-rows: 1
   :widths: 8 22 70

   * - step
     - script
     - what it does
   * - 1
     - ``scripts/01_select.py``
     - Applies the selection to each sweep. Writes a small selection file per sweep: IDs, sweep
       row, ``BEST_Z``, star flag and a cut-flow table (:doc:`selection`).
   * - 2
     - ``scripts/02_randoms.py``
     - Applies the same footprint cuts to the DR11 randoms. Measures the area and writes the
       HEALPix area-fraction and systematics template maps.
   * - 4
     - ``scripts/04_stellar_mass.py``
     - Fits each selected galaxy at fixed :math:`z=` ``BEST_Z`` with every SED code of
       ``sed.codes``. Writes one small file per sweep and code: ``LOGMSTAR``, ``MABS_R`` and their
       errors (:doc:`stellar_masses`).
   * - 5
     - ``scripts/05_vlim.py``
     - Measures the completeness limits, then defines the :math:`M_r` and :math:`M_\star`
       volume-limited samples and writes them with their randoms (:doc:`volume_limited`).
   * - 6
     - ``scripts/06_export.py``
     - Checks the samples against the ``sys_mapping`` / ``sum_stat`` input formats. Completes the
       manifest.

There is no step 3 script. The K-corrections come from the kcorrect backend of step 4.
``scripts/run_all.sh`` runs the steps in order on one machine (:doc:`pipeline`).
``cc_in2p3/submit_all.sh`` submits them as Slurm jobs with dependencies (:doc:`cc_in2p3`).

Design
------

**Nothing is copied from the sweeps.**
  The per-sweep products sit next to the sweeps, one small file per sweep and product. Positions
  and photometry are never copied: they are read back from the sweep at ``SWEEP_ROW``
  (:mod:`ls11samples.catalog`). A selection file costs 22 bytes per object; a stellar-mass file
  costs 20.

**One footprint for galaxies and randoms.**
  The footprint cuts use only columns that the randoms also carry. The same function,
  :func:`ls11samples.selection.footprint_cuts`, is applied to both.

**One redshift.**
  ``BEST_Z`` is chosen once, in step 1 (:ref:`best-z`). The K-corrections, masses, absolute
  magnitudes, volume limits and random redshifts all use it.

**Bits by name.**
  The configuration names the MASKBITS and FITBITS. Their numbers are checked against the
  ``MBIT_n`` / ``FBIT_n`` keywords of each sweep header, so a change of bit layout fails loudly.

**Restartable.**
  Steps 1, 2 and 4 skip outputs that are already up to date; steps 5 and 6 are cheap and always
  rewrite.
  Resubmitting therefore only processes what is missing, for example sweeps downloaded since the
  last run. Files are written under a temporary name and then renamed, so an interrupted job never
  leaves a truncated file.

**Traceable.**
  Every output header records hashes of the configuration it was made with (:ref:`hashes`). These
  hashes are checked before anything is mixed:

  * step 1 refuses to keep a selection file made with another configuration;
  * step 2 refuses selection files made with another footprint, and redoes a random file made
    from other sweeps or another footprint;
  * step 4 refuses to keep mass files made for another selection;
  * step 5 refuses a random file made from other sweeps than the selection it reads.

  ``--overwrite`` replaces such files.

Directory layout
----------------

With the default environment (:doc:`environment`), inputs and products sit under the DR11 region
directory:

.. code-block:: text

   $LS11_DIR/$LS11_REGION/
   ├── sweep/                                    ($LS11_SWEEP_OUT)
   │   ├── 11.0/sweep-<box>.fits                 DR11 sweeps                        input
   │   ├── 11.0-photo-z/sweep-<box>-pz.fits      row-matched photo-z files          input
   │   ├── 11.0-<tag>/sweep-<box>-<tag>.fits     selection                          step 1
   │   └── 11.0-<code>/sweep-<box>-<code>.fits   M*, Mr of one SED code             step 4
   ├── randoms/randoms-<region>-1-<i>.fits       DR11 randoms                       input
   └── <tag>/                                    ($LS11_OUT/<tag>)
       ├── LS11_<tag>_RAND.fits                  footprint randoms                  step 2
       ├── footprint/                            area-fraction maps                 step 2
       ├── systematics/<nside>/                  sys_mapping template maps          step 2
       └── vlim/                                 volume-limited samples, manifest   steps 5-6

``<tag>`` is the name of the selection (``tag`` in the configuration, e.g. ``bgsl`` or
``bgsr21``). ``<code>`` is an SED code (``lephare``, ``cigale``, ``kcorrect``, ...). Several
selections can therefore share the same sweeps. :doc:`data_products` lists the columns and header
keywords of every file.

Conventions
-----------

* **Magnitudes** are AB and corrected for Milky-Way extinction:
  :math:`m = 22.5 - 2.5\log_{10}(F / T_\mathrm{MW})`, with :math:`F` the sweep flux in nanomaggies
  and :math:`T_\mathrm{MW}` = ``MW_TRANSMISSION``. The one exception is the Gaia star test, which
  uses the raw r magnitude, as DESI and DR10 do.
* **Cosmology**: flat ΛCDM, :math:`H_0 = 67.74`, :math:`\Omega_m = 0.3089`, as for the DR10
  samples. Volumes are in :math:`\mathrm{Mpc}^3` (not :math:`h^{-1}`).
* **Stellar masses** are :math:`\log_{10} M_\star/M_\odot` for a Chabrier IMF.
* **Absolute magnitudes** ``MABS_R`` are in the rest-frame DECam r band (no band shift).
* **Errors** are 1σ at fixed redshift. The photo-z term is not included; it can be derived from
  ``BEST_Z_ERR``.
* **In memory**, tables are plain ``dict[str, numpy.ndarray]`` in native byte order
  (:mod:`ls11samples.io`).
