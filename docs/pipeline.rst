Running the pipeline
====================

Quick start
-----------

.. code-block:: bash

   export LS11_CONFIG=config/default.yaml      # or --config on each script
   python scripts/01_select.py                 # selection, per sweep
   python scripts/02_randoms.py                # randoms + maps
   python scripts/04_stellar_mass.py --prepare # LePhare libraries, once
   python scripts/04_stellar_mass.py           # M*, Mr for every code of sed.codes
   python scripts/05_vlim.py                   # volume-limited samples
   python scripts/06_export.py                 # format checks + manifest
   python scripts/validate_run.py              # completeness and sanity of the whole run

``scripts/run_all.sh [options]`` runs the same sequence. Its options are passed to steps 1–5. It
uses the Python interpreter given by ``PY``; the default is the ``dev-full`` environment, which
has the SED codes.

Common options
--------------

Every script accepts these options (:func:`ls11samples.cli.parser`):

.. list-table::
   :widths: 25 75

   * - ``--config FILE``
     - selection YAML (default: ``$LS11_CONFIG``, else ``config/default.yaml``)
   * - ``--nproc N``
     - worker processes (default: ``$LS11_NPROC``, else ``min(8, ncpu)``)
   * - ``--overwrite``
     - redo existing outputs, including ones made with another configuration
   * - ``-v``, ``--verbose``
     - debug logging

Steps 1 and 4 can also be split into parts (:func:`ls11samples.cli.add_slicing`). Part *i* of *n*
processes sweeps *i*, *i + n*, *i + 2n*, … of the sorted sweep list.

.. list-table::
   :widths: 25 75

   * - ``--part I``
     - default: ``$SLURM_ARRAY_TASK_ID``, else 0
   * - ``--nparts N``
     - default: ``$SLURM_ARRAY_TASK_COUNT``, else 1

A Slurm job array therefore splits the sweeps without any extra argument.

Step 1: selection
-----------------

``scripts/01_select.py`` runs :func:`ls11samples.bgsl.run_sweep` on each sweep, in parallel over
``--nproc`` processes. For each sweep it:

1. checks the MASKBITS / FITBITS definitions of the sweep header (:func:`ls11samples.bits.check_header`);
2. reads ``FLUX_R`` and ``MW_TRANSMISSION_R`` of all rows and keeps the rows inside ``galaxy.r_range``;
3. reads the columns the cuts need for those rows only, and applies the footprint and galaxy cuts;
4. reads the photo-z rows of the selected objects and checks their IDs against the sweep;
5. chooses ``BEST_Z``, computes ``STAR_FLAG``, and writes the selection file with its ``CUTFLOW`` HDU.

An existing selection file made with the same configuration is kept. One made with another
configuration raises, unless ``--overwrite`` is given. See :doc:`selection` for the cuts.

Step 2: randoms and maps
------------------------

``scripts/02_randoms.py`` first checks that every processed sweep has its selection file, made
with the same footprint (``FPHASH``). It then reads the DR11 random files matching
``$LS11_RANDOMS``, in chunks of ``randoms.chunk_rows`` rows (:func:`ls11samples.randoms.select_randoms`).
A random is kept when both hold:

* it lies inside the box of a processed sweep (:class:`ls11samples.randoms.BoxSet`);
* it passes :func:`ls11samples.selection.footprint_mask`, the footprint cuts of the galaxies.

The footprint area is

.. math::

   A = \frac{N_\mathrm{kept}}{\sum_\mathrm{files} \mathtt{DENSITY}} \quad [\mathrm{deg}^2],

where ``DENSITY`` (2500 deg\ :sup:`-2` per DR11 file) is read from each file header. While the
files are read, :class:`ls11samples.maps.MapAccumulator` builds two kinds of HEALPix maps from the
kept randoms:

* the area fraction of each pixel, :math:`N_\mathrm{pix} / (\sum \mathtt{DENSITY} \times \Omega_\mathrm{pix})`;
* the mean of each ``maps.quantities`` column per pixel. Depths are converted to 5σ AB magnitudes,
  :math:`22.5 - 2.5\log_{10}(5/\sqrt{\mathrm{ivar}})`.

The Gaia star-density maps of ``maps.gaia_templates`` are copied next to the other templates.

Option: ``--no-maps`` writes the randoms only.

An existing random file is kept only if it is up to date (:func:`ls11samples.randoms.stale_reason`):
it must have been made with the same footprint (``FPHASH``), from the same sweeps (``SWPHASH``)
and from the same random files. Otherwise it is redone, for example after step 1 has selected
sweeps downloaded since the last run. ``--overwrite`` redoes it in any case.

Step 4: stellar masses and absolute magnitudes
----------------------------------------------

``scripts/04_stellar_mass.py`` fits every selected galaxy at fixed :math:`z=` ``BEST_Z`` with each
code of ``sed.codes`` (or ``--code``). The photometry is read from the sweep at ``SWEEP_ROW``. The
output is one file per sweep and code, row-aligned with the selection file. The codes parallelise
internally; sweeps are processed one after the other, in chunks of ``--chunk`` galaxies.

.. list-table::
   :widths: 25 75

   * - ``--code a,b``
     - codes to run (default: ``sed.codes``)
   * - ``--chunk N``
     - galaxies per fit call (default 20000)
   * - ``--fetch``
     - only download what the codes need: LePhare auxiliary data, CIGALE filter registration.
       Needs internet access.
   * - ``--prepare``
     - ``--fetch`` and also build the LePhare libraries

Run ``--fetch`` / ``--prepare`` once before the job arrays, so that parallel tasks never build the
same files. An existing output made from the same selection (``SELHASH``) is kept.

The kcorrect Monte Carlo errors are reproducible: their seed is the CRC32 of the sweep name. See
:doc:`stellar_masses` for what each code returns.

Step 5: volume-limited samples
------------------------------

``scripts/05_vlim.py`` works through the following:

1. **Assembly.** It checks that the random file of step 2 covers the same sweeps (``SWPHASH``),
   and stops otherwise: rerun step 2. It then loads every selection file and the code files, then reads ``RA``, ``DEC``,
   ``EBV`` and the g, r, z fluxes back from the sweeps (:func:`ls11samples.catalog.load`).
2. **Stars.** It leaves out objects with ``STAR_FLAG & vlim.exclude_star_flag``.
3. **Columns.** ``LOGMSTAR`` comes from ``sed.primary`` (or ``--code``) and ``MABS_R`` from
   ``vlim.mr_code``.
4. **Limits and samples.** It measures the :math:`M_r` and :math:`M_\star` completeness limits and
   defines one sample per threshold (:doc:`volume_limited`).
5. **Output.** It writes the ``_DATA``, ``_RAND`` and ``_COLOUR`` files of each sample, together
   with ``vlim_boundaries.fits``, ``vlim_planes.png`` and ``manifest.yaml``.

.. list-table::
   :widths: 25 75

   * - ``--code C``
     - stellar-mass code defining the :math:`M_\star` samples (default: ``sed.primary``). Another
       code writes to ``vlim_<code>/`` instead of ``vlim/``.
   * - ``--seed N``
     - seed of the random subsampling and redshift shuffling (default 1)

Step 6: export checks
---------------------

``scripts/06_export.py`` checks every sample listed in ``vlim/manifest.yaml``:

* the required columns and their types;
* finite positions and redshifts;
* ``COLOUR`` row-aligned with ``DATA``.

It then adds to the manifest the ``sys_mapping`` directories (catalogues, templates per nside,
area-fraction maps), the ``sum_stat`` column names and the list of problems. It exits with status 1
when a check fails.

Using the samples downstream:

.. code-block:: bash

   # sys_mapping
   run_ls10_analysis.py --catalog-dir $LS11_OUT/<tag>/vlim --sample <name> \
       --template-dir $LS11_OUT/<tag>/systematics/<nside:04d> --nside <nside>
   # sum_stat
   measure_joint_sumstat.py --survey custom \
       --data-file <name>_DATA.fits --rand-file <name>_RAND.fits

``sum_stat`` reads ``BEST_Z``, ``LPH_MASS_BEST`` and ``RAND.Z``.

Validating a run
----------------

``scripts/validate_run.py`` (:func:`ls11samples.validate.validate`) checks a whole run made with
the current configuration, over the same sweeps as the jobs (``LS11_SWEEP_LIST`` when it is set):

* every sweep has its selection file and one file per code of ``sed.codes``, each row-aligned
  with the selection (same ``SWEEP_ROW``);
* the fraction of missing values (NaN) per code, over all objects and over the galaxies,
  i.e. ``STAR_FLAG & vlim.exclude_star_flag == 0`` and ``BEST_Z`` inside the redshift range of
  the SED grids;
* the 5th, 50th and 95th percentiles of ``LOGMSTAR``, ``LOGMSTAR_ERR``, ``MABS_R`` and
  ``MABS_R_ERR`` per code, on a random subsample of the galaxies (``--frac``, default 2%);
* the median offset and NMAD of each code against ``sed.primary`` (``--ref``);
* the randoms (number, area, sweeps) and the disk use of every product folder.

It prints a report, writes ``<LS11_OUT>/<tag>/validation.yaml`` and exits with status 1 when a
file is missing or misaligned.

.. code-block:: bash

   LS11_CONFIG=config/bgs_r21_dr10bits.yaml python scripts/validate_run.py

Logging
-------

All scripts log through the ``ls11samples`` logger hierarchy, with timestamps. Step 1 logs one
line per sweep (rows, rows inside the r range, selected, time). Step 2 logs one line per random
file. Step 4 logs one line per sweep and code.
