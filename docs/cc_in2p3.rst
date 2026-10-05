Running at CC-IN2P3
===================

The full DR11 south footprint is processed at CC-IN2P3 with Slurm. The CC-IN2P3 documentation on
job submission, storage and software is at https://doc.cc.in2p3.fr.

Locations
---------

The home directory is small, so the code, environment, caches and SED-code data all live on
``/sps``:

.. list-table::
   :widths: 30 70

   * - code
     - ``/sps/lsst/users/$USER/software/make-sample-ls11``, a clone of the bare repository
       ``/sps/lsst/users/$USER/git/make-sample-ls11.git``
   * - environment
     - ``/sps/lsst/users/$USER/envs/ls11`` (``$LS11_ENV``)
   * - DR11
     - ``/sps/lsst/datasets/desi/legacysurveys/dr11/south``: ``sweep/11.0``,
       ``sweep/11.0-photo-z`` and ``randoms/``
   * - LePhare
     - ``/sps/lsst/users/$USER/lephare/{data,work}``
   * - products
     - next to the sweeps (``sweep/11.0-<tag|code>/``) and in ``dr11/south/<tag>/``
   * - logs
     - ``dr11/south/<config name>_logs/``

To deploy a new version, push to the bare repository, then run ``git pull`` in the clone. The
working directory is never copied.

Environment
-----------

``cc_in2p3/env_ccin2p3.sh`` sets every variable of :doc:`environment` for CC-IN2P3. It is sourced
by hand and by every job. Its defaults are:

* ``LS11_CONFIG`` = ``config/bgs_r21_dr10bits.yaml`` (tag ``bgsr21``);
* ``LS11_RANDOMS`` = ``randoms-south-1-[01].fits``: two random files, 5000 deg\ :sup:`-2`;
* ``LS11_ENV`` and ``LS11_ACTIVATE``: where the Python environment is and how to activate it;
* ``CONDA_PKGS_DIRS`` and ``PIP_CACHE_DIR``, on ``/sps``.

Every variable can be overridden before sourcing the script.

One-off set-up
--------------

Run this on a login node, which has internet access:

.. code-block:: bash

   source cc_in2p3/env_ccin2p3.sh && bash cc_in2p3/setup_env.sh

``setup_env.sh`` does four things:

1. creates ``$LS11_ENV`` from ``environment.yml``;
2. installs CIGALE v2025.1 from a git-lfs clone (its data are in git-lfs; ``np.trapz`` is patched
   to ``np.trapezoid`` for numpy ≥ 2.4);
3. fetches the LePhare data and registers the CIGALE filters
   (``04_stellar_mass.py --fetch --code lephare,cigale,kcorrect``);
4. checks that ``kcorrect``, ``lephare``, ``pcigale`` and ``ls11samples`` import.

Submission
----------

.. code-block:: bash

   source cc_in2p3/env_ccin2p3.sh && bash cc_in2p3/submit_all.sh            # submit
   source cc_in2p3/env_ccin2p3.sh && bash cc_in2p3/submit_all.sh --dry-run  # print the sbatch lines

``submit_all.sh`` submits the steps to the ``htc`` partition with the ``sps`` licence, chained by
``afterok`` dependencies. Each job runs ``cc_in2p3/job_step.sh <script> [options]``, which:

* sources the environment and activates the Python environment;
* sets ``LS11_NPROC`` and ``OMP_NUM_THREADS`` to ``SLURM_CPUS_PER_TASK``;
* runs the script from the repository.

.. list-table::
   :header-rows: 1
   :widths: 14 30 22 34

   * - step
     - script
     - after
     - resources
   * - ``select``
     - ``01_select.py``, array
     - —
     - 2 h, 4 cores, 16 GB
   * - ``randoms``
     - ``02_randoms.py``
     - select
     - 1 day, 4 cores, 32 GB
   * - ``prepare``
     - ``04_stellar_mass.py --prepare``
     - —
     - 6 h, 16 cores, 32 GB
   * - ``fit``
     - ``04_stellar_mass.py --code <code>``, one array per code
     - select, prepare
     - LePhare, CIGALE: 1 day, 16 cores, 32 GB; kcorrect: 6 h, 2 cores, 8 GB
   * - ``vlim``
     - ``05_vlim.py``
     - randoms, all fits
     - 12 h, 4 cores, 128 GB
   * - ``export``
     - ``06_export.py``
     - vlim
     - 2 h, 1 core, 16 GB

**Frozen sweep list.** The DR11 sweeps may still be downloading at submission time.
``submit_all.sh`` therefore writes the list of sweeps that already have their photo-z file to
``<logs>/sweeps_<timestamp>.txt`` and exports it as ``LS11_SWEEP_LIST``. Every array task then
slices the same list (:ref:`sweep-list`). To reuse an earlier list, set ``LS11_SWEEP_LIST``
before submitting.

Four environment variables control what is submitted:

``STEPS``
   the steps to submit. The default is ``select randoms prepare fit``; add ``vlim export`` once the
   masses have been checked.
``CODES``
   the SED codes (default ``lephare cigale kcorrect``).
``NARRAY``
   the number of array tasks (default 100). Task *i* processes sweeps *i*, *i* + ``NARRAY``, …
``PREPARE_JOB``
   the job ID of a running or finished ``prepare`` job. The fit arrays then depend on it, and no
   new ``prepare`` job is submitted.

Restarts
--------

Existing outputs are skipped (:doc:`pipeline`), so resubmitting picks up only what is missing:

* sweeps whose photo-z file has appeared since the last submission: each submission writes a new
  sweep list;
* tasks that failed or ran out of time.

Step 2 redoes the randoms and maps whenever the sweep list has grown, since the footprint area
depends on it. It keeps them when nothing has changed.
