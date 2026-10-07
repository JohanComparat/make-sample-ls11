Installation
============

The pipeline needs a small Python core and, for the stellar masses, external SED-fitting codes
with their own data. The SED codes are not installed with the package: several are not on PyPI,
need data files of 60 MB to 10 GB, or build libraries once. Install only the codes you use, by
following the recipes below. Each recipe ends with a command that checks the result.

.. list-table:: What each part is needed for
   :header-rows: 1
   :widths: 14 22 44 20

   * - part
     - tested version
     - needed for
     - data on disk
   * - core
     - Python 3.11, 3.12
     - steps 1, 2, 5, 6, ``validate_run.py``, the tests
     - none
   * - kcorrect
     - 5.1.9
     - **every SED step**: kcorrect masses and :math:`M_r`; the DECam / WISE filter curves used by
       all codes; the :math:`M_r` completeness limit of step 5
     - shipped in the wheel
   * - LePhare
     - 1.0.0
     - production masses (``sed.primary``)
     - 82 MB + 2.2 GB of libraries
   * - CIGALE
     - 2025.1
     - production masses
     - about 10 GB
   * - eazy
     - 0.8.7
     - benchmark only
     - 200 MB
   * - DSPS
     - 0.4.8
     - benchmark only
     - 60 MB

1. The core environment
-----------------------

Create a conda environment from conda-forge (fitsio and healpy are compiled packages; conda-forge
ships them as binaries), then install the package itself without dependencies:

.. code-block:: bash

   conda create -n ls11 -c conda-forge python=3.12 numpy scipy astropy fitsio healpy pyyaml \
       matplotlib pytest pytest-cov configobj
   conda activate ls11
   git clone https://github.com/JohanComparat/make-sample-ls11.git
   cd make-sample-ls11
   python -m pip install --no-deps -e .
   pytest                                # about 20 s, no data needed

``matplotlib`` is optional (diagnostic figures); ``configobj`` and ``pytest-cov`` are only used by
the tests. ``environment.yml`` holds the same core plus kcorrect and LePhare
(``conda env create -f environment.yml``); it is what ``cc_in2p3/setup_env.sh`` uses on CC-IN2P3.

.. note::

   Maintainer set-up: on the laptop, use the shared ``dev`` / ``dev-full`` environments managed in
   ``~/software/dev_env`` (``dev-full`` has all the SED codes). Do not create a new environment for
   this repository.

2. kcorrect (needed by every SED step)
--------------------------------------

Pure Python; its templates and the DECam / WISE response curves come with the wheel.

.. code-block:: bash

   python -m pip install kcorrect==5.1.9
   python -c "from ls11samples.sed.filters import pivot_wavelength as p; print(round(p('R')))"   # 6421

3. LePhare (production)
-----------------------

The ``lephare`` wheel contains the compiled C++ core. Its data (SEDs, extinction laws, opacities)
and the libraries built from them live in two directories given by environment variables:

.. code-block:: bash

   python -m pip install lephare==1.0.0
   export LEPHAREDIR=/path/to/lephare/data       # 82 MB of auxiliary files
   export LEPHAREWORK=/path/to/lephare/work      # 2.2 GB of libraries
   python scripts/04_stellar_mass.py --fetch --code lephare      # needs internet access
   python scripts/04_stellar_mass.py --prepare --code lephare    # builds the libraries, once

* ``--fetch`` writes the shared DECam / WISE filters into ``$LEPHAREDIR/filt/ls11`` and downloads
  the 393 auxiliary files that ``config/lephare/LS11_zFIX.para`` needs.
* ``--prepare`` builds the galaxy, star and QSO libraries: 12 BC03 Chabrier star-formation
  histories × 221 ages × 13 :math:`E(B-V)` values × 201 redshifts. This takes about 25 min on 16
  to 23 threads (``OMP_NUM_THREADS``). An interrupted build is detected (empty ``.doc`` file) and
  redone; a finished one is reused.
* The output columns come from ``config/lephare/output_ls11.para``.

Check that ``$LEPHAREWORK/lib_mag/BC03_LS11.doc`` exists and is not empty.

4. CIGALE (production)
----------------------

CIGALE is not on PyPI, and its data are stored with git-lfs. A plain ``pip install`` of the git
repository gets pointer files instead of data, and the build fails. Install it as follows:

.. code-block:: bash

   # git-lfs: from your package manager (apt install git-lfs), conda-forge, or the static binary
   # of https://github.com/git-lfs/git-lfs/releases (no installation needed)
   git -c filter.lfs.smudge="git-lfs smudge -- %f" -c filter.lfs.process="git-lfs filter-process" \
       -c filter.lfs.required=true clone --depth 1 --branch v2025.1 \
       https://gitlab.lam.fr/cigale/cigale.git cigale
   # numpy >= 2.4 has no np.trapz: patch the clone once (same signature)
   grep -rl "np\.trapz(" --include=*.py cigale | xargs sed -i "s/np\.trapz(/np.trapezoid(/g"
   # the install builds the SSP database (about 3 min) into a 3 GB wheel: give it room
   TMPDIR=/a/large/disk python -m pip install ./cigale
   python scripts/04_stellar_mass.py --fetch --code cigale    # registers the ls11.* filters

* Check that ``python -c "import pcigale"`` works.
* The clone with its data and build takes about 10 GB.
* The ``ls11.decam_*`` and ``ls11.wise_*`` filters are added to the pcigale database on first use.
  Reinstalling CIGALE removes them, and the next run adds them back.
* Each fit runs in a temporary directory under ``$TMPDIR``, which is removed after the results are
  read.
* ``pcigale``'s ``additionalerror`` is set to 0.
* pcigale floors the errors of extensive properties at 5% (:doc:`stellar_masses`).

5. eazy (benchmark only)
------------------------

.. code-block:: bash

   python -m pip install eazy==0.8.7      # eazy-py; asks for scipy < 1.18
   export EAZY_DATA=/path/to/eazy-photoz
   git clone --depth 1 https://github.com/gbrammer/eazy-photoz.git $EAZY_DATA    # 200 MB

* If ``$EAZY_DATA`` is missing, the backend clones the repository itself.
* eazy warns that ``dust_attenuation`` is missing at import. That package is not needed here.

6. DSPS (benchmark only)
------------------------

.. code-block:: bash

   python -m pip install dsps==0.4.8      # brings jax (CPU is enough)
   export DSPS_DRN=/path/to/dsps
   wget -P $DSPS_DRN https://portal.nersc.gov/project/hacc/aphearin/DSPS_data/ssp_data_fsps_v3.2_lgmet_age.h5

The SSP file is 60 MB: FSPS v3.2 with default settings, a Kroupa IMF and nebular emission. The
backend converts the masses to a Chabrier IMF.

Checking the SED installation
-----------------------------

``pytest tests/test_sed_wrappers.py`` exercises every backend. LePhare and CIGALE run through
stand-ins, so this checks the wrappers, not the codes. The DSPS test runs when its SSP file is
present. The eazy test runs with ``LS11_SLOW_TESTS=1``.

To fit a few real galaxies with the codes you installed, run the benchmark on a small sample (local
sweeps needed):

.. code-block:: bash

   python benchmarks/sed_benchmark.py --codes kcorrect,lephare,cigale --n 200 --n-zsens 0

Building this documentation
---------------------------

.. code-block:: bash

   python -m pip install -r docs/requirements.txt
   make -C docs html            # -> docs/_build/html/index.html
   make -C docs html O=-W       # warnings as errors

The API pages mock fitsio, healpy and the SED codes. The documentation therefore builds with
Sphinx, numpy, scipy, astropy and pyyaml only (``docs/requirements.txt``). Read the Docs installs
exactly this file.
