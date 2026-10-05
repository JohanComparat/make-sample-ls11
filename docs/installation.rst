Installation
============

Requirements
------------

Python ≥ 3.10. There are two levels of dependencies.

**Core** (steps 1, 2, 5, 6, and the tests): numpy, scipy, astropy, fitsio, healpy and pyyaml.
matplotlib is optional; step 5 uses it for its diagnostic figure.

**SED codes** (step 4, the :math:`M_r` completeness curve of step 5, and the benchmark):

.. list-table::
   :header-rows: 1
   :widths: 18 30 52

   * - code
     - package
     - used for
   * - kcorrect
     - ``kcorrect`` ≥ 5.1
     - production; filter curves shared by all codes; :math:`M_r` completeness limit
   * - LePhare
     - ``lephare`` ≥ 1.0
     - production (primary stellar mass)
   * - CIGALE
     - ``pcigale`` ≥ 2025.0
     - production
   * - eazy
     - ``eazy``
     - benchmark only
   * - DSPS
     - ``dsps`` (and ``jax``)
     - benchmark only

All SED codes read the DECam and WISE response curves shipped with kcorrect
(:mod:`ls11samples.sed.filters`). kcorrect is therefore needed whenever another code runs.

Installing
----------

.. code-block:: bash

   pip install -e .          # core
   pip install -e .[sed]     # + kcorrect, lephare, dsps, eazy, pcigale
   pip install -e .[test]    # + pytest
   pip install -e .[docs]    # + sphinx, to build this documentation

On a cluster, ``environment.yml`` creates a conda environment with the core dependencies,
kcorrect, LePhare and the package itself:

.. code-block:: bash

   conda env create -f environment.yml

CIGALE is distributed through git-lfs and is installed separately. ``cc_in2p3/setup_env.sh``
shows how (:doc:`cc_in2p3`).

.. note::

   Maintainer set-up: on the laptop, use the shared ``dev`` environment, which is managed in
   ``~/software/dev_env``. Do not create a new environment for this repository.

Data the SED codes need
-----------------------

Each code fetches or builds its data once, in a cache directory set by an environment variable
(:doc:`environment`):

* **LePhare**: auxiliary files (SEDs, extinction laws) in ``$LEPHAREDIR``, downloaded on first
  use; libraries built in ``$LEPHAREWORK``.
* **CIGALE**: the shared filters, registered in the pcigale database on first use.
* **eazy**: the ``eazy-photoz`` templates and filters, cloned into ``$EAZY_DATA`` when missing.
* **DSPS**: the SSP file ``ssp_data_fsps_v3.2_lgmet_age.h5``, which must be placed in
  ``$DSPS_DRN``.

``scripts/04_stellar_mass.py --fetch`` downloads what the production codes need; it must run on a
node with internet access. ``--prepare`` also builds the LePhare libraries. Run both once before
any job array, so that parallel tasks never build the same files.

Building this documentation
---------------------------

.. code-block:: bash

   make -C docs html            # -> docs/_build/html/index.html
   make -C docs html O=-W       # warnings as errors

The API pages mock fitsio, healpy and the SED codes. The documentation therefore builds with
Sphinx, numpy, astropy and pyyaml only (``docs/requirements.txt``).
