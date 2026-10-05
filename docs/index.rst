make-sample-ls11
================

``make-sample-ls11`` (package ``ls11samples``) builds galaxy samples and matching randoms from the
`Legacy Surveys DR11 <https://www.legacysurvey.org/>`_ for two downstream codes:

* ``sum_stat``, for clustering;
* ``sys_mapping``, for imaging systematics.

It succeeds the DR10 BGS-like / ``LS10_VLIM_*`` samples. It produces:

* a magnitude-limited selection of the DR11 sweeps, with one redshift per object (``BEST_Z``);
* randoms with exactly the same footprint, with HEALPix area-fraction and systematics template
  maps;
* stellar masses and rest-frame r-band absolute magnitudes from three SED codes (LePhare, CIGALE,
  kcorrect), all at fixed redshift;
* volume-limited samples in :math:`M_r` and in :math:`M_\star`, each with its own randoms, in the
  file formats ``sys_mapping`` and ``sum_stat`` read.

The same code runs on a laptop (a few local sweeps) and on the full footprint at CC-IN2P3. Every
path comes from ``LS11_*`` environment variables; every selection choice comes from a YAML file.

.. toctree::
   :maxdepth: 2
   :caption: User guide

   overview
   installation
   environment
   configuration
   pipeline
   cc_in2p3

.. toctree::
   :maxdepth: 2
   :caption: Methods

   selection
   stellar_masses
   volume_limited

.. toctree::
   :maxdepth: 2
   :caption: Reference

   data_products
   benchmark
   testing
   api/index
