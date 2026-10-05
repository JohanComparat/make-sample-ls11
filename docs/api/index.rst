API reference
=============

The package is ``ls11samples``. Tables are plain ``dict[str, numpy.ndarray]`` throughout
(:mod:`ls11samples.io`), configurations are nested dicts read from YAML
(:mod:`ls11samples.config`) and every path comes from :func:`ls11samples.env.get_paths`.

.. toctree::
   :maxdepth: 2

   paths_config
   io
   selection
   randoms_maps
   sed
   vlim
