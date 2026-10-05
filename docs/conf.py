"""Sphinx configuration of the make-sample-ls11 documentation.

Build:  make -C docs html   (or  sphinx-build -b html docs docs/_build/html)
"""

import os
import sys

sys.path.insert(0, os.path.abspath(".."))

from ls11samples import __version__  # noqa: E402 - needs the path above

project = "make-sample-ls11"
author = "Johan Comparat"
copyright = "2026, Johan Comparat"
release = __version__
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "sphinx.ext.intersphinx",
    "sphinx_copybutton",
    "myst_parser",          # the benchmark report is Markdown, included as written
]

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
master_doc = "index"
exclude_patterns = ["_build"]

myst_enable_extensions = ["dollarmath", "colon_fence"]
myst_heading_anchors = 3

# The compiled readers and the SED codes are mocked: the API pages build with numpy, astropy and
# pyyaml only (the SED codes are imported inside the backends anyway).
autodoc_mock_imports = ["fitsio", "healpy", "kcorrect", "lephare", "pcigale", "pcigale_filters",
                        "eazy", "dsps", "jax", "configobj"]
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "member-order": "bysource",
}
autodoc_typehints = "signature"
autodoc_preserve_defaults = True
# return annotations with no entry in the intersphinx inventories
nitpick_ignore = [("py:class", "fitsio.FITSHDR"), ("py:class", "numpy.int64")]
napoleon_google_docstring = False
napoleon_numpy_docstring = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "astropy": ("https://docs.astropy.org/en/stable", None),
}

html_theme = "furo"
html_title = f"make-sample-ls11 {release}"
html_theme_options = {
    "source_repository": "https://github.com/JohanComparat/make-sample-ls11/",
    "source_branch": "main",
    "source_directory": "docs/",
}
