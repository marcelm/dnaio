# Sphinx configuration file

import os
import sys

sys.path.insert(0, os.path.abspath(".."))

project = "dnaio"
copyright = "2021 dnaio authors"
author = "Marcel Martin"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
]

templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

html_theme = "furo"
html_show_sphinx = False

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ["_static"]

# default_role = "obj" (or "any")
