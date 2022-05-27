# Sphinx configuration file

import os
import sys
import time

sys.path.insert(0, os.path.abspath(".."))

project = "dnaio"
copyright = f"{time.gmtime().tm_year} dnaio authors"
author = "Marcel Martin"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_issues",
]

templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

html_theme = "furo"
html_show_sphinx = False
html_title = "dnaio"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ["_static"]

default_role = "obj"  # (or "any")

issues_uri = "https://github.com/marcelm/dnaio/issues/{issue}"
issues_pr_uri = "https://github.com/marcelm/dnaio/pull/{pr}"

autodoc_typehints = "description"
python_use_unqualified_type_names = True
