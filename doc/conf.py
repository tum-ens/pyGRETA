# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.append(os.path.abspath("../code/"))
import sphinx_rtd_theme
import sphinxcontrib.bibtex


# -- Project information -----------------------------------------------------

project = "pyGRETA"
copyright = "ENS 2019"
author = "Kais Siala, Houssame Houmy"


# The full version, including alpha/beta/rc tags
release = "1.1.0"
version = release

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ["sphinx.ext.autodoc", "sphinx_rtd_theme", "sphinxcontrib.bibtex"]

master_doc = "index"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]


# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

html_theme_options = {
    "canonical_url": "",
    # 'analytics_id': 'UA-XXXXXXX-1',  #  Provided by Google in your dashboard
    "logo_only": False,
    "display_version": True,
    "prev_next_buttons_location": "both",
    "style_external_links": False,
    # Toc options
    "collapse_navigation": False,
    "sticky_navigation": True,
    "navigation_depth": 4,
    "includehidden": True,
    "titles_only": False,
}

# -- Options for LatexPDF output ----------------------------------------------

# 'startdocname': '',  # (path) Start file for the documentation to be included in PDF, can be left empty to use default index.rst
# 'targetname': project,  # (str) Output name of the Latex file generated
# 'title':'',  # (str) Title of the Latex/pdf file, can be left empty to use the title of startdocname
# 'author': 'Kais Siala, \\Sergio Alejandro Huezo Rodriguez, \\and Houssame Houmy. ',  # (str) Authors, use \\ to separate authors (e.i. 'John \\and Sarah')
# 'documentclass': '',  # not clear
# 'toctree_only': True  # (bool) Include startdocname in the latex/pdf ? can be used to have different first pages. The first toctree entry in startdocname will be used.

# latex_documents = [(master_doc, project+'.tex', project, 'Kais Siala, Houssame Houmy and Sergio Alejandro Huezo Rodriguez', 'manual', True)]
# Kais Siala \\ Houssame Houmy \\ Sergio Alejandro Huezo Rodriguez \vspace{1cm} \\ Version 1.0.0
latex_documents = [
    (
        master_doc,  # startdocname
        project + ".tex",  # targetname
        project,  # title
        r""" Kais Siala \\ Houssame Houmy \\ Sergio Alejandro Huezo Rodriguez \vspace{1cm} \\ Version 1.0.1""",  # author
        "manual",  # documentclass
        True,
    )
]  # toctree_only

# Remove redundant white pages
latex_elements = {
    "classoptions": "twoside",
    "papersize": "a4paper",
    "pointsize": "11pt",
    "passoptionstopackages": r"""
        \usepackage{charter}
        \usepackage[T1]{fontenc}
        \usepackage{tabulary}
        \usepackage{fancyvrb}
        \usepackage{upquote}
        \usepackage{capt-of}
        \usepackage{needspace}
        \usepackage{inconsolata}
        \makeatletter
        \fancypagestyle{normal}{
        \fancyhf{}
        \fancyfoot[LE,RO]{{\py@HeaderFamily\thepage}}
        \fancyfoot[LO]{{\py@HeaderFamily\nouppercase{\rightmark}}}
        \fancyfoot[RE]{{\py@HeaderFamily\nouppercase{\leftmark}}}
        \fancyhead[LE,RO]{{\py@HeaderFamily \@title, \py@release}}
        \renewcommand{\headrulewidth}{0.4pt}
        \renewcommand{\footrulewidth}{0.4pt}
        }
        \makeatother
    """,
    "fncychap": "",
    "maketitle": "\\maketitle",
}

latex_toplevel_sectioning = "chapter"
# This value determines the topmost sectioning unit. It should be chosen from 'part', 'chapter' or 'section'.
# The default is None; the topmost sectioning unit is switched by documentclass: section is used if documentclass will be howto, otherwise chapter will be used.
# Note that if LaTeX uses \part command, then the numbering of sectioning units one level deep gets off-sync with HTML numbering, because LaTeX numbers continuously \chapter


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []

# Modules to be ignored
autodoc_mock_imports = [
    "gdal",
    "osr",
    "osgeo",
    "numpy",
    "os",
    "glob",
    "psutil",
    "datetime",
    "inspect",
    "sys",
    "math",
    "rasterio",
    "pandas",
    "scipy",
    "geopandas",
    "shapely",
    "fiona",
    "hdf5storage",
    "multiprocessing",
    "itertools",
    "h5netcdf",
    "cProfile",
    "pstats",
    "shutil",
    "pyomo",
    "pdb",
    "logging",
    "networkx",
]
