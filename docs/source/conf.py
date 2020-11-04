# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../../abin_launcher'))
sys.path.insert(0, os.path.abspath('../../check_scripts'))
sys.path.insert(0, os.path.abspath('../../control_launcher'))
sys.path.insert(0, os.path.abspath('../../crontab_scripts'))
sys.path.insert(0, os.path.abspath('../../results_treatment'))
sys.path.insert(0, os.path.abspath('../../'))


# -- Project information -----------------------------------------------------

project = 'CHAINS'
copyright = '2020, Iacobellis Nicolas'
author = 'Iacobellis Nicolas'
master_doc = 'index'

# The full version, including alpha/beta/rc tags
release = 'Beta'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = []
extensions.append('sphinx.ext.autodoc')
extensions.append('sphinx.ext.napoleon')

extensions.append('sphinx.ext.todo')
todo_include_todos=True

# Get argparse informations (https://github.com/alex-rudakov/sphinx-argparse)
extensions.append('sphinxarg.ext')

# View source code (https://www.sphinx-doc.org/en/master/usage/extensions/viewcode.html#module-sphinx.ext.viewcode)
extensions.append('sphinx.ext.viewcode')

# Include source code (https://sphinx-code-include.readthedocs.io/en/latest/)
extensions.append('code_include.extension')

# Embed youtube videos (https://github.com/divi255/sphinxcontrib.youtube)
# extensions.append('sphinxcontrib.yt')

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# Theme official page : https://sphinx-rtd-theme.readthedocs.io/en/stable/index.html
html_theme = 'sphinx_rtd_theme'
html_theme_path = ['_themes', ]

html_theme_options = {
    'navigation_depth': 2
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Define a CSS file for specific modifications to the Sphinx RTD theme (see https://stackoverflow.com/questions/23211695/modifying-content-width-of-the-sphinx-theme-read-the-docs for reference)
def setup(app):
    app.add_css_file('my_theme.css')

# Other options

html_show_copyright = False

# -- Options for LaTeX output ------------------------------------------------

latex_elements = {

    # Paper size ('letterpaper' or 'a4paper')
    'papersize': 'a4paper',

    # Font size ('10pt', '11pt' or '12pt')
    'pointsize': '11pt',

}