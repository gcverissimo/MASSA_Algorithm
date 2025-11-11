# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'MASSA Algorithm'
copyright = '2025, Gabriel C. Veríssimo'
author = 'Gabriel C. Veríssimo'
release = 'v. 2.2.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",               # for Markdown (.md) files
    "sphinx.ext.autodoc",        # for automatic API documentation
    "sphinx.ext.napoleon",       # for Google/Numpy-style docstrings
    "sphinx.ext.viewcode",       # add links to highlighted source code
    "sphinx.ext.githubpages",    # adds .nojekyll for GitHub Pages
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# The master toctree document.
master_doc = "index"

# Language
language = "en"

# Patterns to ignore when looking for source files.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "furo"

# Optional: theme customization
html_theme_options = {
    "sidebar_hide_name": False,
}

# Optional: add your logo or favicon
html_logo = "_static/logo-light.png"
html_favicon = "_static/favicon.ico"

# Custom static files (such as style sheets)
html_static_path = ["_static"]

# Example: Add custom CSS if you want tweaks
# Create docs/source/_static/custom.css and uncomment:
# html_css_files = [
#     "custom.css",
# ]

# -- MyST (Markdown) options -------------------------------------------------
myst_enable_extensions = [
    "colon_fence",       # ::: fenced blocks
    "deflist",           # definition lists
    "linkify",           # autolink URLs
    "substitution",      # {{ variables }}
    "replacements",      # smart quotes, dashes, etc.
]

# -- HTML context for GitHub integration -------------------------------------
# These settings help with "Edit on GitHub" links (optional)
html_context = {
    "display_github": True,
    "github_user": "gcverissimo",
    "github_repo": "Gabriel C. Veríssimo",
    "github_version": "main",
    "conf_py_path": "./",
}

# Prevent broken pages on GitHub Pages (ensures root index.html works)
html_baseurl = "https://gcverissimo.github.io/MASSA_Algorithm/"

