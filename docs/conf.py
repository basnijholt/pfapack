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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

import os
from pathlib import Path

import pfapack

package_path = Path("..").resolve()
docs_path = Path().resolve()


# -- Project information -----------------------------------------------------

project = "pfapack"
copyright = "2019, Bas Nijholt"
author = "Bas Nijholt"
version = pfapack.__version__
release = pfapack.__version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.viewcode",  # View source code from documentation pages
    "sphinx.ext.napoleon",  # numpy-style docstrings
    "sphinxcontrib.apidoc",  # Run sphinx-apidoc when building
    "myst_parser",
]
# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

default_role = "autolink"

# -- Extension configuration -------------------------------------------------

napoleon_google_docstring = False

autosummary_generate = True

apidoc_module_dir = os.path.dirname(pfapack.__file__)
apidoc_output_dir = "reference"
apidoc_separate_modules = True
apidoc_toc_file = False
apidoc_module_first = True
apidoc_extra_args = ["--force"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# The license file has no extension, so Sphinx ignores it by default, so we
# must add it here
html_extra_path = ["LICENSE"]


def _change_alerts_to_admonitions(input_text: str) -> str:
    # Splitting the text into lines
    lines = input_text.split("\n")

    # Placeholder for the edited text
    edited_text = []

    # Mapping of markdown markers to their new format
    mapping = {
        "IMPORTANT": "important",
        "NOTE": "note",
        "TIP": "tip",
        "WARNING": "caution",
    }

    # Variable to keep track of the current block type
    current_block_type = None

    for line in lines:
        # Check if the line starts with any of the markers
        if any(line.strip().startswith(f"> [!{marker}]") for marker in mapping):
            # Find the marker and set the current block type
            current_block_type = next(
                marker for marker in mapping if f"> [!{marker}]" in line
            )
            # Start of a new block
            edited_text.append("```{" + mapping[current_block_type] + "}")
        elif current_block_type and line.strip() == ">":
            # Empty line within the block, skip it
            continue
        elif current_block_type and not line.strip().startswith(">"):
            # End of the current block
            edited_text.append("```")
            edited_text.append(line)  # Add the current line as it is
            current_block_type = None  # Reset the block type
        elif current_block_type:
            # Inside the block, so remove '>' and add the line
            edited_text.append(line.lstrip("> ").rstrip())
        else:
            # Outside any block, add the line as it is
            edited_text.append(line)

    # Join the edited lines back into a single string
    return "\n".join(edited_text)


def change_alerts_to_admonitions(input_file: Path, output_file: Path) -> None:
    """Change markdown alerts to admonitions.

    For example, changes
    > [!NOTE]
    > This is a note.
    to
    ```{note}
    This is a note.
    ```
    """
    with input_file.open("r") as infile:
        content = infile.read()
    new_content = _change_alerts_to_admonitions(content)

    with output_file.open("w") as outfile:
        outfile.write(new_content)


readme_path = package_path / "README.md"
docs_readme_path = docs_path / "README.md"
change_alerts_to_admonitions(readme_path, docs_readme_path)
