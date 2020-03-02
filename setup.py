#!/usr/bin/env python

import sys

from setuptools import find_packages, setup

if sys.version_info < (3, 6):
    print("pfapack requires Python 3.6 or above.")
    sys.exit(1)

with open("README.md") as f:
    readme = f.read()

extras_require = dict(
    docs=[
        "sphinx",
        "sphinx-rtd-theme",
        "m2r",  # markdown support
        "sphinxcontrib.apidoc",  # run sphinx-apidoc when building docs
    ],
    dev=["pre-commit", "bump2version"],
    testing=["pytest", "pytest-cov", "pytest-mypy", "tox"],
)

install_requires = ["scipy", "numpy"]


setup(
    name="pfapack",
    python_requires=">=3.6",
    version="0.1.1",
    classifiers=[
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    description="Efficient numerical computation of the Pfaffian for dense and banded skew-symmetric matrices.",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/basnijholt/pfapack",
    author="Bas Nijholt (package) and M. Wimmer (code)",
    author_email="basnijholt@gmail.com",
    license="MIT",
    packages=find_packages("."),
    install_requires=install_requires,
    extras_require=extras_require,
    zip_safe=False,
)
