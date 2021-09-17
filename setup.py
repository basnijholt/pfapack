#!/usr/bin/env python

import sys

from setuptools import find_packages, setup

if sys.version_info < (3, 7):
    print("pfapack requires Python 3.7 or above.")
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
    dev=["pre-commit"],
    testing=["pytest", "pytest-cov", "pytest-mypy", "tox"],
)

install_requires = ["scipy", "numpy"]


def get_version_and_cmdclass(package_name):
    import os
    from importlib.util import module_from_spec, spec_from_file_location

    spec = spec_from_file_location("version", os.path.join(package_name, "_version.py"))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.__version__, module.cmdclass


version, cmdclass = get_version_and_cmdclass("pfapack")


setup(
    name="pfapack",
    python_requires=">=3.7",
    version=version,
    cmdclass=cmdclass,
    classifiers=[
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
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
