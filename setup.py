#!/usr/bin/env python

import sys
import os
import subprocess
import numpy
from setuptools import find_packages, setup, Extension
from setuptools.command.build_ext import build_ext

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
    from importlib.util import module_from_spec, spec_from_file_location

    spec = spec_from_file_location("version", os.path.join(package_name, "_version.py"))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.__version__, module.cmdclass


version, cmdclass = get_version_and_cmdclass("pfapack")

# Custom build command to compile the C and Fortran libraries
class CustomBuildExtCommand(build_ext):
    """Custom build command to compile the C library."""
    def run(self):
        # Compile the Fortran library
        subprocess.check_call(['make', '-C', 'external/fortran'])
        # Compile the C interface library
        subprocess.check_call(['make', '-C', 'external/c_interface'])
        # Move the shared object file to the correct location
        os.makedirs(os.path.join(self.build_lib, 'pfapack'), exist_ok=True)
        subprocess.check_call(['cp', 'external/c_interface/libcpfapack.so', os.path.join(self.build_lib, 'pfapack')])
        super().run()

setup(
    name="pfapack",
    python_requires=">=3.7",
    version=version,
    cmdclass={'build_ext': CustomBuildExtCommand},
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
    packages=find_packages(),
    install_requires=install_requires,
    extras_require=extras_require,
    zip_safe=False,
    include_dirs=[numpy.get_include()],
)
