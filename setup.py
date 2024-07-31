#!/usr/bin/env python

import sys
import os
import subprocess
from setuptools import find_packages, setup, Command
from setuptools.command.build_py import build_py
from setuptools.command.build_ext import build_ext

if sys.version_info < (3, 9):
    print("pfapack requires Python 3.9 or above.")
    sys.exit(1)

with open("README.md") as f:
    readme = f.read()

extras_require = {
    "docs": ["sphinx", "sphinx-rtd-theme", "m2r", "sphinxcontrib.apidoc"],
    "dev": ["pre-commit"],
    "testing": ["pytest", "pytest-cov", "pytest-mypy", "tox"],
}

install_requires = ["scipy", "numpy", "Cython", "setuptools"]

# Custom command to ensure the C and Fortran libraries are built before the Python package
class CustomBuildExtCommand(build_ext):
    def run(self):
        print("Building external libraries...")
        subprocess.check_call(['make', '-C', 'external/fortran'])
        subprocess.check_call(['make', '-C', 'external/c_interface'])
        lib_path = os.path.abspath("pfapack/libcpfapack.so")
        print(f"Expected location of libcpfapack.so: {lib_path}")
        if os.path.exists(lib_path):
            print("libcpfapack.so found.")
        else:
            print("libcpfapack.so NOT found.")
        super().run()

class CustomBuildPyCommand(build_py):
    """Custom build command to ensure native libraries are built first."""
    def run(self):
        self.run_command('build_ext')
        super().run()

class CustomCleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        os.system('rm -rf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')
        os.system('make clean -C external/c_interface')
        os.system('make clean -C external/fortran')

setup(
    name="pfapack",
    python_requires=">=3.9",
    version="0.0.1",
    cmdclass={
        'build_ext': CustomBuildExtCommand,
        'build_py': CustomBuildPyCommand,
        'clean': CustomCleanCommand,
    },
    classifiers=[
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.9",
    ],
    description="Efficient numerical computation of the Pfaffian for dense and banded skew-symmetric matrices.",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/jjgoings/pfapack",
    author="Bas Nijholt (package) and M. Wimmer (code)",
    author_email="basnijholt@gmail.com",
    license="MIT",
    packages=find_packages(),
    package_data={
        'pfapack': ['libcpfapack.so'],
    },
    include_package_data=True,
    install_requires=install_requires,
    extras_require=extras_require,
    zip_safe=False,
)
