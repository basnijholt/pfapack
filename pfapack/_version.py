#!/usr/bin/env python
try:
    # setuptools_scm is not a runtime dependency. This is used when pip
    # installing from a git repository.
    from setuptools_scm import get_version

    __version__ = get_version(root="..", relative_to=__file__)
except (ImportError, LookupError):
    from importlib.metadata import version, PackageNotFoundError

    try:
        __version__ = version("pfapack")
    except PackageNotFoundError:
        __version__ = "0.0.0+unknown"


if __name__ == "__main__":
    print(__version__)
