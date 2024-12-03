# pfapack Python package
import os
import sys
from pathlib import Path

from pfapack._version import __version__


def _enable_sharedlib_loading():
    """Enable loading of shared libraries from the package directory."""
    if sys.platform == "win32":
        # Get the directory containing the shared libraries
        pkg_dir = Path(__file__).parent
        if hasattr(os, "add_dll_directory"):
            os.add_dll_directory(str(pkg_dir))


# Call this before loading the shared libraries
_enable_sharedlib_loading()

__all__ = ["__version__"]
