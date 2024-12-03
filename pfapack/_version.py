#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

# Is set during `onbuild` if `pip install pfapack` is used
__version__ = ""

if not __version__:
    try:
        import versioningit
    except ImportError:  # pragma: no cover
        import importlib.metadata

        __version__ = importlib.metadata.version("pfapack")
    else:
        PROJECT_DIR = Path(__file__).parent.parent
        __version__ = versioningit.get_version(project_dir=PROJECT_DIR)


if __name__ == "__main__":
    print(__version__)
