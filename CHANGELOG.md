# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.1] - 2026-06-10

### Changed
- Windows wheels are now repaired with `delvewheel`, bundling OpenBLAS and the
  MinGW runtime DLLs. Installing from PyPI no longer requires MSYS2; only
  building from source on Windows does. CI verifies the wheel is
  self-contained by running the test suite with MSYS2 hidden.
  ([#55](https://github.com/basnijholt/pfapack/pull/55))

## [1.1.0] - 2026-06-10

### Fixed
- The skew-symmetry input check in `pfapack.pfaffian` was broken: it took
  `.max()` before `abs()`, so matrices violating skew-symmetry with a purely
  negative deviation passed silently and produced wrong Pfaffians, and for
  complex matrices the comparison was lexicographic rather than by magnitude.
  All public functions now share a correct check with a tolerance relative to
  the matrix scale; non-skew-symmetric input raises `AssertionError`.
  Identified in the [jjgoings/pfapack](https://github.com/jjgoings/pfapack)
  fork by Joshua Goings. ([#52](https://github.com/basnijholt/pfapack/pull/52))

### Changed
- Linux wheels are now built on the `manylinux_2_28` image; the minimum glibc
  rises from 2.17 to 2.28 so that current SciPy releases install alongside
  pfapack. ([#53](https://github.com/basnijholt/pfapack/pull/53))
- Intel macOS wheels are now built on the `macos-15-intel` runner (GitHub
  retired the macOS 13 runners), raising the minimum macOS version for Intel
  wheels accordingly. ([#53](https://github.com/basnijholt/pfapack/pull/53))
- Windows builds use the MSYS2 MinGW64 toolchain; the sdist, wheel, and test
  pipelines were repaired after an extended period of CI breakage.
  ([#53](https://github.com/basnijholt/pfapack/pull/53),
  [#34](https://github.com/basnijholt/pfapack/pull/34))

## [1.0.2] - 2024-12-04

### Fixed
- Version string of 1.0.1 was not set properly.

### Removed
- Dropped Python 3.9 support; pfapack now requires Python ≥ 3.10.
  ([#33](https://github.com/basnijholt/pfapack/pull/33))

### Changed
- Use `numpy.random` instead of the deprecated `numpy.matlib`.
  ([#32](https://github.com/basnijholt/pfapack/pull/32))

## [1.0.1] - 2024-12-04

### Fixed
- Broken macOS and Windows builds of 1.0.0.

## [1.0.0] - 2024-12-03

### Changed
- Switched to `meson-python` for building the C and Fortran code, with
  `pyproject.toml` replacing `setup.py`.
- Wheels are built with `pypa/cibuildwheel` for Linux, macOS, and Windows,
  with OpenBLAS on Linux/Windows and the Accelerate framework on macOS.

### Added
- The original PFAPACK Fortran/C source from Michael Wimmer's website is
  included in the repository for transparency and reproducibility.
- Support for Python up to 3.13.

### Removed
- Dropped support for Python 3.8 and below, and for 32-bit builds.

## [0.3.1] - 2021-09-17

### Fixed
- Packaging fix for the 0.3.0 release.

## [0.3.0] - 2021-09-17

### Added
- Windows support and Python 3.9 testing.

### Removed
- Dropped Python 3.6 support.

## [0.2.2]

### Fixed
- Fix for integer valued arrays.

[1.1.1]: https://github.com/basnijholt/pfapack/releases/tag/v1.1.1
[1.1.0]: https://github.com/basnijholt/pfapack/releases/tag/v1.1.0
[1.0.2]: https://github.com/basnijholt/pfapack/releases/tag/v1.0.2
[1.0.1]: https://github.com/basnijholt/pfapack/releases/tag/v1.0.1
[1.0.0]: https://github.com/basnijholt/pfapack/releases/tag/v1.0.0
[0.3.1]: https://github.com/basnijholt/pfapack/releases/tag/v0.3.1
[0.3.0]: https://github.com/basnijholt/pfapack/releases/tag/v0.3.0
