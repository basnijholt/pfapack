# `pfapack`: Efficient numerical computation of the Pfaffian for dense and banded skew-symmetric matrices

Code and algorithms are taken from [arXiv:1102.3440](https://arxiv.org/abs/1102.3440) which is authored by [Michael Wimmer](https://michaelwimmer.org/).

[![license](https://img.shields.io/github/license/basnijholt/pfapack)](https://github.com/basnijholt/pfapack/blob/master/LICENSE)
[![tests](https://github.com/basnijholt/pfapack/workflows/tests/badge.svg)](https://github.com/basnijholt/pfapack/actions?query=workflow%3Atests)
[![codecov](https://img.shields.io/codecov/c/github/basnijholt/pfapack)](https://codecov.io/gh/basnijholt/pfapack)
[![docs](https://img.shields.io/readthedocs/pfapack)](https://pfapack.readthedocs.io)
[![version](https://img.shields.io/pypi/v/pfapack)](https://pypi.org/project/pfapack/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/pfapack)](https://pypi.org/project/pfapack/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

### Install
Recommended way (because it includes faster C/FORTRAN bindings)
```bash
conda install -c conda-forge pfapack
```

Alternatively use
```bash
pip install pfapack
```

## Usage
```python
from pfapack import pfaffian as pf
import numpy.matlib

A = numpy.matlib.rand(100, 100)
A = A - A.T
pfa1 = pf.pfaffian(A)
pfa2 = pf.pfaffian(A, method="H")
pfa3 = pf.pfaffian_schur(A)

print(pfa1, pfa2, pfa3)
```

If installed with `conda`, C/FORTRAN code is included with Python bindings, use it like:
```python
from pfapack.ctypes import pfaffian as cpf

pfa1 = cpf(A)
pfa2 = cpf(A, method="H")

print(pfa1, pfa2)
```

## Citing
If you have used `pfapack` in your research, please cite it using the following `bib` entry:
```
@article{wimmer2012algorithm,
  title={Efficient numerical computation of the pfaffian for dense and banded skew-symmetric matrices},
  author={Michael Wimmer},
  journal={ACM Transactions on Mathematical Software (TOMS)},
  volume={38},
  number={4},
  pages={1--17},
  year={2012},
  publisher={ACM New York, NY, USA}
}
```

## License
MIT License

## Contributions
- Bas Nijholt
- [Michael Wimmer (author of the algorithms)](https://arxiv.org/abs/1102.3440)
