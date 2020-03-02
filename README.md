# `pfapack`: Efficient numerical computation of the Pfaffian for dense and banded skew-symmetric matrices

Code and algorithms are taken from [arXiv:1102.3440](https://arxiv.org/abs/1102.3440) which is authored by [Michael Wimmer](https://michaelwimmer.org/).

[![license](https://img.shields.io/github/license/basnijholt/pfapack)](https://github.com/basnijholt/pfapack/blob/master/LICENSE)
[![tests](https://github.com/basnijholt/pfapack/workflows/tests/badge.svg)](https://github.com/basnijholt/pfapack/actions?query=workflow%3Atests)
[![codecov](https://img.shields.io/codecov/c/github/basnijholt/pfapack)](https://codecov.io/gh/basnijholt/pfapack)
[![docs](https://img.shields.io/readthedocs/pfapack)](https://pfapack.readthedocs.io)
[![version](https://img.shields.io/pypi/v/pfapack)](https://pypi.org/project/pfapack/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/pfapack)](https://pypi.org/project/pfapack/)

### Install
```bash
pip install pfapack
```

## Usage
```python
from pfapack import pfaffian as pf
import numpy.matlib

# first real matrices
A = numpy.matlib.rand(100, 100)
A = A - A.T
pfa1 = pf.pfaffian(A)
pfa2 = pf.pfaffian(A, method="H")
pfa3 = pf.pfaffian_schur(A)

print(pfa1, pfa2, pfa3)
```

## License
MIT License

## Contributions
- Bas Nijholt
- [Michael Wimmer (author of the algorithms)](https://arxiv.org/abs/1102.3440)
