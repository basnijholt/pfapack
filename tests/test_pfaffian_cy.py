import numpy as np
import numpy.linalg
import numpy.matlib
import pytest

import sys  # isort:skip

sys.path.append("..")

from pfapack import pfaffian_cy as pf  # noqa isort:skip

try:
    from pfapack.ctypes import pfaffian as cpfaffian

    with_ctypes = True
except OSError:
    with_ctypes = False


EPS = 1e-12


def test_pfaffian():
    # Compare the output of the different Pfaffian routines
    # and compare to the determinant

    # first real matrices
    A = numpy.matlib.rand(100, 100)
    A = A - A.T

    pfa1 = pf.pfaffian(A)
    deta = numpy.linalg.det(A)

    assert abs((pfa1 ** 2 - deta) / deta) < EPS

    # then complex matrices
    A = numpy.matlib.rand(100, 100) + 1.0j * numpy.matlib.rand(100, 100)
    A = A - A.T

    pfa1 = pf.pfaffian(A)
    deta = numpy.linalg.det(A)

    assert abs((pfa1 ** 2 - deta) / deta) < EPS
