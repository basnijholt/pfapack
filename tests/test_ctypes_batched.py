"""Tests for the batched Pfaffian interface (pfapack.ctypes.pfaffian_batched)."""

import numpy as np
import pytest

from pfapack import pfaffian as pf
from pfapack.exceptions import InvalidDimensionError, InvalidParameterError

try:
    from pfapack.ctypes import pfaffian as cpfaffian
    from pfapack.ctypes import pfaffian_batched

    with_ctypes = True
except OSError:
    with_ctypes = False

pytestmark = pytest.mark.skipif(not with_ctypes, reason="C library not available")

EPS = 1e-11


def make_skew(batch_shape, n, complex_, seed=0):
    rng = np.random.default_rng(seed)
    a = rng.standard_normal(batch_shape + (n, n))
    if complex_:
        a = a + 1j * rng.standard_normal(batch_shape + (n, n))
    return a - np.swapaxes(a, -1, -2)


# N=6 and N=10 have odd N/2, exercising the (-1)^(N/2) sign correction
# of the internal transpose-elimination trick.
@pytest.mark.parametrize("n", [2, 4, 6, 8, 10])
@pytest.mark.parametrize("complex_", [False, True])
@pytest.mark.parametrize("method", ["P", "H"])
def test_matches_single_matrix(n, complex_, method):
    matrices = make_skew((7,), n, complex_, seed=n)
    result = pfaffian_batched(matrices, method=method)

    expected_c = np.array([cpfaffian(m, method=method) for m in matrices])
    expected_py = np.array([pf.pfaffian(m, method=method) for m in matrices])

    np.testing.assert_allclose(result, expected_c, rtol=EPS, atol=EPS)
    np.testing.assert_allclose(result, expected_py, rtol=EPS, atol=EPS)


@pytest.mark.parametrize("complex_", [False, True])
def test_uplo_lower(complex_):
    matrices = make_skew((5,), 6, complex_, seed=42)
    np.testing.assert_allclose(
        pfaffian_batched(matrices, uplo="L"),
        pfaffian_batched(matrices, uplo="U"),
        rtol=EPS,
        atol=EPS,
    )


def test_known_values():
    block = np.array([[0.0, 1.0], [-1.0, 0.0]])
    batch = np.array([block, 2 * block, -3 * block])
    np.testing.assert_allclose(
        pfaffian_batched(batch), [1.0, 2.0, -3.0], rtol=EPS, atol=EPS
    )


@pytest.mark.parametrize("n", [1, 3, 5])
@pytest.mark.parametrize("complex_", [False, True])
def test_odd_dimension_gives_zero(n, complex_):
    matrices = make_skew((4,), n, complex_, seed=n)
    result = pfaffian_batched(matrices)
    np.testing.assert_array_equal(result, np.zeros(4, dtype=result.dtype))
    expected_py = np.array([pf.pfaffian(m) for m in matrices])
    np.testing.assert_array_equal(result, expected_py)


def test_batch_shape_is_preserved():
    matrices = make_skew((2, 3), 4, complex_=False, seed=1)
    result = pfaffian_batched(matrices)
    assert result.shape == (2, 3)
    expected = np.array([[cpfaffian(m) for m in row] for row in matrices])
    np.testing.assert_allclose(result, expected, rtol=EPS, atol=EPS)


def test_single_matrix_input_gives_scalar():
    matrix = make_skew((), 4, complex_=False, seed=2)
    result = pfaffian_batched(matrix)
    assert np.ndim(result) == 0
    np.testing.assert_allclose(result, cpfaffian(matrix), rtol=EPS, atol=EPS)


def test_empty_batch():
    matrices = np.empty((0, 4, 4))
    result = pfaffian_batched(matrices)
    assert result.shape == (0,)
    assert result.dtype == np.float64


def test_zero_dimensional_matrices():
    matrices = np.empty((3, 0, 0))
    np.testing.assert_array_equal(pfaffian_batched(matrices), np.ones(3))


@pytest.mark.parametrize("complex_", [False, True])
def test_layout_resilience(complex_):
    matrices = make_skew((7,), 8, complex_, seed=3)
    reference = pfaffian_batched(matrices)

    # Fortran-ordered input
    matrices_f = np.asfortranarray(matrices)
    assert not matrices_f.flags["C_CONTIGUOUS"]
    np.testing.assert_allclose(
        pfaffian_batched(matrices_f), reference, rtol=EPS, atol=EPS
    )

    # Non-contiguous strided view
    padded = np.empty(matrices.shape + (2,), dtype=matrices.dtype)
    padded[..., 0] = matrices
    padded[..., 1] = 0
    view = padded[..., 0]
    assert not view.flags["C_CONTIGUOUS"]
    np.testing.assert_allclose(pfaffian_batched(view), reference, rtol=EPS, atol=EPS)


def test_input_is_not_modified():
    matrices = make_skew((4,), 6, complex_=True, seed=4)
    original = matrices.copy()
    pfaffian_batched(matrices)
    np.testing.assert_array_equal(matrices, original)


def test_dtype_upcasting():
    matrices64 = make_skew((5,), 4, complex_=False, seed=5)
    reference = pfaffian_batched(matrices64)

    result32 = pfaffian_batched(matrices64.astype(np.float32))
    assert result32.dtype == np.float64
    np.testing.assert_allclose(result32, reference, rtol=1e-5, atol=1e-5)

    ints = np.array([[[0, 1], [-1, 0]], [[0, -2], [2, 0]]])
    result_int = pfaffian_batched(ints)
    assert result_int.dtype == np.float64
    np.testing.assert_allclose(result_int, [1.0, -2.0], rtol=EPS, atol=EPS)

    cmplx = make_skew((5,), 4, complex_=True, seed=6)
    result64c = pfaffian_batched(cmplx.astype(np.complex64))
    assert result64c.dtype == np.complex128
    np.testing.assert_allclose(result64c, pfaffian_batched(cmplx), rtol=1e-4, atol=1e-4)

    result_list = pfaffian_batched([[[0.0, 1.0], [-1.0, 0.0]]])
    np.testing.assert_allclose(result_list, [1.0], rtol=EPS, atol=EPS)


def test_invalid_shapes_raise():
    with pytest.raises(InvalidDimensionError):
        pfaffian_batched(np.zeros(4))  # 1D
    with pytest.raises(InvalidDimensionError):
        pfaffian_batched(np.zeros((3, 4, 5)))  # not square


def test_invalid_parameters_raise():
    matrices = make_skew((2,), 4, complex_=False, seed=7)
    with pytest.raises(InvalidParameterError):
        pfaffian_batched(matrices, uplo="X")
    with pytest.raises(InvalidParameterError):
        pfaffian_batched(matrices, method="X")
