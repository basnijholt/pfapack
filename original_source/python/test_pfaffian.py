import pfaffian as pf
import numpy.linalg
import numpy.matlib

def test_pfaffian():
    #Compare the output of the different Pfaffian routines
    #and compare to the determinant

    #first real matrices
    A = numpy.matlib.rand(10,10)
    A = A-A.T

    pfa1 = pf.pfaffian(A)
    pfa2 = pf.pfaffian(A, method='H')
    pfa3 = pf.pfaffian_schur(A)
    deta = numpy.linalg.det(A)

    assert abs((pfa1-pfa2)/pfa1) < 1e-13
    assert abs((pfa1-pfa3)/pfa1) < 1e-13
    assert abs((pfa1**2-deta)/deta) < 1e-13

    #then complex matrices
    A = numpy.matlib.rand(10,10)+1.j*numpy.matlib.rand(10,10)
    A = A-A.T

    pfa1 = pf.pfaffian(A)
    pfa2 = pf.pfaffian(A, method='H')

    deta = numpy.linalg.det(A)

    assert abs((pfa1-pfa2)/pfa1) < 1e-13
    assert abs((pfa1**2-deta)/deta) < 1e-13


def test_pfaffian_integer():
    # Run pfapack on the test matrix reported by Amirhossein Ghodrati

    A=numpy.array([[ 0, -1,  1, -1,  1, -1],
                   [ 1,  0, -1, -1, -1, -1],
                   [-1,  1,  0, -1, -1, -1],
                   [ 1,  1,  1,  0, -1, -1],
                   [-1,  1,  1,  1,  0,  1],
                   [ 1,  1,  1,  1, -1,  0]])

    assert abs(pf.pfaffian_schur(A) - 1.0) < 1e-13
    assert abs(pf.pfaffian(A) - 1.0) < 1e-13
    assert abs(pf.pfaffian(A, method='H') - 1.0) < 1e-13


def test_decompositions():
    #Test the LTL^T and Householder decompositions

    #first real matrices
    A = numpy.matlib.rand(100,100)
    A = A-A.T

    T, L, P = pf.skew_LTL(A)

    assert numpy.linalg.norm(P*A*P.T-L*T*L.T)/numpy.linalg.norm(A) < 1e-13

    T, Q = pf.skew_tridiagonalize(A)

    assert numpy.linalg.norm(A-Q*T*Q.T)/numpy.linalg.norm(A) < 1e-13

    #then complex matrices
    A = numpy.matlib.rand(100,100)+1.0j*numpy.matlib.rand(100,100)
    A = A-A.T

    T, L, P = pf.skew_LTL(A)

    assert numpy.linalg.norm(P*A*P.T-L*T*L.T)/numpy.linalg.norm(A) < 1e-13

    T, Q = pf.skew_tridiagonalize(A)

    assert numpy.linalg.norm(A-Q*T*Q.T)/numpy.linalg.norm(A) < 1e-13
