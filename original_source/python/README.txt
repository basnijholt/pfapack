--------
Overview
--------

This directory contains a Python implementation of the skew-symmetric
LTL^T decomposition, the Householder tridiagonalization and the
computation of the Pfaffian based on these algorithms. In contrast to
the Fortran version, this implementation always acts on the full
matrix, rather than only the upper or lower triangle. The Python
implementation is not optimized for performance, but rather serves as
a demonstration of the algorithms.

There are several ways of interfacing Fortran code with Python. If
this is not an option or speed is not critical, the native Python
implementation is a viable alternative.

For real matrices, the Pfaffian can also be computed using the
Hessenberg or real Schur decomposition which do not make use of the
skew-symmetry of the problem. However, the Python implementation of
the more efficient algorithms here is faster than those approaches.

The routines are documented in place.

-----
Files
-----

pfaffian.py:
 Implementation of the skew-symmetric LTL^T decomposition, the
 skew-symmetric Householder tridiagonalization, the corresponding
 Pfaffian routines (based on a partial tridiagonal form) as well
 as a Pfaffian routine based on the real Schur decomposition.

test_pfaffian.py:
 Tests for the decomposition and Pfaffian routines. Can be called
 using 'nosetests'.
