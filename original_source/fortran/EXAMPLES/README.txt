--------
Overview
--------

This folder contains a number of examplary programs demonstrating the
use of the PFAPACK library. The examples only work with fairly small
matrices (4x4) to remain on a tutorial level, but are easily extended
to arbitrarily sized matrices.

example1.f and example1_ws.f use the Fortran77 interface of PFAPACK,
the remaining examples the Fortran95 interface. Examples 1-5 deal with
dense matrices, 6-7 with banded matrices.

-----------
Compilation
-----------

Modify the makefile to use your favorite compiler and LAPACK library
and compile all examples using "make all"

---------------------------------
Overview of the example programs:
---------------------------------

example1.f:
 Compute the Pfaffian of two different 4x4 matrices using the Fortran77
 interface of PFAPACK

example1_ws.f:
 Demonstration of how a workspace query for the computation of
 example1.f should look like

example2.f90:
Identical to example1.f, but using the Fortran95 interface

example3.f90:
 Compute the tridiagonal form of a skew-symmetric matrix under
 orthogonal similarity. This example uses the Householder-based
 algorithm for trdiagonalizing a matrix

example4.f90:
 Compute the LTL^T decomposition of a skew-symmetric matrix using the
 Parlett-Reid algorithm and the lower triangle of the matrix

example4_upper.f90:
 Similar to example4.f90, but using the upper triangle.  Consequently,
 an UTU^T decomposition is computed.

example5.f90:
 Compute the Pfaffian manually from a LTL^T decomposition. Doing so
 might be useful if only the sign of the Pfaffian is needed, i.e. when
 the values as returned by SKPFA and SKPF10 are not needed.

example6.f90:
 Compute the Pfaffian of banded skew-symmetric matrices

example7.f90:
 Compute the tridiagonal form of a banded skew-symmetric matrix under
 orthogonal similarity.
