--------
Overview
--------

This folder contains a number of examplary programs demonstrating the
use of the PFAPACK library from C(++). The examples only work with
fairly small matrices (4x4) to remain on a tutorial level, but are
easily extended to arbitrarily sized matrices.

There are examples accessing the raw Fortran77 interface (example*_f77.*),
as well as examples using the more convenient C interface
(example*_c.c, using the struct definition of complex numbers,
example*_c99.c, using C99 complex numbers, and example*_c++.cc
using C++ complex numbers),

-----------
Compilation
-----------

Modify the makefile to use your favorite compiler and LAPACK library
and compile all examples using "make all".

Note that CC must refer to a ANSI-C++ compatible compiler, ANSIC to
a ANSI-C (C89) compatible compiler, and C99 to a C99 compatible compiler.
If one of these compilers is not available, remove the corresponding
examples from the makefile.

---------------------------------
Overview of the example programs:
---------------------------------

example1_*:
 Compute the pfaffian of a 4x4 matrix

example2_c.c:
 Compute the pfaffian of a 4x4 matrix using the Pfaffian routine that avoids
 floating point overflow

example3_*:
 Compute the Pfaffian of a band matrix

example5_c.c:
 Compute the L T L^T decomposition of a dense matrix

example5_c.c:
 Compute the tridiagonal form of a dense matrix under orthogonal similarity

example6_c.c:
 Compute the tridiagonal form of a band matrix under orthogonal similarity
