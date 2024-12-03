--------
Overview
--------

This directory contains a Matlab implementation of the skew-symmetric
LTL^T decomposition, the Householder tridiagonalization and the
computation of the Pfaffian based on these algorithms. In contrast to
the Fortran version, this implementation always acts on the full
matrix, rather than only the upper or lower triangle. The Matlab
implementation is not optimized for performance, but rather serves as
a demonstration of the algorithms.

With some effort, the fully optimized Fortran implementation can also
be accessed from Matlab. If this is not possible, or speed is not
absolutely critical, the stand-alone Matlab implementation is a viable
alternative. Compiling the matlab functions may increase the performance.

For real matrices, the Pfaffian can also be computed using the
Hessenberg decomposition. Since the Hessenberg decomposition is
implemented efficiently in Matlab (in contrast to the direct Matlab
implementation here), this approach is usually fastest for real
matrices (compilation may change this assessment). For complex
matrices there is no alternative and the algorithms here are the only
option.

The routines are documented in place.

-----
Files
-----

skew_LTL.m:
 Computes the LTL^T decomposition of a skew-symmetric matrix using the
 Parlett-Reid algorithm.

skew_tridiagonalize.m:
 Computes the tridiagonal form of a skew-symmetric matrix under 
 unitary congruence (Householder tridiagonalization).

pfaffian_LTL.m:
 Compute the Pfaffian of a skew-symmetric matrix. Based on the Parlett-Reid
 algorithm.

pfaffian_householder.m:
 Compute the Pfaffian of a skew-symmetric matrix. Based on the Householder
 tridiagonalization.

pfaffian_hessenberg.m:
 Compute the Pfaffian of a real skew-symmetric matrix. Based on the
 real Hessenberg decomposition.

test_pfaffian_real.m, test_pfaffian_complex.m:
 Application of the algorithms to example random matrix.
