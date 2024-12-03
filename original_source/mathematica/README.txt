--------
Overview
--------

This directory contains a Mathematica implementation of the
skew-symmetric LTL^T decomposition, the Householder tridiagonalization
and the computation of the Pfaffian based on these algorithms. In
contrast to the Fortran version, this implementation always acts on
the full matrix, rather than only the upper or lower triangle. The
Mathematica implementation is not optimized for performance, but rather
serves as a demonstration of the algorithms. If Mathematica is the
only option or speed not absolutely critical, this implementation
is nevertheless a viable option.

For real matrices, the Pfaffian can also be computed using the
Hessenberg decomposition. Since the Hessenberg decomposition is
implemented efficiently in Mathematica (in contrast to the direct
Mathematica implementation here), this approach is usually fastest for
real matrices. For complex matrices there is no alternative and the
algorithms here are the only option.

Documentation is in the Mathematica notebook.

-----
Files
-----

pfaffian.nb:
 Implementation of the algorithms together with examples.
