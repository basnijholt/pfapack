--------
CONTENTS
--------

This package contains routines to tridiagonalize skew-symmetric
matrices and compute their Pfaffian. It contains routines written in
FORTRAN that are optimized for speed and numerical stability, and can
deal with dense and banded matrices. It also contains code for python
matlab, and mathematica that is less optimized, but more intuitive.

The root of the package contains:

 README.txt 		this file
 LapackLicense		some notes regarding copyright
 fortran/	     	the implementation in FORTRAN
 fortran/TESTING/	extensive test suite for the FORTRAN
                        implementation
 fortran/EXAMPLES/	example programs in FORTRAN
 c_interface/		a convenient C interface for the FORTRAN routines
 c_interface/TESTING/	test suite for the C interface
 c_interface/EXAMPLES/	example programs in C and C++
 mathematica/		Mathematica implementation
 matlab/                Matlab implementation
 python/		Python implementation + examples

------
CITING
------

If you use this library, please cite the accompanying paper:

M. Wimmer.
"Algorithm 923: Efficient Numerical Computation of the Pfaffian for
Dense and Banded Skew-Symmetric Matrices".
ACM Trans. Math. Software 38, 30 (2012).

------------
AVAILABILITY
------------

The most recent version of this library, including possible bug fixes
and enhancements, will be available at

http://www.michaelwimmer.org/downloads.html

------------
REQUIREMENTS
------------

-FORTRAN implementation:
 In order to use the FORTRAN code, one also needs a LAPACK and BLAS
 library (linked, e.g. as -llapack -lblas). Using an optimized LAPACK
 and BLAS will also benefit the speed of this code.

 The FORTRAN implementation uses the LAPACK/BLAS functions LSAME and
 XERBLA. They are sometimes not included in some machine-specific
 precompiled LAPACK and BLAS libraries. In this case you need to
 obtain them from www.netlib.org/blas or www.netlib.org/lapack (files
 lsame.f and xerbla.f).  For more instructions see "Installation"
 below.

-C interface:
 The C interface requires the Fortran implementation to be built
 successfully. It can be compiled with either a C or C++-compiler.

-Mathematica implementation:
 No extra requirements except Mathematica (tested on versions 7 and 8,
 but should also run on earlier versions).

-Matlab implementation:
 No extra requirements except Matlab (tested on version 2010b, but
 should also work on previous versions)

-python implementation:
 The python implementation requires numpy and scipy.

-----------
INSTALLATON
-----------

-FORTRAN implementation:
 In order to compile the FORTRAN library, change fortran/makefile to use
 your favorite compiler and options, and then type "make". Copy
 libpfapack.a to wherever you need it, along with the
 Fortran95 interface modules, if needed.

 If your LAPACK and BLAS library does not contain LSAME and XERBLA
 (this is the case if the example programs fail to compile with the
 linker complaining about missing symbols such as 'lsame_' and
 'xerbla_), obtain lsame.f and xerbla.f
 (e.g. www.netlib.org/blas/lsame.f and www.netlib.org/blas/xerbla.f)
 and put them in fortran/. Finally, add "lsame.o xerbla.o" to the
 OBJECTS, and compile with "make".

-C interface:
 In order to compile the C interface, change c_interface/makefile
 to use your favorite compiler and options, and then type "make".
 Programs using the C interface need pfapack.h, fortran_pfapack.h
 and fortran.h, as well as libcpfapack.a along with the Fortran
 library libpfapack.a.

 Interfacing Fortran from C usually works out of the box, but may require
 some additional input from the user, especially with regard
 to Fortran function names. Please refer to the corresponding
 c_interface/README.txt.

-Mathematica, Matlab and python implementation:
 No particular installation steps required.

----------------
EXAMPLES AND USE
----------------

All implementations are documented extensively in-place in the
respective files.

-FORTRAN implementation:
 The directory fortran/ contains an extensive test-suite checking the
 implementation in the directory TESTING/. Examples are in the
 directory EXAMPLES/. In all cases modify makefile to use your
 favorite compiler and options, as well as your LAPACK and BLAS
 library.

 For dense matrices, the Householder transformations used in the
 skew-symmetric tridiagonalization are stored in the same format as
 used by the symmetric and Hermitian LAPACK tridiagonalization
 routines. Hence, the LAPACK routines xORGTR/xUNGTR can be used to
 form the transformation matrix explicitely, and xORMTR/xUNMTR can be
 used to multiply another matrix with the transformation, without
 forming it explicitely.

-C interface:
 The directory c_interface contains an extensive test-suite that
 checks the C interface for consistency against the Fortran interface.
 Example programs are included in EXAMPLES/, both showing how to
 access the raw FORTRAN77 interface, and the more comfortable
 C interface. In all cases modify makefile to use your
 favorite compiler and options, as well as your LAPACK and BLAS
 library.

 The C interface expects matrices to be given in FORTRAN format. For
 details see c_interface/README.txt.

-Mathematica implementation:
 The file pfaffian.nb contains both the definitions of the routines as
 well as examples. Details are given in mathematica/README.txt

-Matlab implementation:
 The function files must be in a path where Matlab searches for
 external functions (for example, the working directory). Details
 are given in matlab/README.txt

-python implementation:
 Import the implementation to your python program to use the routines.
 A test program is included, the tests can be run using nosetests.
 Details are given in python/README.txt
