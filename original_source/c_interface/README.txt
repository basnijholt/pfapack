--------
OVERVIEW
--------

The C interface allows you to use the Fortran library and takes care of
allocating the necessary workspace. It is similar to the Fortran95
interface in terms of ease of use.

------------
INSTALLATION
------------

Compiling the C interface requires a ANSI C or C++ compliant compiler.
It may require adapting to the routine naming conventions of your
Fortran compiler (see below in "Mixing Fortran and C"); this is done
by editing the macro fortran_name() in fortran.h.

Change the compiler and options in the makefile and type "make". This will
build the library libcpfapack.a.

--------------------
MIXING FORTRAN AND C
--------------------

The C interface is based on the FORTRAN implementation. Hence it is
necessary to call Fortran routines from C(++). In most cases, this works
without large problems, especially if the C and the Fortran compiler are from
the same collection (e.g. the GNU compiler collection). There are a few
particular things to be aware of:

-Mixing C with Fortran77 vs. Fortran95
 From C, it is only possible to use the Fortran77 routines, as the
 Fortran95 interface uses assumed-shape arrays.

-Fortran routine naming
 Fortran is case-insensitive. This leaves some ambiguity how the
 routine names are translated when they are linked against the calling
 C code, because for C, it does matter. Typically, the routine names
 are all lower-case, for example on Linux systems (although all
 uppercase is in principle also possible). In addition, the Fortran
 compiler usually adds an underscore to the routine name (or even two,
 or something different - this is system-dependent).

 In order to allow for different naming schemes, fortran.h defines
 a macro fortran_name(lcname,UCNAME) which takes the routine name
 both in all lowercase (lcname) and all uppercase (UCNAME), and
 from that input form the correct FORTRAN name. By default, the
 lower-case routine name is chosen and an underscore appended
 (which is valid for most systems). If your system uses a different
 convention, you have to modify that macro.

 fortran_pfapack.h then defines macros such as PFAPACK_sskpfa that make
 use of the appropriate macro fortran_name() to allow for a unfied
 way of calling the Fortran routines.

-Complex numbers
 Complex numbers in Fortran are, as far as their memory layout is
 concerned, a pair of real numbers stored conecutively in memory.
 This can be achieved from C in several ways, and to this end
 fortran.h (which is alos included from fortran_pfapack.h and pfapack.h)
 defines the types floatcmplx and doublecmplx. These types are
 defined differently:

  * when CPLUSPLUS_COMPLEX is defined, floatcmplx and doublecmplx are
    std::complex<float> and std::complex<double>

  * when C99_COMPLEX is defined, floatcmplx and doublecmplx are
    defined as float complex and double complex (the C99 types)

  * when STRUCT_COMPLEX is defined, complex numbers are defined as
    	    typedef struct {
	    	double re, im;
	    } doublecmplx; and
    equivalently for floatcmplx

 If none of these macros are defined, fortran.h defines
 CPLUSPLUS_COMPLEX if compiled with a C++ compiler, C99_COMPLEX if
 compiled with a C99 compiler (which is not also a C++ compiler), and
 STRUCT_COMPLEX else.

-------------
MATRIX FORMAT
-------------

Since the C interface is based on the Fortran implementation,
it expects its input, especially matrices, to be in Fortran
format. "Fortran format" here means that data must be stored in
a contiguous block of memory, which is interpreted for matrices
in column-major order.

Vectors
*******

For vectors this is not a big problem, all of the examples

 double vec1[10];
 double *vec2 = (double *)malloc(sizeof(double) * 10);
 double *vec3 = new double[10]; //only C++
 std::vector<double> vec4(10); //only C++

are valid. The only difference is that Fortran array numbering starts
with 1, wheras C numbering starts with 0. Hence, V(i) in Fortran
corresponds to V[i-1] in C, if the same number i is used in both cases.

Matrices
********

For matrices, the order poses some subtle issues. Fortran matrices are
stored in column-major order, i.e. a 2x3 matrix A is stored as a
consecutive list of numbers { A(1,1), A(2,1), A(1,2), A(2,2), A(1,3),
A(2,3) }. In contrast, C arrays are stored in row-major ordering,
such that for a 2D array m[2][3] we find { m[0][0], m[0][1], m[0][2],
m[1][0], m[1][1], m[1][2] }. Hence, the role of row and column are
interchanged, so that a Fortran routine would interpret a 2D C array
as the transposed matrix.

One way to take this into account is by defining

 double A[2][3]

and accessing the (Fortran) element A(i,j) by A[j-1][i-1].

Another way for making an MxN matrix is to first allocate a contiguous
piece of memory, e.g. by

 double *mat2 = (double *)malloc(sizeof(double) * M*N);
 double *mat3 = new double[M*N]; //only C++
 std::vector<double> mat4(M*N); //only C++

and then accessing element (i,j) through the helper function/macro
dense_fortran defined in fortran_pfapack.h

 mat[dense_fortran(i, j, M)] = ...

where now i and j are counted from 1, and M is the number of rows
in the matrix.

Banded matrices
***************

The storage format of banded matrices is also conveniently constructed
using helper functions. In order to store a NxN banded matrices
with KD sub- or superdiagonals, we need to allocate a memory block
of size (KD+1)*N, e.g.

 double *mat = (double *)malloc(sizeof(double) * (KD+1)*N);

and then fill in the values. If the upper triangle is stored,

 mat[bandupper_fortran(i, j, N, KD)] = ...

if the lower triangle is stored,

 mat[bandlower_fortran(i, j, N, KD)] = ...

where i and j are again counted from 1.

The bandupper_fortran() and bandlower_fortran() functions check if the
position (i,j) is inside the band if compiled with a C++ or C99 compiler.

-----
USAGE
-----

Passing vectors or matrices
***************************

Vectors and matrices are passed (both in the raw Fortran77 as well as in the
C interface) by a pointer to the contiguous memory block. For example,

 double mat1[2][3];
 double *mat2 = (double *)malloc(sizeof(double) * M*N);
 double *mat3 = new double[M*N]; //only C++
 std::vector<double> mat4(M*N); //only C++

 function_call(&mat1[0][0], mat2, mat3, &mat4[0]);

Using the Fortran77 interface directly from C
*********************************************

The Fortran77 interface can also be used directly. For this, it
is necessary to include fortran_pfapack.h. The PFAPACK routines
are then named PFAPACK_ssktrf, etc. Note that Fortran routines take
parameters by reference, rather than value. This means that instead of
passing, for example, an integer, one needs to pass a pointer to that integer.
(Vectors and matrices are passed as a pointer already, nothing changes there.)

A program must be linked with libpfapack.a, as well as LAPACK and BLAS
libraries, and possibly Fortran related libraries (system-specific).

Example programs are in EXAMPLES/

Using the C interface
*********************

The C interface provides access to the PFAPACK functions
XSKTRF, XSKTRD, XSKBTRD, XSKPFA, XSKPF10, XSKBPFA, and XSKBPF10
(with X=S, D, C, or Z). The corresponding C interface
functions are sktrf_x, sktrd_x, skbtrd_x, skpfa_x, skpf10_x, skbpfa_x,
and skbpf10_x (with x=s, d, c, or z). If compiled with a C++ compiler,
overloading allows to call the functions without trailing _x.

The C interface takes care of allocating all necessary workspace.
In functionality, it is therefore similar to the Fortran95 interface.

Errors are indicated by a return value < 0.

A program must be linked with libcpfapack.a and libpfapack.a, as well
as LAPACK and BLAS libraries, and possibly Fortran related libraries
(system-specific).

All of the C interface functions are documented in-place. Several examples
can be found in EXAMPLES.
