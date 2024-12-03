*     Example1: compute the Pfaffian of two different 4x4 matrix
*               Uses the Fortran77 interface
      PROGRAM Example1

      INTEGER LWORK
*     The parameter LWORK sets the workspace size In principle, in this
*     Fortran77 implementation it should be determined separately
*     (Fortran77 does not allow for dynamic memory allocation, hence the
*      workspace size must be fixed at compile time).
*     The program example1-ws.f shows how to do this. Here, LWORK is set
*     to the minimal work space size.
      PARAMETER ( LWORK = 1 )
      DOUBLE PRECISION A(4,4)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION PFAFF
      INTEGER IWORK(4)
      INTEGER INFO

      WRITE (*,*)
     $     "Example program: compute the Pfaffian of two different 4x4 m
     $atrices"

      A=0D0

*     The two matrices can be stored in one dense matrix,
*     in the upper and lower triangle, respectively

*     Matrix 1: sparse matrix
*     Set only the upper triangle

      A(1,3)=1D0
      A(2,4)=1D0

*     Matrix 2: dense matrix
*     Set only the lower triangle

      A(2,1)=-1D0
      A(3,1)=-2D0
      A(4,1)=-3D0
      A(3,2)=-4D0
      A(4,2)=-5D0
      A(4,3)=-6D0

*     
      CALL DSKPFA('U', 'P', 4, A, 4, PFAFF, IWORK, WORK, LWORK, INFO)

      WRITE (*,*) "Pfaffian of example matrix 1: ", PFAFF
      WRITE (*,*) "(The result should be -1.0)"

      END
