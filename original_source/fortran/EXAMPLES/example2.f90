!Fortran 95 interface example:
!compute the Pfaffian of two different 4x4 matrix
!The computation done here is the same as in example1.f
!but using the Fortran95 interface of PFAPACK

PROGRAM Example2
  USE F95_PFAPACK
  IMPLICIT NONE
  REAL(KIND(1D0)) ::A(4,4)
  REAL(KIND(1D0)) :: PFAFF

  WRITE (*,*) "Example program: compute the Pfaffian of two different 4x4 matrices, using the Fortran95 interface"

  A=0D0

!The two matrices can be stored in one dense matrix,
!in the upper and lower triangle, respectively

!Matrix 1: sparse matrix
!Set only the upper triangle

  A(1,3)=1D0
  A(2,4)=1D0

!Matrix 2: dense matrix
!Set only the lower triangle

  A(2,1)=-1D0
  A(3,1)=-2D0
  A(4,1)=-3D0
  A(3,2)=-4D0
  A(4,2)=-5D0
  A(4,3)=-6D0

!First matrix 1
  CALL SKPFA(A, PFAFF)

  WRITE (*,*) "Pfaffian of example matrix 1: ", PFAFF
  WRITE (*,*) "(The result should be -1.0)"

!Then matrix 2
  CALL SKPFA(A, PFAFF, UPLO = 'L')

  WRITE (*,*) "Pfaffian of example matrix 2: ", PFAFF
  WRITE (*,*) "(The result should be 8.0)"


END PROGRAM Example2
