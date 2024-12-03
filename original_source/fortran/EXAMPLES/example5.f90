!Example5: compute the LTL^T decomposition of a skew-symmetric matrix
!          to compute the Pfaffian manually.
!          uses the Fortran95 interface

PROGRAM Example5
  USE F95_PFAPACK
  IMPLICIT NONE
  REAL(KIND(1D0)) :: A(4,4), PFAFF
  INTEGER :: I, IPIV(4)

  A=0D0

  !Set only the lower triangle of the matrix

  A(2,1)=-1D0
  A(3,1)=-2D0
  A(4,1)=-3D0
  A(3,2)=-4D0
  A(4,2)=-5D0
  A(4,3)=-6D0

  !Factorize the matrix to compute only a partial tridiagonal form
  CALL SKTRF(A, IPIV, UPLO='L', MODE='P')

  !Multiply every other entry on the sub-diagonal for the Pfaffian
  PFAFF=1D0
  DO I=1, 3, 2
     PFAFF=PFAFF*A(I+1,I)
  END DO

  !Compute the sign of the determinant of the permutation matrix
  !Every interchange in IPIV gives one "-"
  DO I=1, 4
     IF( IPIV(I) /= I ) PFAFF=-PFAFF
  END DO

  WRITE (*,*) "Result of the manual Pfaffian computation: ", PFAFF
  WRITE (*,*) "(The result should be 8.0)"

END PROGRAM Example5
