!Example7: Compute the tridagonal form of a banded skew-symmetric
!          matrix under an orthogonal similarity transform
!          Uses the Fortran95 interface

PROGRAM Example7
  USE F95_PFAPACK
  USE print_matrix
  IMPLICIT NONE
  REAL(KIND(1D0)) :: A(4,4), AB(3,4), Q(4,4), T(4,4)
  INTEGER KD, I, J

  WRITE (*,*) "Example program: tridiagonalize a banded skew-symmetric matrix"

!We consider matrices with two off-diagonals
  KD=2

!Matrix: dense banded matrix stored in lower band storage format
!The dense storage in A is only for conveniently printing it to the screen,
!AB is what matters for the computation

!The dense form
!Set only the lower triangle
  A=0D0

  A(2,1)=-1D0
  A(3,1)=-2D0
  A(3,2)=-4D0
  A(4,2)=-5D0
  A(4,3)=-6D0

!The banded storage (also lower triangle)
  AB=0D0
  DO I=1, 4
     DO J=1,I
        IF( I<=J+KD ) AB(1+I-J,J)=A(I,j)
     END DO
  END DO

  WRITE (*,*) "The original matrix A:"
  CALL print_skew_mat(A,'L')

  !Now factorize (including a computation of the transformation Q)
  CALL SKBTRD(AB, UPLO='L', Q=Q)

  !Extract the tridiagonal matrix
  T=0D0
  DO I=1,3
     T(I,I+1)=-AB(2,I)
     T(I+1,I)=AB(2,I)
  END DO

  WRITE (*,*) "The matrix A is factorized as Q * T * Q^T"

  WRITE (*,*) "with T="
  CALL print_mat(T)

  WRITE (*,*) "and Q="
  CALL print_mat(Q)

  !Now demonstrate that the factorization is correct:
  A=MATMUL(Q,T)
  T=MATMUL(A, TRANSPOSE(Q))

  WRITE (*,*) "Sanity check: computing Q * T * Q^T should give the original matrix"
  CALL print_mat(T)


END PROGRAM Example7
