!Example4: compute the UTU^T decomposition of a skew-symmetric matrix
!          uses the Fortran95 interface

PROGRAM Example4
  USE F95_PFAPACK
  USE print_matrix
  IMPLICIT NONE
  COMPLEX(KIND(1D0)) :: A(4,4), T(4,4), U(4,4), tmp(4)
  INTEGER :: I, J, IPIV(4)

  A=0D0

  !Set only the upper triangle of the matrix

  A(1,2)=1D0
  A(1,3)=2D0
  A(1,4)=3D0
  A(2,3)=(0,8D0)
  A(2,4)=5D0
  A(3,4)=6D0

  WRITE (*,*) "The example matrix A"
  CALL print_skew_mat(A, 'U')

  !Factorize the matrix
  CALL SKTRF(A, IPIV, UPLO='U')

  !Extract the tridiagonal part
  T=0D0
  DO I=1,3
     T(I,I+1)=A(I,I+1)
     T(I+1,I)=-A(I,I+1)
  END DO

  WRITE (*,*) "The matrix P A P^T is factorized as U * T * U^T"

  WRITE (*,*) "with T="
  CALL print_mat(T)

  !Extract the lower unit triangular matrix L

  U=0D0
  DO J=1,4
     U(J,J)=1D0
     DO I=1,J-2
        U(I,J-1)=A(I,J)
     END DO
  END DO

  WRITE (*,*) "and U="
  CALL print_mat(U)

  !Now demonstrate that the factorization is correct:
  A=MATMUL(U,T)
  T=MATMUL(A, TRANSPOSE(U))

  !apply the inverse permutation to U*T*U^T
  !This is done by doing the interchanges described in IPIV from 1 to N=4
  DO I=1,4
     !interchange rows first
     tmp=T(I,:)
     T(I,:)=T(IPIV(I),:)
     T(IPIV(I),:)=tmp

     !then interchange columns
     tmp=T(:,I)
     T(:,I)=T(:,IPIV(I))
     T(:,IPIV(I))=tmp
  END DO

  WRITE (*,*) "Sanity check: computing P^T* U * T * U^T * P should give the original matrix"
  CALL print_mat(T)

END PROGRAM Example4
