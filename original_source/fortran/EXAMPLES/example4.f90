!Example4: compute the LTL^T decomposition of a skew-symmetric matrix
!          uses the Fortran95 interface

PROGRAM Example4
  USE F95_PFAPACK
  USE print_matrix
  IMPLICIT NONE
  COMPLEX(KIND(1D0)) :: A(4,4), T(4,4), L(4,4), tmp(4)
  INTEGER :: I, J, IPIV(4)

  A=0D0

  !Set only the lower triangle of the matrix

  A(2,1)=-1D0
  A(3,1)=-2D0
  A(4,1)=-3D0
  A(3,2)=-(0,8D0)
  A(4,2)=-5D0
  A(4,3)=-6D0

  WRITE (*,*) "The example matrix A"
  CALL print_skew_mat(A, 'L')

  !Factorize the matrix
  CALL SKTRF(A, IPIV, UPLO='L')

  !Extract the tridiagonal part
  T=0D0
  DO I=1,3
     T(I,I+1)=-A(I+1,I)
     T(I+1,I)=A(I+1,I)
  END DO

  WRITE (*,*) "The matrix P A P^T is factorized as L * T * L^T"

  WRITE (*,*) "with T="
  CALL print_mat(T)

  !Extract the lower unit triangular matrix L

  L=0D0
  DO J=1,4
     L(J,J)=1D0
     DO I=J+2,4
        L(I,J+1)=A(I,J)
     END DO
  END DO

  WRITE (*,*) "and L="
  CALL print_mat(L)

  !Now demonstrate that the factorization is correct:
  A=MATMUL(L,T)
  T=MATMUL(A, TRANSPOSE(L))

  !apply the inverse permutation to L*T*L^T
  !This is done by doing the interchanges described in IPIV from N=4 to 1
  DO I=4,1,-1
     !interchange rows first
     tmp=T(I,:)
     T(I,:)=T(IPIV(I),:)
     T(IPIV(I),:)=tmp

     !then interchange columns
     tmp=T(:,I)
     T(:,I)=T(:,IPIV(I))
     T(:,IPIV(I))=tmp
  END DO

  WRITE (*,*) "Sanity check: computing P^T* L * T * L^T * P should give the original matrix"
  CALL print_mat(T)

END PROGRAM Example4
