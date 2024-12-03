!Example3: compute the tridiagonal form of a skew-symmetric
!matrix under orthogonal similarity and compute the corresponding
!transformation. Uses the Fortrann95 interface of PFAPACK (but
!the Fortran77 version of LAPACK for DORGTR)

PROGRAM Example3
  USE F95_PFAPACK
  USE print_matrix
  IMPLICIT NONE
  REAL(KIND(1D0)) :: A(4,4), T(4,4), TEMP(4,4), TAU(3)
  REAL(KIND(1D0)) :: WORK(3)
  INTEGER :: I, INFO
  EXTERNAL DORGTR

  A=0D0

  !Set only the upper triangle of the matrix

  A(1,2)=1D0
  A(1,3)=2D0
  A(1,4)=3D0
  A(2,3)=4D0
  A(2,4)=5D0
  A(3,4)=6D0

  WRITE (*,*) "The example matrix A"
  CALL print_skew_mat(A, 'U')

  !Factorize the matrix
  CALL SKTRD(A, TAU)

  !Extract the tridiagonal part
  T=0D0
  DO I=1,3
     T(I,I+1)=A(I,I+1)
     T(I+1,I)=-A(I,I+1)
  END DO

  WRITE (*,*) "The matrix A is factorized as Q * T * Q^T"

  WRITE (*,*) "with T="
  CALL print_mat(T)

  !the transformation matrix is stored in A
  !(Call to LAPACK uses only minumum workspace N-1)
  CALL DORGTR('U', 4, A, 4, TAU, WORK, 3, INFO)

  WRITE (*,*) "and Q="
  CALL print_mat(A)
  !Note: Q is a full matrix

  !Now demonstrate that the factorization is correct:
  TEMP=MATMUL(A,T)
  T=MATMUL(TEMP, TRANSPOSE(A))

  WRITE (*,*) "Sanity check: computing Q * T * Q^T should give the original matrix"
  CALL print_mat(T)

END PROGRAM Example3
