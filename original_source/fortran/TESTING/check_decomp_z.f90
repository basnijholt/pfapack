!Functions for double complex

SUBROUTINE check_LTL_F77_z(UPLO, B, RESIDUAL, RELW)
  USE matrix_tools
  USE F77_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO 
  COMPLEX(KIND(1.0D0)) :: B(:,:)
  REAL(KIND(1.0D0)) :: RESIDUAL, RELW
  
  COMPLEX(KIND(1.0D0)), ALLOCATABLE :: A(:,:), Ac(:,:), U(:,:), T(:,:)
  COMPLEX(KIND(1.0D0)), ALLOCATABLE :: TEMP1(:,:), TEMP2(:,:), WORK(:)
  COMPLEX(KIND(1.0D0)) :: QUERY(1)
  INTEGER, ALLOCATABLE :: IPIV(:)
  INTEGER :: INFO, ISTAT, N, LWORK

  N = SIZE(B,2)

  !Allocate memory
  ALLOCATE(A(N,N), Ac(N,N), U(N,N), T(N,N), TEMP1(N,N), &
       TEMP2(N,N), IPIV(N), STAT = ISTAT)

  IF( ISTAT /= 0 ) THEN
     WRITE (*,*) "Ran out of memory in test!"
     STOP
  END IF

  !copy for later comparison
  A=B
  Ac=B

  !compute the factorization
  ! first workspace query
  CALL SKTRF(UPLO, "N", N, A, MAX(1,N), IPIV, QUERY, -1, INFO)
  LWORK = MAX(1, INT(QUERY(1)*RELW))

  ALLOCATE(WORK(LWORK), STAT=ISTAT)

  IF( ISTAT /= 0 ) THEN
     WRITE (*,*) "Ran out of memory in test!"
     STOP
  END IF

  CALL SKTRF(UPLO, "N", N, A, MAX(1,N), IPIV, WORK, MAX(1,LWORK), INFO)

  !Reconstruct the matrix
  CALL extract_unit_tri(UPLO, A, U)
  CALL extract_trid(UPLO, A, T)

  TEMP1=MATMUL(U,T)
  A=TRANSPOSE(U)
  TEMP2=MATMUL(TEMP1,A)

  CALL rowcol_invperm(UPLO, TEMP2, IPIV)

  RESIDUAL = resid(TEMP2, Ac)

  DEALLOCATE(A, Ac, U, T, TEMP1, TEMP2, IPIV, WORK)

END SUBROUTINE check_LTL_F77_z

SUBROUTINE check_LTL_F95_z(UPLO, B, RESIDUAL)
  USE matrix_tools
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO 
  COMPLEX(KIND(1.0D0)) :: B(:,:)
  REAL(KIND(1.0D0)) :: RESIDUAL
  
  COMPLEX(KIND(1.0D0)), ALLOCATABLE :: A(:,:), Ac(:,:), U(:,:), T(:,:)
  COMPLEX(KIND(1.0D0)), ALLOCATABLE :: TEMP1(:,:), TEMP2(:,:)
  INTEGER, ALLOCATABLE :: IPIV(:)
  INTEGER :: ISTAT, N

  N = SIZE(B,2)

  !Allocate memory
  ALLOCATE(A(N,N), Ac(N,N), U(N,N), T(N,N), TEMP1(N,N), &
       TEMP2(N,N), IPIV(N), STAT = ISTAT)

  IF( ISTAT /= 0 ) THEN
     WRITE (*,*) "Ran out of memory in test!"
     STOP
  END IF

  !copy for later comparison
  A=B
  Ac=B

  !compute the factorization
  CALL SKTRF(A, IPIV, UPLO = UPLO)

  !Reconstruct the matrix
  CALL extract_unit_tri(UPLO, A, U)
  CALL extract_trid(UPLO, A, T)

  TEMP1=MATMUL(U,T)
  A=TRANSPOSE(U)
  TEMP2=MATMUL(TEMP1,A)

  CALL rowcol_invperm(UPLO, TEMP2, IPIV)

  RESIDUAL = resid(TEMP2, Ac)

  DEALLOCATE(A, Ac, U, T, TEMP1, TEMP2, IPIV)

END SUBROUTINE check_LTL_F95_z

SUBROUTINE check_LTL_noblock_F77_z(UPLO, B, RESIDUAL)
  USE matrix_tools
  USE F77_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO 
  COMPLEX(KIND(1.0D0)) :: B(:,:)
  REAL(KIND(1.0D0)) :: RESIDUAL
  
  COMPLEX(KIND(1.0D0)), ALLOCATABLE :: A(:,:), Ac(:,:), U(:,:), T(:,:)
  COMPLEX(KIND(1.0D0)), ALLOCATABLE :: TEMP1(:,:), TEMP2(:,:)
  INTEGER, ALLOCATABLE :: IPIV(:)
  INTEGER :: INFO, ISTAT, N

  N = SIZE(B,2)

  !Allocate memory
  ALLOCATE(A(N,N), Ac(N,N), U(N,N), T(N,N), TEMP1(N,N), &
       TEMP2(N,N), IPIV(N), STAT = ISTAT)

  IF( ISTAT /= 0 ) THEN
     WRITE (*,*) "Ran out of memory in test!"
     STOP
  END IF

  !copy for later comparison
  A=B
  Ac=B

  !compute the factorization
  CALL SKTF2(UPLO, "N", N, A, MAX(1,N), IPIV, INFO)

  !Reconstruct the matrix
  CALL extract_unit_tri(UPLO, A, U)
  CALL extract_trid(UPLO, A, T)

  TEMP1=MATMUL(U,T)
  TEMP2=MATMUL(TEMP1,TRANSPOSE(U))

  CALL rowcol_invperm(UPLO, TEMP2, IPIV)

  RESIDUAL = resid(TEMP2, Ac)

  DEALLOCATE(A, Ac, U, T, TEMP1, TEMP2, IPIV)

END SUBROUTINE check_LTL_noblock_F77_z

SUBROUTINE check_LTL_noblock_F95_z(UPLO, B, RESIDUAL)
  USE matrix_tools
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO 
  COMPLEX(KIND(1.0D0)) :: B(:,:)
  REAL(KIND(1.0D0)) :: RESIDUAL
  
  COMPLEX(KIND(1.0D0)), ALLOCATABLE :: A(:,:), Ac(:,:), U(:,:), T(:,:)
  COMPLEX(KIND(1.0D0)), ALLOCATABLE :: TEMP1(:,:), TEMP2(:,:)
  INTEGER, ALLOCATABLE :: IPIV(:)
  INTEGER :: ISTAT, N

  N = SIZE(B,2)

  !Allocate memory
  ALLOCATE(A(N,N), Ac(N,N), U(N,N), T(N,N), TEMP1(N,N), &
       TEMP2(N,N), IPIV(N), STAT = ISTAT)

  IF( ISTAT /= 0 ) THEN
     WRITE (*,*) "Ran out of memory in test!"
     STOP
  END IF

  !copy for later comparison
  A=B
  Ac=B

  !compute the factorization
  CALL SKTF2(A, IPIV, UPLO = UPLO)

  !Reconstruct the matrix
  CALL extract_unit_tri(UPLO, A, U)
  CALL extract_trid(UPLO, A, T)

  TEMP1=MATMUL(U,T)
  A=TRANSPOSE(U)
  TEMP2=MATMUL(TEMP1,A)

  CALL rowcol_invperm(UPLO, TEMP2, IPIV)

  RESIDUAL = resid(TEMP2, Ac)

  DEALLOCATE(A, Ac, U, T, TEMP1, TEMP2, IPIV)

END SUBROUTINE check_LTL_noblock_F95_z

SUBROUTINE check_QTQ_F77_z(UPLO, B, RESIDUAL, RELW)
  USE matrix_tools
  USE F77_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO 
  COMPLEX(KIND(1.0D0)) :: B(:,:)
  REAL(KIND(1.0D0)) :: RESIDUAL, RELW
  
  COMPLEX(KIND(1.0D0)), ALLOCATABLE :: A(:,:), Ac(:,:), Q(:,:), T(:,:), TAU(:)
  REAL(KIND(1.0D0)), ALLOCATABLE :: E(:)
  COMPLEX(KIND(1.0D0)), ALLOCATABLE :: TEMP1(:,:), TEMP2(:,:), WORK(:)
  COMPLEX(KIND(1.0D0)) :: QUERY(1)
  INTEGER :: INFO, ISTAT, N, LWORK

  N = SIZE(B,2)

  !Allocate memory
  ALLOCATE(A(N,N), Ac(N,N), Q(N,N), T(N,N), TEMP1(N,N), &
       TEMP2(N,N), E(N-1), TAU(N-1), STAT = ISTAT)

  IF( ISTAT /= 0 ) THEN
     WRITE (*,*) "Ran out of memory in test!"
     STOP
  END IF

  !copy for later comparison
  A=B
  Ac=B

  !compute the factorization
  CALL SKTRD(UPLO, "N", N, A, MAX(1,N), E, TAU, QUERY, -1, INFO)
  LWORK = MAX(1, INT(QUERY(1)*RELW))

  ALLOCATE(WORK(LWORK), STAT=ISTAT)

  IF( ISTAT /= 0 ) THEN
     WRITE (*,*) "Ran out of memory in test!"
     STOP
  END IF

  CALL SKTRD(UPLO, "N", N, A, MAX(1,N), E, TAU, WORK, MAX(1,LWORK), INFO)

  !Reconstruct the matrix
  CALL extract_unitary(UPLO, A, TAU, Q)
  CALL extract_trid(UPLO, A, T)

  TEMP1=MATMUL(Q,T)
  A=TRANSPOSE(Q)
  TEMP2=MATMUL(TEMP1,A)

  RESIDUAL = resid(TEMP2, Ac)

  DEALLOCATE(A, Ac, Q, T, TEMP1, TEMP2, E, TAU, WORK)

END SUBROUTINE check_QTQ_F77_z

SUBROUTINE check_QTQ_F95_z(UPLO, B, RESIDUAL)
  USE matrix_tools
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO 
  COMPLEX(KIND(1.0D0)) :: B(:,:)
  REAL(KIND(1.0D0)) :: RESIDUAL
  
  COMPLEX(KIND(1.0D0)), ALLOCATABLE :: A(:,:), Ac(:,:), Q(:,:), T(:,:)
  COMPLEX(KIND(1.0D0)), ALLOCATABLE :: TAU(:)
  COMPLEX(KIND(1.0D0)), ALLOCATABLE :: TEMP1(:,:), TEMP2(:,:)
  INTEGER :: ISTAT, N

  N = SIZE(B,2)

  !Allocate memory
  ALLOCATE(A(N,N), Ac(N,N), Q(N,N), T(N,N), TEMP1(N,N), &
       TEMP2(N,N), TAU(N-1), STAT = ISTAT)

  IF( ISTAT /= 0 ) THEN
     WRITE (*,*) "Ran out of memory in test!"
     STOP
  END IF

  !copy for later comparison
  A=B
  Ac=B

  !compute the factorization
  CALL SKTRD(A, TAU, UPLO = UPLO)

  !Reconstruct the matrix
  CALL extract_unitary(UPLO, A, TAU, Q)
  CALL extract_trid(UPLO, A, T)

  TEMP1=MATMUL(Q,T)
  A=TRANSPOSE(Q)
  TEMP2=MATMUL(TEMP1,A)

  RESIDUAL = resid(TEMP2, Ac)

  DEALLOCATE(A, Ac, Q, T, TEMP1, TEMP2, TAU)

END SUBROUTINE check_QTQ_F95_z

SUBROUTINE check_QTQ_noblock_F77_z(UPLO, B, RESIDUAL)
  USE matrix_tools
  USE F77_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO 
  COMPLEX(KIND(1.0D0)) :: B(:,:)
  REAL(KIND(1.0D0)) :: RESIDUAL
  
  COMPLEX(KIND(1.0D0)), ALLOCATABLE :: A(:,:), Ac(:,:), Q(:,:), T(:,:)
  REAL(KIND(1.0D0)), ALLOCATABLE :: E(:)
  COMPLEX(KIND(1.0D0)), ALLOCATABLE :: TEMP1(:,:), TEMP2(:,:), TAU(:)
  INTEGER :: INFO, ISTAT, N

  N = SIZE(B,2)

  !Allocate memory
  ALLOCATE(A(N,N), Ac(N,N), Q(N,N), T(N,N), TEMP1(N,N), &
       TEMP2(N,N), E(N-1), TAU(N-1), STAT = ISTAT)

  IF( ISTAT /= 0 ) THEN
     WRITE (*,*) "Ran out of memory in test!"
     STOP
  END IF

  !copy for later comparison
  A=B
  Ac=B

  !compute the factorization
  CALL SKTD2(UPLO, "N", N, A, MAX(1,N), E, TAU, INFO)

  !Reconstruct the matrix
  CALL extract_unitary(UPLO, A, TAU, Q)
  CALL extract_trid(UPLO, A, T)

  TEMP1=MATMUL(Q,T)
  A=TRANSPOSE(Q)
  TEMP2=MATMUL(TEMP1,A)

  RESIDUAL = resid(TEMP2, Ac)

  DEALLOCATE(A, Ac, Q, T, TEMP1, TEMP2, E, TAU)

END SUBROUTINE check_QTQ_noblock_F77_z

SUBROUTINE check_QTQ_noblock_F95_z(UPLO, B, RESIDUAL)
  USE matrix_tools
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO 
  COMPLEX(KIND(1.0D0)) :: B(:,:)
  REAL(KIND(1.0D0)) :: RESIDUAL
  
  COMPLEX(KIND(1.0D0)), ALLOCATABLE :: A(:,:), Ac(:,:), Q(:,:), T(:,:)
  REAL(KIND(1.0D0)), ALLOCATABLE :: E(:)
  COMPLEX(KIND(1.0D0)), ALLOCATABLE :: TEMP1(:,:), TEMP2(:,:), TAU(:)
  INTEGER :: ISTAT, N

  N = SIZE(B,2)

  !Allocate memory
  ALLOCATE(A(N,N), Ac(N,N), Q(N,N), T(N,N), TEMP1(N,N), &
       TEMP2(N,N), E(N-1), TAU(N-1), STAT = ISTAT)

  IF( ISTAT /= 0 ) THEN
     WRITE (*,*) "Ran out of memory in test!"
     STOP
  END IF

  !copy for later comparison
  A=B
  Ac=B

  !compute the factorization
  CALL SKTD2(A, TAU, UPLO = UPLO)

  !Reconstruct the matrix
  CALL extract_unitary(UPLO, A, TAU, Q)
  CALL extract_trid(UPLO, A, T)

  TEMP1=MATMUL(Q,T)
  A=TRANSPOSE(Q)
  TEMP2=MATMUL(TEMP1,A)

  RESIDUAL = resid(TEMP2, Ac)

  DEALLOCATE(A, Ac, Q, T, TEMP1, TEMP2, E, TAU)

END SUBROUTINE check_QTQ_noblock_F95_z

SUBROUTINE check_QTQ_band_F77_z(UPLO, A, B, RESIDUAL)
  USE matrix_tools
  USE F77_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO 
  COMPLEX(KIND(1.0D0)) :: A(:,:), B(:,:)
  REAL(KIND(1.0D0)) :: RESIDUAL
  
  COMPLEX(KIND(1.0D0)), ALLOCATABLE :: AB(:,:), Q(:,:), T(:,:)
  REAL(KIND(1.0D0)), ALLOCATABLE :: E(:), RWORK(:)
  COMPLEX(KIND(1.0D0)), ALLOCATABLE :: TEMP1(:,:), TEMP2(:,:), WORK(:)
  COMPLEX(KIND(1.0D0)) :: DETQ
  INTEGER :: INFO, ISTAT, N, KD

  N = SIZE(A,2)
  KD = SIZE(A,1) - 1

  !Allocate memory
  ALLOCATE(AB(KD+1,N), Q(N,N), T(N,N), TEMP1(N,N), &
       TEMP2(N,N), E(N-1), RWORK(N), WORK(N), STAT = ISTAT)

  IF( ISTAT /= 0 ) THEN
     WRITE (*,*) "Ran out of memory in test!"
     STOP
  END IF

  !copy to leave A unchanged
  AB=A

  !compute the factorization
  CALL SKBTRD("V", UPLO, "N", N, KD, AB, MAX(1,SIZE(AB,1)), E, &
              DETQ, Q, MAX(1,N), WORK, RWORK, INFO)

  !Reconstruct the matrix
  CALL extract_trid_band(UPLO, AB, T)

  TEMP1=MATMUL(Q,T)
  TEMP2=MATMUL(TEMP1,TRANSPOSE(Q))

  RESIDUAL = resid(TEMP2, B)

  DEALLOCATE(AB, Q, T, TEMP1, TEMP2, E, WORK)

END SUBROUTINE check_QTQ_band_F77_z

SUBROUTINE check_QTQ_band_F95_z(UPLO, A, B, RESIDUAL)
  USE matrix_tools
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO 
  COMPLEX(KIND(1.0D0)) :: A(:,:), B(:,:)
  REAL(KIND(1.0D0)) :: RESIDUAL
  
  COMPLEX(KIND(1.0D0)), ALLOCATABLE :: AB(:,:), Q(:,:), T(:,:)
  COMPLEX(KIND(1.0D0)), ALLOCATABLE :: TEMP1(:,:), TEMP2(:,:)
  INTEGER :: ISTAT, N, KD

  N = SIZE(A,2)
  KD = SIZE(A,1) - 1

  !Allocate memory
  ALLOCATE(AB(KD+1,N), Q(N,N), T(N,N), TEMP1(N,N), &
       TEMP2(N,N), STAT = ISTAT)

  IF( ISTAT /= 0 ) THEN
     WRITE (*,*) "Ran out of memory in test!"
     STOP
  END IF

  !copy to leave A unchanged
  AB=A

  !compute the factorization
  CALL SKBTRD(AB, UPLO=UPLO, Q=Q)

  !Reconstruct the matrix
  CALL extract_trid_band(UPLO, AB, T)

  TEMP1=MATMUL(Q,T)
  TEMP2=MATMUL(TEMP1,TRANSPOSE(Q))

  RESIDUAL = resid(TEMP2, B)

  DEALLOCATE(AB, Q, T, TEMP1, TEMP2)

END SUBROUTINE check_QTQ_band_F95_z
