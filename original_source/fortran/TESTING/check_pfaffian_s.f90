!Functions for single precision

SUBROUTINE check_pfaff_F77_s(UPLO, METHOD, B, RESIDUAL, RELW)
  USE matrix_tools
  USE F77_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, METHOD
  REAL(KIND(1.0E0)) :: B(:,:)
  REAL(KIND(1.0E0)) :: RESIDUAL, RELW

  REAL(KIND(1.0E0)), ALLOCATABLE :: A(:,:)
  REAL(KIND(1.0E0)), ALLOCATABLE :: WORK(:)
  REAL(KIND(1.0E0)) :: QUERY(1), PFAFF, RPFAFF
  INTEGER, ALLOCATABLE :: IWORK(:)
  INTEGER :: INFO, ISTAT, N, LWORK

  N = SIZE(B,2)

  !Allocate memory
  IF( METHOD == 'P') THEN
     ALLOCATE(A(N,N), IWORK(N), STAT = ISTAT)
  ELSE
     ALLOCATE(A(N,N), IWORK(1), STAT = ISTAT)
  END IF

  IF( ISTAT /= 0 ) THEN
     WRITE (*,*) "Ran out of memory in test!"
     STOP
  END IF

  !copy for later comparison
  A=B

  ! first workspace query
  CALL SKPFA(UPLO, METHOD, N, A, MAX(1,N), PFAFF, IWORK, QUERY, -1, INFO)
  IF( METHOD == 'P' ) THEN
     LWORK = MAX(1, INT(QUERY(1)*RELW))
  ELSE
     LWORK = MAX(1, 2*N-1, INT(QUERY(1)*RELW))
  END IF

  ALLOCATE(WORK(LWORK), STAT=ISTAT)

  IF( ISTAT /= 0 ) THEN
     WRITE (*,*) "Ran out of memory in test!"
     STOP
  END IF

  CALL SKPFA(UPLO, METHOD, N, A, MAX(1,N), PFAFF, IWORK, WORK, LWORK, INFO)

  !Now compare to a reference pfaffian
  CALL reference_pfaffian(B, RPFAFF)

  RESIDUAL = resid_pfaffian(N, PFAFF, RPFAFF)

  DEALLOCATE(A, IWORK, WORK)

END SUBROUTINE check_pfaff_F77_s

SUBROUTINE check_pfaff_F95_s(UPLO, METHOD, B, RESIDUAL)
  USE matrix_tools
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, METHOD
  REAL(KIND(1.0E0)) :: B(:,:)
  REAL(KIND(1.0E0)) :: RESIDUAL

  REAL(KIND(1.0E0)), ALLOCATABLE :: A(:,:)
  REAL(KIND(1.0E0)) :: PFAFF, RPFAFF
  INTEGER :: N, ISTAT

  N = SIZE(B,2)

  ALLOCATE(A(N,N), STAT = ISTAT)

  IF( ISTAT /= 0 ) THEN
     WRITE (*,*) "Ran out of memory in test!"
     STOP
  END IF

  !copy for later comparison
  A=B

  CALL SKPFA(A, PFAFF, UPLO = UPLO, MTHD = METHOD)

  !Now compare to a reference pfaffian
  CALL reference_pfaffian(B, RPFAFF)

  RESIDUAL = resid_pfaffian(N, PFAFF, RPFAFF)

END SUBROUTINE check_pfaff_F95_s

SUBROUTINE check_pfaff10_F77_s(UPLO, METHOD, B, RESIDUAL, RELW)
  USE matrix_tools
  USE F77_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, METHOD
  REAL(KIND(1.0E0)) :: B(:,:)
  REAL(KIND(1.0E0)) :: RESIDUAL, RELW

  REAL(KIND(1.0E0)), ALLOCATABLE :: A(:,:)
  REAL(KIND(1.0E0)), ALLOCATABLE :: WORK(:)
  REAL(KIND(1.0E0)) :: QUERY(1), PFAFF(2), RPFAFF
  INTEGER, ALLOCATABLE :: IWORK(:)
  INTEGER :: INFO, ISTAT, N, LWORK

  N = SIZE(B,2)

  !Allocate memory
  IF( METHOD == 'P') THEN
     ALLOCATE(A(N,N), IWORK(N), STAT = ISTAT)
  ELSE
     ALLOCATE(A(N,N), IWORK(1), STAT = ISTAT)
  END IF

  IF( ISTAT /= 0 ) THEN
     WRITE (*,*) "Ran out of memory in test!"
     STOP
  END IF

  !copy for later comparison
  A=B

  ! first workspace query
  CALL SKPF10(UPLO, METHOD, N, A, MAX(1,N), PFAFF, IWORK, QUERY, -1, INFO)
  IF( METHOD == 'P' ) THEN
     LWORK = MAX(1, INT(QUERY(1)*RELW))
  ELSE
     LWORK = MAX(1, 2*N-1, INT(QUERY(1)*RELW))
  END IF

  ALLOCATE(WORK(LWORK), STAT=ISTAT)

  IF( ISTAT /= 0 ) THEN
     WRITE (*,*) "Ran out of memory in test!"
     STOP
  END IF

  CALL SKPF10(UPLO, METHOD, N, A, MAX(1,N), PFAFF, IWORK, WORK, LWORK, INFO)

  !Now compare to a reference pfaffian
  CALL reference_pfaffian(B, RPFAFF)

  RESIDUAL = resid_pfaffian(N, PFAFF, RPFAFF)

  DEALLOCATE(A, IWORK, WORK)

END SUBROUTINE check_pfaff10_F77_s

SUBROUTINE check_pfaff10_F95_s(UPLO, METHOD, B, RESIDUAL)
  USE matrix_tools
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, METHOD
  REAL(KIND(1.0E0)) :: B(:,:)
  REAL(KIND(1.0E0)) :: RESIDUAL

  REAL(KIND(1.0E0)), ALLOCATABLE :: A(:,:)
  REAL(KIND(1.0E0)) :: PFAFF(2), RPFAFF
  INTEGER :: N, ISTAT

  N = SIZE(B,2)

  ALLOCATE(A(N,N), STAT = ISTAT)

  IF( ISTAT /= 0 ) THEN
     WRITE (*,*) "Ran out of memory in test!"
     STOP
  END IF

  !copy for later comparison
  A=B

  CALL SKPF10(A, PFAFF, UPLO = UPLO, MTHD = METHOD)

  !Now compare to a reference pfaffian
  CALL reference_pfaffian(B, RPFAFF)

  RESIDUAL = resid_pfaffian(N, PFAFF, RPFAFF)

END SUBROUTINE check_pfaff10_F95_s

SUBROUTINE check_pfaff_band_F77_s(UPLO, A, B, RESIDUAL)
  USE matrix_tools
  USE F77_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO
  REAL(KIND(1.0E0)) :: A(:,:), B(:,:)
  REAL(KIND(1.0E0)) :: RESIDUAL

  REAL(KIND(1.0E0)) :: PFAFF, RPFAFF
  REAL(KIND(1.0E0)), ALLOCATABLE :: AB(:,:)
  REAL(KIND(1.0E0)), ALLOCATABLE :: WORK(:)
  INTEGER :: INFO, ISTAT, N, KD

  N = SIZE(A,2)
  KD = SIZE(A,1) - 1

  !Allocate memory
  ALLOCATE(AB(KD+1,N), WORK(3*N-1), STAT = ISTAT)

  IF( ISTAT /= 0 ) THEN
     WRITE (*,*) "Ran out of memory in test!"
     STOP
  END IF

  !copy to leave A unchanged
  AB=A

  !compute the factorization
  CALL SKBPFA(UPLO, N, KD, AB, MAX(1,SIZE(AB,1)), PFAFF, WORK, INFO)

  !Now compare to a reference pfaffian
  CALL reference_pfaffian(B, RPFAFF)

  RESIDUAL = resid_pfaffian(N, PFAFF, RPFAFF)

  DEALLOCATE(AB, WORK)

END SUBROUTINE check_pfaff_band_F77_s

SUBROUTINE check_pfaff_band_F95_s(UPLO, A, B, RESIDUAL)
  USE matrix_tools
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO
  REAL(KIND(1.0E0)) :: A(:,:), B(:,:)
  REAL(KIND(1.0E0)) :: RESIDUAL

  REAL(KIND(1.0E0)) :: PFAFF, RPFAFF
  REAL(KIND(1.0E0)), ALLOCATABLE :: AB(:,:)
  INTEGER :: ISTAT, N, KD

  N = SIZE(A,2)
  KD = SIZE(A,1) - 1

  !Allocate memory
  ALLOCATE(AB(KD+1,N), STAT = ISTAT)

  IF( ISTAT /= 0 ) THEN
     WRITE (*,*) "Ran out of memory in test!"
     STOP
  END IF

  !copy to leave A unchanged
  AB=A

  !compute the factorization
  CALL SKBPFA(AB, PFAFF, UPLO = UPLO)

  !Now compare to a reference pfaffian
  CALL reference_pfaffian(B, RPFAFF)

  RESIDUAL = resid_pfaffian(N, PFAFF, RPFAFF)

  DEALLOCATE(AB)

END SUBROUTINE check_pfaff_band_F95_s

SUBROUTINE check_pfaff10_band_F77_s(UPLO, A, B, RESIDUAL)
  USE matrix_tools
  USE F77_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO
  REAL(KIND(1.0E0)) :: A(:,:), B(:,:)
  REAL(KIND(1.0E0)) :: RESIDUAL

  REAL(KIND(1.0E0)) :: PFAFF(2), RPFAFF
  REAL(KIND(1.0E0)), ALLOCATABLE :: AB(:,:)
  REAL(KIND(1.0E0)), ALLOCATABLE :: WORK(:)
  INTEGER :: INFO, ISTAT, N, KD

  N = SIZE(A,2)
  KD = SIZE(A,1) - 1

  !Allocate memory
  ALLOCATE(AB(KD+1,N), WORK(3*N-1), STAT = ISTAT)

  IF( ISTAT /= 0 ) THEN
     WRITE (*,*) "Ran out of memory in test!"
     STOP
  END IF

  !copy to leave A unchanged
  AB=A

  !compute the factorization
  CALL SKBPF10(UPLO, N, KD, AB, MAX(1,SIZE(AB,1)), PFAFF, WORK, INFO)

  !Now compare to a reference pfaffian
  CALL reference_pfaffian(B, RPFAFF)

  RESIDUAL = resid_pfaffian(N, PFAFF, RPFAFF)

  DEALLOCATE(AB, WORK)

END SUBROUTINE check_pfaff10_band_F77_s

SUBROUTINE check_pfaff10_band_F95_s(UPLO, A, B, RESIDUAL)
  USE matrix_tools
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO
  REAL(KIND(1.0E0)) :: A(:,:), B(:,:)
  REAL(KIND(1.0E0)) :: RESIDUAL

  REAL(KIND(1.0E0)) :: PFAFF(2), RPFAFF
  REAL(KIND(1.0E0)), ALLOCATABLE :: AB(:,:)
  INTEGER :: ISTAT, N, KD

  N = SIZE(A,2)
  KD = SIZE(A,1) - 1

  !Allocate memory
  ALLOCATE(AB(KD+1,N), STAT = ISTAT)

  IF( ISTAT /= 0 ) THEN
     WRITE (*,*) "Ran out of memory in test!"
     STOP
  END IF

  !copy to leave A unchanged
  AB=A

  !compute the factorization
  CALL SKBPF10(AB, PFAFF, UPLO = UPLO)

  !Now compare to a reference pfaffian
  CALL reference_pfaffian(B, RPFAFF)

  RESIDUAL = resid_pfaffian(N, PFAFF, RPFAFF)

  DEALLOCATE(AB)

END SUBROUTINE check_pfaff10_band_F95_s
