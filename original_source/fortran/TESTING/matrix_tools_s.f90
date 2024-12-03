!Functions for single precision

!Extract the unit triangular matrix from the return values
!of xSKTRF, xSKTF2
SUBROUTINE extract_unit_tri_s(UPLO, A, B)
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO
  REAL(KIND(1.0E0)) :: A(:,:), B(:,:)

  INTEGER :: J,N

  N=SIZE(A,1)

  B=0D0

  IF( UPLO == 'U' ) THEN
     DO J=1, N
        B(J,J)=1D0
        IF( J > 2 ) B(1:J-2, J-1)=A(1:J-2, J)
     END DO
  ELSE
     DO J=1, N
        B(J,J)=1D0
        IF( J<N-1 ) B(J+2:N, J+1)=A(J+2:N,J)
     END DO
  END IF

END SUBROUTINE extract_unit_tri_s

!Extract the orthogonal/unitary transformation from the return
!values of xSKTRD/xSKTD2
SUBROUTINE extract_unitary_s(UPLO, A, TAU, Q)
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO
  REAL(KIND(1.0E0)) :: A(:,:), TAU(:), Q(:,:)

  REAL(KIND(1.0E0)), ALLOCATABLE :: WORK(:)
  INTEGER :: N, ISTAT, INFO

  N=SIZE(A,1)

  ALLOCATE(WORK(MAX(1,N)), STAT = ISTAT)

  IF( ISTAT /= 0 ) THEN
     WRITE (*,*) "Ran out of memory during test!"
     STOP
  END IF

  Q=A

  CALL SORGTR( UPLO, N, Q, MAX(1,N), TAU, WORK, MAX(1,N), INFO)

  DEALLOCATE(WORK)

END SUBROUTINE extract_unitary_s


!Extract the tridiagonal output from
!of xSKTRF, xSKTF2, xSKTRD, xSKTD2
SUBROUTINE extract_trid_s(UPLO, A, T)
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO
  REAL(KIND(1.0E0)) :: A(:,:), T(:,:)

  INTEGER :: J,N

  N=SIZE(A,1)

  T=0D0

  IF( UPLO == 'U') THEN
     DO J=2, N
        T(J-1,J) = A(J-1,J)
        T(J,J-1) = -A(J-1,J)
     END DO
  ELSE
     DO J=2, N
        T(J,J-1) = A(J,J-1)
        T(J-1,J) = -A(J,J-1)
     END DO
  END IF

END SUBROUTINE extract_trid_s

!Extract the tridiagonal output from
!of xSKBTRD
SUBROUTINE extract_trid_band_s(UPLO, A, T)
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO
  REAL(KIND(1.0E0)) :: A(:,:), T(:,:)

  INTEGER :: J,N,KD

  N=SIZE(A,2)
  KD=SIZE(A,1)-1
  T=0D0

  IF( KD > 0 ) THEN
     IF( UPLO == 'U') THEN
        DO J=2, N
           T(J-1,J) = A(KD,J)
           T(J,J-1) = -A(KD,J)
        END DO
     ELSE
        DO J=2, N
           T(J,J-1) = A(2,J-1)
           T(J-1,J) = -A(2,J-1)
        END DO
     END IF
  END IF

END SUBROUTINE extract_trid_band_s

!make a random skew-symmetric matrix
SUBROUTINE make_skew_mat_s(A, PROB)
  IMPLICIT NONE
  REAL(KIND(1.0E0)) :: A (:,:)
  REAL(KIND(1.0E0)), OPTIONAL :: PROB

  REAL(KIND(1.0E0)) :: P

  INTEGER :: N, I, J

  P = 2.0
  IF( PRESENT(PROB) ) P = PROB

  N=SIZE(A,2)

  A=0D0

  DO I=1, N
     DO J=1, I-1
        IF(RAND(0) <= P) THEN
           A(I,J)=RAND(0)*2-1
           A(J,I)=-A(I,J)
        END IF
     END DO
  END DO
END SUBROUTINE make_skew_mat_s

SUBROUTINE make_skew_mat_banded_s(Au, Ad, A, PROB)
  IMPLICIT NONE
  REAL(KIND(1.0E0)) :: Au(:,:), Ad(:,:), A(:,:)
  REAL(KIND(1.0E0)), OPTIONAL :: PROB

  REAL(KIND(1.0E0)) :: P

  INTEGER :: N, I, J, KD

  P = 2.0
  IF( PRESENT(PROB) ) P = PROB

  N=SIZE(Au,2)
  KD=SIZE(Au,1)

  !Assume without further checks that Ad, A has the same sizes

  Au=0D0
  Ad=0D0
  A=0D0

  DO I=1,N
     DO J=1, I-1
        IF( I-J < KD ) THEN
           IF(RAND(0) <= P) THEN
              A(I,J)=RAND(0)*2-1
              A(J,I)=-A(I,J)

              Au(KD+J-I,I) = A(J,I)
              Ad(1+I-J,J) = A(I,J)
           END IF

        END IF
     END DO
  END DO
END SUBROUTINE make_skew_mat_banded_s


!compute the norm of a matrix
REAL(KIND(1.0E0)) FUNCTION one_norm_s(A)
  IMPLICIT NONE
  REAL(KIND(1.0E0)) :: A (:,:)

  INTEGER :: N
  REAL(KIND(1.0E0)) :: DUMMY, SLANGE
  EXTERNAL DLANGE

  N = SIZE(A,2)

  one_norm_s = SLANGE( "1", N, N, A, N, DUMMY)
  RETURN
END FUNCTION one_norm_s

!compute the relative residual of the difference of matrices
REAL(KIND(1.0E0)) FUNCTION resid_s(A,B)
  USE matrix_tools, ONLY: one_norm
  IMPLICIT NONE
  REAL(KIND(1.0E0)) :: A (:,:), B(:,:)

  INTEGER :: N
  REAL(KIND(1.0E0)) :: normdiff, normb, SLAMCH
  EXTERNAL DLANGE, DLAMCH

  N = SIZE(A,2)

  normb = one_norm(B)
  normdiff = one_norm(A-B)

  IF( N > 0 .AND. normb > 0 ) THEN
     resid_s = normdiff/(N*normb*SLAMCH("E"))
     RETURN
  ELSE
     resid_s = one_norm(A-B)
     RETURN
  END IF
END FUNCTION resid_s

!print a matrix on the screen
SUBROUTINE print_mat_s(A)
  IMPLICIT NONE
  REAL(KIND(1.0E0)) :: A (:,:)

  INTEGER :: N,M,I,J

  N=SIZE(A,2)
  M=SIZE(A,1)

  DO I=1, M
     DO J=1, N
        WRITE (*,'(F6.3)',Advance='NO') A(I,J)
        WRITE (*,'(A)',Advance='NO') "  "
     END DO
     WRITE (*,*)
  END DO

  WRITE (*,*)
END SUBROUTINE print_mat_s

SUBROUTINE rowcol_invperm_s(UPLO, A, IPIV)
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO
  REAL(KIND(1.0E0)) :: A(:,:)
  INTEGER :: IPIV(:)

  INTEGER :: N, I, N1, N2, DN
  REAL(KIND(1.0E0)), ALLOCATABLE :: TEMP(:)

  N=SIZE(A,2)

  ALLOCATE(TEMP(N))

  IF( UPLO == 'U' ) THEN
     N1 = 1
     N2 = N
     DN = 1
  ELSE
     N1 = N
     N2 = 1
     DN = -1
  END IF

  !First switch columns
  DO I=N1, N2, DN
     IF( IPIV(I) /= I ) THEN
        TEMP=A(1:N,I)
        A(1:N,I)=A(1:N,IPIV(I))
        A(1:N,IPIV(I))=TEMP
     END IF
  END DO

  !then switch rows
  DO I=N1, N2, DN
     IF( IPIV(I) /= I ) THEN
        TEMP=A(I, 1:N)
        A(I,1:N)=A(IPIV(I),1:N)
        A(IPIV(I),1:N)=TEMP
     END IF
  END DO

END SUBROUTINE rowcol_invperm_s

SUBROUTINE reference_pfaff_s(A, PFAFF)
  USE F95_PFAPACK
  IMPLICIT NONE
  REAL(KIND(1.0E0)) :: A (:,:), PFAFF

  INTEGER :: N, ISTAT, I
  INTEGER, ALLOCATABLE :: IPIV(:)
  REAL(KIND(1.0E0)), ALLOCATABLE :: Ac (:,:)

  N = SIZE(A,2)

  IF( MOD(N,2) == 1 ) THEN
     PFAFF = 0.0E0
  ELSE IF( N == 0) THEN
     PFAFF = 1.0E0
  ELSE IF( N == 2) THEN
     PFAFF = A(1,2)
  ELSE IF( N == 4 ) THEN
     PFAFF = A(1,2)*A(3,4) - A(1,3)*A(2,4) + A(2,3)*A(1,4)
  ELSE
     !Use the (full) LTL decomposition which has been tested before
     ALLOCATE(IPIV(N), Ac(N,N), STAT=ISTAT)

     IF( ISTAT /=0 ) THEN
        WRITE (*,*) "Ran out of memory in test!"
        STOP
     END IF

     Ac=A
     CALL SKTRF(Ac, IPIV, UPLO='U')

     PFAFF = 1.0E0
     DO I=1, N
        IF( IPIV(I) /= I ) PFAFF = -PFAFF
        IF( MOD(I,2) == 1) then
           PFAFF = PFAFF*Ac(I,I+1)
        end IF
     END DO
  END IF

END SUBROUTINE reference_pfaff_s

!compute the relative (or absolute) error of two Pfaffians
REAL(KIND(1.0E0)) FUNCTION resid_pfaff_s(N, A,B)
  IMPLICIT NONE
  INTEGER :: N
  REAL(KIND(1.0E0)) :: A, B

  REAL(KIND(1.0E0)) :: SLAMCH
  EXTERNAL SLAMCH

  !Take either the relative or absolute error, depending on what's smaller
  IF( ABS(B) > 0.0E0) THEN
     resid_pfaff_s = MIN( ABS(A-B)/(MAX(N,1)*ABS(B)*SLAMCH("E")), &
                          ABS(A-B)/(MAX(N,1)*SLAMCH("E")) )
  ELSE
     resid_pfaff_s = ABS(A-B)/(MAX(N,1)*SLAMCH("E"))
  END IF

END FUNCTION resid_pfaff_s

REAL(KIND(1.0E0)) FUNCTION resid_pfaff10_s(N, A,B)
  IMPLICIT NONE
  INTEGER :: N
  REAL(KIND(1.0E0)) :: A(2), B

  REAL(KIND(1.0E0)) :: SLAMCH
  EXTERNAL SLAMCH

  !Take either the relative or absolute error, depending on what's smaller
  IF( ABS(B) > 0.0E0 ) THEN
     resid_pfaff10_s = &
          MIN( ABS(A(1)*10**A(2)-B)/(MAX(N,1)*ABS(B)*SLAMCH("E")), &
               ABS(A(1)*10**A(2)-B)/(MAX(N,1)*SLAMCH("E")) )
  ELSE
     resid_pfaff10_s = ABS(A(1)*10**A(2)-B)/(MAX(N,1)*SLAMCH("E"))
  END IF

END FUNCTION resid_pfaff10_s
