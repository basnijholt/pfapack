! The routines in this file check if the output of the C interface
! is identical to the Fortran 95 interface. They should be in fact
! identical, as both interfaces call the underlying Fortran 77 routines
! with the same input.

!SKPFA
SUBROUTINE check_skpfa_s(N, A, A2, PFAFF2, UPLO, MTHD)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, MTHD
  INTEGER :: N
  REAL(KIND(1.0E0)) :: A(N, N), A2(N, N), PFAFF2

  REAL(KIND(1.0E0)) :: PFAFF

  CALL SKPFA(A, PFAFF, UPLO = UPLO, MTHD = MTHD)

  IF( ANY(A /= A2) .OR. PFAFF /= PFAFF2) THEN
     PRINT *, "skpfa_s: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

SUBROUTINE check_skpfa_d(N, A, A2, PFAFF2, UPLO, MTHD)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, MTHD
  INTEGER :: N
  REAL(KIND(1.0D0)) :: A(N, N), A2(N, N), PFAFF2

  REAL(KIND(1.0D0)) :: PFAFF

  CALL SKPFA(A, PFAFF, UPLO = UPLO, MTHD = MTHD)

  IF( ANY(A /= A2) .OR. PFAFF /= PFAFF2) THEN
     PRINT *, "skpfa_d: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

SUBROUTINE check_skpfa_c(N, A, A2, PFAFF2, UPLO, MTHD)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, MTHD
  INTEGER :: N
  COMPLEX(KIND(1.0E0)) :: A(N, N), A2(N, N), PFAFF2

  COMPLEX(KIND(1.0E0)) :: PFAFF

  CALL SKPFA(A, PFAFF, UPLO = UPLO, MTHD = MTHD)

  IF( ANY(A /= A2) .OR. PFAFF /= PFAFF2) THEN
     PRINT *, "skpfa_c: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

SUBROUTINE check_skpfa_z(N, A, A2, PFAFF2, UPLO, MTHD)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, MTHD
  INTEGER :: N
  COMPLEX(KIND(1.0D0)) :: A(N, N), A2(N, N), PFAFF2

  COMPLEX(KIND(1.0D0)) :: PFAFF

  CALL SKPFA(A, PFAFF, UPLO = UPLO, MTHD = MTHD)

  IF( ANY(A /= A2) .OR. PFAFF /= PFAFF2) THEN
     PRINT *, "skpfa_z: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

!SKPF10
SUBROUTINE check_skpf10_s(N, A, A2, PFAFF2, UPLO, MTHD)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, MTHD
  INTEGER :: N
  REAL(KIND(1.0E0)) :: A(N, N), A2(N, N), PFAFF2(2)

  REAL(KIND(1.0E0)) :: PFAFF(2)

  CALL SKPF10(A, PFAFF, UPLO = UPLO, MTHD = MTHD)

  IF( ANY(A /= A2) .OR. ANY(PFAFF /= PFAFF2) ) THEN
     PRINT *, "skpf10_s: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

SUBROUTINE check_skpf10_d(N, A, A2, PFAFF2, UPLO, MTHD)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, MTHD
  INTEGER :: N
  REAL(KIND(1.0D0)) :: A(N, N), A2(N, N), PFAFF2(2)

  REAL(KIND(1.0D0)) :: PFAFF(2)

  CALL SKPF10(A, PFAFF, UPLO = UPLO, MTHD = MTHD)

  IF( ANY(A /= A2) .OR. ANY(PFAFF /= PFAFF2) ) THEN
     PRINT *, "skpf10_d: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

SUBROUTINE check_skpf10_c(N, A, A2, PFAFF2, UPLO, MTHD)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, MTHD
  INTEGER :: N
  COMPLEX(KIND(1.0E0)) :: A(N, N), A2(N, N), PFAFF2(2)

  COMPLEX(KIND(1.0E0)) :: PFAFF(2)

  CALL SKPF10(A, PFAFF, UPLO = UPLO, MTHD = MTHD)

  IF( ANY(A /= A2) .OR. ANY(PFAFF /= PFAFF2) ) THEN
     PRINT *, "skpf10_c: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

SUBROUTINE check_skpf10_z(N, A, A2, PFAFF2, UPLO, MTHD)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, MTHD
  INTEGER :: N
  COMPLEX(KIND(1.0D0)) :: A(N, N), A2(N, N), PFAFF2(2)

  COMPLEX(KIND(1.0D0)) :: PFAFF(2)

  CALL SKPF10(A, PFAFF, UPLO = UPLO, MTHD = MTHD)

  IF( ANY(A /= A2) .OR. ANY(PFAFF /= PFAFF2) ) THEN
     PRINT *, "skpf10_z: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

!SKBPFA
SUBROUTINE check_skbpfa_s(N, KD, A, A2, PFAFF2, UPLO)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO
  INTEGER :: N, KD
  REAL(KIND(1.0E0)) :: A(KD+1, N), A2(KD+1, N), PFAFF2

  REAL(KIND(1.0E0)) :: PFAFF

  CALL SKBPFA(A, PFAFF, UPLO = UPLO)

  IF( ANY(A /= A2) .OR. PFAFF /= PFAFF2) THEN
     PRINT *, "skbpfa_s: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

SUBROUTINE check_skbpfa_d(N, KD, A, A2, PFAFF2, UPLO)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO
  INTEGER :: N, KD
  REAL(KIND(1.0D0)) :: A(KD+1, N), A2(KD+1, N), PFAFF2

  REAL(KIND(1.0D0)) :: PFAFF

  CALL SKBPFA(A, PFAFF, UPLO = UPLO)

  IF( ANY(A /= A2) .OR. PFAFF /= PFAFF2) THEN
     PRINT *, "skbpfa_s: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

SUBROUTINE check_skbpfa_c(N, KD, A, A2, PFAFF2, UPLO)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO
  INTEGER :: N, KD
  COMPLEX(KIND(1.0E0)) :: A(KD+1, N), A2(KD+1, N), PFAFF2

  COMPLEX(KIND(1.0E0)) :: PFAFF

  CALL SKBPFA(A, PFAFF, UPLO = UPLO)

  IF( ANY(A /= A2) .OR. PFAFF /= PFAFF2) THEN
     PRINT *, "skbpfa_c: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

SUBROUTINE check_skbpfa_z(N, KD, A, A2, PFAFF2, UPLO)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO
  INTEGER :: N, KD
  COMPLEX(KIND(1.0D0)) :: A(KD+1, N), A2(KD+1, N), PFAFF2

  COMPLEX(KIND(1.0D0)) :: PFAFF

  CALL SKBPFA(A, PFAFF, UPLO = UPLO)

  IF( ANY(A /= A2) .OR. PFAFF /= PFAFF2) THEN
     PRINT *, "skbpfa_z: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

!SKBPF10
SUBROUTINE check_skbpf10_s(N, KD, A, A2, PFAFF2, UPLO)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO
  INTEGER :: N, KD
  REAL(KIND(1.0E0)) :: A(KD+1, N), A2(KD+1, N), PFAFF2(2)

  REAL(KIND(1.0E0)) :: PFAFF(2)

  CALL SKBPF10(A, PFAFF, UPLO = UPLO)

  IF( ANY(A /= A2) .OR. ANY(PFAFF /= PFAFF2) ) THEN
     PRINT *, "skbpf10_s: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

SUBROUTINE check_skbpf10_d(N, KD, A, A2, PFAFF2, UPLO)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO
  INTEGER :: N, KD
  REAL(KIND(1.0D0)) :: A(KD+1, N), A2(KD+1, N), PFAFF2(2)

  REAL(KIND(1.0D0)) :: PFAFF(2)

  CALL SKBPF10(A, PFAFF, UPLO = UPLO)

  IF( ANY(A /= A2) .OR. ANY(PFAFF /= PFAFF2) ) THEN
     PRINT *, "skbpf10_d: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

SUBROUTINE check_skbpf10_c(N, KD, A, A2, PFAFF2, UPLO)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO
  INTEGER :: N, KD
  COMPLEX(KIND(1.0E0)) :: A(KD+1, N), A2(KD+1, N), PFAFF2(2)

  COMPLEX(KIND(1.0E0)) :: PFAFF(2)

  CALL SKBPF10(A, PFAFF, UPLO = UPLO)

  IF( ANY(A /= A2) .OR. ANY(PFAFF /= PFAFF2) ) THEN
     PRINT *, "skbpf10_c: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

SUBROUTINE check_skbpf10_z(N, KD, A, A2, PFAFF2, UPLO)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO
  INTEGER :: N, KD
  COMPLEX(KIND(1.0D0)) :: A(KD+1, N), A2(KD+1, N), PFAFF2(2)

  COMPLEX(KIND(1.0D0)) :: PFAFF(2)

  CALL SKBPF10(A, PFAFF, UPLO = UPLO)

  IF( ANY(A /= A2) .OR. ANY(PFAFF /= PFAFF2) ) THEN
     PRINT *, "skbpf10_z: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

!SKTRF
SUBROUTINE check_sktrf_s(N, A, A2, IPIV2, UPLO, MODE)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, MODE
  INTEGER :: N, IPIV2(N)
  REAL(KIND(1.0E0)) :: A(N, N), A2(N, N)

  INTEGER :: IPIV(N)

  CALL SKTRF(A, IPIV, UPLO = UPLO, MODE = MODE)

  IF( ANY(A /= A2) .OR. ANY(IPIV /= IPIV2) ) THEN
     PRINT *, "sktrf_s: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

SUBROUTINE check_sktrf_d(N, A, A2, IPIV2, UPLO, MODE)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, MODE
  INTEGER :: N, IPIV2(N)
  REAL(KIND(1.0D0)) :: A(N, N), A2(N, N)

  INTEGER :: IPIV(N)

  CALL SKTRF(A, IPIV, UPLO = UPLO, MODE = MODE)

  IF( ANY(A /= A2) .OR. ANY(IPIV /= IPIV2) ) THEN
     PRINT *, "sktrf_s: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

SUBROUTINE check_sktrf_c(N, A, A2, IPIV2, UPLO, MODE)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, MODE
  INTEGER :: N, IPIV2(N)
  COMPLEX(KIND(1.0E0)) :: A(N, N), A2(N, N)

  INTEGER :: IPIV(N)

  CALL SKTRF(A, IPIV, UPLO = UPLO, MODE = MODE)

  IF( ANY(A /= A2) .OR. ANY(IPIV /= IPIV2) ) THEN
     PRINT *, "sktrf_s: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

SUBROUTINE check_sktrf_z(N, A, A2, IPIV2, UPLO, MODE)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, MODE
  INTEGER :: N, IPIV2(N)
  COMPLEX(KIND(1.0D0)) :: A(N, N), A2(N, N)

  INTEGER :: IPIV(N)

  CALL SKTRF(A, IPIV, UPLO = UPLO, MODE = MODE)

  IF( ANY(A /= A2) .OR. ANY(IPIV /= IPIV2) ) THEN
     PRINT *, "sktrf_s: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

!SKTRD
SUBROUTINE check_sktrd_s(N, A, A2, TAU2, UPLO, MODE)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, MODE
  INTEGER :: N
  REAL(KIND(1.0E0)) :: A(N, N), A2(N, N), TAU2(N-1)

  REAL(KIND(1.0E0)) :: TAU(N-1)

  CALL SKTRD(A, TAU, UPLO = UPLO, MODE = MODE)

  IF( ANY(A /= A2) .OR. ANY(TAU /= TAU2) ) THEN
     PRINT *, "sktrd_s: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

SUBROUTINE check_sktrd_d(N, A, A2, TAU2, UPLO, MODE)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, MODE
  INTEGER :: N
  REAL(KIND(1.0D0)) :: A(N, N), A2(N, N), TAU2(N-1)

  REAL(KIND(1.0D0)) :: TAU(N-1)

  CALL SKTRD(A, TAU, UPLO = UPLO, MODE = MODE)

  IF( ANY(A /= A2) .OR. ANY(TAU /= TAU2) ) THEN
     PRINT *, "sktrd_s: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

SUBROUTINE check_sktrd_c(N, A, A2, TAU2, UPLO, MODE)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, MODE
  INTEGER :: N
  COMPLEX(KIND(1.0E0)) :: A(N, N), A2(N, N), TAU2(N-1)

  COMPLEX(KIND(1.0E0)) :: TAU(N-1)

  CALL SKTRD(A, TAU, UPLO = UPLO, MODE = MODE)

  IF( ANY(A /= A2) .OR. ANY(TAU /= TAU2) ) THEN
     PRINT *, "sktrd_s: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

SUBROUTINE check_sktrd_z(N, A, A2, TAU2, UPLO, MODE)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, MODE
  INTEGER :: N
  COMPLEX(KIND(1.0D0)) :: A(N, N), A2(N, N), TAU2(N-1)

  COMPLEX(KIND(1.0D0)) :: TAU(N-1)

  CALL SKTRD(A, TAU, UPLO = UPLO, MODE = MODE)

  IF( ANY(A /= A2) .OR. ANY(TAU /= TAU2) ) THEN
     PRINT *, "sktrd_s: output of C interface differs from Fortran"
     STOP
  END IF
END SUBROUTINE

!SKBTRD
SUBROUTINE check_skbtrd_s(N, KD, A, A2, DETQ2, Q2, UPLO, MODE, VECT)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, MODE, VECT
  INTEGER :: N, KD
  REAL(KIND(1.0E0)) :: A(KD+1, N), A2(KD+1, N), Q2(N, N), DETQ2

  REAL(KIND(1.0E0)) :: Q(N, N), DETQ

  Q = 1

  CALL SKBTRD(A, UPLO = UPLO, MODE = MODE, DETQ = DETQ, Q = Q, VECT = VECT)

  IF( ANY(A /= A2) .OR. DETQ /= DETQ2 ) THEN
     PRINT *, "skbtrd_s: output of C interface differs from Fortran"
     STOP
  END IF

  IF( VECT /= 'N' ) THEN
     IF( ANY(Q /= Q2) ) THEN
        PRINT *, "skbtrd_s: output of C interface differs from Fortran"
        STOP
     END IF
  END IF

END SUBROUTINE

SUBROUTINE check_skbtrd_d(N, KD, A, A2, DETQ2, Q2, UPLO, MODE, VECT)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, MODE, VECT
  INTEGER :: N, KD
  REAL(KIND(1.0D0)) :: A(KD+1, N), A2(KD+1, N), Q2(N, N), DETQ2

  REAL(KIND(1.0D0)) :: Q(N, N), DETQ

  Q = 1

  CALL SKBTRD(A, UPLO = UPLO, MODE = MODE, DETQ = DETQ, Q = Q, VECT = VECT)

  IF( ANY(A /= A2) .OR. DETQ /= DETQ2 ) THEN
     PRINT *, "skbtrd_s: output of C interface differs from Fortran"
     STOP
  END IF

  IF( VECT /= 'N' ) THEN
     IF( ANY(Q /= Q2) ) THEN
        PRINT *, "skbtrd_s: output of C interface differs from Fortran"
        STOP
     END IF
  END IF

END SUBROUTINE

SUBROUTINE check_skbtrd_c(N, KD, A, A2, DETQ2, Q2, UPLO, MODE, VECT)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, MODE, VECT
  INTEGER :: N, KD
  COMPLEX(KIND(1.0E0)) :: A(KD+1, N), A2(KD+1, N), Q2(N, N), DETQ2

  COMPLEX(KIND(1.0E0)) :: Q(N, N), DETQ

  Q = 1

  CALL SKBTRD(A, UPLO = UPLO, MODE = MODE, DETQ = DETQ, Q = Q, VECT = VECT)

  IF( ANY(A /= A2) .OR. DETQ /= DETQ2 ) THEN
     PRINT *, "skbtrd_s: output of C interface differs from Fortran"
     STOP
  END IF

  IF( VECT /= 'N' ) THEN
     IF( ANY(Q /= Q2) ) THEN
        PRINT *, "skbtrd_s: output of C interface differs from Fortran"
        STOP
     END IF
  END IF

END SUBROUTINE

SUBROUTINE check_skbtrd_z(N, KD, A, A2, DETQ2, Q2, UPLO, MODE, VECT)
  USE F95_PFAPACK
  IMPLICIT NONE
  CHARACTER(LEN=1) :: UPLO, MODE, VECT
  INTEGER :: N, KD
  COMPLEX(KIND(1.0D0)) :: A(KD+1, N), A2(KD+1, N), Q2(N, N), DETQ2

  COMPLEX(KIND(1.0D0)) :: Q(N, N), DETQ

  Q = 1

  CALL SKBTRD(A, UPLO = UPLO, MODE = MODE, DETQ = DETQ, Q = Q, VECT = VECT)

  IF( ANY(A /= A2) .OR. DETQ /= DETQ2 ) THEN
     PRINT *, "skbtrd_s: output of C interface differs from Fortran"
     STOP
  END IF

  IF( VECT /= 'N' ) THEN
     IF( ANY(Q /= Q2) ) THEN
        PRINT *, "skbtrd_s: output of C interface differs from Fortran"
        STOP
     END IF
  END IF

END SUBROUTINE
