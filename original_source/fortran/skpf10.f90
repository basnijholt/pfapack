!  Purpose
!  =======
!
!  SKPF10 computes the Pfaffian of a dense skew-symmetric matrix, taking
!  special care to avoid numerical under- or overflow.
!  (at the cost of possible additional round-off errors)
!
!  =========
!
!     SUBROUTINE SKPF10( A, PFAFF, UPLO, MTHD, INFO)
!         <type>(<prec>), INTENT(INOUT) :: A(:,:)
!         <type>(<prec>), INTENT(OUT) :: PFAFF(2)
!         CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO, MTHD
!         INTEGER, INTENT(OUT), OPTIONAL :: INFO
!
!    where
!       <type> ::= REAL | COMPLEX
!       <prec> ::= KIND(1.0) | KIND(1.0D0)
!
!  Arguments
!  =========
!
!  A       (input/output) REAL or COMPLEX array,
!          shape (:,:), size(A,1) == size(A,2) >= 0.
!          On entry, the skew-symmetric matrix A.
!             If UPLO = 'U', the upper triangular part of A contains
!                the upper triangular part of the matrix A, and the
!                strictly lower triangular part of A is not referenced.
!             If UPLO = 'L', the lower triangular part of A contains
!                the lower triangular part of the matrix A, and the
!                strictly upper triangular part of A is not referenced.
!          If the matrix size is odd, A is not referenced. If the matrix
!          size is even, A is overwritten by values generated during
!          the computation.
!
!  PFAFF   (output) REAL or COMPLEX array, dimension 2
!          The value of the Pfaffian in the form
!          PFAFF(1)*10**PFAFF(2).
!
!  UPLO    Optional (input) CHARACTER*1
!          If UPLO is present, then:
!            = 'U':  Upper triangle of A is stored;
!            = 'L':  Lower triangle of A is stored.
!          otherwise UPLO = 'U' is assumed
!
!  MTHD    Optional (input) CHARACTER*1
!          If MTHD is present, then:
!            = 'P': Compute Pfaffian using Parlett-Reid algorithm (recommended)
!            = 'H': Compute Pfaffian using Householder reflections
!          otherwise MTHD = 'P' is assumed
!
!  INFO    Optional (output) INTEGER
!          If INFO is present:
!              = 0:    successful exit
!              < 0:    if INFO = -i, the i-th argument had an illegal value
!              = -100: failed to allocate enough internal memory
!          Otherwise, if INFO is not present, the program stops with
!          an error messafe
!
!  Further Details
!  ===============
!
!
!  For odd-sized matrices, the Pfaffian is 0 by default, hence no
!  computation is needed in this case. For even-sized matrices,
!  the Pfaffian is computed by bringing the skew-symmetric matrix A into
!  a partial tridiagonal form pT, either by computing a partial L pT L^T
!  decomposition (MTHD = 'P'), or by a by a unitary congruence transformation
!  Q^H * A * Q^* = pT (MTHD = 'H').
!  These transformations are computed by the routines DSKTRF or DSKTRD,
!  respectively (for further details see there).
!

SUBROUTINE SSKPF10_F95( A, PFAFF, UPLO, MTHD, INFO)
  USE PFAPACK_PREC, ONLY: precision => singleprec
  USE PFAPACK_MESSAGE, ONLY: MESSAGE
  USE F77_PFAPACK, ONLY: SKPF10
  IMPLICIT NONE
  CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO, MTHD
  INTEGER, INTENT(OUT), OPTIONAL :: INFO
  REAL(precision), INTENT(OUT) :: PFAFF(2)
  REAL(precision), INTENT(INOUT) :: A(:,:)

  CHARACTER(LEN=1) :: LUPLO, LMTHD
  INTEGER :: N, LDA, LWORK, LINFO, ISTAT
  INTEGER, ALLOCATABLE :: IWORK(:)
  REAL(precision) :: QWORK(1)
  REAL(precision), ALLOCATABLE :: WORK(:)

  LOGICAL LSAME
  EXTERNAL LSAME

  ! Figure out array sizes and set defaults
  LINFO = 0
  N = SIZE(A,1)
  LDA = MAX(1,N)

  IF( PRESENT(UPLO) ) THEN
     LUPLO = UPLO
  ELSE
     LUPLO = 'U'
  END IF

  IF( PRESENT(MTHD) ) THEN
     LMTHD = MTHD
  ELSE
     LMTHD = 'P'
  END IF

  ! Test the arguments
  IF( SIZE(A,2) /= N .OR. N<0 ) THEN
     LINFO = -1
  ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') ) THEN
     LINFO = -3
  ELSE IF( .NOT.LSAME(LMTHD,'P') .AND. .NOT.LSAME(LMTHD,'H') ) THEN
     LINFO = -4
  ELSE IF( N > 0 ) THEN
     !Allocate fixed size arrays
     IF( LSAME(LMTHD, 'P') ) THEN
        ALLOCATE(IWORK(N), STAT = ISTAT)
     ELSE
        ALLOCATE(IWORK(1), STAT = ISTAT)
     END IF

     IF( ISTAT == 0 ) THEN
        !Do a workspace query
        CALL SKPF10( LUPLO, LMTHD, N, A, LDA, PFAFF, IWORK, QWORK, -1, LINFO)
        LWORK = INT(QWORK(1))
        ALLOCATE(WORK(LWORK), STAT=ISTAT)

        IF( ISTAT /= 0 ) THEN
           ! Allocating the workspace failed. Try minimal workspace
           DEALLOCATE(WORK)
           CALL MESSAGE(-200, 'SKPF10') ! print warning about not enough wspace

           IF( LSAME(LMTHD, 'H') ) THEN
              LWORK=2*N-1
           ELSE
              LWORK = 1
           END IF

           ALLOCATE(WORK(LWORK), STAT=ISTAT)
        END IF
     END IF

     IF( ISTAT == 0 ) THEN
        CALL SKPF10( LUPLO, LMTHD, N, A, LDA, PFAFF, IWORK, WORK, LWORK, LINFO)
     ELSE
        LINFO=-100 ! not enough memory
     END IF

     DEALLOCATE(WORK)
  ELSE !N=0
     PFAFF(1) = 1.0E0
     PFAFF(2) = 0.0E0
  END IF

  CALL MESSAGE(LINFO, 'SKPF10', INFO, ISTAT)
END SUBROUTINE SSKPF10_F95

SUBROUTINE DSKPF10_F95( A, PFAFF, UPLO, MTHD, INFO)
  USE PFAPACK_PREC, ONLY: precision => doubleprec
  USE PFAPACK_MESSAGE, ONLY: MESSAGE
  USE F77_PFAPACK, ONLY: SKPF10
  IMPLICIT NONE
  CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO, MTHD
  INTEGER, INTENT(OUT), OPTIONAL :: INFO
  REAL(precision), INTENT(OUT) :: PFAFF(2)
  REAL(precision), INTENT(INOUT) :: A(:,:)

  CHARACTER(LEN=1) :: LUPLO, LMTHD
  INTEGER :: N, LDA, LWORK, LINFO, ISTAT
  INTEGER, ALLOCATABLE :: IWORK(:)
  REAL(precision) :: QWORK(1)
  REAL(precision), ALLOCATABLE :: WORK(:)

  LOGICAL LSAME
  EXTERNAL LSAME

  ! Figure out array sizes and set defaults
  LINFO = 0
  N = SIZE(A,1)
  LDA = MAX(1,N)

  IF( PRESENT(UPLO) ) THEN
     LUPLO = UPLO
  ELSE
     LUPLO = 'U'
  END IF

  IF( PRESENT(MTHD) ) THEN
     LMTHD = MTHD
  ELSE
     LMTHD = 'P'
  END IF

  ! Test the arguments
  IF( SIZE(A,2) /= N .OR. N<0 ) THEN
     LINFO = -1
  ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') ) THEN
     LINFO = -3
  ELSE IF( .NOT.LSAME(LMTHD,'P') .AND. .NOT.LSAME(LMTHD,'H') ) THEN
     LINFO = -4
  ELSE IF( N > 0 ) THEN
     !Allocate fixed size arrays
     IF( LSAME(LMTHD, 'P') ) THEN
        ALLOCATE(IWORK(N), STAT = ISTAT)
     ELSE
        ALLOCATE(IWORK(1), STAT = ISTAT)
     END IF

     IF( ISTAT == 0 ) THEN
        !Do a workspace query
        CALL SKPF10( LUPLO, LMTHD, N, A, LDA, PFAFF, IWORK, QWORK, -1, LINFO)
        LWORK = INT(QWORK(1))
        ALLOCATE(WORK(LWORK), STAT=ISTAT)

        IF( ISTAT /= 0 ) THEN
           ! Allocating the workspace failed. Try minimal workspace
           DEALLOCATE(WORK)
           CALL MESSAGE(-200, 'SKPF10') ! print warning about not enough wspace

           IF( LSAME(LMTHD, 'H') ) THEN
              LWORK=2*N-1
           ELSE
              LWORK = 1
           END IF

           ALLOCATE(WORK(LWORK), STAT=ISTAT)
        END IF
     END IF

     IF( ISTAT == 0 ) THEN
        CALL SKPF10( LUPLO, LMTHD, N, A, LDA, PFAFF, IWORK, WORK, LWORK, LINFO)
     ELSE
        LINFO=-100 ! not enough memory
     END IF

     DEALLOCATE(WORK)
  ELSE !N=0
     PFAFF(1) = 1.0D0
     PFAFF(2) = 0.0D0
  END IF

  CALL MESSAGE(LINFO, 'SKPF10', INFO, ISTAT)
END SUBROUTINE DSKPF10_F95

SUBROUTINE CSKPF10_F95( A, PFAFF, UPLO, MTHD, INFO)
  USE PFAPACK_PREC, ONLY: precision => singleprec
  USE PFAPACK_MESSAGE, ONLY: MESSAGE
  USE F77_PFAPACK, ONLY: SKPF10
  IMPLICIT NONE
  CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO, MTHD
  INTEGER, INTENT(OUT), OPTIONAL :: INFO
  COMPLEX(precision), INTENT(OUT) :: PFAFF(2)
  COMPLEX(precision), INTENT(INOUT) :: A(:,:)

  CHARACTER(LEN=1) :: LUPLO, LMTHD
  INTEGER :: N, LDA, LWORK, LINFO, ISTAT
  INTEGER, ALLOCATABLE :: IWORK(:)
  COMPLEX(precision) :: QWORK(1)
  COMPLEX(precision), ALLOCATABLE :: WORK(:)
  REAL(precision), ALLOCATABLE :: RWORK(:)
  LOGICAL LSAME
  EXTERNAL LSAME

  ! Figure out array sizes and set defaults
  LINFO = 0
  N = SIZE(A,1)
  LDA = MAX(1,N)

  IF( PRESENT(UPLO) ) THEN
     LUPLO = UPLO
  ELSE
     LUPLO = 'U'
  END IF

  IF( PRESENT(MTHD) ) THEN
     LMTHD = MTHD
  ELSE
     LMTHD = 'P'
  END IF

  ! Test the arguments
  IF( SIZE(A,2) /= N .OR. N<0 ) THEN
     LINFO = -1
  ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') ) THEN
     LINFO = -3
  ELSE IF( .NOT.LSAME(LMTHD,'P') .AND. .NOT.LSAME(LMTHD,'H') ) THEN
     LINFO = -4
  ELSE IF( N>0 ) THEN
     !Allocate fixed size arrays
     IF( LSAME(LMTHD, 'P') ) THEN
        ALLOCATE(IWORK(N), RWORK(1), STAT = ISTAT)
     ELSE
        ALLOCATE(IWORK(1), RWORK(N-1), STAT = ISTAT)
     END IF

     IF( ISTAT == 0 ) THEN
        !Do a workspace query
        CALL SKPF10( LUPLO, LMTHD, N, A, LDA, PFAFF, IWORK, QWORK, -1, &
                     RWORK, LINFO)
        LWORK = INT(QWORK(1))
        ALLOCATE(WORK(LWORK), STAT=ISTAT)

        IF( ISTAT /= 0 ) THEN
           ! Allocating the workspace failed. Try minimal workspace
           DEALLOCATE(WORK)
           CALL MESSAGE(-200, 'SKPF10') ! print warning about not enough wspace

           IF( LSAME(LMTHD, 'H') ) THEN
              LWORK = N
           ELSE
              LWORK = 1
           END IF

           ALLOCATE(WORK(LWORK), STAT=ISTAT)
        END IF
     END IF

     IF( ISTAT == 0 ) THEN
        CALL SKPF10( LUPLO, LMTHD, N, A, LDA, PFAFF, IWORK, WORK, LWORK, &
                     RWORK, LINFO )
     ELSE
        LINFO=-100 ! not enough memory
     END IF

     DEALLOCATE(WORK)
  ELSE !N=0
     PFAFF(1) = 1.0E0
     PFAFF(2) = 0.0E0
  END IF

  CALL MESSAGE(LINFO, 'SKPF10', INFO, ISTAT)
END SUBROUTINE CSKPF10_F95

SUBROUTINE ZSKPF10_F95( A, PFAFF, UPLO, MTHD, INFO)
  USE PFAPACK_PREC, ONLY: precision => doubleprec
  USE PFAPACK_MESSAGE, ONLY: MESSAGE
  USE F77_PFAPACK, ONLY: SKPF10
  IMPLICIT NONE
  CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO, MTHD
  INTEGER, INTENT(OUT), OPTIONAL :: INFO
  COMPLEX(precision), INTENT(OUT) :: PFAFF(2)
  COMPLEX(precision), INTENT(INOUT) :: A(:,:)

  CHARACTER(LEN=1) :: LUPLO, LMTHD
  INTEGER :: N, LDA, LWORK, LINFO, ISTAT
  INTEGER, ALLOCATABLE :: IWORK(:)
  COMPLEX(precision) :: QWORK(1)
  COMPLEX(precision), ALLOCATABLE :: WORK(:)
  REAL(precision), ALLOCATABLE :: RWORK(:)
  LOGICAL LSAME
  EXTERNAL LSAME

  ! Figure out array sizes and set defaults
  LINFO = 0
  N = SIZE(A,1)
  LDA = MAX(1,N)

  IF( PRESENT(UPLO) ) THEN
     LUPLO = UPLO
  ELSE
     LUPLO = 'U'
  END IF

  IF( PRESENT(MTHD) ) THEN
     LMTHD = MTHD
  ELSE
     LMTHD = 'P'
  END IF

  ! Test the arguments
  IF( SIZE(A,2) /= N .OR. N<0 ) THEN
     LINFO = -1
  ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') ) THEN
     LINFO = -3
  ELSE IF( .NOT.LSAME(LMTHD,'P') .AND. .NOT.LSAME(LMTHD,'H') ) THEN
     LINFO = -4
  ELSE IF( N>0 ) THEN
     !Allocate fixed size arrays
     IF( LSAME(LMTHD, 'P') ) THEN
        ALLOCATE(IWORK(N), RWORK(1), STAT = ISTAT)
     ELSE
        ALLOCATE(IWORK(1), RWORK(N-1), STAT = ISTAT)
     END IF

     IF( ISTAT == 0 ) THEN
        !Do a workspace query
        CALL SKPF10( LUPLO, LMTHD, N, A, LDA, PFAFF, IWORK, QWORK, -1, &
                     RWORK, LINFO)
        LWORK = INT(QWORK(1))
        ALLOCATE(WORK(LWORK), STAT=ISTAT)

        IF( ISTAT /= 0 ) THEN
           ! Allocating the workspace failed. Try minimal workspace
           DEALLOCATE(WORK)
           CALL MESSAGE(-200, 'SKPF10') ! print warning about not enough wspace

           IF( LSAME(LMTHD, 'H') ) THEN
              LWORK = N
           ELSE
              LWORK = 1
           END IF

           ALLOCATE(WORK(LWORK), STAT=ISTAT)
        END IF
     END IF

     IF( ISTAT == 0 ) THEN
        CALL SKPF10( LUPLO, LMTHD, N, A, LDA, PFAFF, IWORK, WORK, LWORK, &
                     RWORK, LINFO )
     ELSE
        LINFO=-100 ! not enough memory
     END IF

     DEALLOCATE(WORK)
  ELSE !N=0
     PFAFF(1) = 1.0D0
     PFAFF(2) = 0.0D0
  END IF

  CALL MESSAGE(LINFO, 'SKPF10', INFO, ISTAT)
END SUBROUTINE ZSKPF10_F95
