!  Purpose
!  =======
!
!  SKBTRD reduces a real or complex skew-symmetric band matrix A to real
!  skew-symmetric tridiagonal form T by an unitary congruence
!  transformation (which reduces to an orthogonal similarity
!  transfromation in the real case):
!     Q^dagger * A * Q^* = T.
!
!  =========
!
!     SUBROUTINE SKBTRD( AB, UPLO, MODE, DETQ, Q, VECT, INFO)
!         <type>(<prec>), INTENT(INOUT) :: AB(:,:)
!         <type>(<prec>), INTENT(OUT), OPTIONAL :: DETQ, Q(:,:)
!         CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO, MODE, VECT
!         INTEGER, INTENT(OUT), OPTIONAL :: INFO
!
!    where
!       <type> ::= REAL | COMPLEX
!       <prec> ::= KIND(1.0) | KIND(1.0D0)
!
!  Arguments
!  =========
!
!  AB      (input/output) REAL or COMPLEX array,
!          shape (:,:), with size(A,1) = KD+1 and size(A,2) = N,
!          where KD is the number of super- and subdiagonals and
!          N >= 0 the order of the matrix.
!             If MODE = 'P', then size(A,1) and size(A,2) must be even.
!          On entry, the upper or lower triangle of the skew-symmetric
!          band matrix A, stored in the first KD+1 rows of the array.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!            if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!            if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=min(n,j+kd).
!          On exit, the zero diagonal elements of AB are left unchanged,
!          if KD > 0, the elements on the first superdiagonal (if UPLO =
!          'U') or the first subdiagonal (if UPLO = 'L') are overwritten
!          by the off-diagonal elements of T; the rest of AB is
!          overwritten by values generated during the reduction. If
!          MODE = 'P', only the off-diagonal entries in the odd rows
!          (columns) are computed for UPLO = 'U' (UPLO = 'L').
!
!  UPLO    Optional (input) CHARACTER*1
!          If UPLO is present, then:
!             = 'U':  Upper triangle of A is stored;
!             = 'L':  Lower triangle of A is stored.
!          otherwise UPLO = 'U' is assumed
!
!  MODE    Optional (input) CHARACTER*1
!          If MODE is present, then:
!             = 'N':  A is fully tridiagonalized
!             = 'P':  A is partially tridiagonalized for Pfaffian computation
!          otherwise MODE = 'N' is assumed
!
!  DETQ    Optional (output) REAL or COMPLEX
!          The value of the determinant of Q, which is a
!          pure phase factor (in the real case, DETQ=1 always).
!          Computing DETQ does not require to form Q explicitly.
!
!  Q       Optional (input/output) REAL or COMPLEX array,
!          shape(:,:), size(Q,1) == size(Q,2) == size(AB,2)
!          On entry:
!            if VECT = 'U', then Q must contain an N-by-N matrix X;
!            if VECT = 'N' or 'V', then Q need not be set.
!          On exit:
!             if VECT = 'V', Q contains the N-by-N unitary matrix Q;
!             if VECT = 'U', Q contains the product X*Q;
!             if VECT = 'N', the array Q is not referenced.
!
!  VECT    Optional (input) CHARACTER*1
!          If VECT is presen:
!            = 'N':  do not form Q;
!            = 'V':  form Q;
!            = 'U':  update a matrix X, by forming X*Q.
!         otherwise, VECT = 'N' is assumed if Q is not present,
!         and VECT='V' if Q is present.
!
!  INFO    Optional (output) INTEGER
!          If INFO is present:
!              = 0:    successful exit
!              < 0:    if INFO = -i, the i-th argument had an illegal value
!              = -100: failed to allocate enough internal memory
!          Otherwise, if INFO is not present and < 0, the program stops with
!          an error messafe
!
!  Further Details
!  ===============
!
!  The storage scheme for the skew-symmetric matrix is identical to the
!  storage scheme for symmetric or Hermitian band matrices in LAPACK,
!  i.e. the diagonal and the KD super- or subdiagonals are stored in an
!  array with KD+1 rows and N columns. Note that the zero diagonal must
!  also be explicitely stored (this was done to keep the structure of
!  the program identical to the symmetric case)
!
!  In particular this means that if
!  - UPLO = 'U', then AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j
!
!    Example: N=5, KD=2
!
!  (  0     a12   a13              )
!  (  -a12  0     a23   a24        )
!  (  -a13  -a23  0     a34   a35  )
!  (        -a24  -a34  0     a45  )
!  (              -a35  -a45  0    )
!
!  is stored as
!
!   x   x   a13 a24 a35
!   x   a12 a23 a34 a45
!   0   0   0   0   0
!
!  where x denotes an unused entry
!
!  - UPLO = 'L', then AB(1+i-j,j) = A(i,j) for j<=i<=min(n,j+kd)
!
!    Example: N=5, KD=2
!
!  (  0     -a21  -a31              )
!  (  a21   0     -a32  -a42        )
!  (  a31   a32   0     -a43  -a53  )
!  (        a42   a43   0     -a54  )
!  (              a53   a54   0     )
!
!  is stored as
!
!   0   0   0   0   0
!   a21 a32 a43 a54 x
!   a31 a42 a53 x   x
!
!  where x denotes an unused entry
!
!  =====================================================================

SUBROUTINE SSKBTRD_F95( AB, UPLO, MODE, DETQ, Q, VECT, INFO )
  USE PFAPACK_PREC, ONLY: precision => singleprec
  USE F77_PFAPACK, ONLY: SKBTRD
  IMPLICIT NONE
  REAL(precision), INTENT(INOUT) :: AB(:,:)
  CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO, MODE, VECT
  REAL(precision), INTENT(OUT), OPTIONAL :: DETQ
  REAL(precision), INTENT(OUT), OPTIONAL :: Q(:,:)
  INTEGER, INTENT(OUT), OPTIONAL :: INFO

  CHARACTER(LEN=1) :: LUPLO, LMODE, LVECT
  INTEGER :: N, KD, LDAB, LINFO, ISTAT
  REAL(precision) :: TEMPQ(1,1)
  REAL(precision), ALLOCATABLE :: E(:), WORK(:)

  LOGICAL LSAME
  EXTERNAL LSAME

  ! Figure out array sizes and set defaults
  LINFO = 0
  KD = SIZE(AB,1)-1
  LDAB = MAX(1,KD+1)
  N= SIZE(AB,2)

  IF( PRESENT(UPLO) ) THEN
     LUPLO=UPLO
  ELSE
     LUPLO='U'
  END IF

  IF( PRESENT(MODE) ) THEN
     LMODE=MODE
  ELSE
     LMODE='N'
  END IF

  IF( .NOT.PRESENT(Q) .AND. .NOT.PRESENT(VECT) ) THEN
     LVECT='N'
  ELSE IF( .NOT.PRESENT(VECT) ) THEN ! .AND. PRESENT(Q)
     LVECT='V'
  ELSE
     LVECT=VECT
  END IF

  ! Test the arguments
  IF( KD<0 .OR. N<0 ) THEN
     LINFO = -1
  ELSE IF( LSAME(LMODE, 'P') .AND. MOD(N,2) /= 0) THEN
     LINFO = -1
  ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') ) THEN
     LINFO = -2
  ELSE IF( .NOT.LSAME(LMODE,'N') .AND. .NOT.LSAME(LMODE,'P') ) THEN
     LINFO = -3
  ELSE IF( PRESENT(Q) ) THEN
     IF( SIZE(Q,1) /= SIZE(Q,2) .OR. SIZE(Q,1) /= N ) THEN
        LINFO = -5
     ENDIF
  ELSE IF( (LSAME(LVECT,'V') .OR. LSAME(LVECT,'U')) &
           .AND. .NOT.PRESENT(Q) ) THEN
     LINFO = -5
  END IF

  IF( LINFO ==0 .AND. N>0 ) THEN
     !Allocate the fixed size arrays
     ALLOCATE(E(N-1), WORK(2*N), STAT=ISTAT)
     !if the allocation fails here, it's fatal

     IF( ISTAT == 0) THEN
        IF( PRESENT(Q) ) THEN
           CALL SKBTRD( LVECT, LUPLO, LMODE, N, KD, AB, LDAB, E, Q, N, WORK, LINFO)
        ELSE
           !Use a dummy Q
           CALL SKBTRD( LVECT, LUPLO, LMODE, N, KD, AB, LDAB, E, TEMPQ, 1, WORK, LINFO)
        END IF
     ELSE
        LINFO=-100 ! not enough memory
     END IF

     DEALLOCATE(WORK, E)
  END IF

  IF( PRESENT(DETQ) ) THEN
     DETQ=1 !in the real case, det(Q) is always 1
  END IF

  CALL MESSAGE(LINFO, 'SKBTRD', INFO, ISTAT)
END SUBROUTINE SSKBTRD_F95


SUBROUTINE DSKBTRD_F95( AB, UPLO, MODE, DETQ, Q, VECT, INFO )
  USE PFAPACK_PREC, ONLY: precision => doubleprec
  USE F77_PFAPACK, ONLY: SKBTRD
  IMPLICIT NONE
  REAL(precision), INTENT(INOUT) :: AB(:,:)
  CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO, MODE, VECT
  REAL(precision), INTENT(OUT), OPTIONAL :: DETQ
  REAL(precision), INTENT(OUT), OPTIONAL :: Q(:,:)
  INTEGER, INTENT(OUT), OPTIONAL :: INFO

  CHARACTER(LEN=1) :: LUPLO, LMODE, LVECT
  INTEGER :: N, KD, LDAB, LINFO, ISTAT
  REAL(precision) :: TEMPQ(1,1)
  REAL(precision), ALLOCATABLE :: E(:), WORK(:)

  LOGICAL LSAME
  EXTERNAL LSAME

  ! Figure out array sizes and set defaults
  LINFO = 0
  KD = SIZE(AB,1)-1
  LDAB = MAX(1,KD+1)
  N= SIZE(AB,2)

  IF( PRESENT(UPLO) ) THEN
     LUPLO=UPLO
  ELSE
     LUPLO='U'
  END IF

  IF( PRESENT(MODE) ) THEN
     LMODE=MODE
  ELSE
     LMODE='N'
  END IF

  IF( .NOT.PRESENT(Q) .AND. .NOT.PRESENT(VECT) ) THEN
     LVECT='N'
  ELSE IF( .NOT.PRESENT(VECT) ) THEN ! .AND. PRESENT(Q)
     LVECT='V'
  ELSE
     LVECT=VECT
  END IF

  ! Test the arguments
  IF( KD<0 .OR. N<0 ) THEN
     LINFO = -1
  ELSE IF( LSAME(LMODE, 'P') .AND. MOD(N,2) /= 0) THEN
     LINFO = -1
  ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') ) THEN
     LINFO = -2
  ELSE IF( .NOT.LSAME(LMODE,'N') .AND. .NOT.LSAME(LMODE,'P') ) THEN
     LINFO = -3
  ELSE IF( PRESENT(Q) ) THEN
     IF( SIZE(Q,1) /= SIZE(Q,2) .OR. SIZE(Q,1) /= N ) THEN
        LINFO = -5
     ENDIF
  ELSE IF( (LSAME(LVECT,'V') .OR. LSAME(LVECT,'U')) &
           .AND. .NOT.PRESENT(Q) ) THEN
     LINFO = -5
  END IF

  IF( LINFO ==0 .AND. N>0 ) THEN
     !Allocate the fixed size arrays
     ALLOCATE(E(N-1), WORK(2*N), STAT=ISTAT)
     !if the allocation fails here, it's fatal

     IF( ISTAT == 0) THEN
        IF( PRESENT(Q) ) THEN
           CALL SKBTRD( LVECT, LUPLO, LMODE, N, KD, AB, LDAB, E, Q, N, WORK, LINFO)
        ELSE
           !Use a dummy Q
           CALL SKBTRD( LVECT, LUPLO, LMODE, N, KD, AB, LDAB, E, TEMPQ, 1, WORK, LINFO)
        END IF
     ELSE
        LINFO=-100 ! not enough memory
     END IF

     DEALLOCATE(WORK, E)
  END IF

  IF( PRESENT(DETQ) ) THEN
     DETQ=1 !in the real case, det(Q) is always 1
  END IF

  CALL MESSAGE(LINFO, 'SKBTRD', INFO, ISTAT)
END SUBROUTINE DSKBTRD_F95

SUBROUTINE CSKBTRD_F95( AB, UPLO, MODE, DETQ, Q, VECT, INFO )
  USE PFAPACK_PREC, ONLY: precision => singleprec
  USE F77_PFAPACK, ONLY: SKBTRD
  IMPLICIT NONE
  COMPLEX(precision), INTENT(INOUT) :: AB(:,:)
  CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO, VECT, MODE
  COMPLEX(precision), INTENT(OUT), OPTIONAL :: DETQ
  COMPLEX(precision), INTENT(OUT), OPTIONAL :: Q(:,:)
  INTEGER, INTENT(OUT), OPTIONAL :: INFO

  CHARACTER(LEN=1) :: LUPLO, LVECT, LMODE
  INTEGER :: N, KD, LDAB, LINFO, ISTAT
  COMPLEX(precision) :: TEMPQ(1,1), TEMPDETQ
  REAL(precision), ALLOCATABLE :: E(:), RWORK(:)
  COMPLEX(precision), ALLOCATABLE :: WORK(:)

  LOGICAL LSAME
  EXTERNAL LSAME

  ! Figure out array sizes and set defaults
  LINFO = 0
  KD = SIZE(AB,1)-1
  LDAB = MAX(1,KD+1)
  N= SIZE(AB,2)

  IF( PRESENT(UPLO) ) THEN
     LUPLO=UPLO
  ELSE
     LUPLO='U'
  END IF

  IF( PRESENT(MODE) ) THEN
     LMODE=MODE
  ELSE
     LMODE='N'
  END IF

  IF( .NOT.PRESENT(Q) .AND. .NOT.PRESENT(VECT) ) THEN
     LVECT='N'
  ELSE IF( .NOT.PRESENT(VECT) ) THEN ! .AND. PRESENT(Q)
     LVECT='V'
  ELSE
     LVECT=VECT
  END IF

  ! Test the arguments
  IF( KD<0 .OR. N<0 ) THEN
     LINFO = -1
  ELSE IF( LSAME(LMODE, 'P') .AND. MOD(N,2) /= 0) THEN
     LINFO = -1
  ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') ) THEN
     LINFO = -2
  ELSE IF( .NOT.LSAME(LMODE,'N') .AND. .NOT.LSAME(LMODE,'P') ) THEN
     LINFO = -3
  ELSE IF( PRESENT(Q) ) THEN
     IF( SIZE(Q,1) /= SIZE(Q,2) .OR. SIZE(Q,1) /= N ) THEN
        LINFO = -5
     ENDIF
  ELSE IF( (LSAME(LVECT,'V') .OR. LSAME(LVECT,'U')) &
           .AND. .NOT.PRESENT(Q) ) THEN
     LINFO = -5
  END IF

  IF( LINFO ==0 .AND. N>0 ) THEN
     !Allocate the fixed size arrays
     ALLOCATE(E(N-1), WORK(N), RWORK(N), STAT=ISTAT)
     !if the allocation fails here, it's fatal

     IF( ISTAT == 0) THEN
        IF( PRESENT(Q) ) THEN
           IF( PRESENT(DETQ) ) THEN
              CALL SKBTRD( LVECT, LUPLO, LMODE, N, KD, AB, LDAB, E, DETQ, &
                           Q, N, WORK, RWORK, LINFO)
           ELSE
              CALL SKBTRD( LVECT, LUPLO, LMODE, N, KD, AB, LDAB, E, TEMPDETQ, &
                           Q, N, WORK, RWORK, LINFO)
           END IF
        ELSE
           IF( PRESENT(DETQ) ) THEN
              CALL SKBTRD( LVECT, LUPLO, LMODE, N, KD, AB, LDAB, E, DETQ, &
                           TEMPQ, 1, WORK, RWORK, LINFO)
           ELSE
              CALL SKBTRD( LVECT, LUPLO, LMODE, N, KD, AB, LDAB, E, TEMPDETQ, &
                           TEMPQ, 1, WORK, RWORK, LINFO)
           END IF
        END IF
     ELSE
        LINFO=-100 ! not enough memory
     END IF

     DEALLOCATE(RWORK, WORK, E)
  END IF

  CALL MESSAGE(LINFO, 'SKBTRD', INFO, ISTAT)
END SUBROUTINE CSKBTRD_F95

SUBROUTINE ZSKBTRD_F95( AB, UPLO, MODE, DETQ, Q, VECT, INFO )
  USE PFAPACK_PREC, ONLY: precision => doubleprec
  USE F77_PFAPACK, ONLY: SKBTRD
  IMPLICIT NONE
  COMPLEX(precision), INTENT(INOUT) :: AB(:,:)
  CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO, VECT, MODE
  COMPLEX(precision), INTENT(OUT), OPTIONAL :: DETQ
  COMPLEX(precision), INTENT(OUT), OPTIONAL :: Q(:,:)
  INTEGER, INTENT(OUT), OPTIONAL :: INFO

  CHARACTER(LEN=1) :: LUPLO, LVECT, LMODE
  INTEGER :: N, KD, LDAB, LINFO, ISTAT
  COMPLEX(precision) :: TEMPQ(1,1), TEMPDETQ
  REAL(precision), ALLOCATABLE :: E(:), RWORK(:)
  COMPLEX(precision), ALLOCATABLE :: WORK(:)

  LOGICAL LSAME
  EXTERNAL LSAME

  ! Figure out array sizes and set defaults
  LINFO = 0
  KD = SIZE(AB,1)-1
  LDAB = MAX(1,KD+1)
  N= SIZE(AB,2)

  IF( PRESENT(UPLO) ) THEN
     LUPLO=UPLO
  ELSE
     LUPLO='U'
  END IF

  IF( PRESENT(MODE) ) THEN
     LMODE=MODE
  ELSE
     LMODE='N'
  END IF

  IF( .NOT.PRESENT(Q) .AND. .NOT.PRESENT(VECT) ) THEN
     LVECT='N'
  ELSE IF( .NOT.PRESENT(VECT) ) THEN ! .AND. PRESENT(Q)
     LVECT='V'
  ELSE
     LVECT=VECT
  END IF

  ! Test the arguments
  IF( KD<0 .OR. N<0 ) THEN
     LINFO = -1
  ELSE IF( LSAME(LMODE, 'P') .AND. MOD(N,2) /= 0) THEN
     LINFO = -1
  ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') ) THEN
     LINFO = -2
  ELSE IF( .NOT.LSAME(LMODE,'N') .AND. .NOT.LSAME(LMODE,'P') ) THEN
     LINFO = -3
  ELSE IF( PRESENT(Q) ) THEN
     IF( SIZE(Q,1) /= SIZE(Q,2) .OR. SIZE(Q,1) /= N ) THEN
        LINFO = -5
     ENDIF
  ELSE IF( (LSAME(LVECT,'V') .OR. LSAME(LVECT,'U')) &
           .AND. .NOT.PRESENT(Q) ) THEN
     LINFO = -5
  END IF

  IF( LINFO ==0 .AND. N>0 ) THEN
     !Allocate the fixed size arrays
     ALLOCATE(E(N-1), WORK(N), RWORK(N), STAT=ISTAT)
     !if the allocation fails here, it's fatal

     IF( ISTAT == 0) THEN
        IF( PRESENT(Q) ) THEN
           IF( PRESENT(DETQ) ) THEN
              CALL SKBTRD( LVECT, LUPLO, LMODE, N, KD, AB, LDAB, E, DETQ, &
                           Q, N, WORK, RWORK, LINFO)
           ELSE
              CALL SKBTRD( LVECT, LUPLO, LMODE, N, KD, AB, LDAB, E, TEMPDETQ, &
                           Q, N, WORK, RWORK, LINFO)
           END IF
        ELSE
           IF( PRESENT(DETQ) ) THEN
              CALL SKBTRD( LVECT, LUPLO, LMODE, N, KD, AB, LDAB, E, DETQ, &
                           TEMPQ, 1, WORK, RWORK, LINFO)
           ELSE
              CALL SKBTRD( LVECT, LUPLO, LMODE, N, KD, AB, LDAB, E, TEMPDETQ, &
                           TEMPQ, 1, WORK, RWORK, LINFO)
           END IF
        END IF
     ELSE
        LINFO=-100 ! not enough memory
     END IF

     DEALLOCATE(RWORK, WORK, E)
  END IF

  CALL MESSAGE(LINFO, 'SKBTRD', INFO, ISTAT)
END SUBROUTINE ZSKBTRD_F95
