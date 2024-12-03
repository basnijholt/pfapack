!  Purpose
!  =======
!
!  SKBPFA computes the Pfaffian of a banded skew-symmetric matrix.
!
!  =========
!
!     SUBROUTINE SKBPFA( AB, PFAFF, UPLO, INFO)
!         <type>(<prec>), INTENT(INOUT) :: A(:,:)
!         <type>(<prec>), INTENT(OUT) :: PFAFF
!         CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
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
!          If N is odd, AB is not referenced.
!          If N is even:
!          On entry, the upper or lower triangle of the skew-symmetric
!          band matrix A, stored in the first KD+1 rows of the array.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!            if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!            if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=min(n,j+kd).
!          On exit, AB is overwritten with values generated during the
!          computation.
!
!  PFAFF   (output) REAL or COMPLEX
!          The value of the Pfaffian.
!
!  UPLO    Optional (input) CHARACTER*1
!          If UPLO is present, then:
!            = 'U':  Upper triangle of A is stored;
!            = 'L':  Lower triangle of A is stored.
!          otherwise UPLO = 'U' is assumed
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
!  For odd-sized matrices, the Pfaffian is 0 by default, hence no
!  computation is needed in this case. For even-sized matrices,
!  the Pfaffian is computed by bringing the skew-symmetric matrix A into
!  tridiagonal form T by a unitary congruence transformation:
!  Q^H * A * Q^* = T.
!  This transformation is computed by the routine SKBTRD (for further
!  details see there)

SUBROUTINE SSKBPFA_F95( AB, PFAFF, UPLO, INFO )
  USE PFAPACK_PREC, ONLY: precision => singleprec
  USE PFAPACK_MESSAGE, ONLY: MESSAGE
  USE F77_PFAPACK, ONLY: SKBPFA
  IMPLICIT NONE
  CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
  INTEGER, INTENT(OUT), OPTIONAL :: INFO
  REAL(precision), INTENT(OUT) :: PFAFF
  REAL(precision), INTENT(INOUT) :: AB(:,:)

  CHARACTER(LEN=1) :: LUPLO
  INTEGER :: N, KD, LDAB, LINFO, ISTAT
  REAL(precision), ALLOCATABLE :: WORK(:)

  LOGICAL LSAME
  EXTERNAL LSAME

  ! Figure out array sizes and set defaults
  LINFO = 0
  KD = SIZE(AB,1)-1
  N = SIZE(AB,2)
  LDAB = MAX(1,KD+1)

  IF( PRESENT(UPLO) ) THEN
     LUPLO = UPLO
  ELSE
     LUPLO = 'U'
  END IF

  ! Test the arguments
  IF( KD<0 .OR. N<0 ) THEN
     LINFO = -1
  ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') ) THEN
     LINFO = -3
  ELSE IF( N>0 ) THEN
     !Allocate fixed size arrays
     ALLOCATE(WORK(3*N-1), STAT = ISTAT)

     IF( ISTAT == 0) THEN
        CALL SKBPFA( LUPLO, N, KD, AB, LDAB, PFAFF, WORK, LINFO )
     ELSE
        LINFO=-100 ! not enough memory
     END IF

     DEALLOCATE(WORK)
  ELSE !N=0
     PFAFF = 1.0E0
  END IF

  CALL MESSAGE(LINFO, 'SKBPFA', INFO, ISTAT)
END SUBROUTINE SSKBPFA_F95

SUBROUTINE DSKBPFA_F95( AB, PFAFF, UPLO, INFO )
  USE PFAPACK_PREC, ONLY: precision => doubleprec
  USE PFAPACK_MESSAGE, ONLY: MESSAGE
  USE F77_PFAPACK, ONLY: SKBPFA
  IMPLICIT NONE
  CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
  INTEGER, INTENT(OUT), OPTIONAL :: INFO
  REAL(precision), INTENT(OUT) :: PFAFF
  REAL(precision), INTENT(INOUT) :: AB(:,:)

  CHARACTER(LEN=1) :: LUPLO
  INTEGER :: N, KD, LDAB, LINFO, ISTAT
  REAL(precision), ALLOCATABLE :: WORK(:)

  LOGICAL LSAME
  EXTERNAL LSAME

  ! Figure out array sizes and set defaults
  LINFO = 0
  KD = SIZE(AB,1)-1
  N = SIZE(AB,2)
  LDAB = MAX(1,KD+1)

  IF( PRESENT(UPLO) ) THEN
     LUPLO = UPLO
  ELSE
     LUPLO = 'U'
  END IF

  ! Test the arguments
  IF( KD<0 .OR. N<0 ) THEN
     LINFO = -1
  ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') ) THEN
     LINFO = -3
  ELSE IF( N>0 ) THEN
     !Allocate fixed size arrays
     ALLOCATE(WORK(3*N-1), STAT = ISTAT)

     IF( ISTAT == 0) THEN
        CALL SKBPFA( LUPLO, N, KD, AB, LDAB, PFAFF, WORK, LINFO )
     ELSE
        LINFO=-100 ! not enough memory
     END IF

     DEALLOCATE(WORK)
  ELSE !N=0
     PFAFF = 1.0D0
  END IF

  CALL MESSAGE(LINFO, 'SKBPFA', INFO, ISTAT)
END SUBROUTINE DSKBPFA_F95

SUBROUTINE CSKBPFA_F95( AB, PFAFF, UPLO, INFO )
  USE PFAPACK_PREC, ONLY: precision => singleprec
  USE PFAPACK_MESSAGE, ONLY: MESSAGE
  USE F77_PFAPACK, ONLY: SKBPFA
  IMPLICIT NONE
  CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
  INTEGER, INTENT(OUT), OPTIONAL :: INFO
  COMPLEX(precision), INTENT(OUT) :: PFAFF
  COMPLEX(precision), INTENT(INOUT) :: AB(:,:)

  CHARACTER(LEN=1) :: LUPLO
  INTEGER :: N, KD, LDAB, LINFO, ISTAT
  COMPLEX(precision), ALLOCATABLE :: WORK(:)
  REAL(precision), ALLOCATABLE :: RWORK(:)
  LOGICAL LSAME
  EXTERNAL LSAME

  ! Figure out array sizes and set defaults
  LINFO = 0
  KD = SIZE(AB,1)-1
  N = SIZE(AB,2)
  LDAB = MAX(1,KD+1)

  IF( PRESENT(UPLO) ) THEN
     LUPLO = UPLO
  ELSE
     LUPLO = 'U'
  END IF

  ! Test the arguments
  IF( KD<0 .OR. N<0 ) THEN
     LINFO = -1
  ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') ) THEN
     LINFO = -3
  ELSE IF( N>0 ) THEN
     !Allocate fixed size arrays
     ALLOCATE(WORK(N), RWORK(2*N-1), STAT = ISTAT)

     IF( ISTAT == 0) THEN
        CALL SKBPFA( LUPLO, N, KD, AB, LDAB, PFAFF, WORK, RWORK, LINFO )
     ELSE
        LINFO=-100 ! not enough memory
     END IF

     DEALLOCATE(WORK, RWORK)
  ELSE !N=0
     PFAFF = 1.0E0
  END IF

  CALL MESSAGE(LINFO, 'SKBPFA', INFO, ISTAT)
END SUBROUTINE CSKBPFA_F95

SUBROUTINE ZSKBPFA_F95( AB, PFAFF, UPLO, INFO )
  USE PFAPACK_PREC, ONLY: precision => doubleprec
  USE PFAPACK_MESSAGE, ONLY: MESSAGE
  USE F77_PFAPACK, ONLY: SKBPFA
  IMPLICIT NONE
  CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
  INTEGER, INTENT(OUT), OPTIONAL :: INFO
  COMPLEX(precision), INTENT(OUT) :: PFAFF
  COMPLEX(precision), INTENT(INOUT) :: AB(:,:)

  CHARACTER(LEN=1) :: LUPLO
  INTEGER :: N, KD, LDAB, LINFO, ISTAT
  COMPLEX(precision), ALLOCATABLE :: WORK(:)
  REAL(precision), ALLOCATABLE :: RWORK(:)
  LOGICAL LSAME
  EXTERNAL LSAME

  ! Figure out array sizes and set defaults
  LINFO = 0
  KD = SIZE(AB,1)-1
  N = SIZE(AB,2)
  LDAB = MAX(1,KD+1)

  IF( PRESENT(UPLO) ) THEN
     LUPLO = UPLO
  ELSE
     LUPLO = 'U'
  END IF

  ! Test the arguments
  IF( KD<0 .OR. N<0 ) THEN
     LINFO = -1
  ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') ) THEN
     LINFO = -3
  ELSE IF( N>0 ) THEN
     !Allocate fixed size arrays
     ALLOCATE(WORK(N), RWORK(2*N-1), STAT = ISTAT)

     IF( ISTAT == 0) THEN
        CALL SKBPFA( LUPLO, N, KD, AB, LDAB, PFAFF, WORK, RWORK, LINFO )
     ELSE
        LINFO=-100 ! not enough memory
     END IF

     DEALLOCATE(WORK, RWORK)
  ELSE !N=0
     PFAFF = 1.0D0
  END IF

  CALL MESSAGE(LINFO, 'SKBPFA', INFO, ISTAT)
END SUBROUTINE ZSKBPFA_F95
