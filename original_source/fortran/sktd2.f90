!  Purpose
!  =======
!
!  SKTD2 reduces a real or complex skew-symmetric matrix A to real
!  skew-symmetric tridiagonal form T by a unitary congruence
!  transformation (which reduces to an orthogonal similarity
!  transformation in the real case):
!  Q^T * A * Q = T.
!
!  This subroutine uses the unblocked version of the Householder
!  algorithm.
!
!  =========
!
!     SUBROUTINE SKTD2( A, TAU, UPLO, MODE, INFO)
!         <type>(<prec>), INTENT(INOUT) :: A(:,:)
!         <type>(<prec>), INTENT(OUT) :: TAU(:)
!         CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO, MODE
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
!             If MODE = 'P', then size(A,1) and size(A,2) must be even.
!          On entry, the skew-symmetric matrix A.
!            If UPLO = 'U', the leading N-by-N upper triangular part
!              of A contains the upper triangular part of the matrix A,
!              and the strictly lower triangular part of A is not referenced.
!            If UPLO = 'L', the leading N-by-N lower triangular part
!              of A contains the lower triangular part of the matrix A,
!              and the strictly upper triangular part of A is not referenced.
!          On exit, if MODE = 'N' (default):
!            If UPLO = 'U', the diagonal and first superdiagonal
!              of A are overwritten by the corresponding elements of the
!              tridiagonal matrix T, and the elements above the first
!              superdiagonal, with the array TAU, represent the unitary
!              matrix Q as a product of elementary reflectors;
!            If UPLO = 'L', the diagonal and first subdiagonal of A are over-
!              written by the corresponding elements of the tridiagonal
!              matrix T, and the elements below the first subdiagonal, with
!              the array TAU, represent the unitary matrix Q as a product
!              of elementary reflectors.
!            See Further Details, also for information about MODE = 'P'.
!
!  TAU     (output) REAL or COMPLEX array,
!          shape (:), size(TAU) = size(A,1)-1
!          The scalar factors of the elementary reflectors (see Further
!          Details).
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
!  The normal use for SKTD2 is to compute the tridiagonal form of
!  a skew-symmetric matrix under an orthogonal similarity transformation,
!  and chosen by setting MODE = 'N' ("normal" mode). The other
!  use of SKTD2 is the computation the Pfaffian of a skew-symmetric matrix,
!  which only requires a partial tridiagonalization, this mode is chosen
!  by setting MODE = 'P' ("Pfaffian" mode).
!
!  Normal mode (MODE = 'N'):
!  ========================
!
!  The routine computes a tridiagonal matrix T and an orthogonal Q such
!  that A = Q * T * Q^T .
!
!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(n-1) . . . H(2) H(1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
!  A(1:i-1,i+1), and tau in TAU(i).
!
!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(1) H(2) . . . H(n-1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
!  and tau in TAU(i).
!
!  The contents of A on exit are illustrated by the following examples
!  with n = 5:
!
!  if UPLO = 'U':                       if UPLO = 'L':
!
!    (  0   e   v2  v3  v4 )              (  0                  )
!    (      0   e   v3  v4 )              (  e   0              )
!    (          0   e   v4 )              (  v1  e   0          )
!    (              0   e  )              (  v1  v2  e   0      )
!    (                  0  )              (  v1  v2  v3  e   0  )
!
!  where d and e denote diagonal and off-diagonal elements of T, and vi
!  denotes an element of the vector defining H(i).
!
!  The LAPACK routine DORGTR can be used to form the transformation
!  matrix explicitely, and DORMTR can be used to multiply another
!  matrix without forming the transformation.
!
!  Pfaffian mode (MODE = 'P'):
!  ==========================
!
!  For computing the Pfaffian, it is enough to bring A into a partial
!  tridiagonal form. In particular, assuming n even, it is enough to
!  bring A into a form with A(i,j) = A(j,i) = 0 for i > j+1 with j odd
!  (this is computed if UPLO = 'L'), or A(i,j) = A(j,i) = 0 for
!  i > j-1 with j even (this is computed if UPLO = 'U'). Note that
!  only the off-diagonal entries in the odd columns (if UPLO = 'L')
!  or in the even columns (if UPLU = 'U') are properly computed by SKTD2.
!
!  A is brought into this special form pT using an orthogonal matrix Q:
!  A = Q * pT * Q^T
!
!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(n-1) H(n-3) . . . H(3) H(1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v^T
!
!  where tau is a real scalar, and v is a real vector with
!  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
!  A(1:i-1,i+1), and tau in TAU(i).
!
!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(1) H(3) . . . H(n-3) H(n-1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v^T
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
!  and tau in TAU(i).
!
!  The contents of A on exit are illustrated by the following examples
!  with n = 6:
!
!  if UPLO = 'U':                       if UPLO = 'L':
!
!    (  0   e   x   v3  x   v5 )        (  0                      )
!    (      0   x   v3  x   v5 )        (  e   0                  )
!    (          0   e   x   v5 )        (  v1  x   0              )
!    (              0   x   v5 )        (  v1  x   e   0          )
!    (                  0   e  )        (  v1  x   v3  x   0      )
!    (                      0  )        (  v1  x   v3  x   e   0  )
!
!  where d and e denote diagonal and off-diagonal elements of T, vi
!  denotes an element of the vector defining H(i), and x denotes an
!  element not computed by SKTD2.
!
!  =====================================================================
!

SUBROUTINE SSKTD2_F95( A, TAU, UPLO, MODE, INFO )
  USE PFAPACK_PREC, ONLY: precision => singleprec
  USE PFAPACK_MESSAGE, ONLY: MESSAGE
  USE F77_PFAPACK, ONLY: SKTD2
  IMPLICIT NONE
  REAL(precision), INTENT(INOUT) :: A(:,:)
  REAL(precision), INTENT(OUT) :: TAU(:)
  CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO, MODE
  INTEGER, INTENT(OUT), OPTIONAL :: INFO

  CHARACTER(LEN=1) :: LUPLO, LMODE
  INTEGER :: N, LDA, LINFO, ISTAT
  REAL(precision), ALLOCATABLE :: E(:)

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

  IF( PRESENT(MODE) ) THEN
     LMODE = MODE
  ELSE
     LMODE = 'N'
  END IF

  ! Test the arguments
  IF( SIZE(A,2) /= N .OR. N<0 ) THEN
     LINFO = -1
  ELSE IF( N > 0 .AND. SIZE(TAU) /= N-1 ) THEN
     LINFO = -2
  ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') ) THEN
     LINFO = -3
  ELSE IF( .NOT.LSAME(LMODE,'N') .AND. .NOT.LSAME(LMODE,'P') ) THEN
     LINFO = -4
  ELSE IF( N>0 ) THEN
     !Allocate the fixed size arrays first
     ALLOCATE(E(N-1), STAT=ISTAT)
     !if the allocation fails here, it's fatal

     IF( ISTAT == 0 ) THEN
        CALL SKTD2( LUPLO, LMODE, N, A, LDA, E, TAU, LINFO )
     ELSE
        LINFO=-100 ! not enough memory
     END IF

     DEALLOCATE(E)
  END IF

  CALL MESSAGE(LINFO, 'SKTD2', INFO, ISTAT)
END SUBROUTINE SSKTD2_F95


SUBROUTINE DSKTD2_F95( A, TAU, UPLO, MODE, INFO )
  USE PFAPACK_PREC, ONLY: precision => doubleprec
  USE PFAPACK_MESSAGE, ONLY: MESSAGE
  USE F77_PFAPACK, ONLY: SKTD2
  IMPLICIT NONE
  REAL(precision), INTENT(INOUT) :: A(:,:)
  REAL(precision), INTENT(OUT) :: TAU(:)
  CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO, MODE
  INTEGER, INTENT(OUT), OPTIONAL :: INFO

  CHARACTER(LEN=1) :: LUPLO, LMODE
  INTEGER :: N, LDA, LINFO, ISTAT
  REAL(precision), ALLOCATABLE :: E(:)

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

  IF( PRESENT(MODE) ) THEN
     LMODE = MODE
  ELSE
     LMODE = 'N'
  END IF

  ! Test the arguments
  IF( SIZE(A,2) /= N .OR. N<0 ) THEN
     LINFO = -1
  ELSE IF( N > 0 .AND. SIZE(TAU) /= N-1 ) THEN
     LINFO = -2
  ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') ) THEN
     LINFO = -3
  ELSE IF( .NOT.LSAME(LMODE,'N') .AND. .NOT.LSAME(LMODE,'P') ) THEN
     LINFO = -4
  ELSE IF( N>0 ) THEN
     !Allocate the fixed size arrays first
     ALLOCATE(E(N-1), STAT=ISTAT)
     !if the allocation fails here, it's fatal

     IF( ISTAT == 0 ) THEN
        CALL SKTD2( LUPLO, LMODE, N, A, LDA, E, TAU, LINFO )
     ELSE
        LINFO=-100 ! not enough memory
     END IF

     DEALLOCATE(E)
  END IF

  CALL MESSAGE(LINFO, 'SKTD2', INFO, ISTAT)
END SUBROUTINE DSKTD2_F95

SUBROUTINE CSKTD2_F95( A, TAU, UPLO, MODE, INFO )
  USE PFAPACK_PREC, ONLY: precision => singleprec
  USE PFAPACK_MESSAGE, ONLY: MESSAGE
  USE F77_PFAPACK, ONLY: SKTD2
  IMPLICIT NONE
  COMPLEX(precision), INTENT(INOUT) :: A(:,:)
  COMPLEX(precision), INTENT(OUT) :: TAU(:)
  CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO, MODE
  INTEGER, INTENT(OUT), OPTIONAL :: INFO

  CHARACTER(LEN=1) :: LUPLO, LMODE
  INTEGER :: N, LDA, LINFO, ISTAT
  REAL(precision), ALLOCATABLE :: E(:)

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

  IF( PRESENT(MODE) ) THEN
     LMODE = MODE
  ELSE
     LMODE = 'N'
  END IF

  ! Test the arguments
  IF( SIZE(A,2) /= N .OR. N<0 ) THEN
     LINFO = -1
  ELSE IF( N > 0 .AND. SIZE(TAU) /= N-1 ) THEN
     LINFO = -2
  ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') ) THEN
     LINFO = -3
  ELSE IF( .NOT.LSAME(LMODE,'N') .AND. .NOT.LSAME(LMODE,'P') ) THEN
     LINFO = -4
  ELSE IF( N>0 ) THEN
     !Allocate the fixed size arrays first
     ALLOCATE(E(N-1), STAT=ISTAT)
     !if the allocation fails here, it's fatal

     IF( ISTAT == 0 ) THEN
        CALL SKTD2( LUPLO, LMODE, N, A, LDA, E, TAU, LINFO )
     ELSE
        LINFO=-100 ! not enough memory
     END IF

     DEALLOCATE(E)
  END IF

  CALL MESSAGE(LINFO, 'SKTD2', INFO, ISTAT)
END SUBROUTINE CSKTD2_F95


SUBROUTINE ZSKTD2_F95( A, TAU, UPLO, MODE, INFO )
  USE PFAPACK_PREC, ONLY: precision => doubleprec
  USE PFAPACK_MESSAGE, ONLY: MESSAGE
  USE F77_PFAPACK, ONLY: SKTD2
  IMPLICIT NONE
  COMPLEX(precision), INTENT(INOUT) :: A(:,:)
  COMPLEX(precision), INTENT(OUT) :: TAU(:)
  CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO, MODE
  INTEGER, INTENT(OUT), OPTIONAL :: INFO

  CHARACTER(LEN=1) :: LUPLO, LMODE
  INTEGER :: N, LDA, LINFO, ISTAT
  REAL(precision), ALLOCATABLE :: E(:)

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

  IF( PRESENT(MODE) ) THEN
     LMODE = MODE
  ELSE
     LMODE = 'N'
  END IF

  ! Test the arguments
  IF( SIZE(A,2) /= N .OR. N<0 ) THEN
     LINFO = -1
  ELSE IF( N > 0 .AND. SIZE(TAU) /= N-1 ) THEN
     LINFO = -2
  ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') ) THEN
     LINFO = -3
  ELSE IF( .NOT.LSAME(LMODE,'N') .AND. .NOT.LSAME(LMODE,'P') ) THEN
     LINFO = -4
  ELSE IF( N>0 ) THEN
     !Allocate the fixed size arrays first
     ALLOCATE(E(N-1), STAT=ISTAT)
     !if the allocation fails here, it's fatal

     IF( ISTAT == 0 ) THEN
        CALL SKTD2( LUPLO, LMODE, N, A, LDA, E, TAU, LINFO )
     ELSE
        LINFO=-100 ! not enough memory
     END IF

     DEALLOCATE(E)
  END IF

  CALL MESSAGE(LINFO, 'SKTD2', INFO, ISTAT)
END SUBROUTINE ZSKTD2_F95
