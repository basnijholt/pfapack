!  Purpose
!  =======
!
!  SKTRF computes the factorization of a skew-symmetric matrix A
!  using the Parlett-Reid algorithm:
!
!     P*A*P^T = U*T*U^T  or  P*A*P^T = L*T*L^T
!
!  where U (or L) unit upper (lower) triangular matrix (^T denotes
!  the transpose), T is a skew-symmetric tridiagonal matrix and P
!  is a permutation matrix. In addition to being unit triangular,
!  U(1:n-1,n)=0 and L(2:n,1)=0.
!  Instead of a full tridiagonalization, SKTRF can also compute a
!  partial tridiagonal form for computing the Pfaffian.
!
!  This subroutine uses, if enough memory is available, the
!  blocked version of the algorithm.
!
!  =========
!
!     SUBROUTINE SKTRF( A, IPIV, UPLO, MODE, INFO)
!         <type>(<prec>), INTENT(INOUT) :: A(:,:)
!         INTEGER, INTENT(OUT) :: IPIV(:)
!         CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO, MODE
!         INTEGER, INTENT(OUT), OPTIONAL :: INFO
!
!    where
!       <type> ::= REAL | COMPLEX
!       <prec> ::= KIND(1.0) | KIND(1.0D0)
!
!
!  Arguments
!  =========
!
!  A       (input/output) REAL or COMPLEX array,
!          shape (:,:), size(A,1) == size(A,2) >= 0.
!             If MODE = 'P', then size(A,1) and size(A,2) must be even.
!          On entry, the symmetric matrix A.
!             If UPLO = 'U', the leading n-by-n upper triangular part
!                of A contains the upper triangular part of the matrix A,
!                and the strictly lower triangular part of A is not referenced.
!             If UPLO = 'L', the leading n-by-n lower triangular part
!                of A contains the lower triangular part of the matrix A,
!                and the strictly upper triangular part of A is not referenced.
!          On exit, the tridiagonal matrix T and the multipliers used
!          to obtain the factor U or L (see below for further details).
!
!  IPIV    (output) INTEGER array, dimension (N)
!          Information about the permutation matrix P: row and column
!          i are interchanged with IPIV(i). If UPLO = 'U', those
!          interchanges are done in the order i = N ... 1, if UPLO = 'L'
!          in the order i = 1 ... N.
!
!  UPLO    Optional (input) CHARACTER*1
!          If UPLO is present, it specifies whether the upper or lower
!          triangular part of the skew-symmetric matrix A is stored:
!            = 'U':  Upper triangular
!            = 'L':  Lower triangular
!          Otherwise, UPLO = 'U' is assumed.
!
!  MODE    Optional (input) CHARACTER*1
!          If MODE is present, then:
!            = 'N':  A is fully tridiagonalized
!            = 'P':  A is partially tridiagonalized for Pfaffian computation
!                   (details see below)
!          otherwise MODE = 'N' is assumed
!
!  INFO    Optional (output) INTEGER
!          If INFO is present:
!            = 0: successful exit
!            < 0: if INFO = -k, the k-th argument had an illegal value
!            > 0: if INFO = k, the off-diagonal entry in the k-th row
!                              (UPLO = 'U') or k-th column (UPLO = 'L')
!                              is exactly zero.
!            = -100: failed to allocate enough internal memory
!          Otherwise, if INFO is not present and < 0, the program stops with
!          an error message.
!
!  Further Details
!  ===============
!
!  The normal use for SKTRF is to compute the U T U^T or L T L^T
!  decomposition of a skew-symmetric matrix with pivoting. This mode
!  is chosen by setting MODE = 'N' ("normal" mode). The other
!  use of SKTRF is the computation the Pfaffian of a skew-symmetric matrix,
!  which only requires a partial computation of T, this mode is chosen
!  by setting MODE = 'P' ("Pfaffian" mode).
!
!  Normal mode (MODE = 'N'):
!  ========================
!
!  If UPLO = 'U', the U*T*U^T decomposition of A is computed. U is a
!  upper triangular unit matrix with the additional constraint
!  U(1:n-1,n) = 0, and T a tridiagonal matrix. The upper diagonal
!  of T is stored on exit in A(i,i+1) for i = 1 .. n-1. The column
!  U(1:i-1, i) is stored in A(1:i-1,i+1).
!
!  If UPLO = 'L', the L*T*L^T decomposition of A is computed. L is a
!  lower triangular unit matrix with the additional constraint
!  L(2:n,1) = 0, and T a tridiagonal matrix. The lower diagonal
!  of T is stored on exit in A(i+1,i) for i = 1 .. n-1. The column
!  L(i+1:n, i) is stored in A(i+1:n,i-1).
!
!  The contents of A on exit are illustrated by the following examples
!  with n = 5:
!
!  if UPLO = 'U':                       if UPLO = 'L':
!
!    (  0   e   u2  u3  u4 )              (  0                  )
!    (      0   e   u3  u4 )              (  e   0              )
!    (          0   e   u4 )              (  l2  e   0          )
!    (              0   e  )              (  l2  l3  e   0      )
!    (                  0  )              (  l2  l3  l4  e   0  )
!
!  where e denotes the off-diagonal elements of T, and ui (li)
!  denotes an element of the i-th column of U (L).
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
!  or in the even columns (if UPLU = 'U') are properly computed by SKTRF.
!
!  If UPLO = 'U', the U*pT*U^T decomposition of A is computed. U is a
!  upper triangular unit matrix with the additional constraint
!  U(1:i-1,i) = 0 for even i, and pT a partially tridiagonal matrix.
!  The entries in the odd rows of the upper diagonal of pT are stored
!  on exit in A(i,i+1) for i odd. The column U(1:i-1, i) for odd i
!  is stored in A(1:i-1,i+1).
!
!  If UPLO = 'L', the L*pT*L^T decomposition of A is computed. L is a
!  lower triangular unit matrix with the additional constraint
!  L(i+1:n,i) = 0 for odd i, and pT a partially tridiagonal matrix.
!  The entries in odd columns in the lower diagonal of pT are stored
!  on exit in A(i+1,i) for i odd. The column L(i+1:n, i) for i odd
!  is stored in A(i+1:n,i-1).
!
!  The contents of A on exit are illustrated by the following examples
!  with n = 6:
!
!  if UPLO = 'U':                       if UPLO = 'L':
!
!    (  0   e   x   u3  x   u5 )              (  0                    )
!    (      0   x   u3  x   u5 )              (  e   0                )
!    (          0   e   x   u5 )              (  l2  x   0            )
!    (              0   x   u5 )              (  l2  x   e   0        )
!    (                  0   e  )              (  l2  x   l4  x   0    )
!    (                      0  )              (  l2  x   l4  x   e  0 )
!
!  where e denotes the off-diagonal elements of T, ui (li)
!  denotes an element of the i-th column of U (L), and x denotes an
!  element not computed by SKTRF.
!
!  =====================================================================
!

SUBROUTINE SSKTRF_F95( A, IPIV, UPLO, MODE, INFO )
  USE PFAPACK_PREC, ONLY: precision => singleprec
  USE PFAPACK_MESSAGE, ONLY: MESSAGE
  USE F77_PFAPACK, ONLY: SKTRF
  IMPLICIT NONE
  REAL(precision), INTENT(INOUT) :: A(:,:)
  INTEGER, INTENT(OUT) :: IPIV(:)
  CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO, MODE
  INTEGER, INTENT(OUT), OPTIONAL :: INFO

  CHARACTER(LEN=1) :: LUPLO, LMODE
  INTEGER :: N, LDA, LWORK, LINFO, ISTAT
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

  IF( PRESENT(MODE) ) THEN
     LMODE = MODE
  ELSE
     LMODE = 'N'
  END IF

  ! Test the arguments
  IF( SIZE(A,2) /= N .OR. N<0 ) THEN
     LINFO = -1
  ELSE IF( SIZE(IPIV) /= N ) THEN
     LINFO = -2
  ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') ) THEN
     LINFO = -3
  ELSE IF( .NOT.LSAME(LMODE,'N') .AND. .NOT.LSAME(LMODE,'P') ) THEN
     LINFO = -4
  ELSE IF( N>0 ) THEN

     !First do a workspace query
     CALL SKTRF( LUPLO, LMODE, N, A, LDA, IPIV, QWORK, -1, LINFO)
     LWORK = INT(QWORK(1))
     ALLOCATE(WORK(LWORK), STAT=ISTAT)

     IF( ISTAT /= 0 ) THEN
        ! Allocating the workspace failed. Try minimal workspace
        DEALLOCATE(WORK)
        CALL MESSAGE(-200, 'SKTRF') ! print warning about not enough workspace

        LWORK=1
        ALLOCATE(WORK(LWORK), STAT=ISTAT)
     END IF

     IF( ISTAT == 0 ) THEN
        CALL SKTRF( LUPLO, LMODE, N, A, LDA, IPIV, WORK, LWORK, LINFO )
     ELSE
        LINFO=-100 ! not enough memory
     END IF

     DEALLOCATE(WORK)
  END IF

  CALL MESSAGE(LINFO, "SKTRF", INFO, ISTAT)
END SUBROUTINE SSKTRF_F95

SUBROUTINE DSKTRF_F95( A, IPIV, UPLO, MODE, INFO )
  USE PFAPACK_PREC, ONLY: precision => doubleprec
  USE PFAPACK_MESSAGE, ONLY: MESSAGE
  USE F77_PFAPACK, ONLY: SKTRF
  IMPLICIT NONE
  REAL(precision), INTENT(INOUT) :: A(:,:)
  INTEGER, INTENT(OUT) :: IPIV(:)
  CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO, MODE
  INTEGER, INTENT(OUT), OPTIONAL :: INFO

  CHARACTER(LEN=1) :: LUPLO, LMODE
  INTEGER :: N, LDA, LWORK, LINFO, ISTAT
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

  IF( PRESENT(MODE) ) THEN
     LMODE = MODE
  ELSE
     LMODE = 'N'
  END IF

  ! Test the arguments
  IF( SIZE(A,2) /= N .OR. N<0 ) THEN
     LINFO = -1
  ELSE IF( SIZE(IPIV) /= N ) THEN
     LINFO = -2
  ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') ) THEN
     LINFO = -3
  ELSE IF( .NOT.LSAME(LMODE,'N') .AND. .NOT.LSAME(LMODE,'P') ) THEN
     LINFO = -4
  ELSE IF( N>0 ) THEN

     !First do a workspace query
     CALL SKTRF( LUPLO, LMODE, N, A, LDA, IPIV, QWORK, -1, LINFO)
     LWORK = INT(QWORK(1))
     ALLOCATE(WORK(LWORK), STAT=ISTAT)

     IF( ISTAT /= 0 ) THEN
        ! Allocating the workspace failed. Try minimal workspace
        DEALLOCATE(WORK)
        CALL MESSAGE(-200, 'SKTRF') ! print warning about not enough workspace

        LWORK=1
        ALLOCATE(WORK(LWORK), STAT=ISTAT)
     END IF

     IF( ISTAT == 0 ) THEN
        CALL SKTRF( LUPLO, LMODE, N, A, LDA, IPIV, WORK, LWORK, LINFO )
     ELSE
        LINFO=-100 ! not enough memory
     END IF

     DEALLOCATE(WORK)
  END IF

  CALL MESSAGE(LINFO, "SKTRF", INFO, ISTAT)
END SUBROUTINE DSKTRF_F95

SUBROUTINE CSKTRF_F95( A, IPIV, UPLO, MODE, INFO )
  USE PFAPACK_PREC, ONLY: precision => singleprec
  USE PFAPACK_MESSAGE, ONLY: MESSAGE
  USE F77_PFAPACK, ONLY: SKTRF
  IMPLICIT NONE
  COMPLEX(precision), INTENT(INOUT) :: A(:,:)
  INTEGER, INTENT(OUT) :: IPIV(:)
  CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO, MODE
  INTEGER, INTENT(OUT), OPTIONAL :: INFO

  CHARACTER(LEN=1) :: LUPLO, LMODE
  INTEGER :: N, LDA, LWORK, LINFO, ISTAT
  COMPLEX(precision) :: QWORK(1)
  COMPLEX(precision), ALLOCATABLE :: WORK(:)

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
  ELSE IF( SIZE(IPIV) /= N ) THEN
     LINFO = -2
  ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') ) THEN
     LINFO = -3
  ELSE IF( .NOT.LSAME(LMODE,'N') .AND. .NOT.LSAME(LMODE,'P') ) THEN
     LINFO = -4
  ELSE IF( N>0 ) THEN

     !First do a workspace query
     CALL SKTRF( LUPLO, LMODE, N, A, LDA, IPIV, QWORK, -1, LINFO)
     LWORK = INT(QWORK(1))
     ALLOCATE(WORK(LWORK), STAT=ISTAT)

     IF( ISTAT /= 0 ) THEN
        ! Allocating the workspace failed. Try minimal workspace
        DEALLOCATE(WORK)
        CALL MESSAGE(-200, 'SKTRF') ! print warning about not enough workspace

        LWORK=1
        ALLOCATE(WORK(LWORK), STAT=ISTAT)
     END IF

     IF( ISTAT == 0 ) THEN
        CALL SKTRF( LUPLO, LMODE, N, A, LDA, IPIV, WORK, LWORK, LINFO )
     ELSE
        LINFO=-100 ! not enough memory
     END IF

     DEALLOCATE(WORK)
  END IF

  CALL MESSAGE(LINFO, "SKTRF", INFO, ISTAT)
END SUBROUTINE CSKTRF_F95

SUBROUTINE ZSKTRF_F95( A, IPIV, UPLO, MODE, INFO )
  USE PFAPACK_PREC, ONLY: precision => doubleprec
  USE PFAPACK_MESSAGE, ONLY: MESSAGE
  USE F77_PFAPACK, ONLY: SKTRF
  IMPLICIT NONE
  COMPLEX(precision), INTENT(INOUT) :: A(:,:)
  INTEGER, INTENT(OUT) :: IPIV(:)
  CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO, MODE
  INTEGER, INTENT(OUT), OPTIONAL :: INFO

  CHARACTER(LEN=1) :: LUPLO, LMODE
  INTEGER :: N, LDA, LWORK, LINFO, ISTAT
  COMPLEX(precision) :: QWORK(1)
  COMPLEX(precision), ALLOCATABLE :: WORK(:)

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
  ELSE IF( SIZE(IPIV) /= N ) THEN
     LINFO = -2
  ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') ) THEN
     LINFO = -3
  ELSE IF( .NOT.LSAME(LMODE,'N') .AND. .NOT.LSAME(LMODE,'P') ) THEN
     LINFO = -4
  ELSE IF( N>0 ) THEN

     !First do a workspace query
     CALL SKTRF( LUPLO, LMODE, N, A, LDA, IPIV, QWORK, -1, LINFO)
     LWORK = INT(QWORK(1))
     ALLOCATE(WORK(LWORK), STAT=ISTAT)

     IF( ISTAT /= 0 ) THEN
        ! Allocating the workspace failed. Try minimal workspace
        DEALLOCATE(WORK)
        CALL MESSAGE(-200, 'SKTRF') ! print warning about not enough workspace

        LWORK=1
        ALLOCATE(WORK(LWORK), STAT=ISTAT)
     END IF

     IF( ISTAT == 0 ) THEN
        CALL SKTRF( LUPLO, LMODE, N, A, LDA, IPIV, WORK, LWORK, LINFO )
     ELSE
        LINFO=-100 ! not enough memory
     END IF

     DEALLOCATE(WORK)
  END IF

  CALL MESSAGE(LINFO, "SKTRF", INFO, ISTAT)
END SUBROUTINE ZSKTRF_F95
