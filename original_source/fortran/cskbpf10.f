      SUBROUTINE CSKBPF10( UPLO, N, KD, AB, LDAB, PFAFF, WORK, RWORK,
     $                     INFO )
*
*  -- Written on 10/25/2010
*     Michael Wimmer, Universiteit Leiden

*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, N
*     ..
*     .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            PFAFF( 2 )
      COMPLEX            AB( LDAB, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CSKBPF10 computes the Pfaffian of a complex skew-symmetric band
*  matrix, taking special care to avoid numerical under- or overflow.
*  (at the cost of possible additional round-off errors).
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of superdiagonals of the matrix A if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
*
*  AB      (input/output) COMPLEX array, dimension (LDAB,N)
*          If N is odd, AB is not referenced.
*          If N is even:
*          On entry, the upper or lower triangle of the skew-symmetric
*          band matrix A, stored in the first KD+1 rows of the array.
*          The j-th column of A is stored in the j-th column of the
*          array AB as follows:
*            if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*            if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=min(n,j+kd).
*          On exit, AB is overwritten with values generated during the
*          computation.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD+1.
*
*  PFAFF   (output) COMPLEX array, dimension 2
*          The value of the Pfaffian in the form
*          PFAFF(1)*10**PFAFF(2) with PFAFF(2) purely real.
*
*  WORK    (workspace) COMPLEX array, dimension (N)
*
*  RWORK   (workspace) REAL array, dimension (2*N-1)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The Pfaffian is computed by bringing the skew-symmetric matrix A into
*  a partial tridiagonal form pT by a unitary congruence transformation:
*  Q^H * A * Q^* = pT.
*  This transformation is computed by the routine CSKBTRD (for further
*  details see there)
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ),
     $                     ONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I
      COMPLEX            DETQ
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, CSKBTRD, CMUL10
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      UPPER = LSAME( UPLO, 'U' )
*
      INFO = 0
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CSKBPF10', -INFO )
         RETURN
      END IF

      PFAFF( 1 ) = ONE
      PFAFF( 2 ) = ZERO

*     Quick return if possible
      IF( N.EQ.0 ) THEN
         RETURN
      ELSE IF( MOD(N,2).EQ.1 ) THEN
         PFAFF( 1 ) = ZERO
         RETURN
      END IF
*     Reduce to tridiagonal form

      CALL CSKBTRD("N", UPLO, "P", N, KD, AB, LDAB, RWORK(1),
     $             DETQ, WORK, 1, WORK, RWORK(N),
     $             INFO)

      PFAFF( 1 ) = DETQ

      IF( UPPER ) THEN
*        Multiply every other entry on the superdiagonal
         DO 10 I = 1, N-1, 2
            CALL CMUL10( PFAFF, CMPLX(RWORK( I )) )
 10      CONTINUE
      ELSE
*        Multiply every other entry on the superdiagonal
         DO 20 I = 1, N-1, 2
            CALL CMUL10( PFAFF, -CMPLX(RWORK( I )) )
 20      CONTINUE
      END IF

      RETURN

*     end of CSKBPF10

      END
