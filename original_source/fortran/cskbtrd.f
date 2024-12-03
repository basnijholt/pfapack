      SUBROUTINE CSKBTRD( VECT, UPLO, MODE, N, KD, AB, LDAB, E, DETQ,
     $                    Q, LDQ, WORK, RWORK, INFO )
*
*  -- Written on 10/25/2010
*     Michael Wimmer, Universiteit Leiden
*
*  -- Derived from the LAPACK routine ZHBTRD (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, VECT, MODE
      INTEGER            INFO, KD, LDAB, LDQ, N
      COMPLEX            DETQ
*     ..
*     .. Array Arguments ..
      REAL               E( * ), RWORK( * )
      COMPLEX            AB( LDAB, * ), Q( LDQ, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CSKBTRD reduces a complex skew-symmetric band matrix A to real
*  skew-symmetric tridiagonal form T by a unitary congruence
*  transformation: Q^dagger * A * Q^* = T. Alternatively, the routine can
*  also compute a partial tridiagonal form useful for computing the Pfaffian.
*
*  Arguments
*  =========
*
*  VECT    (input) CHARACTER*1
*          = 'N':  do not form Q;
*          = 'V':  form Q;
*          = 'U':  update a matrix X, by forming X*Q.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  MODE    (input) CHARACTER*1
*          = 'N':  A is fully tridiagonalized
*          = 'P':  A is partially tridiagonalized for Pfaffian computation
*                  (details see below)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0. N must be even if MODE = 'P'.
*
*  KD      (input) INTEGER
*          The number of superdiagonals of the matrix A if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
*
*  AB      (input/output) COMPLEX array, dimension (LDAB,N)
*          On entry, the upper or lower triangle of the skew-symmetric
*          band matrix A, stored in the first KD+1 rows of the array.
*          The j-th column of A is stored in the j-th column of the
*          array AB as follows:
*            if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*            if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=min(n,j+kd).
*          On exit, the zero diagonal elements of AB are left unchanged,
*          if KD > 0, the elements on the first superdiagonal (if UPLO =
*          'U') or the first subdiagonal (if UPLO = 'L') are overwritten
*          by the off-diagonal elements of T; the rest of AB is
*          overwritten by values generated during the reduction. If
*          MODE = 'P', only the off-diagonal entries in the odd rows
*          (columns) are computed for UPLO = 'U' (UPLO = 'L').
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD+1.
*
*
*  E       (output) REAL array, dimension (N-1)
*          The off-diagonal elements of the tridiagonal matrix T:
*          E(i) = T(i,i+1) if UPLO = 'U'; E(i) = T(i+1,i) if UPLO = 'L'.
*
*  DETQ    (output) COMPLEX
*          The value of the determinant of Q, which is a
*          pure phase factor. Always computed, even if Q is
*          not explicitely formed.
*
*  Q       (input/output) COMPLEX array, dimension (LDQ,N)
*          On entry, if VECT = 'U', then Q must contain an N-by-N
*          matrix X; if VECT = 'N' or 'V', then Q need not be set.
*
*          On exit:
*          if VECT = 'V', Q contains the N-by-N unitary matrix Q;
*          if VECT = 'U', Q contains the product X*Q;
*          if VECT = 'N', the array Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.
*          LDQ >= 1, and LDQ >= N if VECT = 'V' or 'U'.
*
*  WORK    (workspace) COMPLEX array, dimension (N)
*
*  RWORK   (workspace) REAL array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The storage scheme for the skew-symmetric matrix is identical to the
*  storage scheme for symmetric or Hermitian band matrices in LAPACK,
*  i.e. the diagonal and the KD super- or subdiagonals are stored in an
*  array with KD+1 rows and N columns. Note that the zero diagonal must
*  also be explicitely stored (this was done to keep the structure of
*  the program identical to the symmetric case)
*
*  In particular this means that if
*  - UPLO = 'U', then AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j
*
*    Example: N=5, KD=2
*
*  (  0     a12   a13              )
*  (  -a12  0     a23   a24        )
*  (  -a13  -a23  0     a34   a35  )
*  (        -a24  -a34  0     a45  )
*  (              -a35  -a45  0    )
*
*  is stored as
*
*   x   x   a13 a24 a35
*   x   a12 a23 a34 a45
*   0   0   0   0   0
*
*  where x denotes an unused entry
*
*  - UPLO = 'L', then AB(1+i-j,j) = A(i,j) for j<=i<=min(n,j+kd)
*
*    Example: N=5, KD=2
*
*  (  0     -a21  -a31              )
*  (  a21   0     -a32  -a42        )
*  (  a31   a32   0     -a43  -a53  )
*  (        a42   a43   0     -a54  )
*  (              a53   a54   0     )
*
*  is stored as
*
*   0   0   0   0   0
*   a21 a32 a43 a54 x
*   a31 a42 a53 x   x
*
*  where x denotes an unused entry
*
*  =====================================================================
*
*     .. Parameters ..
      REAL   ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ),
     $                   CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            INITQ, UPPER, WANTQ, NORMAL
      INTEGER            I, I2, IBL, INCA, INCX, IQAEND, IQB, IQEND, J,
     $                   J1, J1END, J1INC, J2, JEND, JIN, JINC, K, KD1,
     $                   KDM1, KDN, L, LAST, LEND, NQ, NR, NRT, STEP
      REAL               ABST
      COMPLEX            T, TEMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, CLARGV, CLARTG, CLARTV,
     $                   CLASET, CROT, CSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CONJG, MAX, MIN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INITQ = LSAME( VECT, 'V' )
      WANTQ = INITQ .OR. LSAME( VECT, 'U' )
      UPPER = LSAME( UPLO, 'U' )
      NORMAL = LSAME( MODE, 'N' )
      KD1 = KD + 1
      KDM1 = KD - 1
      INCX = LDAB - 1
      IQEND = 1
*
      INFO = 0
      IF( .NOT.WANTQ .AND. .NOT.LSAME( VECT, 'N' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NORMAL .AND. .NOT.LSAME( MODE, 'P' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( .NOT.NORMAL .AND. MOD(N,2).EQ.1 ) THEN
*     If STEP == 2, we need an even-dimensional matrix
         INFO = -4
      ELSE IF( KD.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.KD1 ) THEN
         INFO = -7
      ELSE IF( LDQ.LT.MAX( 1, N ) .AND. WANTQ ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CSKBTRD', -INFO )
         RETURN
      END IF

      IF( NORMAL ) THEN
         STEP = 1
      ELSE
         STEP = 2
      END IF

*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Initialize Q to the unit matrix, if needed
*
      IF( INITQ )
     $   CALL CLASET( 'Full', N, N, CZERO, CONE, Q, LDQ )
*
*     Wherever possible, plane rotations are generated and applied in
*     vector operations of length NR over the index set J1:J2:KD1.
*
*     The real cosines and complex sines of the plane rotations are
*     stored in the arrays RWORK and WORK.
*
      INCA = KD1*LDAB
      KDN = MIN( N-1, KD )
      IF( UPPER ) THEN
*
         IF( KD.GT.1 ) THEN
*
*           Reduce to real skew-symmetric tridiagonal form, working with
*           the upper triangle
*
            NR = 0
            J1 = KDN + 2
            J2 = 1
*
            DO 90 I = 1, N - 2, STEP
*
*              Reduce i-th row of matrix to tridiagonal form
*
               DO 80 K = KDN + 1, 2, -1

                  IF( STEP.EQ.2 .AND. K.EQ.2 ) THEN

*     Skip the entry that was generated in the even row I+1
                     J1 = J1 + KDN + 1
                     NR = NR - 1

*     Skip the loop with K = 2
                     GOTO 80
                  END IF

                  J1 = J1 + KDN
                  J2 = J2 + KDN
*
                  IF( NR.GT.0 ) THEN
*
*                    generate plane rotations to annihilate nonzero
*                    elements which have been created outside the band
*
                     CALL CLARGV( NR, AB( 1, J1-1 ), INCA, WORK( J1 ),
     $                            KD1, RWORK( J1 ), KD1 )
*
*                    apply rotations from the right
*
*
*                    Dependent on the the number of diagonals either
*                    CLARTV or CROT is used
*
                     IF( NR.GE.2*KD-1 ) THEN
                        DO 10 L = 1, KD - 1
                           CALL CLARTV( NR, AB( L+1, J1-1 ), INCA,
     $                                  AB( L, J1 ), INCA, RWORK( J1 ),
     $                                  WORK( J1 ), KD1 )
   10                   CONTINUE
*
                     ELSE
                        JEND = J1 + ( NR-1 )*KD1
                        DO 20 JINC = J1, JEND, KD1
                           CALL CROT( KDM1, AB( 2, JINC-1 ), 1,
     $                                AB( 1, JINC ), 1, RWORK( JINC ),
     $                                WORK( JINC ) )
   20                   CONTINUE
                     END IF
                  END IF
*
*
                  IF( K.GT.2 ) THEN
                     IF( K.LE.N-I+1 ) THEN
*
*                       generate plane rotation to annihilate a(i,i+k-1)
*                       within the band
*
                        CALL CLARTG( AB( KD-K+3, I+K-2 ),
     $                               AB( KD-K+2, I+K-1 ), RWORK(I+K-1),
     $                               WORK( I+K-1 ), TEMP )
                        AB( KD-K+3, I+K-2 ) = TEMP
*
*                       apply rotation from the right
*
                        CALL CROT( K-3, AB( KD-K+4, I+K-2 ), 1,
     $                             AB( KD-K+3, I+K-1 ), 1, RWORK(I+K-1),
     $                             WORK( I+K-1 ) )
                     END IF
                     NR = NR + 1
                     J1 = J1 - KDN - 1
                  END IF
*
*                 apply plane rotations from both sides to diagonal
*                 blocks
*
*                 Diagonal blocks are invariant in the skew-symmetric case
*
*                 apply plane rotations from the left
*
                  IF( NR.GT.0 ) THEN
*                     No complex conjugation, as we have Q^T.
*                     But needed below for accumulating Q
                     IF( 2*KD-1.LT.NR ) THEN
*
*                    Dependent on the the number of diagonals either
*                    CLARTV or CROT is used
*
                        DO 30 L = 1, KD - 1
                           IF( J2+L.GT.N ) THEN
                              NRT = NR - 1
                           ELSE
                              NRT = NR
                           END IF
                           IF( NRT.GT.0 )
     $                        CALL CLARTV( NRT, AB( KD-L, J1+L ), INCA,
     $                                     AB( KD-L+1, J1+L ), INCA,
     $                                     RWORK( J1 ), WORK( J1 ), KD1)
   30                   CONTINUE
                     ELSE
                        J1END = J1 + KD1*( NR-2 )
                        IF( J1END.GE.J1 ) THEN
                           DO 40 JIN = J1, J1END, KD1
                              CALL CROT( KD-1, AB( KD-1, JIN+1 ), INCX,
     $                                   AB( KD, JIN+1 ), INCX,
     $                                   RWORK( JIN ), WORK( JIN ) )
   40                      CONTINUE
                        END IF
                        LEND = MIN( KDM1, N-J2 )
                        LAST = J1END + KD1
                        IF( LEND.GT.0 )
     $                     CALL CROT( LEND, AB( KD-1, LAST+1 ), INCX,
     $                                AB( KD, LAST+1 ), INCX,
     $                                RWORK( LAST ), WORK( LAST ) )
                     END IF
                  END IF
*
                  IF( WANTQ ) THEN
*
*                    accumulate product of plane rotations in Q
*
                     IF( INITQ ) THEN
*
*                 take advantage of the fact that Q was
*                 initially the Identity matrix
*
                        IQEND = MAX( IQEND, J2 )
                        I2 = MAX( 0, K-3 )
                        IQAEND = 1 + I*KD
                        IF( K.EQ.2 )
     $                     IQAEND = IQAEND + KD
                        IQAEND = MIN( IQAEND, IQEND )
                        DO 50 J = J1, J2, KD1
                           IBL = I - I2 / KDM1
                           I2 = I2 + 1
                           IQB = MAX( 1, J-IBL )
                           NQ = 1 + IQAEND - IQB
                           IQAEND = MIN( IQAEND+KD, IQEND )
                           CALL CROT( NQ, Q( IQB, J-1 ), 1, Q( IQB, J ),
     $                                1, RWORK( J ), CONJG(WORK( J )) )
   50                   CONTINUE
                     ELSE
*
                        DO 60 J = J1, J2, KD1
                           CALL CROT( N, Q( 1, J-1 ), 1, Q( 1, J ), 1,
     $                                RWORK( J ), CONJG(WORK( J ))  )
   60                   CONTINUE
                     END IF
*
                  END IF
*
                  IF( J2+KDN.GT.N ) THEN
*
*                    adjust J2 to keep within the bounds of the matrix
*
                     NR = NR - 1
                     J2 = J2 - KDN - 1
                  END IF
*
                  DO 70 J = J1, J2, KD1
*
*                    create nonzero element a(j-1,j+kd) outside the band
*                    and store it in WORK
*
                     WORK( J+KD ) = WORK( J )*AB( 1, J+KD )
                     AB( 1, J+KD ) = RWORK( J )*AB( 1, J+KD )
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
         END IF
*
         DETQ = CONE

         IF( KD.GT.0 ) THEN
*
*           make off-diagonal elements real and copy them to E
*
            DO 100 I = 1, N - 1
               T = AB( KD, I+1 )
               ABST = ABS( T )
               AB( KD, I+1 ) = ABST
               E( I ) = ABST
               IF( ABST.NE.ZERO ) THEN
                  T = T / ABST
               ELSE
                  T = CONE
               END IF
               IF( I.LT.N-1 )
     $            AB( KD, I+2 ) = AB( KD, I+2 )*CONJG(T)
               IF( WANTQ ) THEN
                  CALL CSCAL( N, T, Q( 1, I+1 ), 1 )
               END IF
               DETQ = DETQ * T
  100       CONTINUE
         ELSE
*
*           set E to zero if original matrix was diagonal
*
            DO 110 I = 1, N - 1
               E( I ) = ZERO
  110       CONTINUE
         END IF
      ELSE
*
         IF( KD.GT.1 ) THEN
*
*           Reduce to real skew-symmetric tridiagonal form, working with
*           the lower triangle
*
            NR = 0
            J1 = KDN + 2
            J2 = 1
*
            DO 210 I = 1, N - 2, STEP
*
*              Reduce i-th column of matrix to tridiagonal form
*
               DO 200 K = KDN + 1, 2, -1

                  IF( STEP.EQ.2 .AND. K.EQ.2 ) THEN

*     Skip the entry that was generated in the even row I+1
                     J1 = J1 + KDN + 1
                     NR = NR - 1

*     Skip the loop with K = 2
                     GOTO 200
                  END IF

                  J1 = J1 + KDN
                  J2 = J2 + KDN
*
                  IF( NR.GT.0 ) THEN
*
*                    generate plane rotations to annihilate nonzero
*                    elements which have been created outside the band
*
                     CALL CLARGV( NR, AB( KD1, J1-KD1 ), INCA,
     $                            WORK( J1 ), KD1, RWORK( J1 ), KD1 )
*
*                    apply plane rotations from one side
*
*
*                    Dependent on the the number of diagonals either
*                    CLARTV or CROT is used
*
                     IF( NR.GT.2*KD-1 ) THEN
                        DO 130 L = 1, KD - 1
                           CALL CLARTV( NR, AB( KD1-L, J1-KD1+L ), INCA,
     $                                  AB( KD1-L+1, J1-KD1+L ), INCA,
     $                                  RWORK( J1 ), WORK( J1 ), KD1 )
  130                   CONTINUE
                     ELSE
                        JEND = J1 + KD1*( NR-1 )
                        DO 140 JINC = J1, JEND, KD1
                           CALL CROT( KDM1, AB( KD, JINC-KD ), INCX,
     $                                AB( KD1, JINC-KD ), INCX,
     $                                RWORK( JINC ), WORK( JINC ) )
  140                   CONTINUE
                     END IF
*
                  END IF
*
                  IF( K.GT.2 ) THEN
                     IF( K.LE.N-I+1 ) THEN
*
*                       generate plane rotation to annihilate a(i+k-1,i)
*                       within the band
*
                        CALL CLARTG( AB( K-1, I ), AB( K, I ),
     $                               RWORK( I+K-1 ), WORK( I+K-1 ),TEMP)
                        AB( K-1, I ) = TEMP
*
*                       apply rotation from the left
*
                        CALL CROT( K-3, AB( K-2, I+1 ), LDAB-1,
     $                             AB( K-1, I+1 ), LDAB-1, RWORK(I+K-1),
     $                             WORK( I+K-1 ) )
                     END IF
                     NR = NR + 1
                     J1 = J1 - KDN - 1
                  END IF
*
*                 apply plane rotations from both sides to diagonal
*                 blocks
*
*                 Not necessary in skew-symmetric case
*
*                 apply plane rotations from the right
*
*
*                    Dependent on the the number of diagonals either
*                    CLARTV or CROT is used
*
                  IF( NR.GT.0 ) THEN
*                     Not needed here, we have Q^T,
*                     but we need it further below when accumulating Q
                     IF( NR.GT.2*KD-1 ) THEN
                        DO 150 L = 1, KD - 1
                           IF( J2+L.GT.N ) THEN
                              NRT = NR - 1
                           ELSE
                              NRT = NR
                           END IF
                           IF( NRT.GT.0 )
     $                        CALL CLARTV( NRT, AB( L+2, J1-1 ), INCA,
     $                                     AB( L+1, J1 ), INCA,
     $                                     RWORK( J1 ), WORK( J1 ), KD1)
  150                   CONTINUE
                     ELSE
                        J1END = J1 + KD1*( NR-2 )
                        IF( J1END.GE.J1 ) THEN
                           DO 160 J1INC = J1, J1END, KD1
                              CALL CROT( KDM1, AB( 3, J1INC-1 ), 1,
     $                                   AB( 2, J1INC ), 1,
     $                                   RWORK( J1INC ), WORK( J1INC ) )
  160                      CONTINUE
                        END IF
                        LEND = MIN( KDM1, N-J2 )
                        LAST = J1END + KD1
                        IF( LEND.GT.0 )
     $                     CALL CROT( LEND, AB( 3, LAST-1 ), 1,
     $                                AB( 2, LAST ), 1, RWORK( LAST ),
     $                                WORK( LAST ) )
                     END IF
                  END IF
*
*
*
                  IF( WANTQ ) THEN
*
*                    accumulate product of plane rotations in Q
*
*                    Note that up to now we multiplied with Q^H from the left
*                    Hence, to accumulate Q, we have to multiply with the
*                    hermitian conjugate of the Givens rotations

                     IF( INITQ ) THEN
*
*                 take advantage of the fact that Q was
*                 initially the Identity matrix
*
                        IQEND = MAX( IQEND, J2 )
                        I2 = MAX( 0, K-3 )
                        IQAEND = 1 + I*KD
                        IF( K.EQ.2 )
     $                     IQAEND = IQAEND + KD
                        IQAEND = MIN( IQAEND, IQEND )
                        DO 170 J = J1, J2, KD1
                           IBL = I - I2 / KDM1
                           I2 = I2 + 1
                           IQB = MAX( 1, J-IBL )
                           NQ = 1 + IQAEND - IQB
                           IQAEND = MIN( IQAEND+KD, IQEND )
                           CALL CROT( NQ, Q( IQB, J-1 ), 1, Q( IQB, J ),
     $                                1, RWORK( J ), CONJG(WORK( J )) )
  170                   CONTINUE
                     ELSE
*
                        DO 180 J = J1, J2, KD1
                           CALL CROT( N, Q( 1, J-1 ), 1, Q( 1, J ), 1,
     $                                RWORK( J ), CONJG(WORK( J )) )
  180                   CONTINUE
                     END IF
                  END IF
*
                  IF( J2+KDN.GT.N ) THEN
*
*                    adjust J2 to keep within the bounds of the matrix
*
                     NR = NR - 1
                     J2 = J2 - KDN - 1
                  END IF
*
                  DO 190 J = J1, J2, KD1
*
*                    create nonzero element a(j+kd,j-1) outside the
*                    band and store it in WORK
*
                     WORK( J+KD ) = WORK( J )*AB( KD1, J )
                     AB( KD1, J ) = RWORK( J )*AB( KD1, J )
  190             CONTINUE
  200          CONTINUE
  210       CONTINUE
         END IF
*
         DETQ = CONE

         IF( KD.GT.0 ) THEN
*
*           make off-diagonal elements real and copy them to E
*
            DO 220 I = 1, N - 1
               T = AB( 2, I )
               ABST = ABS( T )
               AB( 2, I ) = ABST
               E( I ) = ABST
               IF( ABST.NE.ZERO ) THEN
                  T = T / ABST
               ELSE
                  T = CONE
               END IF
               IF( I.LT.N-1 )
* Note: difference here because of Q^T instead of Q^H (*T^* instead of *T)
     $            AB( 2, I+1 ) = AB( 2, I+1 )*CONJG(T)
               IF( WANTQ ) THEN
                  CALL CSCAL( N, T, Q( 1, I+1 ), 1 )
               END IF
               DETQ = DETQ * T
 220       CONTINUE
         ELSE
*
*           set E to zero if original matrix was diagonal
*
            DO 230 I = 1, N - 1
               E( I ) = ZERO
  230       CONTINUE
         END IF

      END IF
*
      RETURN
*
*     End of CSKBTRD
*
      END
