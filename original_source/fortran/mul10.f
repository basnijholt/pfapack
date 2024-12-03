*     Routines to multiply numbers and avoiding overflow
      SUBROUTINE SMUL10(A, B)
      REAL A(2)
      REAL B
      REAL EXPONENT, SLAMCH
      INTEGER IEXPONENT
      EXTERNAL SLAMCH

      A(1) = A(1) * B

      IF( A(1).EQ.0.0E0 ) THEN
         A(1) = 0.0E0
         A(2) = 0.0E0
      ELSE
         EXPONENT = LOG10(ABS(A(1)))

         IF (EXPONENT .GE. 0E0) THEN
            IEXPONENT=INT(EXPONENT)
         ELSE
            IEXPONENT=INT(EXPONENT)-1
         END IF

*     If 10**IEXPONENT is smaller than the safe minimum for inversion,
*     the number is considered 0
         IF( SLAMCH( "S" ) > 10E0**IEXPONENT ) THEN
            A(1) = 0.0E0
            A(2) = 0.0E0
         ELSE
            A(1) = A(1)/10E0**IEXPONENT
            A(2) = A(2) + IEXPONENT
         END IF
      END IF

      END

      SUBROUTINE DMUL10(A, B)
      DOUBLE PRECISION A(2)
      DOUBLE PRECISION B
      DOUBLE PRECISION EXPONENT, DLAMCH
      INTEGER IEXPONENT
      EXTERNAL DLAMCH

      A(1) = A(1) * B

      IF( A(1).EQ.0.0D0 ) THEN
         A(1) = 0.0D0
         A(2) = 0.0D0
      ELSE
         EXPONENT = LOG10(ABS(A(1)))

         IF (EXPONENT .GE. 0D0) THEN
            IEXPONENT=INT(EXPONENT)
         ELSE
            IEXPONENT=INT(EXPONENT)-1
         END IF

*     If 10**IEXPONENT is smaller than the safe minimum for inversion,
*     the number is considered 0
         IF( DLAMCH( "S" ) > 10D0**IEXPONENT ) THEN
            A(1) = 0.0D0
            A(2) = 0.0D0
         ELSE
            A(1) = A(1)/10D0**IEXPONENT
            A(2) = A(2) + IEXPONENT
         END IF
      END IF

      END

      SUBROUTINE CMUL10(A, B)
      COMPLEX A(2)
      COMPLEX B
      REAL EXPONENT, SLAMCH
      INTEGER IEXPONENT
      EXTERNAL SLAMCH
      A(1) = A(1) * B

      IF( A(1).EQ.(0.0E0, 0.0E0) ) THEN
         A(1) = (0.0E0, 0.0E0)
         A(2) = (0.0E0, 0.0E0)
      ELSE
         EXPONENT = LOG10(ABS(A(1)))

         IF (EXPONENT .GE. 0E0) THEN
            IEXPONENT=INT(EXPONENT)
         ELSE
            IEXPONENT=INT(EXPONENT)-1
         END IF

*     If 10**IEXPONENT is smaller than the safe minimum for inversion,
*     the number is considered 0
         IF( SLAMCH( "S" ) > 10E0**IEXPONENT ) THEN
            A(1) = (0.0E0, 0.0E0)
            A(2) = (0.0E0, 0.0E0)
         ELSE
            A(1) = A(1)/10E0**IEXPONENT
            A(2) = A(2) + IEXPONENT
         END IF
      END IF

      END

      SUBROUTINE ZMUL10(A, B)
      DOUBLE COMPLEX A(2)
      DOUBLE COMPLEX B
      DOUBLE PRECISION EXPONENT, DLAMCH
      INTEGER IEXPONENT
      EXTERNAL DLAMCH

      A(1) = A(1) * B

      IF( A(1).EQ.(0.0D0, 0.0D0) ) THEN
         A(1) = (0.0D0, 0.0D0)
         A(2) = (0.0D0, 0.0D0)
      ELSE
         EXPONENT = LOG10(ABS(A(1)))

         IF (EXPONENT .GE. 0D0) THEN
            IEXPONENT=INT(EXPONENT)
         ELSE
            IEXPONENT=INT(EXPONENT)-1
         END IF

*     If 10**IEXPONENT is smaller than the safe minimum for inversion,
*     the number is considered 0
         IF( DLAMCH( "S" ) > 10D0**IEXPONENT ) THEN
            A(1) = (0.0D0, 0.0D0)
            A(2) = (0.0D0, 0.0D0)
         ELSE
            A(1) = A(1)/10D0**IEXPONENT
            A(2) = A(2) + IEXPONENT
         END IF
      END IF

      END
