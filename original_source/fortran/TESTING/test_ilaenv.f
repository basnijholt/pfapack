      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
*     Replacement for ILAENV to check different blocking parameters
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4

      INTEGER            IPARMS( 10 )
      COMMON             / LAENV / IPARMS
      SAVE               / LAENV /

      ILAENV = IPARMS( ISPEC )

      RETURN

      END

      SUBROUTINE SLAENV( ISPEC, VALUE )
      INTEGER ISPEC, VALUE

      INTEGER            IPARMS( 10 )
      COMMON             / LAENV / IPARMS
      SAVE               / LAENV /

      IPARMS( ISPEC ) = VALUE

      END
