  MODULE F77_PFAPACK

    !Interfaces for the Pfaffian routines

    INTERFACE SKPFA
       SUBROUTINE SSKPFA( UPLO, MTHD, N, A, LDA, PFAFF, &
                          IWORK, WORK, LWORK, INFO)
         USE PFAPACK_PREC, ONLY: singleprec
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, MTHD
         INTEGER, INTENT(IN) :: LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK( * )
         REAL(singleprec), INTENT(OUT) :: PFAFF
         REAL(singleprec), INTENT(INOUT) :: A( LDA, * )
         REAL(singleprec), INTENT(OUT) :: WORK( * )
       END SUBROUTINE SSKPFA

       SUBROUTINE DSKPFA( UPLO, MTHD, N, A, LDA, PFAFF, &
                          IWORK, WORK, LWORK, INFO)
         USE PFAPACK_PREC, ONLY: doubleprec
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, MTHD
         INTEGER, INTENT(IN) :: LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK( * )
         REAL(doubleprec), INTENT(OUT) :: PFAFF
         REAL(doubleprec), INTENT(INOUT) :: A( LDA, * )
         REAL(doubleprec), INTENT(OUT) :: WORK( * )
       END SUBROUTINE DSKPFA

       SUBROUTINE CSKPFA( UPLO, MTHD, N, A, LDA, PFAFF, &
                          IWORK, WORK, LWORK, RWORK, INFO)
         USE PFAPACK_PREC, ONLY: singleprec
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, MTHD
         INTEGER, INTENT(IN) :: LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK( * )
         COMPLEX(singleprec), INTENT(OUT) :: PFAFF
         REAL(singleprec), INTENT(OUT) :: RWORK( * )
         COMPLEX(singleprec), INTENT(INOUT) :: A( LDA, * )
         COMPLEX(singleprec), INTENT(OUT) :: WORK( * )
       END SUBROUTINE CSKPFA

       SUBROUTINE ZSKPFA( UPLO, MTHD, N, A, LDA, PFAFF, &
                          IWORK, WORK, LWORK, RWORK, INFO)
         USE PFAPACK_PREC, ONLY: doubleprec
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, MTHD
         INTEGER, INTENT(IN) :: LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK( * )
         COMPLEX(doubleprec), INTENT(OUT) :: PFAFF
         REAL(doubleprec), INTENT(OUT) :: RWORK( * )
         COMPLEX(doubleprec), INTENT(INOUT) :: A( LDA, * )
         COMPLEX(doubleprec), INTENT(OUT) :: WORK( * )
       END SUBROUTINE ZSKPFA
    END INTERFACE

    INTERFACE SKPF10
       SUBROUTINE SSKPF10( UPLO, MTHD, N, A, LDA, PFAFF, &
                           IWORK, WORK, LWORK, INFO)
         USE PFAPACK_PREC, ONLY: singleprec
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, MTHD
         INTEGER, INTENT(IN) :: LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK( * )
         REAL(singleprec), INTENT(OUT) :: PFAFF( 2 )
         REAL(singleprec), INTENT(INOUT) :: A( LDA, * )
         REAL(singleprec), INTENT(OUT) :: WORK( * )
       END SUBROUTINE SSKPF10

       SUBROUTINE DSKPF10( UPLO, MTHD, N, A, LDA, PFAFF, &
                           IWORK, WORK, LWORK, INFO)
         USE PFAPACK_PREC, ONLY: doubleprec
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, MTHD
         INTEGER, INTENT(IN) :: LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK( * )
         REAL(doubleprec), INTENT(OUT) :: PFAFF( 2 )
         REAL(doubleprec), INTENT(INOUT) :: A( LDA, * )
         REAL(doubleprec), INTENT(OUT) :: WORK( * )
       END SUBROUTINE DSKPF10

       SUBROUTINE CSKPF10( UPLO, MTHD, N, A, LDA, PFAFF, &
                           IWORK, WORK, LWORK, RWORK, INFO)
         USE PFAPACK_PREC, ONLY: singleprec
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, MTHD
         INTEGER, INTENT(IN) :: LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK( * )
         COMPLEX(singleprec), INTENT(OUT) :: PFAFF( 2 )
         REAL(singleprec), INTENT(OUT) :: RWORK( * )
         COMPLEX(singleprec), INTENT(INOUT) :: A( LDA, * )
         COMPLEX(singleprec), INTENT(OUT) :: WORK( * )
       END SUBROUTINE CSKPF10

       SUBROUTINE ZSKPF10( UPLO, MTHD, N, A, LDA, PFAFF, &
                           IWORK, WORK, LWORK, RWORK, INFO)
         USE PFAPACK_PREC, ONLY: doubleprec
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, MTHD
         INTEGER, INTENT(IN) :: LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK( * )
         COMPLEX(doubleprec), INTENT(OUT) :: PFAFF( 2 )
         REAL(doubleprec), INTENT(OUT) :: RWORK( * )
         COMPLEX(doubleprec), INTENT(INOUT) :: A( LDA, * )
         COMPLEX(doubleprec), INTENT(OUT) :: WORK( * )
       END SUBROUTINE ZSKPF10
    END INTERFACE

    INTERFACE SKBPFA
      SUBROUTINE SSKBPFA( UPLO, N, KD, AB, LDAB, PFAFF, WORK, INFO )
        USE PFAPACK_PREC, ONLY: singleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO
        INTEGER, INTENT(IN) :: KD, LDAB, N
        INTEGER, INTENT(OUT) :: INFO
        REAL(singleprec), INTENT(OUT) :: PFAFF
        REAL(singleprec), INTENT(INOUT) :: AB( LDAB, * )
        REAL(singleprec), INTENT(OUT) :: WORK( * )
      END SUBROUTINE SSKBPFA

      SUBROUTINE DSKBPFA( UPLO, N, KD, AB, LDAB, PFAFF, WORK, INFO )
         USE PFAPACK_PREC, ONLY: doubleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO
        INTEGER, INTENT(IN) :: KD, LDAB, N
        INTEGER, INTENT(OUT) :: INFO
        REAL(doubleprec), INTENT(OUT) :: PFAFF
        REAL(doubleprec), INTENT(INOUT) :: AB( LDAB, * )
        REAL(doubleprec), INTENT(OUT) :: WORK( * )
      END SUBROUTINE DSKBPFA

      SUBROUTINE CSKBPFA( UPLO, N, KD, AB, LDAB, PFAFF, WORK, &
                          RWORK, INFO )
        USE PFAPACK_PREC, ONLY: singleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO
        INTEGER, INTENT(IN) :: KD, LDAB, N
        INTEGER, INTENT(OUT) :: INFO
        COMPLEX(singleprec), INTENT(OUT) :: PFAFF
        REAL(singleprec), INTENT(OUT) :: RWORK( * )
        COMPLEX(singleprec), INTENT(INOUT) :: AB( LDAB, * )
        COMPLEX(singleprec), INTENT(OUT) :: WORK( * )
      END SUBROUTINE CSKBPFA

      SUBROUTINE ZSKBPFA( UPLO, N, KD, AB, LDAB, PFAFF, WORK, &
                          RWORK, INFO )
        USE PFAPACK_PREC, ONLY: doubleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO
        INTEGER, INTENT(IN) :: KD, LDAB, N
        INTEGER, INTENT(OUT) :: INFO
        COMPLEX(doubleprec), INTENT(OUT) :: PFAFF
        REAL(doubleprec), INTENT(OUT) :: RWORK( * )
        COMPLEX(doubleprec), INTENT(INOUT) :: AB( LDAB, * )
        COMPLEX(doubleprec), INTENT(OUT) :: WORK( * )
      END SUBROUTINE ZSKBPFA
   END INTERFACE

   INTERFACE SKBPF10
      SUBROUTINE SSKBPF10( UPLO, N, KD, AB, LDAB, PFAFF, WORK, INFO )
        USE PFAPACK_PREC, ONLY: singleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO
        INTEGER, INTENT(IN) :: KD, LDAB, N
        INTEGER, INTENT(OUT) :: INFO
        REAL(singleprec), INTENT(OUT) :: PFAFF( 2 )
        REAL(singleprec), INTENT(INOUT) :: AB( LDAB, * )
        REAL(singleprec), INTENT(OUT) :: WORK( * )
      END SUBROUTINE SSKBPF10

      SUBROUTINE DSKBPF10( UPLO, N, KD, AB, LDAB, PFAFF, WORK, INFO )
        USE PFAPACK_PREC, ONLY: doubleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO
        INTEGER, INTENT(IN) :: KD, LDAB, N
        INTEGER, INTENT(OUT) :: INFO
        REAL(doubleprec), INTENT(OUT) :: PFAFF( 2 )
        REAL(doubleprec), INTENT(INOUT) :: AB( LDAB, * )
        REAL(doubleprec), INTENT(OUT) :: WORK( * )
      END SUBROUTINE DSKBPF10

      SUBROUTINE CSKBPF10( UPLO, N, KD, AB, LDAB, PFAFF, WORK, &
                          RWORK, INFO )
        USE PFAPACK_PREC, ONLY: singleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO
        INTEGER, INTENT(IN) :: KD, LDAB, N
        INTEGER, INTENT(OUT) :: INFO
        COMPLEX(singleprec), INTENT(OUT) :: PFAFF( 2 )
        REAL(singleprec), INTENT(OUT) :: RWORK( * )
        COMPLEX(singleprec), INTENT(INOUT) :: AB( LDAB, * )
        COMPLEX(singleprec), INTENT(OUT) :: WORK( * )
      END SUBROUTINE CSKBPF10

      SUBROUTINE ZSKBPF10( UPLO, N, KD, AB, LDAB, PFAFF, WORK, &
                          RWORK, INFO )
        USE PFAPACK_PREC, ONLY: doubleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO
        INTEGER, INTENT(IN) :: KD, LDAB, N
        INTEGER, INTENT(OUT) :: INFO
        COMPLEX(doubleprec), INTENT(OUT) :: PFAFF( 2 )
        REAL(doubleprec), INTENT(OUT) :: RWORK( * )
        COMPLEX(doubleprec), INTENT(INOUT) :: AB( LDAB, * )
        COMPLEX(doubleprec), INTENT(OUT) :: WORK( * )
      END SUBROUTINE ZSKBPF10
   END INTERFACE

   !Interfaces for the Parlett-Reid tridiagonalization routines

   INTERFACE SKTRF
      SUBROUTINE SSKTRF( UPLO, MODE, N, A, LDA, IPIV, &
                         WORK, LWORK, INFO)
        USE PFAPACK_PREC, ONLY: singleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO, MODE
        INTEGER, INTENT(IN) :: LDA, LWORK, N
        INTEGER, INTENT(OUT) :: INFO
        INTEGER, INTENT(OUT) :: IPIV( * )
        REAL(singleprec), INTENT(INOUT) :: A( LDA, * )
        REAL(singleprec), INTENT(OUT) :: WORK( * )
      END SUBROUTINE SSKTRF

      SUBROUTINE DSKTRF( UPLO, MODE, N, A, LDA, IPIV, &
                         WORK, LWORK, INFO)
        USE PFAPACK_PREC, ONLY: doubleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO, MODE
        INTEGER, INTENT(IN) :: LDA, LWORK, N
        INTEGER, INTENT(OUT) :: INFO
        INTEGER, INTENT(OUT) :: IPIV( * )
        REAL(doubleprec), INTENT(INOUT) :: A( LDA, * )
        REAL(doubleprec), INTENT(OUT) :: WORK( * )
      END SUBROUTINE DSKTRF

      SUBROUTINE CSKTRF( UPLO, MODE, N, A, LDA, IPIV, &
                         WORK, LWORK, INFO)
        USE PFAPACK_PREC, ONLY: singleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO, MODE
        INTEGER, INTENT(IN) :: LDA, LWORK, N
        INTEGER, INTENT(OUT) :: INFO
        INTEGER, INTENT(OUT) :: IPIV( * )
        COMPLEX(singleprec), INTENT(INOUT) :: A( LDA, * )
        COMPLEX(singleprec), INTENT(OUT) :: WORK( * )
      END SUBROUTINE CSKTRF

      SUBROUTINE ZSKTRF( UPLO, MODE, N, A, LDA, IPIV, &
                         WORK, LWORK, INFO)
        USE PFAPACK_PREC, ONLY: doubleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO, MODE
        INTEGER, INTENT(IN) :: LDA, LWORK, N
        INTEGER, INTENT(OUT) :: INFO
        INTEGER, INTENT(OUT) :: IPIV( * )
        COMPLEX(doubleprec), INTENT(INOUT) :: A( LDA, * )
        COMPLEX(doubleprec), INTENT(OUT) :: WORK( * )
      END SUBROUTINE ZSKTRF
   END INTERFACE

   INTERFACE SKTF2
      SUBROUTINE SSKTF2( UPLO, MODE, N, A, LDA, IPIV, INFO )
        USE PFAPACK_PREC, ONLY: singleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO, MODE
        INTEGER, INTENT(IN) :: LDA, N
        INTEGER, INTENT(OUT) :: INFO
        INTEGER, INTENT(OUT) :: IPIV( * )
        REAL(singleprec), INTENT(INOUT) :: A( LDA, * )
      END SUBROUTINE SSKTF2

      SUBROUTINE DSKTF2( UPLO, MODE, N, A, LDA, IPIV, INFO )
        USE PFAPACK_PREC, ONLY: doubleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO, MODE
        INTEGER, INTENT(IN) :: LDA, N
        INTEGER, INTENT(OUT) :: INFO
        INTEGER, INTENT(OUT) :: IPIV( * )
        REAL(doubleprec), INTENT(INOUT) :: A( LDA, * )
      END SUBROUTINE DSKTF2

      SUBROUTINE CSKTF2( UPLO, MODE, N, A, LDA, IPIV, INFO )
        USE PFAPACK_PREC, ONLY: singleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO, MODE
        INTEGER, INTENT(IN) :: LDA, N
        INTEGER, INTENT(OUT) :: INFO
        INTEGER, INTENT(OUT) :: IPIV( * )
        COMPLEX(singleprec), INTENT(INOUT) :: A( LDA, * )
      END SUBROUTINE CSKTF2

      SUBROUTINE ZSKTF2( UPLO, MODE, N, A, LDA, IPIV, INFO )
        USE PFAPACK_PREC, ONLY: doubleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO, MODE
        INTEGER, INTENT(IN) :: LDA, N
        INTEGER, INTENT(OUT) :: INFO
        INTEGER, INTENT(OUT) :: IPIV( * )
        COMPLEX(doubleprec), INTENT(INOUT) :: A( LDA, * )
      END SUBROUTINE ZSKTF2
   END INTERFACE

   !Interfaces for the Householder tridiagonalization routines

   INTERFACE SKTRD
      SUBROUTINE SSKTRD( UPLO, MODE, N, A, LDA, E, TAU, WORK, LWORK, INFO )
        USE PFAPACK_PREC, ONLY: singleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO, MODE
        INTEGER, INTENT(IN) :: LDA, LWORK, N
        INTEGER, INTENT(OUT) :: INFO
        REAL(singleprec), INTENT(OUT) :: E( * ), TAU( * ), WORK( * )
        REAL(singleprec), INTENT(INOUT) :: A( LDA, * )
      END SUBROUTINE SSKTRD

      SUBROUTINE DSKTRD( UPLO, MODE, N, A, LDA, E, TAU, WORK, LWORK, INFO )
        USE PFAPACK_PREC, ONLY: doubleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO, MODE
        INTEGER, INTENT(IN) :: LDA, LWORK, N
        INTEGER, INTENT(OUT) :: INFO
        REAL(doubleprec), INTENT(OUT) :: E( * ), TAU( * ), WORK( * )
        REAL(doubleprec), INTENT(INOUT) :: A( LDA, * )
      END SUBROUTINE DSKTRD

      SUBROUTINE CSKTRD( UPLO, MODE, N, A, LDA, E, TAU, WORK, LWORK, INFO )
        USE PFAPACK_PREC, ONLY: singleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO, MODE
        INTEGER, INTENT(IN) :: LDA, LWORK, N
        INTEGER, INTENT(OUT) :: INFO
        REAL(singleprec), INTENT(OUT) :: E( * )
        COMPLEX(singleprec), INTENT(OUT) :: TAU( * ), WORK( * )
        COMPLEX(singleprec), INTENT(INOUT) :: A( LDA, * )
      END SUBROUTINE CSKTRD

      SUBROUTINE ZSKTRD( UPLO, MODE, N, A, LDA, E, TAU, WORK, LWORK, INFO )
        USE PFAPACK_PREC, ONLY: doubleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO, MODE
        INTEGER, INTENT(IN) :: LDA, LWORK, N
        INTEGER, INTENT(OUT) :: INFO
        REAL(doubleprec), INTENT(OUT) :: E( * )
        COMPLEX(doubleprec), INTENT(OUT) :: TAU( * ), WORK( * )
        COMPLEX(doubleprec), INTENT(INOUT) :: A( LDA, * )
      END SUBROUTINE ZSKTRD
   END INTERFACE

   INTERFACE SKTD2
      SUBROUTINE SSKTD2( UPLO, MODE, N, A, LDA, E, TAU, INFO )
        USE PFAPACK_PREC, ONLY: singleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO, MODE
        INTEGER, INTENT(IN) :: LDA, N
        INTEGER, INTENT(OUT) :: INFO
        REAL(singleprec), INTENT(OUT) :: E( * ), TAU( * )
        REAL(singleprec), INTENT(INOUT) :: A( LDA, * )
      END SUBROUTINE SSKTD2

      SUBROUTINE DSKTD2( UPLO, MODE, N, A, LDA, E, TAU, INFO )
        USE PFAPACK_PREC, ONLY: doubleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO, MODE
        INTEGER, INTENT(IN) :: LDA, N
        INTEGER, INTENT(OUT) :: INFO
        REAL(doubleprec), INTENT(OUT) :: E( * ), TAU( * )
        REAL(doubleprec), INTENT(INOUT) :: A( LDA, * )
      END SUBROUTINE DSKTD2

      SUBROUTINE CSKTD2( UPLO, MODE, N, A, LDA, E, TAU, INFO )
        USE PFAPACK_PREC, ONLY: singleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO, MODE
        INTEGER, INTENT(IN) :: LDA, N
        INTEGER, INTENT(OUT) :: INFO
        REAL(singleprec), INTENT(OUT) :: E( * )
        COMPLEX(singleprec), INTENT(OUT) :: TAU( * )
        COMPLEX(singleprec), INTENT(INOUT) :: A( LDA, * )
      END SUBROUTINE CSKTD2

      SUBROUTINE ZSKTD2( UPLO, MODE, N, A, LDA, E, TAU, INFO )
        USE PFAPACK_PREC, ONLY: doubleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO, MODE
        INTEGER, INTENT(IN) :: LDA, N
        INTEGER, INTENT(OUT) :: INFO
        REAL(doubleprec), INTENT(OUT) :: E( * )
        COMPLEX(doubleprec), INTENT(OUT) :: TAU( * )
        COMPLEX(doubleprec), INTENT(INOUT) :: A( LDA, * )
      END SUBROUTINE ZSKTD2
   END INTERFACE

   INTERFACE SKBTRD
      SUBROUTINE SSKBTRD( VECT, UPLO, MODE, N, KD, AB, LDAB, E, &
                          Q, LDQ, WORK, INFO )
        USE PFAPACK_PREC, ONLY: singleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO, VECT, MODE
        INTEGER, INTENT(IN) :: KD, LDAB, LDQ, N
        INTEGER, INTENT(OUT) :: INFO
        REAL(singleprec), INTENT(OUT) :: E( * ), Q( LDQ, * ), WORK( * )
        REAL(singleprec), INTENT(INOUT) :: AB( LDAB, * )
      END SUBROUTINE SSKBTRD

      SUBROUTINE DSKBTRD( VECT, UPLO, MODE, N, KD, AB, LDAB, E, &
                          Q, LDQ, WORK, INFO )
        USE PFAPACK_PREC, ONLY: doubleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO, VECT, MODE
        INTEGER, INTENT(IN) :: KD, LDAB, LDQ, N
        INTEGER, INTENT(OUT) :: INFO
        REAL(doubleprec), INTENT(OUT) :: E( * ), Q( LDQ, * ), WORK( * )
        REAL(doubleprec), INTENT(INOUT) :: AB( LDAB, * )
      END SUBROUTINE DSKBTRD

      SUBROUTINE CSKBTRD( VECT, UPLO, MODE, N, KD, AB, LDAB, E, DETQ, &
                          Q, LDQ, WORK, RWORK, INFO )
        USE PFAPACK_PREC, ONLY: singleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO, VECT, MODE
        INTEGER, INTENT(IN) :: KD, LDAB, LDQ, N
        INTEGER, INTENT(OUT) :: INFO
        COMPLEX(singleprec), INTENT(OUT) :: DETQ
        REAL(singleprec), INTENT(OUT) :: E( * ), RWORK( * )
        COMPLEX(singleprec), INTENT(OUT) :: Q( LDQ, * ), WORK( * )
        COMPLEX(singleprec), INTENT(INOUT) :: AB( LDAB, * )
      END SUBROUTINE CSKBTRD

      SUBROUTINE ZSKBTRD( VECT, UPLO, MODE, N, KD, AB, LDAB, E, DETQ, &
                          Q, LDQ, WORK, RWORK, INFO )
        USE PFAPACK_PREC, ONLY: doubleprec
        CHARACTER(LEN=1), INTENT(IN) :: UPLO, VECT, MODE
        INTEGER, INTENT(IN) :: KD, LDAB, LDQ, N
        INTEGER, INTENT(OUT) :: INFO
        COMPLEX(doubleprec), INTENT(OUT) :: DETQ
        REAL(doubleprec), INTENT(OUT) :: E( * ), RWORK( * )
        COMPLEX(doubleprec), INTENT(OUT) :: Q( LDQ, * ), WORK( * )
        COMPLEX(doubleprec), INTENT(INOUT) :: AB( LDAB, * )
      END SUBROUTINE ZSKBTRD
   END INTERFACE

 END MODULE F77_PFAPACK
