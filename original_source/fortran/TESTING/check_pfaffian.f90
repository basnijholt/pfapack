MODULE check_pfaffian
  INTERFACE check_pfaff_F77
     SUBROUTINE check_pfaff_F77_s(UPLO, METHOD, B, RESIDUAL, RELW)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO, METHOD
       REAL(KIND(1.0E0)) :: B(:,:)
       REAL(KIND(1.0E0)) :: RESIDUAL, RELW
     END SUBROUTINE check_pfaff_F77_s

     SUBROUTINE check_pfaff_F77_d(UPLO, METHOD, B, RESIDUAL, RELW)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO, METHOD
       REAL(KIND(1.0D0)) :: B(:,:)
       REAL(KIND(1.0D0)) :: RESIDUAL, RELW
     END SUBROUTINE check_pfaff_F77_d

     SUBROUTINE check_pfaff_F77_c(UPLO, METHOD, B, RESIDUAL, RELW)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO, METHOD
       COMPLEX(KIND(1.0E0)) :: B(:,:)
       REAL(KIND(1.0E0)) :: RESIDUAL, RELW
     END SUBROUTINE check_pfaff_F77_c

     SUBROUTINE check_pfaff_F77_z(UPLO, METHOD, B, RESIDUAL, RELW)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO, METHOD
       COMPLEX(KIND(1.0D0)) :: B(:,:)
       REAL(KIND(1.0D0)) :: RESIDUAL, RELW
     END SUBROUTINE check_pfaff_F77_z
  END INTERFACE

  INTERFACE check_pfaff_F95
     SUBROUTINE check_pfaff_F95_s(UPLO, METHOD, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO, METHOD
       REAL(KIND(1.0E0)) :: B(:,:)
       REAL(KIND(1.0E0)) :: RESIDUAL
     END SUBROUTINE check_pfaff_F95_s

     SUBROUTINE check_pfaff_F95_d(UPLO, METHOD, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO, METHOD
       REAL(KIND(1.0D0)) :: B(:,:)
       REAL(KIND(1.0D0)) :: RESIDUAL
     END SUBROUTINE check_pfaff_F95_d

     SUBROUTINE check_pfaff_F95_c(UPLO, METHOD, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO, METHOD
       COMPLEX(KIND(1.0E0)) :: B(:,:)
       REAL(KIND(1.0E0)) :: RESIDUAL
     END SUBROUTINE check_pfaff_F95_c

     SUBROUTINE check_pfaff_F95_z(UPLO, METHOD, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO, METHOD
       COMPLEX(KIND(1.0D0)) :: B(:,:)
       REAL(KIND(1.0D0)) :: RESIDUAL
     END SUBROUTINE check_pfaff_F95_z
  END INTERFACE

  INTERFACE check_pfaff10_F77
     SUBROUTINE check_pfaff10_F77_s(UPLO, METHOD, B, RESIDUAL, RELW)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO, METHOD
       REAL(KIND(1.0E0)) :: B(:,:)
       REAL(KIND(1.0E0)) :: RESIDUAL, RELW
     END SUBROUTINE check_pfaff10_F77_s

     SUBROUTINE check_pfaff10_F77_d(UPLO, METHOD, B, RESIDUAL, RELW)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO, METHOD
       REAL(KIND(1.0D0)) :: B(:,:)
       REAL(KIND(1.0D0)) :: RESIDUAL, RELW
     END SUBROUTINE check_pfaff10_F77_d

     SUBROUTINE check_pfaff10_F77_c(UPLO, METHOD, B, RESIDUAL, RELW)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO, METHOD
       COMPLEX(KIND(1.0E0)) :: B(:,:)
       REAL(KIND(1.0E0)) :: RESIDUAL, RELW
     END SUBROUTINE check_pfaff10_F77_c

     SUBROUTINE check_pfaff10_F77_z(UPLO, METHOD, B, RESIDUAL, RELW)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO, METHOD
       COMPLEX(KIND(1.0D0)) :: B(:,:)
       REAL(KIND(1.0D0)) :: RESIDUAL, RELW
     END SUBROUTINE check_pfaff10_F77_z
  END INTERFACE

  INTERFACE check_pfaff10_F95
     SUBROUTINE check_pfaff10_F95_s(UPLO, METHOD, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO, METHOD
       REAL(KIND(1.0E0)) :: B(:,:)
       REAL(KIND(1.0E0)) :: RESIDUAL
     END SUBROUTINE check_pfaff10_F95_s

     SUBROUTINE check_pfaff10_F95_d(UPLO, METHOD, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO, METHOD
       REAL(KIND(1.0D0)) :: B(:,:)
       REAL(KIND(1.0D0)) :: RESIDUAL
     END SUBROUTINE check_pfaff10_F95_d

     SUBROUTINE check_pfaff10_F95_c(UPLO, METHOD, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO, METHOD
       COMPLEX(KIND(1.0E0)) :: B(:,:)
       REAL(KIND(1.0E0)) :: RESIDUAL
     END SUBROUTINE check_pfaff10_F95_c

     SUBROUTINE check_pfaff10_F95_z(UPLO, METHOD, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO, METHOD
       COMPLEX(KIND(1.0D0)) :: B(:,:)
       REAL(KIND(1.0D0)) :: RESIDUAL
     END SUBROUTINE check_pfaff10_F95_z
  END INTERFACE

  INTERFACE check_pfaff_band_F77
     SUBROUTINE check_pfaff_band_F77_s(UPLO, A, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO
       REAL(KIND(1.0E0)) :: A(:,:), B(:,:)
       REAL(KIND(1.0E0)) :: RESIDUAL
     END SUBROUTINE check_pfaff_band_F77_s

     SUBROUTINE check_pfaff_band_F77_d(UPLO, A, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO
       REAL(KIND(1.0D0)) :: A(:,:), B(:,:)
       REAL(KIND(1.0D0)) :: RESIDUAL
     END SUBROUTINE check_pfaff_band_F77_d

     SUBROUTINE check_pfaff_band_F77_c(UPLO, A, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO
       COMPLEX(KIND(1.0E0)) :: A(:,:), B(:,:)
       REAL(KIND(1.0E0)) :: RESIDUAL
     END SUBROUTINE check_pfaff_band_F77_c

     SUBROUTINE check_pfaff_band_F77_z(UPLO, A, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO
       COMPLEX(KIND(1.0D0)) :: A(:,:), B(:,:)
       REAL(KIND(1.0D0)) :: RESIDUAL
     END SUBROUTINE check_pfaff_band_F77_z
  END INTERFACE

  INTERFACE check_pfaff_band_F95
     SUBROUTINE check_pfaff_band_F95_s(UPLO, A, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO
       REAL(KIND(1.0E0)) :: A(:,:), B(:,:)
       REAL(KIND(1.0E0)) :: RESIDUAL
     END SUBROUTINE check_pfaff_band_F95_s

     SUBROUTINE check_pfaff_band_F95_d(UPLO, A, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO
       REAL(KIND(1.0D0)) :: A(:,:), B(:,:)
       REAL(KIND(1.0D0)) :: RESIDUAL
     END SUBROUTINE check_pfaff_band_F95_d

     SUBROUTINE check_pfaff_band_F95_c(UPLO, A, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO
       COMPLEX(KIND(1.0E0)) :: A(:,:), B(:,:)
       REAL(KIND(1.0E0)) :: RESIDUAL
     END SUBROUTINE check_pfaff_band_F95_c

     SUBROUTINE check_pfaff_band_F95_z(UPLO, A, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO
       COMPLEX(KIND(1.0D0)) :: A(:,:), B(:,:)
       REAL(KIND(1.0D0)) :: RESIDUAL
     END SUBROUTINE check_pfaff_band_F95_z
  END INTERFACE

  INTERFACE check_pfaff10_band_F77
     SUBROUTINE check_pfaff10_band_F77_s(UPLO, A, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO
       REAL(KIND(1.0E0)) :: A(:,:), B(:,:)
       REAL(KIND(1.0E0)) :: RESIDUAL
     END SUBROUTINE check_pfaff10_band_F77_s

     SUBROUTINE check_pfaff10_band_F77_d(UPLO, A, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO
       REAL(KIND(1.0D0)) :: A(:,:), B(:,:)
       REAL(KIND(1.0D0)) :: RESIDUAL
     END SUBROUTINE check_pfaff10_band_F77_d

     SUBROUTINE check_pfaff10_band_F77_c(UPLO, A, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO
       COMPLEX(KIND(1.0E0)) :: A(:,:), B(:,:)
       REAL(KIND(1.0E0)) :: RESIDUAL
     END SUBROUTINE check_pfaff10_band_F77_c

     SUBROUTINE check_pfaff10_band_F77_z(UPLO, A, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO
       COMPLEX(KIND(1.0D0)) :: A(:,:), B(:,:)
       REAL(KIND(1.0D0)) :: RESIDUAL
     END SUBROUTINE check_pfaff10_band_F77_z
  END INTERFACE

  INTERFACE check_pfaff10_band_F95
     SUBROUTINE check_pfaff10_band_F95_s(UPLO, A, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO
       REAL(KIND(1.0E0)) :: A(:,:), B(:,:)
       REAL(KIND(1.0E0)) :: RESIDUAL
     END SUBROUTINE check_pfaff10_band_F95_s

     SUBROUTINE check_pfaff10_band_F95_d(UPLO, A, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO
       REAL(KIND(1.0D0)) :: A(:,:), B(:,:)
       REAL(KIND(1.0D0)) :: RESIDUAL
     END SUBROUTINE check_pfaff10_band_F95_d

     SUBROUTINE check_pfaff10_band_F95_c(UPLO, A, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO
       COMPLEX(KIND(1.0E0)) :: A(:,:), B(:,:)
       REAL(KIND(1.0E0)) :: RESIDUAL
     END SUBROUTINE check_pfaff10_band_F95_c

     SUBROUTINE check_pfaff10_band_F95_z(UPLO, A, B, RESIDUAL)
       IMPLICIT NONE
       CHARACTER(LEN=1) :: UPLO
       COMPLEX(KIND(1.0D0)) :: A(:,:), B(:,:)
       REAL(KIND(1.0D0)) :: RESIDUAL
     END SUBROUTINE check_pfaff10_band_F95_z
  END INTERFACE

END MODULE check_pfaffian
