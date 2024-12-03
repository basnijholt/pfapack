  MODULE matrix_tools

    !print a matrix
    INTERFACE print_mat
       SUBROUTINE print_mat_s(A)
         IMPLICIT NONE
         REAL(KIND(1.0E0)) :: A (:,:)
       END SUBROUTINE print_mat_s

       SUBROUTINE print_mat_d(A)
         IMPLICIT NONE
         REAL(KIND(1.0D0)) :: A (:,:)
       END SUBROUTINE print_mat_d

       SUBROUTINE print_mat_c(A)
         IMPLICIT NONE
         COMPLEX(KIND(1.0E0)) :: A (:,:)
       END SUBROUTINE print_mat_c

       SUBROUTINE print_mat_z(A)
         IMPLICIT NONE
         COMPLEX(KIND(1.0D0)) :: A (:,:)
       END SUBROUTINE print_mat_z
    END INTERFACE

    !compute the norm of a matrix
    INTERFACE one_norm
       REAL(KIND(1.0E0)) FUNCTION one_norm_s(A)
         IMPLICIT NONE
         REAL(KIND(1.0E0)) :: A (:,:)
       END FUNCTION one_norm_s

       REAL(KIND(1.0D0)) FUNCTION one_norm_d(A)
         IMPLICIT NONE
         REAL(KIND(1.0D0)) :: A (:,:)
       END FUNCTION one_norm_d

       REAL(KIND(1.0E0)) FUNCTION one_norm_c(A)
         IMPLICIT NONE
         COMPLEX(KIND(1.0E0)) :: A (:,:)
       END FUNCTION one_norm_c

       REAL(KIND(1.0D0)) FUNCTION one_norm_z(A)
         IMPLICIT NONE
         COMPLEX(KIND(1.0D0)) :: A (:,:)
       END FUNCTION one_norm_z
    END INTERFACE

    !compute the relative residual of the difference of two matrices
    INTERFACE resid
       REAL(KIND(1.0E0)) FUNCTION resid_s(A,B)
         IMPLICIT NONE
         REAL(KIND(1.0E0)) :: A (:,:), B(:,:)
       END FUNCTION resid_s

       REAL(KIND(1.0D0)) FUNCTION resid_d(A,B)
         IMPLICIT NONE
         REAL(KIND(1.0D0)) :: A (:,:), B(:,:)
       END FUNCTION resid_d

       REAL(KIND(1.0E0)) FUNCTION resid_c(A,B)
         IMPLICIT NONE
         COMPLEX(KIND(1.0E0)) :: A (:,:), B(:,:)
       END FUNCTION resid_c

       REAL(KIND(1.0D0)) FUNCTION resid_z(A,B)
         IMPLICIT NONE
         COMPLEX(KIND(1.0D0)) :: A (:,:), B(:,:)
       END FUNCTION resid_z
    END INTERFACE

    !make a random, skew-symmetric matrix
    INTERFACE make_skew_mat
       SUBROUTINE make_skew_mat_s(A, PROB)
         IMPLICIT NONE
         REAL(KIND(1.0E0)) :: A (:,:)
         REAL(KIND(1.0E0)), OPTIONAL :: PROB
       END SUBROUTINE make_skew_mat_s

       SUBROUTINE make_skew_mat_d(A, PROB)
         IMPLICIT NONE
         REAL(KIND(1.0D0)) :: A (:,:)
         REAL(KIND(1.0D0)), OPTIONAL :: PROB
       END SUBROUTINE make_skew_mat_d

       SUBROUTINE make_skew_mat_c(A, PROB)
         IMPLICIT NONE
         COMPLEX(KIND(1.0E0)) :: A (:,:)
         REAL(KIND(1.0E0)), OPTIONAL :: PROB
       END SUBROUTINE make_skew_mat_c

       SUBROUTINE make_skew_mat_z(A, PROB)
         IMPLICIT NONE
         COMPLEX(KIND(1.0D0)) :: A (:,:)
         REAL(KIND(1.0D0)), OPTIONAL :: PROB
       END SUBROUTINE make_skew_mat_z
    END INTERFACE

    INTERFACE make_skew_mat_banded
       SUBROUTINE make_skew_mat_banded_s(Au, Ad, A, PROB)
         IMPLICIT NONE
         REAL(KIND(1.0E0)) :: Au(:,:), Ad(:,:), A(:,:)
         REAL(KIND(1.0E0)), OPTIONAL :: PROB
       END SUBROUTINE make_skew_mat_banded_s

       SUBROUTINE make_skew_mat_banded_d(Au, Ad, A, PROB)
         IMPLICIT NONE
         REAL(KIND(1.0D0)) :: Au(:,:), Ad(:,:), A(:,:)
         REAL(KIND(1.0D0)), OPTIONAL :: PROB
       END SUBROUTINE make_skew_mat_banded_d

       SUBROUTINE make_skew_mat_banded_c(Au, Ad, A, PROB)
         IMPLICIT NONE
         COMPLEX(KIND(1.0E0)) :: Au(:,:), Ad(:,:), A(:,:)
         REAL(KIND(1.0E0)), OPTIONAL :: PROB
       END SUBROUTINE make_skew_mat_banded_c

       SUBROUTINE make_skew_mat_banded_z(Au, Ad, A, PROB)
         IMPLICIT NONE
         COMPLEX(KIND(1.0D0)) :: Au(:,:), Ad(:,:), A(:,:)
         REAL(KIND(1.0D0)), OPTIONAL :: PROB
       END SUBROUTINE make_skew_mat_banded_z
    END INTERFACE

    !extract lower and upper unit triangular matrices
    INTERFACE extract_unit_tri
       SUBROUTINE extract_unit_tri_s(UPLO, A, B)
         IMPLICIT NONE
         CHARACTER(LEN=1) :: UPLO
         REAL(KIND(1.0E0)) :: A(:,:), B(:,:)
       END SUBROUTINE extract_unit_tri_s

       SUBROUTINE extract_unit_tri_d(UPLO, A, B)
         IMPLICIT NONE
         CHARACTER(LEN=1) :: UPLO
         REAL(KIND(1.0D0)) :: A(:,:), B(:,:)
       END SUBROUTINE extract_unit_tri_d

       SUBROUTINE extract_unit_tri_c(UPLO, A, B)
         IMPLICIT NONE
         CHARACTER(LEN=1) :: UPLO
         COMPLEX(KIND(1.0E0)) :: A(:,:), B(:,:)
       END SUBROUTINE extract_unit_tri_c

       SUBROUTINE extract_unit_tri_z(UPLO, A, B)
         IMPLICIT NONE
         CHARACTER(LEN=1) :: UPLO
         COMPLEX(KIND(1.0D0)) :: A(:,:), B(:,:)
       END SUBROUTINE extract_unit_tri_z
    END INTERFACE

    INTERFACE extract_trid
       SUBROUTINE extract_trid_s(UPLO, A, T)
         IMPLICIT NONE
         CHARACTER(LEN=1) :: UPLO
         REAL(KIND(1.0E0)) :: A(:,:), T(:,:)
       END SUBROUTINE extract_trid_s

       SUBROUTINE extract_trid_d(UPLO, A, T)
         IMPLICIT NONE
         CHARACTER(LEN=1) :: UPLO
         REAL(KIND(1.0D0)) :: A(:,:), T(:,:)
       END SUBROUTINE extract_trid_d

       SUBROUTINE extract_trid_c(UPLO, A, T)
         IMPLICIT NONE
         CHARACTER(LEN=1) :: UPLO
         COMPLEX(KIND(1.0E0)) :: A(:,:), T(:,:)
       END SUBROUTINE extract_trid_c

       SUBROUTINE extract_trid_z(UPLO, A, T)
         IMPLICIT NONE
         CHARACTER(LEN=1) :: UPLO
         COMPLEX(KIND(1.0D0)) :: A(:,:), T(:,:)
       END SUBROUTINE extract_trid_z
    END INTERFACE

    INTERFACE extract_trid_band
       SUBROUTINE extract_trid_band_s(UPLO, A, T)
         IMPLICIT NONE
         CHARACTER(LEN=1) :: UPLO
         REAL(KIND(1.0E0)) :: A(:,:), T(:,:)
       END SUBROUTINE extract_trid_band_s

       SUBROUTINE extract_trid_band_d(UPLO, A, T)
         IMPLICIT NONE
         CHARACTER(LEN=1) :: UPLO
         REAL(KIND(1.0D0)) :: A(:,:), T(:,:)
       END SUBROUTINE extract_trid_band_d

       SUBROUTINE extract_trid_band_c(UPLO, A, T)
         IMPLICIT NONE
         CHARACTER(LEN=1) :: UPLO
         COMPLEX(KIND(1.0E0)) :: A(:,:), T(:,:)
       END SUBROUTINE extract_trid_band_c

       SUBROUTINE extract_trid_band_z(UPLO, A, T)
         IMPLICIT NONE
         CHARACTER(LEN=1) :: UPLO
         COMPLEX(KIND(1.0D0)) :: A(:,:), T(:,:)
       END SUBROUTINE extract_trid_band_z
    END INTERFACE

    !build the orthogonal/unitary transformation from individual Householder
    INTERFACE extract_unitary
       SUBROUTINE extract_unitary_s(UPLO, A, TAU, Q)
         IMPLICIT NONE
         CHARACTER(LEN=1) :: UPLO
         REAL(KIND(1.0E0)) :: A(:,:), Q(:,:), TAU(:)
       END SUBROUTINE extract_unitary_s

       SUBROUTINE extract_unitary_d(UPLO, A, TAU, Q)
         IMPLICIT NONE
         CHARACTER(LEN=1) :: UPLO
         REAL(KIND(1.0D0)) :: A(:,:), Q(:,:), TAU(:)
       END SUBROUTINE extract_unitary_d

       SUBROUTINE extract_unitary_c(UPLO, A, TAU, Q)
         IMPLICIT NONE
         CHARACTER(LEN=1) :: UPLO
         COMPLEX(KIND(1.0E0)) :: A(:,:), Q(:,:), TAU(:)
       END SUBROUTINE extract_unitary_c

       SUBROUTINE extract_unitary_z(UPLO, A, TAU, Q)
         IMPLICIT NONE
         CHARACTER(LEN=1) :: UPLO
         COMPLEX(KIND(1.0D0)) :: A(:,:), Q(:,:), TAU(:)
       END SUBROUTINE extract_unitary_z
    END INTERFACE

    INTERFACE rowcol_invperm
       SUBROUTINE rowcol_invperm_s(UPLO, A, IPIV)
         IMPLICIT NONE
         CHARACTER(LEN=1) :: UPLO
         REAL(KIND(1.0E0)) :: A(:,:)
         INTEGER :: IPIV(:)
       END SUBROUTINE rowcol_invperm_s

       SUBROUTINE rowcol_invperm_d(UPLO, A, IPIV)
         IMPLICIT NONE
         CHARACTER(LEN=1) :: UPLO
         REAL(KIND(1.0D0)) :: A(:,:)
         INTEGER :: IPIV(:)
       END SUBROUTINE rowcol_invperm_d

       SUBROUTINE rowcol_invperm_c(UPLO, A, IPIV)
         IMPLICIT NONE
         CHARACTER(LEN=1) :: UPLO
         COMPLEX(KIND(1.0E0)) :: A(:,:)
         INTEGER :: IPIV(:)
       END SUBROUTINE rowcol_invperm_c

       SUBROUTINE rowcol_invperm_z(UPLO, A, IPIV)
         IMPLICIT NONE
         CHARACTER(LEN=1) :: UPLO
         COMPLEX(KIND(1.0D0)) :: A(:,:)
         INTEGER :: IPIV(:)
       END SUBROUTINE rowcol_invperm_z
    END INTERFACE

    !Reference Pfaffian calculation
    INTERFACE reference_pfaffian
       SUBROUTINE reference_pfaff_s(A, PFAFF)
         IMPLICIT NONE
         REAL(KIND(1.0E0)) :: A (:,:), PFAFF
       END SUBROUTINE reference_pfaff_s

       SUBROUTINE reference_pfaff_d(A, PFAFF)
         IMPLICIT NONE
         REAL(KIND(1.0D0)) :: A (:,:), PFAFF
       END SUBROUTINE reference_pfaff_d

       SUBROUTINE reference_pfaff_c(A, PFAFF)
         IMPLICIT NONE
         COMPLEX(KIND(1.0E0)) :: A (:,:), PFAFF
       END SUBROUTINE reference_pfaff_c

       SUBROUTINE reference_pfaff_z(A, PFAFF)
         IMPLICIT NONE
         COMPLEX(KIND(1.0D0)) :: A (:,:), PFAFF
       END SUBROUTINE reference_pfaff_z
    END INTERFACE

    INTERFACE resid_pfaffian
       REAL(KIND(1.0E0)) FUNCTION resid_pfaff_s(N, A,B)
         IMPLICIT NONE
         INTEGER :: N
         REAL(KIND(1.0E0)) :: A, B
       END FUNCTION resid_pfaff_s

       REAL(KIND(1.0E0)) FUNCTION resid_pfaff10_s(N, A,B)
         IMPLICIT NONE
         INTEGER :: N
         REAL(KIND(1.0E0)) :: A(2), B
       END FUNCTION resid_pfaff10_s

       REAL(KIND(1.0D0)) FUNCTION resid_pfaff_d(N, A,B)
         IMPLICIT NONE
         INTEGER :: N
         REAL(KIND(1.0D0)) :: A, B
       END FUNCTION resid_pfaff_d

       REAL(KIND(1.0D0)) FUNCTION resid_pfaff10_d(N, A,B)
         IMPLICIT NONE
         INTEGER :: N
         REAL(KIND(1.0D0)) :: A(2), B
       END FUNCTION resid_pfaff10_d

       REAL(KIND(1.0E0)) FUNCTION resid_pfaff_c(N, A,B)
         IMPLICIT NONE
         INTEGER :: N
         COMPLEX(KIND(1.0E0)) :: A, B
       END FUNCTION resid_pfaff_c

       REAL(KIND(1.0E0)) FUNCTION resid_pfaff10_c(N, A,B)
         IMPLICIT NONE
         INTEGER :: N
         COMPLEX(KIND(1.0E0)) :: A(2), B
       END FUNCTION resid_pfaff10_c

       REAL(KIND(1.0D0)) FUNCTION resid_pfaff_z(N, A,B)
         IMPLICIT NONE
         INTEGER :: N
         COMPLEX(KIND(1.0D0)) :: A, B
       END FUNCTION resid_pfaff_z

       REAL(KIND(1.0D0)) FUNCTION resid_pfaff10_z(N, A,B)
         IMPLICIT NONE
         INTEGER :: N
         COMPLEX(KIND(1.0D0)) :: A(2), B
       END FUNCTION resid_pfaff10_z
    END INTERFACE

  END MODULE matrix_tools
