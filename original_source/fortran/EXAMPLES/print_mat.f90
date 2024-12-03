!Little helper function to print a matrix on the screen
MODULE print_matrix
  INTERFACE print_mat
     SUBROUTINE print_mat_d(A)
       IMPLICIT NONE
       REAL(KIND(1.0D0)) :: A (:,:)
     END SUBROUTINE print_mat_d

     SUBROUTINE print_mat_z(A)
       IMPLICIT NONE
       COMPLEX(KIND(1.0D0)) :: A (:,:)
     END SUBROUTINE print_mat_z
  END INTERFACE

  INTERFACE print_skew_mat
     SUBROUTINE print_skew_mat_d(A, UPLO)
       IMPLICIT NONE
       REAL(KIND(1.0D0)) :: A (:,:)
       CHARACTER(LEN=1) :: UPLO
     END SUBROUTINE print_skew_mat_d

     SUBROUTINE print_skew_mat_z(A, UPLO)
       IMPLICIT NONE
       COMPLEX(KIND(1.0D0)) :: A (:,:)
       CHARACTER(LEN=1) :: UPLO
     END SUBROUTINE print_skew_mat_z
  END INTERFACE
END MODULE print_matrix
   
!print a full matrix on the screen
SUBROUTINE print_mat_d(A)
  IMPLICIT NONE
  REAL(KIND(1.0D0)) :: A (:,:)

  INTEGER :: N,M,I,J

  N=SIZE(A,2)
  M=SIZE(A,1)

  DO I=1, M
     DO J=1, N
        WRITE (*,'(F6.3)',Advance='NO') A(I,J)
        WRITE (*,'(A)',Advance='NO') "  "
     END DO
     WRITE (*,*)
  END DO

  WRITE (*,*)
END SUBROUTINE print_mat_d

SUBROUTINE print_mat_z(A)
  IMPLICIT NONE
  COMPLEX(KIND(1.0D0)) :: A (:,:)

  INTEGER :: N,M,I,J

  N=SIZE(A,2)
  M=SIZE(A,1)

  DO I=1, M
     DO J=1, N
        WRITE (*,'(A)',Advance='NO') "("
        WRITE (*,'(F6.2, F6.2)',Advance='NO') A(I,J)
        WRITE (*,'(A)',Advance='NO') " )  "
     END DO
     WRITE (*,*)
  END DO

  WRITE (*,*)
END SUBROUTINE print_mat_z

!print a skew-symmetric matrix on the screen
!only accesses either the upper or lower triangle
SUBROUTINE print_skew_mat_d(A, UPLO)
  IMPLICIT NONE
  REAL(KIND(1.0D0)) :: A (:,:)
  CHARACTER(LEN=1) :: UPLO

  INTEGER :: N,M,I,J

  N=SIZE(A,2)
  M=SIZE(A,1)

  DO I=1, M
     DO J=1, N
        IF( UPLO == 'U' ) THEN
           IF( I<=J ) THEN
              WRITE (*,'(F6.3)',Advance='NO') A(I,J)
           ELSE 
              WRITE (*,'(F6.3)',Advance='NO') -A(J,I)
           END IF
        ELSE
           IF( I>=J ) THEN
              WRITE (*,'(F6.3)',Advance='NO') A(I,J)
           ELSE 
              WRITE (*,'(F6.3)',Advance='NO') -A(J,I)
           END IF
        END IF
        WRITE (*,'(A)',Advance='NO') "  "
     END DO
     WRITE (*,*)
  END DO

  WRITE (*,*)
END SUBROUTINE print_skew_mat_d

SUBROUTINE print_skew_mat_z(A, UPLO)
  IMPLICIT NONE
  COMPLEX(KIND(1.0D0)) :: A (:,:)
  CHARACTER(LEN=1) :: UPLO

  INTEGER :: N,M,I,J

  N=SIZE(A,2)
  M=SIZE(A,1)

  DO I=1, M
     DO J=1, N
        WRITE (*,'(A)',Advance='NO') "("
        IF( UPLO == 'U' ) THEN
           IF( I<=J ) THEN
              WRITE (*,'(F6.2, F6.2)',Advance='NO') A(I,J)
           ELSE 
              WRITE (*,'(F6.2, F6.2)',Advance='NO') -A(J,I)
           END IF
        ELSE
           IF( I>=J ) THEN
              WRITE (*,'(F6.2, F6.2)',Advance='NO') A(I,J)
           ELSE 
              WRITE (*,'(F6.2, F6.2)',Advance='NO') -A(J,I)
           END IF
        END IF
        WRITE (*,'(A)',Advance='NO') " )  "
     END DO
     WRITE (*,*)
  END DO

  WRITE (*,*)
END SUBROUTINE print_skew_mat_z
