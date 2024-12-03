PROGRAM test_decompositions

  WRITE (*,*) "Checking single precision decompositions ..."
  CALL test_decomp_s()
  WRITE (*,*) "ok"

  WRITE (*,*) "Checking double precision decompositions ..."
  CALL test_decomp_d()
  WRITE (*,*) "ok"

  WRITE (*,*) "Checking single complex decompositions ..."
  CALL test_decomp_c()
  WRITE (*,*) "ok"

  WRITE (*,*) "Checking double complex decompositions ..."
  CALL test_decomp_z()
  WRITE (*,*) "ok"

  WRITE(*,*) "********************************************"
  WRITE(*,*) "All decomposition tests passed successfully!"
  WRITE(*,*) "********************************************"

END PROGRAM test_decompositions

SUBROUTINE test_decomp_s()
  USE check_decomp
  USE matrix_tools
  IMPLICIT NONE
  INTEGER, PARAMETER :: MAX_SIZE = 20
  INTEGER, DIMENSION(5) :: NB, NBMIN, NX, WORKSP
  INTEGER :: N, ISTAT, IPARAM, K, IDENS
  REAL(KIND(1.0E0)), ALLOCATABLE :: A(:,:), Au(:,:), Ad(:,:)
  REAL(KIND(1.0E0)) :: RES, DENSITY(5)

  !First check the dense decomposition interfaces with blocking
  !This includes different blocking parameters

  NB = (/1, 2, 5, 8, 25/)
  NBMIN = (/1, 1, 5, 1, 3/)
  NX = (/1, 1, 3, 10, 10/)
  WORKSP = (/100, 100, 100, 50, 100/)

  DENSITY = (/0.1, 0.2, 0.5, 0.75, 1.0/)

  DO IPARAM = 1, 5
     CALL SLAENV( 1, NB(IPARAM) )
     CALL SLAENV( 2, NBMIN(IPARAM) )
     CALL SLAENV( 3, NX(IPARAM) )

     DO IDENS = 1 ,5
        DO N = 0, MAX_SIZE

           ALLOCATE(A(N,N), STAT = ISTAT)

           IF( ISTAT /= 0 ) THEN
              WRITE (*,*) "Ran out of memory in test!", N, ISTAT
              STOP
           END IF

           CALL make_skew_mat(A, DENSITY(IDENS))

           !First check the LTL decomposition
           CALL check_LTL_F77("U", A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F77", "SKTRF", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_LTL_F77("L", A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F77","SKTRF", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_LTL_F95("U", A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F95", "SKTRF", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_LTL_F95("L", A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F95","SKTRF", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           !Then check the Householder tridiagonalization

           CALL check_QTQ_F77("U", A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F77", "SKTRD", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_QTQ_F77("L", A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F77", "SKTRD", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_QTQ_F95("U", A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F95", "SKTRD", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_QTQ_F95("L", A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F95", "SKTRD", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           DEALLOCATE(A)
        END DO
     END DO
  END DO

  !Check the unblocked versions SKTF2, SKTD2
  DO IDENS = 1 ,5
     DO N = 0, MAX_SIZE

        ALLOCATE(A(N,N), STAT = ISTAT)

        IF( ISTAT /= 0 ) THEN
           WRITE (*,*) "Ran out of memory in test!", N, ISTAT
           STOP
        END IF

        CALL make_skew_mat(A, DENSITY(IDENS))

        !First check the LTL decomposition
        CALL check_LTL_noblock_F77("U", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTF2", RES)
        END IF

        CALL check_LTL_noblock_F77("L", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTF2", RES)
        END IF

        CALL check_LTL_noblock_F95("U", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTF2", RES)
        END IF

        CALL check_LTL_noblock_F95("L", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTF2", RES)
        END IF

        !Then check the Householder tridiagonalization

        CALL check_QTQ_noblock_F77("U", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTD2", RES)
        END IF

        CALL check_QTQ_noblock_F77("L", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTD2", RES)
        END IF

        CALL check_QTQ_noblock_F95("U", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTD2", RES)
        END IF

        CALL check_QTQ_noblock_F95("L", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTD2", RES)
        END IF

        DEALLOCATE(A)
     END DO
  END DO

  DO IDENS = 1, 5
     DO N=0, MAX_SIZE
        DO K=0, N+2 !beyond the reasonable band width, but should still work!
           ALLOCATE(A(N,N), Au(K+1, N), Ad(K+1, N), STAT=ISTAT)

           IF( ISTAT /= 0 ) THEN
              WRITE (*,*) "Ran out of memory in test!"
              STOP
           END IF

           CALL make_skew_mat_banded(Au, Ad, A, DENSITY(IDENS))

           CALL check_QTQ_band_F77("U", Au, A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_unblocked("F77", "SKBTRD", RES)
           END IF

           CALL check_QTQ_band_F77("L", Ad, A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_unblocked("F77", "SKBTRD", RES)
           END IF

           CALL check_QTQ_band_F95("U", Au, A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_unblocked("F95", "SKBTRD", RES)
           END IF

           CALL check_QTQ_band_F95("L", Ad, A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_unblocked("F95", "SKBTRD", RES)
           END IF

           DEALLOCATE(A, Ad, Au)
        END DO
     END DO
  END DO

END SUBROUTINE test_decomp_s

SUBROUTINE test_decomp_d()
  USE check_decomp
  USE matrix_tools
  IMPLICIT NONE
  INTEGER, PARAMETER :: MAX_SIZE = 20
  INTEGER, DIMENSION(5) :: NB, NBMIN, NX, WORKSP
  INTEGER :: N, ISTAT, IPARAM, K, IDENS
  REAL(KIND(1.0D0)), ALLOCATABLE :: A(:,:), Au(:,:), Ad(:,:)
  REAL(KIND(1.0D0)) :: RES, DENSITY(5)

  !First check the dense decomposition interfaces with blocking
  !This includes different blocking parameters

  NB = (/1, 2, 5, 8, 25/)
  NBMIN = (/1, 1, 5, 1, 3/)
  NX = (/1, 1, 3, 10, 10/)
  WORKSP = (/100, 100, 100, 50, 100/)

  DENSITY = (/0.1, 0.2, 0.5, 0.75, 1.0/)

  DO IPARAM = 1, 5
     CALL SLAENV( 1, NB(IPARAM) )
     CALL SLAENV( 2, NBMIN(IPARAM) )
     CALL SLAENV( 3, NX(IPARAM) )

     DO IDENS = 1, 5
        DO N = 0, MAX_SIZE

           ALLOCATE(A(N,N), STAT = ISTAT)

           IF( ISTAT /= 0 ) THEN
              WRITE (*,*) "Ran out of memory in test!"
              STOP
           END IF

           CALL make_skew_mat(A, DENSITY(IDENS))

           !First check the LTL decomposition

           CALL check_LTL_F77("U", A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F77", "SKTRF", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_LTL_F77("L", A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F77","SKTRF", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_LTL_F95("U", A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F95", "SKTRF", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_LTL_F95("L", A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F95","SKTRF", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           !Then check the Householder tridiagonalization

           CALL check_QTQ_F77("U", A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F77", "SKTRD", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_QTQ_F77("L", A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F77", "SKTRD", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_QTQ_F95("U", A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F95", "SKTRD", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_QTQ_F95("L", A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F95", "SKTRD", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           DEALLOCATE(A)
        END DO
     END DO
  END DO

  !Check the unblocked versions SKTF2, SKTD2
  DO IDENS = 1 ,5
     DO N = 0, MAX_SIZE

        ALLOCATE(A(N,N), STAT = ISTAT)

        IF( ISTAT /= 0 ) THEN
           WRITE (*,*) "Ran out of memory in test!", N, ISTAT
           STOP
        END IF

        CALL make_skew_mat(A, DENSITY(IDENS))

        !First check the LTL decomposition
        CALL check_LTL_noblock_F77("U", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTF2", RES)
        END IF

        CALL check_LTL_noblock_F77("L", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTF2", RES)
        END IF

        CALL check_LTL_noblock_F95("U", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTF2", RES)
        END IF

        CALL check_LTL_noblock_F95("L", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTF2", RES)
        END IF

        !Then check the Householder tridiagonalization

        CALL check_QTQ_noblock_F77("U", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTD2", RES)
        END IF

        CALL check_QTQ_noblock_F77("L", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTD2", RES)
        END IF

        CALL check_QTQ_noblock_F95("U", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTD2", RES)
        END IF

        CALL check_QTQ_noblock_F95("L", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTD2", RES)
        END IF

        DEALLOCATE(A)
     END DO
  END DO

  !Now check the band decompositions
  !(no blocking here)

  DO IDENS = 1, 5
     DO N=0, MAX_SIZE
        DO K=0, N+2 !beyond the reasonable band width, but should still work!
           ALLOCATE(A(N,N), Au(K+1, N), Ad(K+1, N), STAT=ISTAT)

           IF( ISTAT /= 0 ) THEN
              WRITE (*,*) "Ran out of memory in test!"
              STOP
           END IF

           CALL make_skew_mat_banded(Au, Ad, A, DENSITY(IDENS))

           CALL check_QTQ_band_F77("U", Au, A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_unblocked("F77", "SKBTRD", RES)
           END IF

           CALL check_QTQ_band_F77("L", Ad, A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_unblocked("F77", "SKBTRD", RES)
           END IF

           CALL check_QTQ_band_F95("U", Au, A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_unblocked("F95", "SKBTRD", RES)
           END IF

           CALL check_QTQ_band_F95("L", Ad, A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_unblocked("F95", "SKBTRD", RES)
           END IF

           DEALLOCATE(A, Ad, Au)
        END DO
     END DO
  END DO

END SUBROUTINE test_decomp_d

SUBROUTINE test_decomp_c()
  USE check_decomp
  USE matrix_tools
  IMPLICIT NONE
  INTEGER, PARAMETER :: MAX_SIZE = 20
  INTEGER, DIMENSION(5) :: NB, NBMIN, NX, WORKSP
  INTEGER :: N, ISTAT, IPARAM, K, IDENS
  COMPLEX(KIND(1.0E0)), ALLOCATABLE :: A(:,:), Au(:,:), Ad(:,:)
  REAL(KIND(1.0E0)) :: RES, DENSITY(5)

  !First check the dense decomposition interfaces with blocking
  !This includes different blocking parameters

  NB = (/1, 2, 5, 8, 25/)
  NBMIN = (/1, 1, 5, 1, 3/)
  NX = (/1, 1, 3, 10, 10/)
  WORKSP = (/100, 100, 100, 50, 100/)

  DENSITY = (/0.1, 0.2, 0.5, 0.75, 1.0/)

  DO IPARAM = 1, 5
     CALL SLAENV( 1, NB(IPARAM) )
     CALL SLAENV( 2, NBMIN(IPARAM) )
     CALL SLAENV( 3, NX(IPARAM) )

     DO IDENS = 1, 5
        DO N = 0, MAX_SIZE

           ALLOCATE(A(N,N), STAT = ISTAT)

           IF( ISTAT /= 0 ) THEN
              WRITE (*,*) "Ran out of memory in test!"
              STOP
           END IF

           CALL make_skew_mat(A, DENSITY(IDENS))

           !First check the LTL decomposition

           CALL check_LTL_F77("U", A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F77", "SKTRF", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_LTL_F77("L", A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F77","SKTRF", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_LTL_F95("U", A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F95", "SKTRF", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_LTL_F95("L", A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F95","SKTRF", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           !Then check the Householder tridiagonalization

           CALL check_QTQ_F77("U", A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F77", "SKTRD", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_QTQ_F77("L", A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F77", "SKTRD", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_QTQ_F95("U", A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F95", "SKTRD", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_QTQ_F95("L", A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F95", "SKTRD", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           DEALLOCATE(A)
        END DO
     END DO
  END DO

  !Check the unblocked versions SKTF2, SKTD2
  DO IDENS = 1 ,5
     DO N = 0, MAX_SIZE

        ALLOCATE(A(N,N), STAT = ISTAT)

        IF( ISTAT /= 0 ) THEN
           WRITE (*,*) "Ran out of memory in test!", N, ISTAT
           STOP
        END IF

        CALL make_skew_mat(A, DENSITY(IDENS))

        !First check the LTL decomposition
        CALL check_LTL_noblock_F77("U", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTF2", RES)
        END IF

        CALL check_LTL_noblock_F77("L", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTF2", RES)
        END IF

        CALL check_LTL_noblock_F95("U", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTF2", RES)
        END IF

        CALL check_LTL_noblock_F95("L", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTF2", RES)
        END IF

        !Then check the Householder tridiagonalization

        CALL check_QTQ_noblock_F77("U", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTD2", RES)
        END IF

        CALL check_QTQ_noblock_F77("L", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTD2", RES)
        END IF

        CALL check_QTQ_noblock_F95("U", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTD2", RES)
        END IF

        CALL check_QTQ_noblock_F95("L", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTD2", RES)
        END IF

        DEALLOCATE(A)
     END DO
  END DO

  !Now check the band decompositions
  !(no blocking here)

  DO IDENS = 1, 5
     DO N=0, MAX_SIZE
        DO K=0, N+2 !beyond the reasonable band width, but should still work!
           ALLOCATE(A(N,N), Au(K+1, N), Ad(K+1, N), STAT=ISTAT)

           IF( ISTAT /= 0 ) THEN
              WRITE (*,*) "Ran out of memory in test!"
              STOP
           END IF

           CALL make_skew_mat_banded(Au, Ad, A, DENSITY(IDENS))

           CALL check_QTQ_band_F77("U", Au, A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_unblocked("F77", "SKBTRD", RES)
           END IF

           CALL check_QTQ_band_F77("L", Ad, A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_unblocked("F77", "SKBTRD", RES)
           END IF

           CALL check_QTQ_band_F95("U", Au, A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_unblocked("F95", "SKBTRD", RES)
           END IF

           CALL check_QTQ_band_F95("L", Ad, A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_unblocked("F95", "SKBTRD", RES)
           END IF

           DEALLOCATE(A, Ad, Au)
        END DO
     END DO
  END DO
END SUBROUTINE test_decomp_c

SUBROUTINE test_decomp_z()
  USE check_decomp
  USE matrix_tools
  IMPLICIT NONE
  INTEGER, PARAMETER :: MAX_SIZE = 20
  INTEGER, DIMENSION(5) :: NB, NBMIN, NX, WORKSP
  INTEGER :: N, ISTAT, IPARAM, K, IDENS
  COMPLEX(KIND(1.0D0)), ALLOCATABLE :: A(:,:), Au(:,:), Ad(:,:)
  REAL(KIND(1.0D0)) :: RES, DENSITY(5)

  !First check the dense decomposition interfaces with blocking
  !This includes different blocking parameters

  NB = (/1, 2, 5, 8, 25/)
  NBMIN = (/1, 1, 5, 1, 3/)
  NX = (/1, 1, 3, 10, 10/)
  WORKSP = (/100, 100, 100, 50, 100/)

  DENSITY = (/0.1, 0.2, 0.5, 0.75, 1.0/)

  DO IPARAM = 1, 5
     CALL SLAENV( 1, NB(IPARAM) )
     CALL SLAENV( 2, NBMIN(IPARAM) )
     CALL SLAENV( 3, NX(IPARAM) )

     DO IDENS = 1, 5
        DO N = 0, MAX_SIZE

           ALLOCATE(A(N,N), STAT = ISTAT)

           IF( ISTAT /= 0 ) THEN
              WRITE (*,*) "Ran out of memory in test!"
              STOP
           END IF

           CALL make_skew_mat(A, DENSITY(IDENS))

           !First check the LTL decomposition

           CALL check_LTL_F77("U", A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F77", "SKTRF", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_LTL_F77("L", A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F77","SKTRF", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_LTL_F95("U", A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F95", "SKTRF", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_LTL_F95("L", A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F95","SKTRF", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           !Then check the Householder tridiagonalization

           CALL check_QTQ_F77("U", A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F77", "SKTRD", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_QTQ_F77("L", A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F77", "SKTRD", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_QTQ_F95("U", A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F95", "SKTRD", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_QTQ_F95("L", A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_blocked("F95", "SKTRD", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           DEALLOCATE(A)
        END DO
     END DO
  END DO

  !Check the unblocked versions SKTF2, SKTD2
  DO IDENS = 1 ,5
     DO N = 0, MAX_SIZE

        ALLOCATE(A(N,N), STAT = ISTAT)

        IF( ISTAT /= 0 ) THEN
           WRITE (*,*) "Ran out of memory in test!", N, ISTAT
           STOP
        END IF

        CALL make_skew_mat(A, DENSITY(IDENS))

        !First check the LTL decomposition
        CALL check_LTL_noblock_F77("U", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTF2", RES)
        END IF

        CALL check_LTL_noblock_F77("L", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTF2", RES)
        END IF

        CALL check_LTL_noblock_F95("U", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTF2", RES)
        END IF

        CALL check_LTL_noblock_F95("L", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTF2", RES)
        END IF

        !Then check the Householder tridiagonalization

        CALL check_QTQ_noblock_F77("U", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTD2", RES)
        END IF

        CALL check_QTQ_noblock_F77("L", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTD2", RES)
        END IF

        CALL check_QTQ_noblock_F95("U", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTD2", RES)
        END IF

        CALL check_QTQ_noblock_F95("L", A, RES)

        IF( RES > 20.0 ) THEN
           CALL error_unblocked("F95", "SKTD2", RES)
        END IF

        DEALLOCATE(A)
     END DO
  END DO

  !Now check the band decompositions
  !(no blocking here)

  DO IDENS = 1, 5
     DO N=0, MAX_SIZE
        DO K=0, N+2 !beyond the reasonable band width, but should still work!
           ALLOCATE(A(N,N), Au(K+1, N), Ad(K+1, N), STAT=ISTAT)

           IF( ISTAT /= 0 ) THEN
              WRITE (*,*) "Ran out of memory in test!"
              STOP
           END IF

           CALL make_skew_mat_banded(Au, Ad, A, DENSITY(IDENS))

           CALL check_QTQ_band_F77("U", Au, A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_unblocked("F77", "SKBTRD", RES)
           END IF

           CALL check_QTQ_band_F77("L", Ad, A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_unblocked("F77", "SKBTRD", RES)
           END IF

           CALL check_QTQ_band_F95("U", Au, A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_unblocked("F95", "SKBTRD", RES)
           END IF

           CALL check_QTQ_band_F95("L", Ad, A, RES)

           IF( RES > 20.0 ) THEN
              CALL error_unblocked("F95", "SKBTRD", RES)
           END IF

           DEALLOCATE(A, Ad, Au)
        END DO
     END DO
  END DO
END SUBROUTINE test_decomp_z

SUBROUTINE error_blocked(INTF, ROUT, RES, NB, NBMIN, NX, REALWS)
  CHARACTER(LEN=*) :: INTF, ROUT
  REAL(KIND(1.0D0)) :: RES
  INTEGER :: NB, NBMIN, NX, REALWS

  WRITE (*,*) "Error in the ",INTF, " interface of routine ",ROUT
  WRITE (*,*) "Residual is  ", RES
  WRITE (*,*) "NB=",NB, "NBMIN=", NBMIN, "NX=", NX, REALWS, "% of workspace"

  STOP
END SUBROUTINE error_blocked

SUBROUTINE error_unblocked(INTF, ROUT, RES)
  CHARACTER(LEN=*) :: INTF, ROUT
  REAL(KIND(1.0D0)) :: RES

  WRITE (*,*) "Error in the ",INTF, " interface of routine ",ROUT
  WRITE (*,*) "Residual is  ", RES

  STOP
END SUBROUTINE error_unblocked
