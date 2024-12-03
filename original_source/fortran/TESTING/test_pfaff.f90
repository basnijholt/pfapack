PROGRAM test_pfaff
  IMPLICIT NONE

  WRITE (*,*) "Checking single precision Pfaffian ..."
  CALL test_pfaff_s()
  WRITE (*,*) "ok"

  WRITE (*,*) "Checking double precision Pfaffian ..."
  CALL test_pfaff_d()
  WRITE (*,*) "ok"

  WRITE (*,*) "Checking single complex Pfaffian ..."
  CALL test_pfaff_c()
  WRITE (*,*) "ok"

  WRITE (*,*) "Checking double complex Pfaffian ..."
  CALL test_pfaff_z()
  WRITE (*,*) "ok"

  WRITE (*,*) "***************************************"
  WRITE (*,*) "All Pfaffian tests passed successfully!"
  WRITE (*,*) "***************************************"

END PROGRAM test_pfaff

SUBROUTINE test_pfaff_s()
  USE check_pfaffian
  USE matrix_tools
  USE test_error
  IMPLICIT NONE
  INTEGER, PARAMETER :: MAX_SIZE = 20
  INTEGER, DIMENSION(5) :: NB, NBMIN, NX, WORKSP
  INTEGER :: N, ISTAT, IPARAM, K, IDENS
  REAL(KIND(1.0E0)), ALLOCATABLE :: A(:,:), Au(:,:), Ad(:,:)
  REAL(KIND(1.0E0)) :: RES, DENSITY(5)

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

           !First check the Pfaffian from the LTL decomposition

           CALL check_pfaff_F77("U", 'P', A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPFA",RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F77("L", "P", A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F95("U", 'P', A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F95("L", "P", A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           !Then check the Pfaffian from the Householder decomposition

           CALL check_pfaff_F77("U", "H", A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F77("L", "H", A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F95("U", 'H', A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F95("L", "H", A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           !Now the version avoiding overflow

           !First check the Pfaffian from the LTL decomposition

           CALL check_pfaff10_F77("U", 'P', A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F77("L", "P", A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F95("U", 'P', A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F95("L", "P", A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           !Then check the Pfaffian from the Householder decomposition

           CALL check_pfaff10_F77("U", "H", A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F77("L", "H", A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F95("U", 'H', A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F95("L", "H", A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           DEALLOCATE(A)
        END DO
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

           !First the regular Pfaffian

           CALL check_pfaff_band_F77("U", Au, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F77", "SKBPFA", RES)
           END IF

           CALL check_pfaff_band_F77("L", Ad, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F77", "SKBPFA", RES)
           END IF

           CALL check_pfaff_band_F95("U", Au, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F95", "SKBPFA", RES)
           END IF

           CALL check_pfaff_band_F95("L", Ad, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F95", "SKBPFA", RES)
           END IF

           !Then the version avoiding overflow

           CALL check_pfaff10_band_F77("U", Au, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F77", "SKBPF10", RES)
           END IF

           CALL check_pfaff10_band_F77("L", Ad, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F77", "SKBPF10", RES)
           END IF

           CALL check_pfaff10_band_F95("U", Au, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F95", "SKBPF10", RES)
           END IF

           CALL check_pfaff10_band_F95("L", Ad, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F95", "SKBPF10", RES)
           END IF

           DEALLOCATE(A, Ad, Au)
        END DO
     END DO
  END DO

END SUBROUTINE test_pfaff_s

SUBROUTINE test_pfaff_d()
  USE check_pfaffian
  USE matrix_tools
  USE test_error
  IMPLICIT NONE
  INTEGER, PARAMETER :: MAX_SIZE = 20
  INTEGER, DIMENSION(5) :: NB, NBMIN, NX, WORKSP
  INTEGER :: N, ISTAT, IPARAM, K, IDENS
  REAL(KIND(1.0D0)), ALLOCATABLE :: A(:,:), Au(:,:), Ad(:,:)
  REAL(KIND(1.0D0)) :: RES, DENSITY(5)

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

           !First check the Pfaffian from the LTL decomposition

           CALL check_pfaff_F77("U", 'P', A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPFA",RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F77("L", "P", A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F95("U", 'P', A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F95("L", "P", A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           !Then check the Pfaffian from the Householder decomposition

           CALL check_pfaff_F77("U", "H", A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F77("L", "H", A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F95("U", 'H', A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F95("L", "H", A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           !Now the version avoiding overflow

           !First check the Pfaffian from the LTL decomposition

           CALL check_pfaff10_F77("U", 'P', A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F77("L", "P", A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F95("U", 'P', A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F95("L", "P", A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           !Then check the Pfaffian from the Householder decomposition

           CALL check_pfaff10_F77("U", "H", A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F77("L", "H", A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F95("U", 'H', A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F95("L", "H", A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           DEALLOCATE(A)
        END DO
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

           !First the regular Pfaffian

           CALL check_pfaff_band_F77("U", Au, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F77", "SKBPFA", RES)
           END IF

           CALL check_pfaff_band_F77("L", Ad, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F77", "SKBPFA", RES)
           END IF

           CALL check_pfaff_band_F95("U", Au, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F95", "SKBPFA", RES)
           END IF

           CALL check_pfaff_band_F95("L", Ad, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F95", "SKBPFA", RES)
           END IF

           !Then the version avoiding overflow

           CALL check_pfaff10_band_F77("U", Au, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F77", "SKBPF10", RES)
           END IF

           CALL check_pfaff10_band_F77("L", Ad, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F77", "SKBPF10", RES)
           END IF

           CALL check_pfaff10_band_F95("U", Au, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F95", "SKBPF10", RES)
           END IF

           CALL check_pfaff10_band_F95("L", Ad, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F95", "SKBPF10", RES)
           END IF

           DEALLOCATE(A, Ad, Au)
        END DO
     END DO
  END DO

END SUBROUTINE test_pfaff_d

SUBROUTINE test_pfaff_c()
  USE check_pfaffian
  USE matrix_tools
  USE test_error
  IMPLICIT NONE
  INTEGER, PARAMETER :: MAX_SIZE = 20
  INTEGER, DIMENSION(5) :: NB, NBMIN, NX, WORKSP
  INTEGER :: N, ISTAT, IPARAM, K, IDENS
  COMPLEX(KIND(1.0E0)), ALLOCATABLE :: A(:,:), Au(:,:), Ad(:,:)
  REAL(KIND(1.0E0)) :: RES, DENSITY(5)

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

           !First check the Pfaffian from the LTL decomposition

           CALL check_pfaff_F77("U", 'P', A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPFA",RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F77("L", "P", A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F95("U", 'P', A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F95("L", "P", A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           !Then check the Pfaffian from the Householder decomposition

           CALL check_pfaff_F77("U", "H", A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F77("L", "H", A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F95("U", "H", A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F95("L", "H", A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           !Now the version avoiding overflow

           !First check the Pfaffian from the LTL decomposition

           CALL check_pfaff10_F77("U", 'P', A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F77("L", "P", A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F95("U", 'P', A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F95("L", "P", A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           !Then check the Pfaffian from the Householder decomposition

           CALL check_pfaff10_F77("U", "H", A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F77("L", "H", A, RES, WORKSP(IPARAM)/100.0E0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F95("U", 'H', A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F95("L", "H", A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           DEALLOCATE(A)
        END DO
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

           !First the regular Pfaffian

           CALL check_pfaff_band_F77("U", Au, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F77", "SKBPFA", RES)
           END IF

           CALL check_pfaff_band_F77("L", Ad, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F77", "SKBPFA", RES)
           END IF

           CALL check_pfaff_band_F95("U", Au, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F95", "SKBPFA", RES)
           END IF

           CALL check_pfaff_band_F95("L", Ad, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F95", "SKBPFA", RES)
           END IF

           !Then the version avoiding overflow

           CALL check_pfaff10_band_F77("U", Au, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F77", "SKBPF10", RES)
           END IF

           CALL check_pfaff10_band_F77("L", Ad, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F77", "SKBPF10", RES)
           END IF

           CALL check_pfaff10_band_F95("U", Au, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F95", "SKBPF10", RES)
           END IF

           CALL check_pfaff10_band_F95("L", Ad, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F95", "SKBPF10", RES)
           END IF

           DEALLOCATE(A, Ad, Au)
        END DO
     END DO
  END DO

END SUBROUTINE test_pfaff_c

SUBROUTINE test_pfaff_z()
  USE check_pfaffian
  USE matrix_tools
  USE test_error
  IMPLICIT NONE
  INTEGER, PARAMETER :: MAX_SIZE = 20
  INTEGER, DIMENSION(5) :: NB, NBMIN, NX, WORKSP
  INTEGER :: N, ISTAT, IPARAM, K, IDENS
  COMPLEX(KIND(1.0D0)), ALLOCATABLE :: A(:,:), Au(:,:), Ad(:,:)
  REAL(KIND(1.0D0)) :: RES, DENSITY(5)

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

           !First check the Pfaffian from the LTL decomposition

           CALL check_pfaff_F77("U", 'P', A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPFA",RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F77("L", "P", A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F95("U", 'P', A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F95("L", "P", A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           !Then check the Pfaffian from the Householder decomposition

           CALL check_pfaff_F77("U", "H", A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F77("L", "H", A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F95("U", "H", A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff_F95("L", "H", A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPFA", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           !Now the version avoiding overflow

           !First check the Pfaffian from the LTL decomposition

           CALL check_pfaff10_F77("U", 'P', A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F77("L", "P", A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F95("U", 'P', A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F95("L", "P", A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           !Then check the Pfaffian from the Householder decomposition

           CALL check_pfaff10_F77("U", "H", A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F77("L", "H", A, RES, WORKSP(IPARAM)/100.0D0)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F77", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F95("U", 'H', A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           CALL check_pfaff10_F95("L", "H", A, RES)

           IF( RES > 200.0  ) THEN
              CALL error_blocked("F95", "SKPF10", RES, NB(IPARAM), NBMIN(IPARAM), &
                   NX(IPARAM), WORKSP(IPARAM))
           END IF

           DEALLOCATE(A)
        END DO
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

           !First the regular Pfaffian

           CALL check_pfaff_band_F77("U", Au, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F77", "SKBPFA", RES)
           END IF

           CALL check_pfaff_band_F77("L", Ad, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F77", "SKBPFA", RES)
           END IF

           CALL check_pfaff_band_F95("U", Au, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F95", "SKBPFA", RES)
           END IF

           CALL check_pfaff_band_F95("L", Ad, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F95", "SKBPFA", RES)
           END IF

           !Then the version avoiding overflow

           CALL check_pfaff10_band_F77("U", Au, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F77", "SKBPF10", RES)
           END IF

           CALL check_pfaff10_band_F77("L", Ad, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F77", "SKBPF10", RES)
           END IF

           CALL check_pfaff10_band_F95("U", Au, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F95", "SKBPF10", RES)
           END IF

           CALL check_pfaff10_band_F95("L", Ad, A, RES)

           IF( RES > 2000.0 ) THEN
              CALL error_unblocked("F95", "SKBPF10", RES)
           END IF

           DEALLOCATE(A, Ad, Au)
        END DO
     END DO
  END DO

END SUBROUTINE test_pfaff_z
