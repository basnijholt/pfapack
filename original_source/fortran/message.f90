MODULE PFAPACK_MESSAGE
  INTERFACE
     SUBROUTINE MESSAGE( LINFO, NAME, INFO, ISTAT )
       INTEGER, INTENT(IN) :: LINFO
       CHARACTER(LEN = *), INTENT(IN) :: NAME
       INTEGER, OPTIONAL, INTENT(OUT) :: INFO
       INTEGER, OPTIONAL, INTENT(IN) :: ISTAT
     END SUBROUTINE MESSAGE
  END INTERFACE
END MODULE PFAPACK_MESSAGE

SUBROUTINE MESSAGE( LINFO, NAME, INFO, ISTAT )
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: LINFO
  CHARACTER(LEN = *), INTENT(IN) :: NAME
  INTEGER, INTENT(OUT), OPTIONAL :: INFO
  INTEGER, INTENT(IN), OPTIONAL :: ISTAT

  ! special error codes:
  ! LINFO = -100 -> not enough memory
  ! LINFO = -200 -> not optimal workspace

  IF( LINFO == -200 ) THEN
     WRITE (*,*) 'Failed to allocate optimal workspace for ', NAME
     IF( PRESENT(ISTAT) ) THEN
        WRITE(*,*) 'ALLOCATE failed with status ', ISTAT
     END IF
     WRITE (*,*) 'trying minimal workspace (routine may run less efficient)'
     IF( PRESENT(INFO) ) THEN
        INFO = LINFO
     END IF
     RETURN
  ELSE IF( LINFO == -100 ) THEN
     WRITE (*,*) 'Failed to allocate enough memory for ', NAME
     IF( PRESENT(ISTAT) ) THEN
        WRITE(*,*) 'ALLOCATE failed with status ', ISTAT
     END IF
  END IF

  IF( LINFO < 0 .AND. .NOT.PRESENT(INFO) ) THEN
     WRITE (*,*) 'Program failed in subroutine ', NAME, &
          ' with INFO = ', LINFO
     STOP
  END IF

  IF( PRESENT(INFO) ) THEN
     INFO=LINFO
  END IF

END SUBROUTINE MESSAGE
