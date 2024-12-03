MODULE test_error
  INTERFACE error_blocked
     SUBROUTINE error_blocked_s(INTF, ROUT, RES, NB, NBMIN, NX, REALWS)
       CHARACTER(LEN=*) :: INTF, ROUT
       REAL(KIND(1.0E0)) :: RES
       INTEGER :: NB, NBMIN, NX, REALWS
     END SUBROUTINE error_blocked_s

     SUBROUTINE error_blocked_d(INTF, ROUT, RES, NB, NBMIN, NX, REALWS)
       CHARACTER(LEN=*) :: INTF, ROUT
       REAL(KIND(1.0D0)) :: RES
       INTEGER :: NB, NBMIN, NX, REALWS
     END SUBROUTINE error_blocked_d
  END INTERFACE

  INTERFACE error_unblocked
     SUBROUTINE error_unblocked_s(INTF, ROUT, RES)
       CHARACTER(LEN=*) :: INTF, ROUT
       REAL(KIND(1.0E0)) :: RES
     END SUBROUTINE error_unblocked_s

     SUBROUTINE error_unblocked_d(INTF, ROUT, RES)
       CHARACTER(LEN=*) :: INTF, ROUT
       REAL(KIND(1.0D0)) :: RES
     END SUBROUTINE error_unblocked_d
  END INTERFACE
END MODULE test_error

SUBROUTINE error_blocked_s(INTF, ROUT, RES, NB, NBMIN, NX, REALWS)
  CHARACTER(LEN=*) :: INTF, ROUT
  REAL(KIND(1.0E0)) :: RES
  INTEGER :: NB, NBMIN, NX, REALWS

  WRITE (*,*) "Error in the ",INTF, " interface of routine ",ROUT
  WRITE (*,*) "Residual is  ", RES
  WRITE (*,*) "NB=",NB, "NBMIN=", NBMIN, "NX=", NX, REALWS, "% of workspace"

  STOP
END SUBROUTINE error_blocked_s

SUBROUTINE error_blocked_d(INTF, ROUT, RES, NB, NBMIN, NX, REALWS)
  CHARACTER(LEN=*) :: INTF, ROUT
  REAL(KIND(1.0D0)) :: RES
  INTEGER :: NB, NBMIN, NX, REALWS

  WRITE (*,*) "Error in the ",INTF, " interface of routine ",ROUT
  WRITE (*,*) "Residual is  ", RES
  WRITE (*,*) "NB=",NB, "NBMIN=", NBMIN, "NX=", NX, REALWS, "% of workspace"

  STOP
END SUBROUTINE error_blocked_d

SUBROUTINE error_unblocked_s(INTF, ROUT, RES)
  CHARACTER(LEN=*) :: INTF, ROUT
  REAL(KIND(1.0E0)) :: RES

  WRITE (*,*) "Error in the ",INTF, " interface of routine ",ROUT
  WRITE (*,*) "Residual is  ", RES

  STOP
END SUBROUTINE error_unblocked_s

SUBROUTINE error_unblocked_d(INTF, ROUT, RES)
  CHARACTER(LEN=*) :: INTF, ROUT
  REAL(KIND(1.0D0)) :: RES

  WRITE (*,*) "Error in the ",INTF, " interface of routine ",ROUT
  WRITE (*,*) "Residual is  ", RES

  STOP
END SUBROUTINE error_unblocked_d
