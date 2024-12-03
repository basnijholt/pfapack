*     Fortran 77 interface example: 
*     workspace query for program example1.f

      PROGRAM Example1ws
      
      DOUBLE PRECISION A(4,4)
      DOUBLE PRECISION WORK(1)
      DOUBLE PRECISION PFAFF
      INTEGER IWORK(4)
      INTEGER INFO

      WRITE (*,*) "Workspace query for example1.f"

*     Perform workspace query
      CALL DSKPFA('U', 'P', 4, A, 4, PFAFF, IWORK, WORK, -1, INFO)

      WRITE (*,*) "Optimal workspace for example1.f: ", INT(WORK(1))
      WRITE (*,*) "(optimal workspace is individual for each computer)"
      END
