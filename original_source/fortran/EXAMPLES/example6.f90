!Example6: Compute the Pfaffian of a banded skew-symmetric matrix
!          Uses the Fortran95 interface

PROGRAM Example6
  USE F95_PFAPACK
  IMPLICIT NONE
  REAL(KIND(1D0)) :: Aupper(3,4), Alower(3,4)
  REAL(KIND(1D0)) :: PFAFF
  INTEGER KD

  WRITE (*,*) "Example program: compute the Pfaffian of two banded 4x4 matrices, using the Fortran95 interface"

!We consider matrices with two off-diagonals
  KD=2

  Aupper=0D0
  Alower=0D0

!Matrix1: sparse banded matrix stored in upper band storage format
!entries: A(1,3)=1, A(2,4)=1, all other are zero

  Aupper(KD+1+1-3, 3)=1D0
  Aupper(KD+1+2-4, 4)=1D0

!Matrix 2: dense banded matrix stored in lower band storage format
!entries: A(2,1)=-1, A(3,1)=-2, A(3,2)=-4, A(4,2)=-5, A(4,3)=-6

  Alower(1+2-1, 1)=-1D0
  Alower(1+3-1, 1)=-2D0
  Alower(1+3-2, 2)=-4D0
  Alower(1+4-2, 2)=-5D0
  Alower(1+4-3, 3)=-6D0

!First matrix 1
  CALL SKBPFA(Aupper, PFAFF, UPLO='U')

  WRITE (*,*) "Pfaffian of example matrix 1: ", PFAFF
  WRITE (*,*) "(The result should be approx. -1.0)"

!Then matrix 2
  CALL SKBPFA(Alower, PFAFF, UPLO = 'L')

  WRITE (*,*) "Pfaffian of example matrix 2: ", PFAFF
  WRITE (*,*) "(The result should be approx. -4.0)"

END PROGRAM Example6
