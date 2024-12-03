#include <iostream>
#include <vector>

#include "fortran_pfapack.h"

using namespace std;

int main()
{
  //dense complex example
  int N=4;

  vector<complex<double> > A(N*N);

  //build up a skewsymmetric matrix

  fill(A.begin(), A.end(), 0.0);

  //LAPACK uses FORTRAN arrays
  //A(i,j)=A[i-1+(j-1)*N]

  A[dense_fortran(1,2,N)]=1.0;
  A[dense_fortran(2,1,N)]=-1.0;

  A[dense_fortran(1,3,N)]=2.0;
  A[dense_fortran(3,1,N)]=-2.0;

  A[dense_fortran(1,4,N)]=3.0;
  A[dense_fortran(4,1,N)]=-3.0;

  A[dense_fortran(2,3,N)]=complex<double>(0, 4.0);
  A[dense_fortran(3,2,N)]=-complex<double>(0, 4.0);

  A[dense_fortran(2,4,N)]=5.0;
  A[dense_fortran(4,2,N)]=-5.0;

  A[dense_fortran(3,4,N)]=6.0;
  A[dense_fortran(4,3,N)]=-6.0;

  //now compute the pfaffian working with the lower triangle

  complex<double> pfaffian;
  int info=0;
  vector<complex<double> > WORK(N);
  vector<double> RWORK(N-1);
  vector<int> IWORK(N);
  int LWORK;

  //workspace query for optimal workspace
  LWORK=-1;
  PFAPACK_zskpfa("L", "P", &N, &A[0], &N, &pfaffian, &IWORK[0],
		 &WORK[0], &LWORK, &RWORK[0], &info);

  LWORK=static_cast<int>(real(WORK[0]));
  WORK.resize(LWORK);

  //now do the real calculation
  PFAPACK_zskpfa("L", "P", &N, &A[0], &N, &pfaffian, &IWORK[0],
		 &WORK[0], &LWORK, &RWORK[0], &info);

  cout << "The pfaffian is " << pfaffian << endl;

  //Now repeat the calculation using the upper triangle
  //(which is still untouched)

  //just for fun, use only the minimal WORK size,
  //the function should still work (though maybe slower)

  WORK.resize(N);
  PFAPACK_zskpfa("U", "P", &N, &A[0], &N, &pfaffian, &IWORK[0],
		 &WORK[0], &N, &RWORK[0], &info);

  cout << "The pfaffian is " << pfaffian << endl;

  cout << "Those two numbers should be the same and equal to approx. (-4, 12)"
       << endl;

  return 0;
}
