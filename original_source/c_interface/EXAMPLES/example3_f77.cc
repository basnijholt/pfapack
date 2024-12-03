#include <iostream>
#include <vector>
#include <cstdlib>

#include "fortran_pfapack.h"

using namespace std;

int main()
{
  //banded real example

  int N=4;
  int KD=2;
  int KD1=KD+1;

  vector<double> A((KD1)*N);

  //build up a skewsymmetric matrix

  fill(A.begin(), A.end(), 0.0);

  //LAPACK uses FORTRAN arrays
  //A(i,j)=A[i-1+(j-1)*N]

  A[bandlower_fortran(3,1,N,KD)]=1.0;

  A[bandlower_fortran(4,2,N,KD)]=1.0;

  //now compute the pfaffian working with the lower triangle

  double pfaffian;
  int info=0;
  vector<double> WORK(3*N-1);

  //now do the real calculation
  PFAPACK_dskbpfa("L", &N, &KD, &A[0], &KD1, &pfaffian, &WORK[0],
		  &info);

  cout << "The pfaffian is " << pfaffian << endl;
  cout << "(result should be approx. -1.0)" << endl;

  return 0;
}
