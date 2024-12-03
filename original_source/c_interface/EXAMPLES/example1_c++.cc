#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>

using namespace std;

/* It would be OK to add a

#define CPLUSPLUS_COMPLEX

   here. However, since the code only compiles with a C++ compiler
   anyways, this is chosen automatically when including pfapack.h
*/

#include "pfapack.h"

int main()
{
  /* dense complex example */
  int N=4;

  vector<doublecmplx> A(N*N);

  /* build up a skewsymmetric matrix */
  fill(A.begin(), A.end(), 0.0);

  A[dense_fortran(1,2,N)]=1.0;
  A[dense_fortran(2,1,N)]=-1.0;

  A[dense_fortran(1,3,N)]=2.0;
  A[dense_fortran(3,1,N)]=-2.0;

  A[dense_fortran(1,4,N)]=3.0;
  A[dense_fortran(4,1,N)]=-3.0;

  A[dense_fortran(2,3,N)]=doublecmplx(0.0, 4.0);
  A[dense_fortran(3,2,N)]=-doublecmplx(0.0, 4.0);

  A[dense_fortran(2,4,N)]=5.0;
  A[dense_fortran(4,2,N)]=-5.0;

  A[dense_fortran(3,4,N)]=6.0;
  A[dense_fortran(4,3,N)]=-6.0;

  complex<double> pfaffian;
  int info;

  /* Compute the pfaffian using the lower triangle and the Parlett-Reid
     algorithm */

  info = skpfa(N, &A[0], &pfaffian, "L", "P");
  assert(info == 0);

  cout << "The pfaffian is " << pfaffian << endl;

  /* Compute the pfaffian using the upper triangle (which is untouched)
     and the Householder algorithm */

  info = skpfa(N, &A[0], &pfaffian, "U", "H");
  assert(info == 0);

  cout << "The pfaffian is " << pfaffian << endl;

  cout << "Those two numbers should be equal and be approx. (-4, 12)"
       << endl;

  return 0;
}
