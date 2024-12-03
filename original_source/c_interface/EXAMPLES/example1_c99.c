#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define C99_COMPLEX

#include "pfapack.h"

int main()
{
  /* dense complex example */
  int N=4;
  int info;
  double complex pfaffian;

  doublecmplx *A=(doublecmplx *)malloc(sizeof(doublecmplx)*N*N);

  /* build up a skewsymmetric matrix */
  memset(A, 0, sizeof(doublecmplx)*N*N);

  A[dense_fortran(1,2,N)]=1.0;
  A[dense_fortran(2,1,N)]=-1.0;

  A[dense_fortran(1,3,N)]=2.0;
  A[dense_fortran(3,1,N)]=-2.0;

  A[dense_fortran(1,4,N)]=3.0;
  A[dense_fortran(4,1,N)]=-3.0;

  A[dense_fortran(2,3,N)]=4.0*_Complex_I;
  A[dense_fortran(3,2,N)]=-4.0*_Complex_I;

  A[dense_fortran(2,4,N)]=5.0;
  A[dense_fortran(4,2,N)]=-5.0;

  A[dense_fortran(3,4,N)]=6.0;
  A[dense_fortran(4,3,N)]=-6.0;

  /* Compute the pfaffian using the lower triangle and the Parlett-Reid
     algorithm */

  info = skpfa_z(N, A, &pfaffian, "L", "P");
  assert(info == 0);

  printf("The pfaffian is %f + i %f\n", creal(pfaffian), cimag(pfaffian));

  /* Compute the pfaffian using the upper triangle (which is untouched)
     and the Householder algorithm */

  info = skpfa_z(N, A, &pfaffian, "U", "H");
  assert(info == 0);

  printf("The pfaffian is %f + i %f\n", creal(pfaffian), cimag(pfaffian));

  printf("Those two numbers should be equal and be approx. -4 + i 12\n");

  return 0;
}
