#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define STRUCT_COMPLEX

#include "pfapack.h"

int main()
{
  /* dense complex example */
  int N=4;
  int info;
  doublecmplx pfaffian;

  doublecmplx *A=(doublecmplx *)malloc(sizeof(doublecmplx)*N*N);

  /* build up a skewsymmetric matrix */
  memset(A, 0, sizeof(doublecmplx)*N*N);

  A[dense_fortran(1,2,N)].re=1.0;
  A[dense_fortran(2,1,N)].re=-1.0;

  A[dense_fortran(1,3,N)].re=2.0;
  A[dense_fortran(3,1,N)].re=-2.0;

  A[dense_fortran(1,4,N)].re=3.0;
  A[dense_fortran(4,1,N)].re=-3.0;

  A[dense_fortran(2,3,N)].im=4.0;
  A[dense_fortran(3,2,N)].im=-4.0;

  A[dense_fortran(2,4,N)].re=5.0;
  A[dense_fortran(4,2,N)].re=-5.0;

  A[dense_fortran(3,4,N)].re=6.0;
  A[dense_fortran(4,3,N)].re=-6.0;

  /* Compute the pfaffian using the lower triangle and the Parlett-Reid
     algorithm */

  info = skpfa_z(N, A, &pfaffian, "L", "P");
  assert(info == 0);

  printf("The pfaffian is %f + i %f\n", pfaffian.re, pfaffian.im);

  /* Compute the pfaffian using the upper triangle (which is untouched)
     and the Householder algorithm */

  info = skpfa_z(N, A, &pfaffian, "U", "H");
  assert(info == 0);

  printf("The pfaffian is %f + i %f\n", pfaffian.re, pfaffian.im);

  printf("Those two numbers should be equal and be approx. -4 + i 12\n");

  return 0;
}
