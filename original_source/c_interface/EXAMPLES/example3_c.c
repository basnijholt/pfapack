#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "pfapack.h"

int main()
{
  /* banded real example */
  int N=4;
  int KD=2;
  int KD1=KD+1;
  int info;
  double pfaffian;

  double * Au = (double *)malloc(sizeof(double)* KD1*N);
  double * Ad = (double *)malloc(sizeof(double)* KD1*N);

  /* build up a skewsymmetric band matrix */

  memset(Au, 0, sizeof(double)* KD1*N);
  memset(Ad, 0, sizeof(double)* KD1*N);

  Au[bandupper_fortran(1,3,N,KD)]=-1.0;

  Au[bandupper_fortran(2,4,N,KD)]=-1.0;

  Ad[bandlower_fortran(3,1,N,KD)]=1.0;

  Ad[bandlower_fortran(4,2,N,KD)]=1.0;

  /* Pfaffian from the upper triangle */

  info = skbpfa_d(N, KD, Au, &pfaffian, "U");
  assert(info == 0);

  printf("The pfaffian is %f\n", pfaffian);

  /* Pfaffian from the lower triangle */

  info = skbpfa_d(N, KD, Ad, &pfaffian, "L");
  assert(info == 0);

  printf("The pfaffian is %f\n", pfaffian);

  printf("Those two numbers should be equal and be approx. -1\n");

  return 0;
}
