#define _XOPEN_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <pari/pari.h>

#include "mont.h"
#include "constants.h"


int main() {
  srand48(1);
  
  // Only implementing for p=3 mod 4, for the moment
  assert(class_mod_4 == 3);

  proj A = { fp2_0, fp2_1 };

  proj pt;
  proj2 P, Q,R,S;
  uintbig k;
  for (int i = 0; i < 200; i++) {
    do {
      fp2_random(&pt.x); fp2_random(&pt.z);
    } while (!is_on_curve(&pt,&A));
    P.x = pt.x;
    P.z = pt.z;
    xLIFT(&P.y, &A, &pt);

    assert(xy_is_on_curve(&A, &P));

    xyNEG(&Q, &P);
    xyADD(&R, &A, &P, &Q);

    assert(xy_is_on_curve(&A, &R));
    assert(xy_is_zero(&R));

    xyDBL(&R, &A, &P);
    xyDBL(&S, &A, &R);
    xyADD(&R, &A, &P, &Q);

    assert(xy_is_on_curve(&A, &R));
    assert(xy_is_zero(&R));

    do {
      fp2_random(&pt.x); fp2_random(&pt.z);
    } while (!is_on_curve(&pt,&A));
    Q.x = pt.x;
    Q.z = pt.z;
    xLIFT(&Q.y, &A, &pt);

    assert(xy_is_on_curve(&A, &Q));

    xyADD(&R, &A, &P, &Q);


    xyDBL(&R, &A, &P);
    xyADD(&R, &A, &R, &P);
    xyADD(&R, &A, &R, &P);
    xyADD(&R, &A, &R, &P);
    xyADD(&R, &A, &R, &P);
    xyADD(&R, &A, &R, &P);


    uintbig_set(&k, 7);
    xyMUL(&S, &A, &P, &k);

    assert(xy_equal(&R,&S));


  }
  
  printf("    \033[1;32mAll tests passed\033[0m\n");
  exit(0);
}
