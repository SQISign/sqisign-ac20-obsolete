#define _XOPEN_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "mont.h"
#include "constants.h"

int main() {
  srand48(1);

  // Only implementing for p=3 mod 4, for the moment
  assert(class_mod_4 == 3);

  proj A = { fp2_0, fp2_1 };

  proj P, Q, Qk;
  for (int i = 0; i < 200; i++) {
    fp2_random(&P.x); fp2_random(&P.z);
    fp2_random(&Q.x); fp2_random(&Q.z);
    uintbig k, scal, ord;
    uintbig_set(&k, i+3);
    xMUL(&Qk, &A, &Q, &k);



    const long *factors;
    long len;
    if (is_on_curve(&P, &A)) {
      len = p_plus_len;
      factors = p_plus_fact;
      uintbig_add3(&ord, &p, &uintbig_1);
    } else {
      len = p_minus_len;
      factors = p_minus_fact;
      uintbig_sub3(&ord, &p, &uintbig_1);
    }
    long fact = factors[i % len];

    // Cofactor multiplication
    long rem = uintbig_div3_64(&ord, &ord, fact);
    assert(rem == 0);
    xMUL(&P, &A, &P, &ord);

    proj Z;
    uintbig_set(&scal, fact);
    xMUL(&Z, &A, &P, &scal);
    assert(fp2_iszero(&Z.z));

    if (!fp2_iszero(&P.z)) {
      proj AA = A;
      xISOG(&A, &Q, &P, fact);
      xISOG(&AA, &Qk, &P, fact);

      fp2_sub2(&AA.x, &A.x);
      fp2_sub2(&AA.z, &A.z);
      assert(fp2_iszero(&AA.x) && fp2_iszero(&AA.z));

      xMUL(&Q, &A, &Q, &k);
      fp2_mul2(&Q.x, &Qk.z);
      fp2_mul2(&Q.z, &Qk.x);
      fp2_sub2(&Q.x, &Q.z);
      assert(fp2_iszero(&Q.x));
    }
  }

  printf("    \033[1;32mAll tests passed\033[0m\n");
  exit(0);
}
