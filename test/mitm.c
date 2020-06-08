#define _XOPEN_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "mont.h"
#include "constants.h"
#include "two_walks.h"

int main() {
  srand48(1);
  
  // Only implementing for p=3 mod 4, for the moment
  assert(class_mod_4 == 3);

  two_walk phi;
  phi.A = (proj){ fp2_0, fp2_1 };
  proj B, P4, tmp;

  for (int i = 0; i < 100; i++) {
  phi.A = (proj){ fp2_0, fp2_1 };
    phi.len = 4 + i % 7;

    // sample point on curve
    while (true) {
      fp2_random(&phi.ker.x); fp2_random(&phi.ker.z);
      if (!is_on_curve(&phi.ker, &phi.A))
	continue;
      // multiply by cofactor
      xMUL(&phi.ker, &phi.A, &phi.ker, &p_even_cofactor);
      for (int i = 0; i < two_tors_height - phi.len - 2; i++)
	xDBL(&phi.ker, &phi.A, &phi.ker);
      P4 = phi.ker;
      xDBL(&phi.ker, &phi.A, &phi.ker);
      xDBL(&phi.ker, &phi.A, &phi.ker);
      // check order
      tmp = phi.ker;
      for (int i = 1; i < phi.len; i++)
	xDBL(&tmp, &phi.A, &tmp);
      if (!fp2_iszero(&tmp.x) && !mont_iszero(&tmp))
	break;
    }
    
    // Evaluate isogeny
    eval_walk(&phi, &B, &P4);

    // Change orientation of B
    fp2 a, b;
    proj j1, j2, inf;
    jinv256(&j1, &B);
    xDBL(&tmp, &B, &P4);
    xDBL(&inf, &B, &tmp);
    assert(mont_iszero(&inf));
    fp2_add3(&a, &tmp.x, &tmp.x);
    fp2_add2(&a, &tmp.x);
    fp2_mul2(&a, &B.z);
    fp2_mul2(&B.x, &tmp.z);
    fp2_add2(&B.x, &a);
    fp2_mul2(&B.x, &P4.z);
    fp2_mul3(&a, &P4.x, &tmp.z);
    fp2_mul3(&b, &P4.z, &tmp.x);
    fp2_sub2(&a, &b);
    fp2_mul2(&B.z, &a);
    jinv256(&j2, &B);
    assert(mont_equal(&j1, &j2));

    // MITM
    assert(MITM(&phi, &phi.A, &B, phi.len));
    // check
    eval_walk(&phi, &B, &P4);
    jinv256(&j2, &B);
    assert(mont_equal(&j1, &j2));
    
    phi.A = B;
  }
  
  printf("    \033[1;32mAll tests passed\033[0m\n");
  exit(0);
}
