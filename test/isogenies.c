#define _XOPEN_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "mont.h"
#include "constants.h"
#include "isogenies.h"

isog_degree order(const proj *P, const proj *A, const long *fact, const long *mult, long len) {
  proj tmp;
  uintbig cof;
  isog_degree deg = degree_co((isog_degree){ 0 }, mult, len);
  for (int j = 0; j < len; j++) {
    degree_unset(&deg, j);
    degree_to_uint(&cof, deg, fact, len);
    xMUL(&tmp, A, P, &cof);
    uintbig_set(&cof, fact[j]);
    uint8_t v = 0;
    for ( ; !mont_iszero(&tmp); v++) {
      xMUL(&tmp, A, &tmp, &cof);
    }
    degree_set(&deg, j, v);
  }
  return deg;
}

int main() {
  srand48(1);
  // Only implementing for p=3 mod 4, for the moment
  assert(class_mod_4 == 3);

  proj A = { fp2_0, fp2_1 }, Abak;

  odd_isogeny phi;
  proj P, Pbak;
  for (int oncurve = 0; oncurve < 2; oncurve++) {
    proj *kernel, *nkernel;
    isog_degree *deg_k, *deg_nk;
    const uintbig *cofactor;
    const long *fact, *mult;
    long len;

    if (oncurve) {
      kernel = &phi.kernel_plus; nkernel =  &phi.kernel_minus;
      deg_k = &phi.deg_plus; deg_nk = &phi.deg_minus;
      cofactor = &p_plus_odd_cofactor;
      fact = p_plus_fact; mult = p_plus_mult;
      len = p_plus_len;
    } else {
      kernel = &phi.kernel_minus; nkernel =  &phi.kernel_plus;
      deg_k = &phi.deg_minus; deg_nk = &phi.deg_plus;
      cofactor = &p_minus_odd_cofactor;
      fact = p_minus_fact; mult = p_minus_mult;
      len = p_minus_len;
    }

    for (int i = 0; i < 10; i++) {
      // sample point on curve
      do {
	fp2_random(&kernel->x); fp2_random(&kernel->z);
      } while (is_on_curve(kernel, &A) != oncurve);
      // multiply by cofactor
      xMUL(kernel, &A, kernel, cofactor);
      // computer order
      *deg_k = order(kernel, &A, fact, mult, len);

      // trivial kernel on twist
      nkernel->x = fp2_1;
      nkernel->z = fp2_0;
      degree_one(deg_nk);

      // Sample point on twist
      do {
	fp2_random(&P.x); fp2_random(&P.z);
      } while (is_on_curve(&P, &A) == oncurve);

      // Evaluate isogeny
      Abak = A; Pbak = P;
      eval(&A, &phi, &P);
      // Compute and evaluate dual, check image equality
      A = Abak;
      dual(&A, &phi);
      eval(&A, &phi, &P);
      assert(mont_equal(&A, &Abak));
      // Check dual isogeny equation
      uintbig ord;
      degree_to_uint(&ord, *deg_k, fact, len);
      xMUL(&Pbak, &A, &Pbak, &ord);
      assert(mont_equal(&Pbak, &P));
    }
  }

  printf("    \033[1;32mAll tests passed\033[0m\n");
  exit(0);
}
