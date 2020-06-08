#define _XOPEN_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "mont.h"
#include "constants.h"
#include "isogenies.h"
#include "toolbox.h"

static isog_degree mont_order(const proj *P, const proj *A, const long *fact, const long *mult, long len) {
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


  float accumulated_time_ms = 0.;
  int repetitions = 50;
  clock_t t;
  
  proj A = { fp2_0, fp2_1 };

  odd_isogeny phi;
  proj P[2];

    
    for (int i = 0; i < repetitions; i++) {
      // sample point on curve
      do {
       fp2_random(&phi.kernel_plus.x); fp2_random(&phi.kernel_plus.z);
      } while (!is_on_curve(&phi.kernel_plus, &A));
      // multiply by cofactor
      xMUL(&phi.kernel_plus, &A, &phi.kernel_plus, &p_plus_odd_cofactor);
      // computer order
      phi.deg_plus = mont_order(&phi.kernel_plus, &A, p_plus_fact, p_plus_mult, p_plus_len);
    

      do {
       fp2_random(&phi.kernel_minus.x); fp2_random(&phi.kernel_minus.z);
      } while (is_on_curve(&phi.kernel_minus, &A));
      // multiply by cofactor
      xMUL(&phi.kernel_minus, &A, &phi.kernel_minus, &p_minus_odd_cofactor);
      // computer order
      phi.deg_minus = mont_order(&phi.kernel_minus, &A, p_minus_fact, p_minus_mult, p_minus_len);
    


      // Sample point
     fp2_random(&P[0].x); fp2_random(&P[0].z);
     fp2_random(&P[1].x); fp2_random(&P[1].z);


      t = tic();
      eval_mult(&A, &phi, P, 2);
      //dual(&A, &phi);
      //fp2_sqrt(&P[0].x);
      accumulated_time_ms += toc(t);
    }
    printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));
    accumulated_time_ms = 0.;




  
  printf("    \033[1;32mAll tests passed\033[0m\n");
  exit(0);
}
