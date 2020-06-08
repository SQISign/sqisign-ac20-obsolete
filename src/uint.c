#include "uintbig.h"
#include <gmp.h>

#define N_LIMBS (4 * 64 / GMP_LIMB_BITS)

const uintbig uintbig_1 = { 1, 0, 0, 0 };

void uintbig_set(uintbig *x, uint64_t y) {
  x->c[0] = y;
  x->c[1] = x->c[2] = x->c[3] = 0;
}

bool uintbig_bit(uintbig const *x, uint64_t k) {
  return x->c[k / 64] >> (k % 64) & 1;
}

bool uintbig_add3(uintbig *x, uintbig const *y, uintbig const *z) {
  return mpn_add_n(x->c, y->c, z->c, N_LIMBS);
}
bool uintbig_sub3(uintbig *x, uintbig const *y, uintbig const *z) {
  return mpn_sub_n(x->c, y->c, z->c, N_LIMBS);
}

void uintbig_mul3_64(uintbig *x, uintbig const *y, uint64_t z) {
  mpn_mul_1(x->c, y->c, N_LIMBS, z);
}
uint64_t uintbig_div3_64(uintbig *x, uintbig const *y, uint64_t z) {
  return mpn_divmod_1(x->c, y->c, N_LIMBS, z);
}
