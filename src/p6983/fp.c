#include "fp.h"
#include "rng.h"
#include <gmp.h>
#include <pari/pari.h>

#define FP_LIMBS (4 * 64 / GMP_LIMB_BITS)

const fp fp_0 = {0, 0, 0, 0};
const fp fp_1 = {1, 0, 0, 0};
const uintbig p = { 0x50ca4291ffffffff, 0xd1b004f94a5952c9,
		    0xb25c76a437728f3b, 0xa3091565b678a990 };
const uintbig mp = { 0xaf35bd6e00000001, 0x2e4ffb06b5a6ad36, 
		     0x4da3895bc88d70c4, 0x5cf6ea9a4987566f };

void fp_set(fp *x, uint64_t y) {
  x->x.c[0] = y;
  x->x.c[1] = x->x.c[2] = x->x.c[3] = 0;
}

void fp_cswap(fp *x, fp *y, bool c) {
  uint64_t tmp;
  for (int i = 0; i < 4*c; i++) {
    tmp = y->x.c[i];
    y->x.c[i] = x->x.c[i];
    x->x.c[i] = tmp;
  }
}

void fp_enc(fp *x, uintbig const *y) {
  x->x.c[0] = y->c[0]; x->x.c[1] = y->c[1];
  x->x.c[2] = y->c[2]; x->x.c[3] = y->c[3];
}
void fp_dec(uintbig *x, fp const *y) {
  x->c[0] = y->x.c[0]; x->c[1] = y->x.c[1];
  x->c[2] = y->x.c[2]; x->c[3] = y->x.c[3];
}

void fp_add2(fp *x, fp const *y) { fp_add3(x, x, y); }
void fp_sub2(fp *x, fp const *y) { fp_sub3(x, x, y); }
void fp_mul2(fp *x, fp const *y) { fp_mul3(x, x, y); }

void fp_add3(fp *x, fp const *y, fp const *z) {
  mp_limb_t carry = mpn_add_n(x->x.c, y->x.c, z->x.c, FP_LIMBS);
  if (carry) {
    mpn_add_n(x->x.c, x->x.c, mp.c, FP_LIMBS);
  } else if (x->x.c[3] > p.c[3] ||
	     (x->x.c[3] == p.c[3] && x->x.c[2] > p.c[2]) ||
	     (x->x.c[3] == p.c[3] && x->x.c[2] == p.c[2] && x->x.c[1] > p.c[1]) ||
	     (x->x.c[3] == p.c[3] && x->x.c[2] == p.c[2] && x->x.c[1] == p.c[1] && x->x.c[0] >= p.c[0])
	     ) {
    mpn_sub_n(x->x.c, x->x.c, p.c, FP_LIMBS);
  }
}

void fp_sub3(fp *x, fp const *y, fp const *z) {
  mp_limb_t borrow = mpn_sub_n(x->x.c, y->x.c, z->x.c, FP_LIMBS);
  if (borrow) {
    mpn_sub_n(x->x.c, x->x.c, mp.c, FP_LIMBS);
  }
}

void fp_mul3(fp *x, fp const *y, fp const *z) {
  uint64_t tmp[8], thrash[5];
  mpn_mul_n(tmp, y->x.c, z->x.c, FP_LIMBS);
  mpn_tdiv_qr(thrash, x->x.c, 0, tmp, 2 * FP_LIMBS, p.c, FP_LIMBS);
}

void fp_sq1(fp *x) { fp_sq2(x, x); }
void fp_sq2(fp *x, fp const *y) { fp_mul3(x, y, y); }
void fp_inv(fp *x) {
  mpz_t res, mpzx, mpzp;
  mpz_init(res);
  mpz_roinit_n(mpzx, x->x.c, FP_LIMBS);
  mpz_roinit_n(mpzp, p.c, FP_LIMBS);
  mpz_invert(res, mpzx, mpzp);
  int i = 0;
  for (; i < res->_mp_size; ++i) {
    x->x.c[i] = ((uint64_t*)res->_mp_d)[i];
  }
  for (; i < 4; ++i) {
    x->x.c[i] = 0;
  }
  mpz_clear(res);
}
bool fp_issquare(fp *x) {
  mpz_t mpzx, mpzp;
  mpz_roinit_n(mpzx, x->x.c, FP_LIMBS);
  mpz_roinit_n(mpzp, p.c, FP_LIMBS);
  int s = mpz_legendre(mpzx, mpzp);
  return s+1;
}

void fp_sqrt(fp *x) {
    mpz_t mpzx, mpzp, mpzsqrt, mpzp14;
    mpz_init(mpzsqrt);
    mpz_init(mpzp14);

    mpz_roinit_n(mpzx, x->x.c, FP_LIMBS);
    mpz_roinit_n(mpzp, p.c, FP_LIMBS);
    mpz_add_ui(mpzp14, mpzp, 1);
    mpz_div_ui(mpzp14, mpzp14, 4);
    mpz_powm(mpzsqrt,mpzx,mpzp14,mpzp);

    int i = 0;
    for (; i < mpzsqrt->_mp_size; ++i) {
        x->x.c[i] = ((uint64_t*)mpzsqrt->_mp_d)[i];
    }
    for (; i < 4; ++i) {
        x->x.c[i] = 0;
    }

    mpz_clear(mpzsqrt);
    mpz_clear(mpzp14);
}

void fp_random(fp *x) {
  uint64_t thrash;
  randombytes(x->x.c + 0, 32);
  mpn_tdiv_qr(&thrash, x->x.c, 0, x->x.c, FP_LIMBS, p.c, FP_LIMBS);
}
