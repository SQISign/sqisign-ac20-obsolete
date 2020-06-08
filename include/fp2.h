#ifndef FP2_H
#define FP2_H

#include "fp.h"
#include "rng.h"

typedef struct fp2 {
  fp re, im;
} fp2;

extern const fp2 fp2_0;
#define fp2_1 ((fp2){ fp_1, fp_0 })
#define fp2_i ((fp2){ fp_0, fp_1 })

// Return 3 + 2i
static inline fp2 fp2_non_residue() {
  fp2 res;
  fp_add3(&res.im, &fp_1, &fp_1);
  fp_add3(&res.re, &res.im, &fp_1);
  return res;
}

static inline bool fp2_iszero(const fp2 *a) {
  return fp_iszero(&a->re) && fp_iszero(&a->im);
}

static inline bool fp2_equal(const fp2 *a, const fp2 *b) {
  return (uintbig_equal(&a->re.x,&b->re.x)) && (uintbig_equal(&a->im.x,&b->im.x));
}

void fp2_frob2(fp2 *x, const fp2 *y);
static inline void fp2_frob1(fp2 *x) { fp2_frob2(x, x); }

static inline void fp2_set(fp2 *x, uint64_t y) {
  fp_set(&x->re, y);
  fp_set(&x->im, 0);
}
static inline void fp2_cswap(fp2 *x, fp2 *y, bool c) {
  fp_cswap(&x->re, &y->re, c);
  fp_cswap(&x->im, &y->im, c);
}

static inline void fp2_add2(fp2 *x, fp2 const *y) {
  fp_add2(&x->re, &y->re);
  fp_add2(&x->im, &y->im);
}
static inline void fp2_sub2(fp2 *x, fp2 const *y) {
  fp_sub2(&x->re, &y->re);
  fp_sub2(&x->im, &y->im);
}

static inline void fp2_add3(fp2 *x, fp2 const *y, fp2 const *z) {
  fp_add3(&x->re, &y->re, &z->re);
  fp_add3(&x->im, &y->im, &z->im);
}
static inline void fp2_sub3(fp2 *x, fp2 const *y, fp2 const *z) {
  fp_sub3(&x->re, &y->re, &z->re);
  fp_sub3(&x->im, &y->im, &z->im);
}

void fp2_mul3(fp2 *x, fp2 const *y, fp2 const *z);
static inline void fp2_mul2(fp2 *x, fp2 const *y) { fp2_mul3(x, x, y); }

void fp2_sq2(fp2 *x, fp2 const *y);
static inline void fp2_sq1(fp2 *x) { fp2_sq2(x, x); }

void fp2_inv(fp2 *x);
bool fp2_issquare(fp2 *x);

static inline void fp2_random(fp2 *x) {
  fp_random(&x->re);
  fp_random(&x->im);
}


static inline void fp2_neg1(fp2 *x)
{
  fp_neg1(&x->re);
  fp_neg1(&x->im);
}

static inline void fp2_neg2(fp2 *x, const fp2 *y)
{
  fp_neg2(&x->re, &y->re);
  fp_neg2(&x->im, &y->im);
}

void fp2_exp(fp2 *res, fp2 const *x, uintbig const *k);

void fp2_sqrt(fp2 *x);

bool fp2_dlp_naive(long *res, const fp2 *h, const fp2 *g, long ell);

#endif
