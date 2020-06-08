#ifndef FP_H
#define FP_H

#include "uintbig.h"

/* fp is in the Montgomery domain, so interpreting that
   as an integer should never make sense.
   enable compiler warnings when mixing up uintbig and fp. */
typedef struct fp {
    uintbig x;
} fp;

extern const uintbig p;
extern const fp fp_0;
extern const fp fp_1;

void fp_set(fp *x, uint64_t y);
void fp_cswap(fp *x, fp *y, bool c);

static inline bool fp_iszero(const fp *a) { return uintbig_iszero(&a->x); }

void fp_enc(fp *x, uintbig const *y); /* encode to Montgomery representation */
void fp_dec(uintbig *x, fp const *y); /* decode from Montgomery representation */

void fp_add2(fp *x, fp const *y);
void fp_sub2(fp *x, fp const *y);
void fp_mul2(fp *x, fp const *y);

void fp_add3(fp *x, fp const *y, fp const *z);
void fp_sub3(fp *x, fp const *y, fp const *z);
void fp_mul3(fp *x, fp const *y, fp const *z);

void fp_sq1(fp *x);
void fp_sq2(fp *x, fp const *y);
void fp_inv(fp *x);
bool fp_issquare(fp *x);
void fp_sqrt(fp *x);

void fp_random(fp *x);

extern long long fp_mul_count;

static inline void fp_neg1(fp *x)
{
  fp_sub3(x,&fp_0,x);
}

static inline void fp_neg2(fp *x,const fp *y)
{
  fp_sub3(x,&fp_0,y);
}

#endif
