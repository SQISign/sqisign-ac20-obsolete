#ifndef UINT_H
#define UINT_H

#include <stdbool.h>
#include <stdint.h>

#define BITS 256

// 256 bits unsigned integers
typedef struct uintbig {
    uint64_t c[4];
} uintbig;

extern const uintbig uintbig_1;

void uintbig_set(uintbig *x, uint64_t y);

bool uintbig_bit(uintbig const *x, uint64_t k);

static inline bool uintbig_iszero(const uintbig *a) {
  return !(a->c[0] || a->c[1] || a->c[2] || a->c[3]);
}

static inline bool uintbig_equal(const uintbig *a, const uintbig *b) {
  return ((a->c[0]==b->c[0]) && (a->c[1]==b->c[1]) && (a->c[2]==b->c[2]) && (a->c[3]==b->c[3]));
}

bool uintbig_add3(uintbig *x, uintbig const *y, uintbig const *z); /* returns carry */
bool uintbig_sub3(uintbig *x, uintbig const *y, uintbig const *z); /* returns borrow */

void uintbig_mul3_64(uintbig *x, uintbig const *y, uint64_t z);
uint64_t uintbig_div3_64(uintbig *x, uintbig const *y, uint64_t z); /* returns remainder */

#endif
