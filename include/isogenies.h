#ifndef ISOGENIES_H
#define ISOGENIES_H

#include "constants.h"
#include "uintbig.h"
#include "mont.h"

// A packed struct representing the degree of an isogeny as a bitfield
// of valuations corresponding to p_minus_fact or p_plus_fact (see
// constants.c)
typedef struct isog_degree {
  uint64_t val;    // The first byte is used for the valuation of 3 or
                   // 5; the rest is used for all other valuations, 2
                   // bits per factor.

                   // It is obviously limited to 29 factors
} isog_degree;

/* Getters and setters for isog_degree */

static inline uint8_t degree_get(isog_degree deg, int i) {
  return i == 0 ? deg.val & 0xff : (deg.val >> (i*2 + 6)) & 0x3;
}
// Set degree to 1
static inline void degree_one(isog_degree *deg) { deg->val = 0; }
// Set valuation of i-th factor to 0
static inline void degree_unset(isog_degree *deg, int i) {
  deg->val &= ~(i == 0 ? 0xff : (3ll << (i*2 + 6)));
}
// Set valuation of i-th factor. Do call degree_one or degree_unset
// before this!
static inline void degree_set(isog_degree *deg, int i, uint64_t val) {
  deg->val |= val << (i*2 + !!i * 6);
}

/*
  Operations on the degree that depend on parameters.
  
  These expect parameters from constants.c, such as p_minus_fact, etc.
 */

// Return the complement of the degree, i.e. (p ± 1)/deg
static inline isog_degree degree_co(isog_degree deg, const long *mult, long len) {
  isog_degree res = { 0 };
  for (int i = 0; i < len; i++) {
    degree_set(&res, i, mult[i] - degree_get(deg, i));
  }
  return res;
}
static inline void degree_to_uint(uintbig *res, isog_degree deg, const long *fact, long len) {
  *res = uintbig_1;
  for (int i = 0; i < len; i++) {
    uint8_t val = degree_get(deg, i);
    for (int j = 0; j < val; j++)
      uintbig_mul3_64(res, res, fact[i]);
  }
}

// An isogeny of odd degree dividing (p²-1)
typedef struct odd_isogeny {
  proj kernel_plus, kernel_minus;
  isog_degree deg_plus, deg_minus;
} odd_isogeny;

// Evaluate isogeny phi : A -> ?? at point P.
// A is set to the image curve, and P to the image point.
void eval(proj *A, const odd_isogeny *phi, proj *P);

void eval_mult(proj *A, const odd_isogeny *phi, proj *P, int n);


// Compute the dual of phi : A -> ??.
// phi is set to the dual and A is set to the new domain.
void dual(proj *A, odd_isogeny *phi);



#endif
