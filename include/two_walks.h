#ifndef TWO_WALKS_H
#define TWO_WALKS_H

#include "mont.h"
#include "isomorphism.h"
#include "constants.h"

// An isogeny of degree 2^len.
//
// A is the domain curve
// ker is the kernel of order 2^len
//
// Condition: [2^(len-1)]ker must not be the (0,0) point.
typedef struct two_walk {
  proj A;
  proj ker;
  long len;
} two_walk;

// Evaluate 2-isogeny of kernel K at P
void two_isog(const proj *K, proj *P);
// The dual of the above
void two_isog_dual(const proj *K, proj *P);

// Evaluate 2-isogeny walk phi : A -> B at point P.
// B is set to the image curve, and P to the image point.
void eval_walk(const two_walk *phi, proj *B, proj *P);

// same but P is an array of length cardinality
void eval_walk_mult(const two_walk *phi, proj *B, proj *P, long cardinality);

// Internal function. You probably don't want to call this.
void eval_walk_rec(proj *A, proj *K, long len, bool advance, proj *P, long stacklen);

// Compute the dual walk to phi : A -> ?? ≃ B.
//
// P ∈ B is first converted to an element of the image of phi using an
// isomorphism, and then pushed through the dual
void eval_dual(const two_walk *phi, const proj *B, proj *P);

// Compute the dual of phi
void dual_walk(two_walk* phi);

// finds isom such that phi*(isom inverse) can be evaluated
// set phi_new phi*(isom inverse)
// then set P to phi_new(isom(P)), and A to the target curve
void eval_walk_isom(isomorphism *isom, two_walk *phi_new, proj *B, proj *R, const two_walk *phi, const proj *P);

void eval_walk_isom_mult(isomorphism *isom, two_walk *phi_new, proj *B, const two_walk *phi, proj *P, long cardinality);

// Find a walk phi : from -> two of length len by MITM.
// It uses a hash table to accomodate for 2^tab_size entries.
//
// It assumes that the walk from both ends does *not* start with the
// isogeny of kernel (0,0).
//
// returns true if walk is found, false otherwise
//
// Do not call with tab_size = 0
bool MITM_cutoff(two_walk *phi, const proj *from, const proj *to, long len, long tab_size);
static inline bool MITM(two_walk *phi, const proj *from, const proj *to, long len) {
  return MITM_cutoff(phi, from, to, len, len / 2);
}
bool MITM2(two_walk *eta, const proj *from, const proj *to, long length);


#endif
