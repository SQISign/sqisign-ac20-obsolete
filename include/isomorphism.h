#ifndef ISOMORPHISM_H
#define ISOMORPHISM_H

#include "mont.h"
#include <stdio.h>

// The j-invariant divided by 256
void jinv256(proj *j, const proj *A);

// Represents the isomorphism that maps (X:Z) ↦ ( (Nx X - Nz Z) : (D Z) )
typedef struct isomorphism {
  fp2 Nx, Nz, D;
} isomorphism;

// Given curves A and B, computes an isomorphism A -> B
//
// It works for curves j-invariant 0 or 1728, however this is probably
// not the function you're looking for.
void mont_isom(isomorphism *isom, const proj *A, const proj *B);

// Change A to an equivalent A-invariant, and produce associated
// isomorphism
void rand_isom(isomorphism *isom, proj *A);

static inline void trivial_isom(isomorphism *isom) {
	isomorphism id = {fp2_1,fp2_0,fp2_1}; *isom = id;
}

// Apply isomorphism to point P
static inline void mont_isom_apply(const isomorphism *isom, proj *P) {
  fp2 tmp;
  fp2_mul2(&P->x, &isom->Nx);
  fp2_mul3(&tmp, &P->z, &isom->Nz);
  fp2_sub2(&P->x, &tmp);
  fp2_mul2(&P->z, &isom->D);
}


#endif
