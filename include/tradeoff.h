#ifndef TRADEOFF_H
#define TRADEOFF_H

#include <pari/pari.h>
#include "isogenies.h"
#include "two_walks.h"

//Given that F is at most T, the easiest way to include F in the norm of the output is in the strong approximation step
//
//same as klpt_general_power but for the norm of the output must be F 2^e where F is the integer whose factorization matrix is fm. In practice F is a divisor of T.
GEN klpt_general_power_T(GEN I, GEN K, GEN l,GEN fm)


// this is to translate the output of klpt to the signing isogeny, the slight modification is that the norm is a power of 2 time a divisor of T. This can be accomplished by applying ideal_to_isogeny_two on I1,
// and then translating I2 to an isogeny (which is already a suboperation performed in ideal_to_isogeny_two)
//
// T = global_setup.gen_odd_torsion
// I = I1 * I2 is a left O0-ideal of norm dividing T^3 2^e  for some positive integer e where I1 has norm dividing T^2 2^e and is a compatible input for ideal_to_isogeny_two, I_2 has norm dividing T
// J = I1 + O0*T^2
// K is a left O0-ideal equivalent to J of norm a power of 2
// Finds phi such that phi_I = phi_2 * phi_1 * phi_J
// Finds L equivalent to I1 of norm dividing T^2
void ideal_to_isogeny_two_T(two_walk_long *phi_res2,special_isogeny *phi_resT, GEN *L, special_isogeny *phi_L, GEN I, GEN J, GEN K, const special_isogeny *phi_J, const two_walk_long *phi_K);


//Then, we need a way to compress the signing isogeny
//
//compress a 2^e T isogeny to a sequence of integers
void compress_T(uint64_t *zip, const two_walk *walk2, const special_isogeny *walkT, long len);

//Inverse of the above
// A is the starting curve
void decompress_T(special_isogeny *walkT,two_walk *walk2, proj *A, const uint64_t *zip, long len);
