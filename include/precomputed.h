#ifndef PRECOMPUTED_H
#define PRECOMPUTED_H

#include <pari/pari.h>
#include <stdbool.h>
#include "constants.h"
#include "mont.h"
#include "tedwards.h"


// each basis entry is a triple of the form P,Q,P+Q
extern proj torsion_basis[][3];
extern proj torsion_basis_sum[3];
extern point torsion_basis_ted_sum[3];
extern proj torsion_basis_twist[][3];
extern proj torsion_basis_twist_sum[3];
extern point torsion_basis_twist_ted_sum[3];
extern proj torsion_basis_two[3];


struct precomp_struct {
    // quaternion data

    GEN p; // the prime
    GEN B; // the quaternion algebra
    GEN qf; // the quadratic form defined by the reduced norm with respect to the standard basis
    GEN O0; // the cannonical maximal order
    GEN one;
    GEN i;
    GEN j;
    GEN ji;
    GEN torsion_fm; // factorisation matrix of the available torsion

    GEN O0_b1; // 1
    GEN O0_b2; // i
    GEN O0_b3; // (1-ji)/2
    GEN O0_b4; // (i+j)/2

    GEN O0_to_standard;
    GEN standard_to_O0;


    // elliptic curve data

    proj E0;

    GEN *action_2, *action_3, *action_4;
    GEN *action_twist_2, *action_twist_3, *action_twist_4;
    GEN action_two_2, action_two_3, action_two_4;

    GEN gen_p_plus_fact, gen_p_minus_fact; // factorisation of p+1 and p-1
    GEN gen_p_plus_primary, gen_p_minus_primary; // primary decomposition (list of prime powers)
    GEN gen_odd_torsion;
};


extern struct precomp_struct global_setup;

void init_precomputations();


long ell_to_index(long ell, bool *twist);
long ell_to_e(long ell);

#endif
