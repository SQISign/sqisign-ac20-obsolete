
#ifndef IDISO_H
#define IDISO_H

#include <pari/pari.h>
#include "isogenies.h"
#include "two_walks.h"


typedef struct special_isogeny {
  // phi1 : source -> something isomorphic to middle
  // phi2 : middle -> target
  proj source;
  proj target;

  proj middle;

  odd_isogeny phi1;

  odd_isogeny phi2;
  bool phi2_set;

  odd_isogeny phi2_dual;
  bool phi2_dual_set;

} special_isogeny;


// A chain of two_walk
typedef struct two_walk_long {
    two_walk *phi;
    long len;
} two_walk_long;


void action_from_elle(GEN *m1, GEN *m2, GEN *m3, GEN *m4, long ell, long e);


// Stuff to translate between isogenies and ideals
GEN kernel_to_ideal_gen_action(GEN v, GEN m1, GEN m2, GEN m3, GEN m4, long ell, long e);
GEN endo_to_kernel_action(GEN endo, GEN m1, GEN m2, GEN m3, GEN m4, long ell, long e);

GEN kernel_to_ideal_action_O0(GEN v, GEN m1, GEN m2, GEN m3, GEN m4, long ell, long e);
GEN ideal_to_kernel_action_O0(GEN I, GEN m1, GEN m2, GEN m3, GEN m4, long ell, long e);

GEN kernel_to_ideal_O0_ell(GEN v, long ell);
GEN ideal_to_kernel_O0_ell(GEN I, long ell);

GEN ideal_to_kernel_O0_T(GEN I, GEN factorisation_norm);
GEN kernel_to_ideal_O0_T(GEN coeff);

odd_isogeny ideal_to_isogeny_O0_T(GEN I, GEN factorisation_norm);
two_walk ideal_to_isogeny_O0_two(GEN I);

GEN kernel_to_ideal_gen_O0_ell(GEN v, long ell, long *e);
GEN endo_to_kernel_O0_ell(GEN alpha, long ell, long e);

GEN torsion_crt_decompose (GEN v, bool twist);
GEN torsion_crt_compose (GEN coeff, bool twist);

void famat_to_degree(isog_degree *deg, isog_degree *deg_twist, GEN f);
proj coeff_to_E0(GEN coeff, bool twist);
void gentobig(uintbig *res, GEN a);
void two_walk_stol(two_walk_long *phil, const two_walk *phi);

// Evaluate special isogeny phi : A -> ?? at point P.
// sets P to the image point
void init_trivial_two_walk_long(two_walk_long *phi);
void free_two_walk_long(two_walk_long *phi);
void copy_two_walk_long(two_walk_long *copy, const two_walk_long *phi);
proj eval_special(proj *A, special_isogeny *phi, const proj *P);
void two_walk_composition_ls(two_walk_long *phi, const two_walk_long *phi2, const two_walk *phi1);
void two_walk_composition_sl(two_walk_long *phi, const two_walk *phi2, const two_walk_long *phi1);
void two_walk_composition_ss(two_walk_long *phi, const two_walk *phi2, const two_walk *phi1);
void two_walk_composition_ll(two_walk_long *phi, const two_walk_long *phi2, const two_walk_long *phi1);
odd_isogeny push_odd_isogeny_through_two_walk(const odd_isogeny *phi_odd, proj *phi_odd_source, const two_walk *phi_two);
odd_isogeny push_odd_isogeny_through_two_walk_long(const odd_isogeny *phi_odd, proj *phi_odd_source, const two_walk_long *phi_two);
two_walk push_two_walk_through_odd_isogeny(const two_walk *phi_two, const odd_isogeny *phi_odd, const proj *phi_odd_source);
void push_two_walk_long_through_odd_isogeny(two_walk_long *phi, const two_walk_long *phi_two, const odd_isogeny *phi_odd, const proj *phi_odd_source);
two_walk push_two_walk_through_special_isogeny(const two_walk *phi_two, special_isogeny *phi_special);
//two_walk_long push_two_walk_long_through_special_isogeny(const two_walk_long *phi_two, const special_isogeny *phi_special);


// J is an ideal of norm dividing (global_setup.gen_odd_torsion)^2
// I is an equivalent ideal of norm a power of 2
special_isogeny special_ideal_to_isogeny(GEN J, GEN I, const two_walk_long *phi_I);

// T = global_setup.gen_odd_torsion
// f = two_tors_height
// I is a left O0-ideal of norm dividing T^2 2^{2f+delta}
// J is a left O0-ideal containing I of norm gcd(T^2,n(I))
// K is a left O0-ideal equivalent to J of norm a power of 2
// Finds phi such that phi_I = phi * phi_J
// Finds L equivalent to I of norm dividing T^2
void ideal_to_isogeny_two_2f_delta(two_walk_long *phi, GEN *L, special_isogeny *phi_L, GEN I, GEN J, GEN K, proj *phi_K_basis, proj *phi_K_target, int delta, GEN I_long);

// T = global_setup.gen_odd_torsion
// I is a left O0-ideal of norm dividing T^2 2^e for some positive integer e
// J = I + O0*T^2
// K is a left O0-ideal equivalent to J of norm a power of 2
// Finds phi such that phi_I = phi * phi_J
// Finds L equivalent to I of norm dividing T^2
void ideal_to_isogeny_two(two_walk_long *phi_res, GEN *L, special_isogeny *phi_L, GEN I, GEN J, GEN K, const special_isogeny *phi_J, const two_walk_long *phi_K, bool endpoint_close_to_E0);

void ideal_to_isogeny_O0_two_long(two_walk_long *phi, GEN *L, special_isogeny *phi_L, GEN I, bool endpoint_close_to_E0);

#endif
