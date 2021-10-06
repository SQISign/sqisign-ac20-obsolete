#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "uintbig.h"
#include <assert.h>

extern const long class_mod_4;
extern const long two_tors_height;
// The cofactor of 2^two_tors_height in p±1
extern const uintbig p_even_cofactor;
extern const long security_level;

// the signing isogeny has degree 2^signing_length
extern const long signing_length;
// we have the equality signin_length = two_tors_height * (signing_length_two_tors_height_step -1 ) + last_step_length
extern const long signing_length_two_tors_height_step;
extern const long last_step_length;


// The useful odd factors in p-1
extern const long p_minus_len;
extern const long p_minus_fact[];
extern const long p_minus_mult[];
// The cofactor of the useful odd torsion in p-1
extern const uintbig p_minus_odd_cofactor;

// The useful odd factors in p+1
extern const long p_plus_len;
extern const long p_plus_fact[];
extern const long p_plus_mult[];
// The cofactor of the useful odd torsion in p+1
extern const uintbig p_plus_odd_cofactor;

// the multiplicities to take to obtain log2(p) bits of torsion (for commitment)
extern const long p_minus_mult_com[];
extern const long p_plus_mult_com[];

// the multiplicities to take to obtain log2(p)/2 bits of torsion (for challenge)
extern const long p_minus_mult_cha[];
extern const long p_plus_mult_cha[];

// inverse mapping of p_plus_fact and p_minus_fact
// Warning: unsafe if ell is not in the factors!
static inline long ell_to_index(long ell, bool *twist) {
  *twist = false;
  for (const long *f = p_plus_fact; *f <= ell; f++)
    if (*f == ell)
      return f - p_plus_fact;
  *twist = true;
  for (const long *f = p_minus_fact; *f <= ell; f++)
    if (*f == ell)
      return f - p_minus_fact;
  assert(0);
  return(0);
}
static inline long ell_to_e(long ell) {
  if (ell == 2)
    return two_tors_height;
  bool twist;
  int index = ell_to_index(ell, &twist);
  return (twist ? p_minus_mult : p_plus_mult)[index];
}

#endif
