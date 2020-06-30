#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "uintbig.h"

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

#endif
