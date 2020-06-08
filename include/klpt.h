
#ifndef KLPT_H
#define KLPT_H

#include <pari/pari.h>

// runs KLPT for the left ideal I in the special order of the quaternion algebra A
// the result is an ideal equivalent to I of norm dividing the integer whose factorisation matrix is fm
// Assumes the basis of A is 1, i, j, j*i, where i^2 = -1 and j^2 = -p
GEN klpt_special_smooth(GEN I, GEN fm);

// same as above, when I is of norm a small power of two (in which case one cannot find an equivalent prime ideal of small norm)
GEN klpt_special_smooth_small_2e_input(GEN I, GEN fm);

GEN klpt_general_power(GEN I, GEN K, GEN l);


#endif




