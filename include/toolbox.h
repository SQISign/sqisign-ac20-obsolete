
#ifndef TOOLBOX_H
#define TOOLBOX_H

#include <pari/pari.h>
#include <stdbool.h>

clock_t tic();
float tac(); /* time in ms since last tic */
float TAC(const char *str); /* same, but prints it with label 'str' */
float toc(const clock_t t); /* time in ms since t */
float TOC(const clock_t t, const char *str); /* same, but prints it with label 'str' */

// computes the factorisation matrix of x1 by x2 given their factorisation matrices f1 and f2
GEN famat_div(GEN f1, GEN f2);

// returns the first divisor in f1 (first prime in the list with non-zero exponent)
// if f2 is given, it is set to the updated factorisation, where the returned factor has been removed
GEN famat_pop(GEN f1, GEN* f2);

// returns a random divisor of the factorisation matrix f1, with product at least B
// if f2 is given, it is set to f1 divided by the random divisor
GEN famat_random(GEN f1, GEN B);

// returns the product
GEN famat_prod(GEN f);

// returns the n-th prime divisor in f, where primes are counted with multiplicity (the 3rd prime of 2*3^2*5 is 3 because 2,3,3,5)
GEN famat_get_ith(GEN f, GEN n);

int cornacchia_extended(GEN N, GEN *x, GEN *y);

// solve x^2 + y^2 + p(u^2 + v^2) = M, with (u,v) != (0,0)
// when parity != 0, ensures that (x+v) and (y+u) are not both even
// (this means that x + y*i + u*j + v*ji is not is 2*Order(1,i,(1-ji)/2, (i+j)/2)
GEN norm_equation_special(GEN p, GEN M, long parity, bool randomized);

GEN lattice_nearest_plane(GEN lat, GEN target, long flag);

#endif




