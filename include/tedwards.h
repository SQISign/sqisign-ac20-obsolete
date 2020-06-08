#ifndef TEDWARDS_H
#define TEDWARDS_H

#include <pari/pari.h>
#include "uintbig.h"
#include "fp2.h"
#include "mont.h"

// a*x^2+y^2=1+d*x^2*y^2

typedef struct point {
    fp2 x;
    fp2 y;
    fp2 z;
    fp2 t; // t = x*y/z
} point;

extern const point ted_0;

bool ted_is_on_curve(point const *P, proj const *E);
bool ted_equal(point const *P, point const *Q);

void ted_double(point *Q, proj const *E, point const *P);
void ted_add(point *S, proj const *E, point const *P, point const *Q);
void ted_mul(point *res, point const *P, proj const *E, uintbig const *k);

void ted_miller_phi(fp2 *phi, proj const *E, point const *P1, point const *P2, point const *Q, bool dou);
void ted_miller(fp2 *res, fp2 *res2, proj const *E, point const *P, point const *Q, point const *Q2, uintbig const *k);
void ted_weil(fp2 *res, proj const *E, point const *P, point const *Q, uintbig const *k);

void ted_neg(point *Q, point const *P);
bool ted_iszero(point const *P);
void mont_to_ted(proj *E, proj const *A, bool twist);
void mont_to_ted_point(point *Q, proj const *A, proj const *P);
void ted_to_mont_point(proj *Q, point const *P);

bool ted_bidim_log_weil(long *a, long *b, const proj *E, const point *Q, const point *P1, const point *P2, long ell);
bool ted_bidim_log(GEN *a, GEN *b, const proj *E, const point *Q, const point *P1, const point *P2, long ell, long e);

#endif
