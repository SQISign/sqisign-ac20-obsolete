#ifndef MONT_H
#define MONT_H

#include "uintbig.h"
#include "fp2.h"

// curve of the form y^2 = x^3 + (a/c)*x^2 + x
// twist of the form B*y^2 = x^3 + (a/c)*x^2 + x where B = Fp2_inv(fp2_non_residue());

/* P^1 over fp2. */
typedef struct proj {
    fp2 x;
    fp2 z;
} proj;

/* P^2 over fp2. */
typedef struct proj2 {
    fp2 x;
    fp2 y;
    fp2 z;
} proj2;

void xDBL(proj *Q, proj const *A, proj const *P);
void xADD(proj *S, proj const *P, proj const *Q, proj const *PQ);
void xDBLADD(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A);
void xMUL(proj *Q, proj const *A, proj const *P, uintbig const *k);
void xBIDIM(proj *S, proj const *A, proj const *P, uintbig const *k, proj const *Q, uintbig const *l, proj const *PQ);
void xISOG_many(proj *A, proj *P, int n, proj const *K, long long k);
static inline void xISOG(proj *A, proj *P, proj const *K, long long k) { xISOG_many(A, P, 1, K, k); }
void xISOG_old(proj *A, proj *P, proj const *K, long long k);

bool is_on_curve(const proj *P, const proj *A);
bool mont_equal(proj const *P1, proj const *P2);
bool mont_iszero(proj const *P);


void normalize_proj(proj *A);

// returns false of it is on the curve, true if it is on the twist
bool xLIFT(fp2 *y, const proj *A, const proj *P);

// Given x(P) and x(Q) both in A or both not in A, computes x(P±Q)
void xBILIFT(proj *PQ1, proj *PQ2, const proj *P, const proj *Q, const proj *A);

// computes P1+P2 or P1-P2
// slow, should not be used for fast arithmetic, only when xADD cannot be used
void mont_add(proj *Q, proj const *A, proj const *P1, proj const *P2);


void xyADD(proj2 *Q, proj const *A, proj2 const *P1, proj2 const *P2);
void xyDBL(proj2 *Q, proj const *A, proj2 const *P1);
void xyNEG(proj2 *Q, proj2 const *P1);
void xyMUL(proj2 *Q, proj const *A, proj2 const *P, uintbig const *k);
bool xy_is_on_curve(const proj *A, const proj2 *P);
bool xy_is_zero(const proj2 *P);
bool xy_equal(const proj2 *P1, const proj2 *P2);
void xtoxy(proj2 *Q, const proj *A, const proj *P);
void xytox(proj *Q, const proj2 *P);

//perform DLP on points of order 2^e, not very optimized but works
bool mont_two_DLP(uint64_t *a,const proj *A, const proj *Q, const proj *P,const proj *PQ,long e);

#endif
