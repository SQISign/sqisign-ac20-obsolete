#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "steps.h"
#include "mont.h"
#include "uintbig.h"
#include "poly.h"

void xDBLADD(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A24)
{
    fp2 tmp0, tmp1, tmp2;        //requires precomputation of A24=(A+2C:4C)
    proj P_copy = *P, Q_copy = *Q;

    fp2_add3(&tmp0, &P->x, &P->z);
    fp2_sub3(&tmp1, &P->x, &P->z);
    fp2_sq2(&R->x, &tmp0);
    fp2_sub3(&tmp2, &Q->x, &Q->z);
    fp2_add3(&S->x, &Q->x, &Q->z);
    fp2_mul2(&tmp0, &tmp2);
    fp2_sq2(&R->z, &tmp1);
    fp2_mul2(&tmp1, &S->x);
    fp2_sub3(&tmp2, &R->x, &R->z);
    fp2_mul2(&R->z, &A24->z);
    fp2_mul2(&R->x, &R->z);
    fp2_mul3(&S->x, &A24->x, &tmp2);
    fp2_sub3(&S->z, &tmp0, &tmp1);
    fp2_add2(&R->z, &S->x);
    fp2_add3(&S->x, &tmp0, &tmp1);
    fp2_mul2(&R->z, &tmp2);
    fp2_sq1(&S->z);
    fp2_sq1(&S->x);
    fp2_mul2(&S->z, &PQ->x);
    fp2_mul2(&S->x, &PQ->z);

    if (mont_iszero(&Q_copy)) { *S = P_copy; } // doesn't work without this check
}

void xDBL(proj *Q, proj const *A, proj const *P)
{
    fp2 a, b, c;
    fp2_add3(&a, &P->x, &P->z);
    fp2_sq1(&a);
    fp2_sub3(&b, &P->x, &P->z);
    fp2_sq1(&b);
    fp2_sub3(&c, &a, &b);
    fp2_add2(&b, &b); fp2_add2(&b, &b); /* multiplication by 4 */
    fp2_mul2(&b, &A->z);
    fp2_mul3(&Q->x, &a, &b);
    fp2_add3(&a, &A->z, &A->z); /* multiplication by 2 */
    fp2_add2(&a, &A->x);
    fp2_mul2(&a, &c);
    fp2_add2(&a, &b);
    fp2_mul3(&Q->z, &a, &c);
}

void xADD(proj *S, proj const *P, proj const *Q, proj const *PQ)
{
    fp2 a, b, c, d;
    fp2_add3(&a, &P->x, &P->z);
    fp2_sub3(&b, &P->x, &P->z);
    fp2_add3(&c, &Q->x, &Q->z);
    fp2_sub3(&d, &Q->x, &Q->z);
    fp2_mul2(&a, &d);
    fp2_mul2(&b, &c);
    fp2_add3(&c, &a, &b);
    fp2_sub3(&d, &a, &b);
    fp2_sq1(&c);
    fp2_sq1(&d);
    fp2_mul3(&S->x, &PQ->z, &c);
    fp2_mul3(&S->z, &PQ->x, &d);
}

/* Montgomery ladder. */
/* P must not be the unique point of order 2. */
/* not constant-time! */
void xMUL(proj *Q, proj const *A, proj const *P, uintbig const *k)
{
    proj R = *P;
    proj A24;
    const proj Pcopy = *P; /* in case Q = P */

    Q->x = fp2_1;
    Q->z = fp2_0;

    fp2_add3(&A24.x, &A->z, &A->z);    //precomputation of A24=(A+2C:4C)
    fp2_add3(&A24.z, &A24.x, &A24.x);
    fp2_add2(&A24.x, &A->x);

    unsigned long i = BITS;
    while (--i && !uintbig_bit(k, i));

    do {

        bool bit = uintbig_bit(k, i);

        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp2_cswap(&Q->x, &R.x, bit);
        //fp2_cswap(&Q->z, &R.z, bit);

        xDBLADD(Q, &R, Q, &R, &Pcopy, &A24);

        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp2_cswap(&Q->x, &R.x, bit);
        //fp2_cswap(&Q->z, &R.z, bit);

    } while (i--);
}

// returns false of it is on the curve, true if it is on the twist
bool xLIFT(fp2 *y, const proj *A, const proj *P) {

    if (fp2_iszero(&P->z)) return false;

    fp2 z2, tmp1, tmp2, y2;

    // (X^2 + Z^2) C
    fp2_sq2(&tmp1, &P->x);
    fp2_sq2(&z2, &P->z);
    fp2_add2(&tmp1, &z2);
    fp2_mul2(&tmp1, &A->z);

    // X^2C + AXZ + Z^2C
    fp2_mul3(&tmp2, &P->x, &P->z);
    fp2_mul2(&tmp2, &A->x);
    fp2_add2(&tmp1, &tmp2);

    // X^3C + AX^2Z + XZ^2C = Z^3(Cx^3 + Ax^2 + Cx) = Z^3 C (B*y^2) = Z C (B*Y^2) // x = X/Z
    fp2_mul2(&tmp1, &P->x);

    // (ZC)^(-1)
    fp2_mul3(&tmp2, &A->z, &P->z);

    assert(!fp2_iszero(&tmp2));

    fp2_inv(&tmp2);

    // (B*Y^2)
    fp2_mul3(&y2, &tmp1, &tmp2);

    *y = y2;

    if (fp2_issquare(&y2)) { // on the curve
        fp2_sqrt(y);
        return false;
    }
    else { // on the twist
        fp2 tmp = fp2_non_residue();
        fp2_mul2(y, &tmp);
        fp2_sqrt(y);
        return true;
    }

}

// Given x(P) and x(Q) both in A or both not in A, computes x(P±Q)
void xBILIFT(proj *PQ1, proj *PQ2, const proj *P, const proj *Q, const proj *A) {
  fp2 xPxQ, xPzQ, zPxQ, zPzQ, disc;
  fp2_mul3(&xPxQ, &P->x, &Q->x);
  fp2_mul3(&xPzQ, &P->x, &Q->z);
  fp2_mul3(&zPxQ, &P->z, &Q->x);
  fp2_mul3(&zPzQ, &P->z, &Q->z);

  // denominator
  fp2_sub3(&PQ1->z, &xPzQ, &zPxQ);
  fp2_sq1(&PQ1->z);
  fp2_mul2(&PQ1->z, &A->z);
  PQ2->z = PQ1->z;

  // trace
  fp2_add3(&PQ1->x, &xPzQ, &zPxQ);
  fp2_add3(&PQ2->x, &xPxQ, &zPzQ);
  fp2_mul2(&PQ1->x, &PQ2->x);
  fp2_mul2(&PQ1->x, &A->z);
  fp2_mul3(&PQ2->x, &xPxQ, &zPzQ);
  fp2_mul2(&PQ2->x, &A->x);
  fp2_add2(&PQ2->x, &PQ2->x);
  fp2_add2(&PQ1->x, &PQ2->x);

  // discriminant
  fp2_sub3(&disc, &xPxQ, &zPzQ);
  fp2_sq1(&disc);
  fp2_mul2(&disc, &A->z);
  fp2_mul2(&disc, &PQ1->z);
  fp2_sq2(&PQ2->x, &PQ1->x);
  fp2_sub3(&disc, &PQ2->x, &disc);
  fp2_sqrt(&disc);

  // finish off
  fp2_add3(&PQ2->x, &PQ1->x, &disc);
  fp2_sub2(&PQ1->x, &disc);
}

// compute S = k*P + l*Q, with PQ = P+Q
void xBIDIM(proj *S, proj const *A, proj const *P, uintbig const *k, proj const *Q, uintbig const *l, proj const *PQ) {
    proj A24;

    fp2_add3(&A24.x, &A->z, &A->z);    //precomputation of A24=(A+2C:4C)
    fp2_add3(&A24.z, &A24.x, &A24.x);
    fp2_add2(&A24.x, &A->x);

    const proj Pcopy = *P;
    const proj Qcopy = *Q;
    const proj PQcopy = *PQ;

    unsigned long i = BITS;
    while (--i && !uintbig_bit(k, i) && !uintbig_bit(l, i));

    proj a_b, a1_b, a_b1, a1_b1;
    proj T;

    a_b.x = fp2_1;
    a_b.z = fp2_0;

    a1_b = *P;
    a_b1 = *Q;
    a1_b1 = *PQ;

    do {
        bool bitk = uintbig_bit(k, i);
        bool bitl = uintbig_bit(l, i);

        if (!bitk && !bitl) {
            xDBLADD(&T, &a1_b, &a_b, &a1_b, &Pcopy, &A24);
            xADD(&a_b1, &a_b, &a_b1, &Qcopy);
            xADD(&a1_b1, &a_b, &a1_b1, &PQcopy);
            a_b = T;
        }
        else if (bitk && !bitl) {
            xADD(&a_b1, &a_b, &a1_b1, &PQcopy);
            xDBLADD(&T, &a_b, &a1_b, &a_b, &Pcopy, &A24);
            xADD(&a1_b1, &a1_b, &a1_b1, &Qcopy);
            a1_b = T;
        }
        else if (!bitk && bitl) {
            xADD(&a1_b, &a_b, &a1_b1, &PQcopy);
            xDBLADD(&T, &a_b, &a_b1, &a_b, &Qcopy, &A24);
            xADD(&a1_b1, &a_b1, &a1_b1, &Pcopy);
            a_b1 = T;
        }
        else {
            xDBLADD(&T, &a1_b, &a1_b1, &a1_b, &Qcopy, &A24);
            xADD(&a_b, &a_b, &a1_b1, &PQcopy);
            xADD(&a_b1, &a1_b1, &a_b1, &Pcopy);
            a1_b1 = T;
        }

    } while (i--);
    *S = a_b;
}


//simultaneous square-and-multiply, computes x^exp and y^exp
void exp_by_squaring_(fp2* x, fp2* y, uint64_t exp)
{
        fp2 result1, result2;
        fp2_set(&result1, 1);
        fp2_set(&result2, 1);

    while (exp)
    {
        if (exp & 1){
          fp2_mul2(&result1, x);
          fp2_mul2(&result2, y);
        }

        fp2_sq1(x);
        fp2_sq1(y);
        exp >>= 1;
    }

    fp2_cswap(&result1, x, 1);
    fp2_cswap(&result2, y, 1);

}

/* goal: compute coeffs of Z^0,Z^1,Z^2 in the poly */
/* F0(Z,P)Q^2 + F1(Z,P)Q + F2(Z,P) */
/* = (Z-P)^2*Q^2 - 2*((Z*P+1)*(Z+P)+2*A*Z*P)*Q + (Z*P-1)^2 */
/* = (P-Q)^2*Z^2 - 2*((Q*P+1)*(Q+P)+2*A*Q*P)*Z + (Q*P-1)^2 */
/* but multiply by denominators to avoid divisions: */
/* out[2] = Az*Pz^2*Qz^2*(P-Q)^2 */
/* out[1] = Az*Pz^2*Qz^2*-2*((Q*P+1)*(Q+P)+2*A*Q*P) */
/* out[0] = Az*Pz^2*Qz^2*(Q*P-1)^2 */
/* i.e.: */
/* out[2] = Az*(Px*Qz-Qx*Pz)^2 */
/* out[1] = -2*(Az*(Qx*Px+Qz*Pz)*(Qx*Pz+Px*Qz)+2*Ax*Qx*Qz*Px*Pz) */
/* out[0] = Az*(Qx*Px-Qz*Pz)^2 */
void biquad(fp2 *out,const proj *P,const proj *Q,const proj *A)
{
  fp2 PxQx; fp2_mul3(&PxQx,&P->x,&Q->x);
  fp2 PxQz; fp2_mul3(&PxQz,&P->x,&Q->z);
  fp2 PzQx; fp2_mul3(&PzQx,&P->z,&Q->x);
  fp2 PzQz; fp2_mul3(&PzQz,&P->z,&Q->z);

  fp2 s; fp2_add3(&s,&PxQx,&PzQz);
  fp2 t; fp2_add3(&t,&PzQx,&PxQz);
  fp2_mul3(&out[1],&s,&t);
  fp2_mul2(&out[1],&A->z);
  fp2 u; fp2_mul3(&u,&PxQx,&PzQz);
  fp2_mul2(&u,&A->x);
  fp2_add2(&u,&u);
  fp2_add2(&out[1],&u);
  fp2_add2(&out[1],&out[1]);
  fp2_neg1(&out[1]); /* XXX: push through other computations? */

  fp2_sub3(&out[2],&PxQz,&PzQx);
  fp2_sq1(&out[2]);
  fp2_mul2(&out[2],&A->z);

  fp2_sub3(&out[0],&PxQx,&PzQz);
  fp2_sq1(&out[0]);
  fp2_mul2(&out[0],&A->z);
}

/* same as biquad but with Q inverted */
void biquad_inv(fp2 *out,const proj *P,const proj *Q,const proj *A)
{
  fp2 PxQz; fp2_mul3(&PxQz,&P->x,&Q->z);
  fp2 PxQx; fp2_mul3(&PxQx,&P->x,&Q->x);
  fp2 PzQz; fp2_mul3(&PzQz,&P->z,&Q->z);
  fp2 PzQx; fp2_mul3(&PzQx,&P->z,&Q->x);

  fp2 s; fp2_add3(&s,&PxQz,&PzQx);
  fp2 t; fp2_add3(&t,&PzQz,&PxQx);
  fp2_mul3(&out[1],&s,&t);
  fp2_mul2(&out[1],&A->z);
  fp2 u; fp2_mul3(&u,&PxQz,&PzQx);
  fp2_mul2(&u,&A->x);
  fp2_add2(&u,&u);
  fp2_add2(&out[1],&u);
  fp2_add2(&out[1],&out[1]);
  fp2_neg1(&out[1]); /* XXX: push through other computations? */

  fp2_sub3(&out[2],&PxQx,&PzQz);
  fp2_sq1(&out[2]);
  fp2_mul2(&out[2],&A->z);

  fp2_sub3(&out[0],&PxQz,&PzQx);
  fp2_sq1(&out[0]);
  fp2_mul2(&out[0],&A->z);
}

/* biquad and biquad_inv */
void biquad_both(fp2 *out,fp2 *outinv,const proj *P,const proj *Q,const proj *A)
{
  fp2 PxQx; fp2_mul3(&PxQx,&P->x,&Q->x);
  fp2 PxQz; fp2_mul3(&PxQz,&P->x,&Q->z);
  fp2 PzQx; fp2_mul3(&PzQx,&P->z,&Q->x);
  fp2 PzQz; fp2_mul3(&PzQz,&P->z,&Q->z);
  fp2 PPQQ; fp2_mul3(&PPQQ,&PxQx,&PzQz);
  fp2_add2(&PPQQ,&PPQQ);
  fp2_mul2(&PPQQ,&A->x);

  fp2 s,t;

  fp2_add3(&s,&PxQx,&PzQz);
  fp2_add3(&t,&PzQx,&PxQz);
  fp2_mul3(&out[1],&s,&t);
  fp2_mul2(&out[1],&A->z);
  fp2_add2(&out[1],&PPQQ);
  fp2_add2(&out[1],&out[1]);
  fp2_neg1(&out[1]); /* XXX: push through other computations? */

  fp2_sub3(&out[2],&PxQz,&PzQx);
  fp2_sq1(&out[2]);
  fp2_mul2(&out[2],&A->z);

  fp2_sub3(&out[0],&PxQx,&PzQz);
  fp2_sq1(&out[0]);
  fp2_mul2(&out[0],&A->z);

  /* ----- */

  fp2_add3(&s,&PxQz,&PzQx);
  fp2_add3(&t,&PzQz,&PxQx);
  fp2_mul3(&outinv[1],&s,&t);
  fp2_mul2(&outinv[1],&A->z);
  fp2_add2(&outinv[1],&PPQQ);
  fp2_add2(&outinv[1],&outinv[1]);
  fp2_neg1(&outinv[1]); /* XXX: push through other computations? */

  outinv[2] = out[0];
  outinv[0] = out[2];
}

/* biquad specifically for Q=1 */
/* i.e.: */
/* out[1] = -2*(Az*(Px+Pz)^2+2*Ax*Px*Pz) */
/* out[2] = out[0] = Az*(Px-Pz)^2 */
void biquad_1(fp2 *out,const proj *P,const proj *A)
{
  fp2 Pplus; fp2_add3(&Pplus,&P->x,&P->z);
  fp2 Pminus; fp2_sub3(&Pminus,&P->x,&P->z);
  fp2_sq1(&Pplus);
  fp2_sq1(&Pminus);

  fp2_mul3(&out[0],&Pminus,&A->z);

  out[2] = out[0];

  fp2_sub3(&out[1],&Pminus,&Pplus);
  fp2_mul2(&out[1],&A->x);

  fp2 t; fp2_mul3(&t,&Pplus,&A->z);
  fp2_add2(&t,&t);
  fp2_sub2(&out[1],&t);

/*
  proj Q;
  Q.x = fp2_1;
  Q.z = fp2_1;
  fp2 outcheck[3];
  biquad(outcheck,P,&Q,A);
  assert(!memcmp(out,&outcheck,3*sizeof(fp2)));
*/
}

/* biquad specifically for Q=-1 */
/* i.e.: */
/* out[1] = 2*(Az*(Px-Pz)^2+2*Ax*Px*Pz) */
/* out[2] = out[0] = Az*(Px+Pz)^2 */
void biquad_minus1(fp2 *out,const proj *P,const proj *A)
{
  fp2 Pplus; fp2_add3(&Pplus,&P->x,&P->z);
  fp2 Pminus; fp2_sub3(&Pminus,&P->x,&P->z);
  fp2_sq1(&Pplus);
  fp2_sq1(&Pminus);

  fp2_mul3(&out[0],&Pplus,&A->z);

  out[2] = out[0];

  fp2_sub3(&out[1],&Pplus,&Pminus);
  fp2_mul2(&out[1],&A->x);

  fp2 t; fp2_mul3(&t,&Pminus,&A->z);
  fp2_add2(&t,&t);
  fp2_add2(&out[1],&t);

/*
  proj Q;
  Q.z = fp2_1;
  fp2_sub3(&Q.x,&fp2_0,&Q.z);
  fp2 outcheck[3];
  biquad(outcheck,P,&Q,A);
  assert(!memcmp(out,&outcheck,3*sizeof(fp2)));
*/
}

/* biquad_1 and biquad_minus1 */
void biquad_pm1(fp2 *outplus,fp2 *outminus,const proj *P,const proj *A)
{
  fp2 Pplus; fp2_add3(&Pplus,&P->x,&P->z);
  fp2 Pminus; fp2_sub3(&Pminus,&P->x,&P->z);
  fp2_sq1(&Pplus);
  fp2_sq1(&Pminus);

  fp2_mul3(&outplus[0],&Pminus,&A->z);
  outplus[2] = outplus[0];
  fp2_mul3(&outminus[0],&Pplus,&A->z);
  outminus[2] = outminus[0];

  fp2 u;
  fp2_sub3(&u,&Pminus,&Pplus);
  fp2_mul2(&u,&A->x);

  fp2 t;
  fp2_add3(&t,&outminus[0],&outminus[0]);
  fp2_sub3(&outplus[1],&u,&t);

  fp2_add3(&t,&outplus[0],&outplus[0]);
  fp2_sub3(&outminus[1],&t,&u);
}

/* computes the isogeny with kernel point K of order k */
/* returns the new curve coefficient A and the image of P */
/* (obviously) not constant time in k */
void xISOG_many(proj *A, proj *P, int n, proj const *K, long long k)
{
    assert (k >= 3);
    assert (k % 2 == 1);

    long long sqrtvelu = 0;
    long long bs = 0;
    long long gs = 0;

    steps(&bs,&gs,k);
    if (bs) {
      sqrtvelu = 1;
      assert(bs > 0);
      assert(gs > 0);
      assert(!(bs&1));
    }

    fp2 tmp0, tmp1, tmp2, tmp3, tmp4;

    proj Aed;
    fp2_add3(&Aed.z, &A->z, &A->z);  //compute twisted Edwards curve coefficients
    fp2_add3(&Aed.x, &A->x, &Aed.z);
    fp2_sub3(&Aed.z, &A->x, &Aed.z);

    fp2 Psum[n], Pdif[n];
    for (int r = 0; r < n; r++) {
      fp2_add3(Psum + r, &(P+r)->x, &(P+r)->z);   //precomputations
      fp2_sub3(Pdif + r, &(P+r)->x, &(P+r)->z);
    }

    int Minit[k];
    proj M[k]; /* XXX: use less space */

    for (long long s = 0;s < k;++s) Minit[s] = 0;

    M[1] = *K; Minit[1] = 1;
    xDBL(&M[2], A, K); Minit[2] = 1;

    if (sqrtvelu) {
      for (long long s = 3;s < k;++s) {
        if (s & 1) {
          long long i = s/2;
          assert(s == 2*i+1);
          if (i < bs) {
            if (s == 3) {
              assert(Minit[1]);
              assert(Minit[2]);
              xADD(&M[s],&M[2],&M[1],&M[1]);
              Minit[s] = 1;
              continue;
            }
            assert(Minit[s-2]);
            assert(Minit[s-4]);
            assert(Minit[2]);
            xADD(&M[s],&M[s-2],&M[2],&M[s-4]);
            Minit[s] = 1;
            continue;
          }
        } else {
          long long i = s/2-1;
          assert(s == 2*i+2);
          if (i < (k-1)/2-2*bs*gs) {
            if (s == 4) {
              assert(Minit[2]);
              xDBL(&M[s],A,&M[2]);
              Minit[s] = 1;
              continue;
            }
            assert(Minit[s-2]);
            assert(Minit[s-4]);
            assert(Minit[2]);
            xADD(&M[s],&M[s-2],&M[2],&M[s-4]);
            Minit[s] = 1;
            continue;
          }
        }

        if (bs > 0) {
          if (s == 2*bs) {
            assert(Minit[bs-1]);
            assert(Minit[bs+1]);
            assert(Minit[2]);
            xADD(&M[s],&M[bs+1],&M[bs-1],&M[2]);
            Minit[s] = 1;
            continue;
          } else if (s == 4*bs) {
            assert(Minit[2*bs]);
            xDBL(&M[s],A,&M[2*bs]);
            Minit[s] = 1;
            continue;
          } else if (s == 6*bs) {
            assert(Minit[2*bs]);
            assert(Minit[4*bs]);
            xADD(&M[s],&M[4*bs],&M[2*bs],&M[2*bs]);
            Minit[s] = 1;
            continue;
          } else if (s%(4*bs) == 2*bs) {
            long long j = s/(4*bs);
            assert(s == 2*bs*(2*j+1));
            if (j < gs) {
              assert(Minit[s-4*bs]);
              assert(Minit[s-8*bs]);
              assert(Minit[4*bs]);
              xADD(&M[s],&M[s-4*bs],&M[4*bs],&M[s-8*bs]);
              Minit[s] = 1;
              continue;
            }
          }
        }
      }
    } else {
      for (long long i = 3;i <= (k-1)/2;++i) {
	Minit[i] = 1;
        xADD(&M[i],&M[i-1],K,&M[i-2]);
      }
    }

    proj Abatch;
    Abatch.x = fp2_1;
    Abatch.z = fp2_1;
    proj Qbatch[n];
    for (int r = 0; r < n; r++) {
      Qbatch[r].x = fp2_1;
      Qbatch[r].z = fp2_1;
    }

    if (sqrtvelu) {
      long long TIlen = 2*bs+poly_tree1size(bs);
      fp2 TI[TIlen];

      for (long long i = 0;i < bs;++i) {
        assert(Minit[2*i+1]);
        fp2_neg2(&TI[2*i],&M[2*i+1].x);
        TI[2*i+1] = M[2*i+1].z;
      }

      poly_tree1(TI+2*bs,TI,bs);

      fp2 TP[3*gs*n];
      fp2 TPinv[3*gs*n];
      fp2 T1[3*gs];
      fp2 Tminus1[3*gs];

      for (long long j = 0;j < gs;++j) {
        assert(Minit[2*bs*(2*j+1)]);
	for (int r = 0; r < n; r++)
	  biquad_both(TP+3*(j+gs*r),TPinv+3*(j+gs*r),&M[2*bs*(2*j+1)],P+r,A);
        biquad_pm1(T1+3*j,Tminus1+3*j,&M[2*bs*(2*j+1)],A);
      }

      for (int r = 0; r < n; r++) {
	poly_multiprod2(TP+3*gs*r,gs);
	poly_multiprod2(TPinv+3*gs*r,gs);
      }
      poly_multiprod2_selfreciprocal(T1,gs);
      poly_multiprod2_selfreciprocal(Tminus1,gs);

      long long precompsize = poly_multieval_precomputesize(bs,2*gs+1);
      fp2 precomp[precompsize];
      poly_multieval_precompute(precomp,bs,2*gs+1,TI,TI+2*bs);

      fp2 v[bs];

      for (int r = 0; r < n; r++) {
	poly_multieval_postcompute(v,bs,TP+3*gs*r,2*gs+1,TI,TI+2*bs,precomp);
	Qbatch[r].z = v[0];
	for (long long i = 1;i < bs;++i) fp2_mul2(&Qbatch[r].z,&v[i]);

	poly_multieval_postcompute(v,bs,TPinv+3*gs*r,2*gs+1,TI,TI+2*bs,precomp);
	Qbatch[r].x = v[0];
	for (long long i = 1;i < bs;++i) fp2_mul2(&Qbatch[r].x,&v[i]);
      }

      poly_multieval_postcompute(v,bs,T1,2*gs+1,TI,TI+2*bs,precomp);
      Abatch.x = v[0];
      for (long long i = 1;i < bs;++i) fp2_mul2(&Abatch.x,&v[i]);

      poly_multieval_postcompute(v,bs,Tminus1,2*gs+1,TI,TI+2*bs,precomp);
      Abatch.z = v[0];
      for (long long i = 1;i < bs;++i) fp2_mul2(&Abatch.z,&v[i]);

      for (long long i = 0;i < (k-1)/2-2*bs*gs;++i) {
        assert(Minit[2*i+2]);
        fp2_sub3(&tmp1, &M[2*i+2].x, &M[2*i+2].z);
        fp2_add3(&tmp0, &M[2*i+2].x, &M[2*i+2].z);
        fp2_mul2(&Abatch.x,&tmp1);
        fp2_mul2(&Abatch.z,&tmp0);
	for (int r = 0; r < n; r++) {
	  fp2_mul3(&tmp4, &tmp1, Psum+r);
	  fp2_mul3(&tmp3, &tmp0, Pdif+r);
	  fp2_add3(&tmp2, &tmp3, &tmp4);
	  fp2_mul2(&Qbatch[r].x, &tmp2);
	  fp2_sub3(&tmp2, &tmp3, &tmp4);
	  fp2_mul2(&Qbatch[r].z, &tmp2);
	}
      }
    } else {
      for (long long i = 1;i <= (k-1)/2;++i) {
        assert(Minit[i]);
        fp2_sub3(&tmp1, &M[i].x, &M[i].z);
        fp2_add3(&tmp0, &M[i].x, &M[i].z);
        if (i > 1) {
          fp2_mul2(&Abatch.x,&tmp1);
          fp2_mul2(&Abatch.z,&tmp0);
        } else {
          Abatch.x = tmp1;
          Abatch.z = tmp0;
        }
	for (int r = 0; r < n; r++) {
	  fp2_mul3(&tmp4, &tmp1, Psum+r);
	  fp2_mul3(&tmp3, &tmp0, Pdif+r);
	  fp2_add3(&tmp2, &tmp3, &tmp4);
	  if (i > 1) {
	    fp2_mul2(&Qbatch[r].x, &tmp2);
	  } else {
	    Qbatch[r].x = tmp2;
	  }
	  fp2_sub3(&tmp2, &tmp3, &tmp4);
	  if (i > 1) {
	    fp2_mul2(&Qbatch[r].z, &tmp2);
	  } else {
	    Qbatch[r].z = tmp2;
	  }
	}
      }
    }


    // point evaluation
    for (int r = 0; r < n; r++) {
      fp2_sq1(&Qbatch[r].x);
      fp2_sq1(&Qbatch[r].z);
      fp2_mul2(&P[r].x, &Qbatch[r].x);
      fp2_mul2(&P[r].z, &Qbatch[r].z);
    }

    //compute Aed.x^k, Aed.z^k
    exp_by_squaring_(&Aed.x, &Aed.z, k);

    //compute prod.x^8, prod.z^8
    fp2_sq1(&Abatch.x);
    fp2_sq1(&Abatch.x);
    fp2_sq1(&Abatch.x);
    fp2_sq1(&Abatch.z);
    fp2_sq1(&Abatch.z);
    fp2_sq1(&Abatch.z);

    //compute image curve parameters
    fp2_mul2(&Aed.z, &Abatch.x);
    fp2_mul2(&Aed.x, &Abatch.z);

    //compute Montgomery params
    fp2_add3(&A->x, &Aed.x, &Aed.z);
    fp2_sub3(&A->z, &Aed.x, &Aed.z);
    fp2_add2(&A->x, &A->x);
}


void xISOG_old(proj *A, proj *P, proj const *K, long long k)
{
    assert (k >= 3);
    assert (k % 2 == 1);

    fp2 tmp0, tmp1, tmp2, Psum, Pdif;
    proj Q, Aed, prod;

    fp2_add3(&Aed.z, &A->z, &A->z);  //compute twisted Edwards curve coefficients
    fp2_add3(&Aed.x, &A->x, &Aed.z);
    fp2_sub3(&Aed.z, &A->x, &Aed.z);

    fp2_add3(&Psum, &P->x, &P->z);   //precomputations
    fp2_sub3(&Pdif, &P->x, &P->z);

    fp2_sub3(&prod.x, &K->x, &K->z);
    fp2_add3(&prod.z, &K->x, &K->z);

    fp2_mul3(&tmp1, &prod.x, &Psum);
    fp2_mul3(&tmp0, &prod.z, &Pdif);
    fp2_add3(&Q.x, &tmp0, &tmp1);
    fp2_sub3(&Q.z, &tmp0, &tmp1);

    proj M[k];

    M[0] = *K;
    xDBL(&M[1], A, K);
    for (long long i = 2; i < k / 2; ++i) {
        xADD(&M[i], &M[(i - 1)], K, &M[(i - 2)]);
    }

    for (long long i = 1; i < k / 2; ++i) {
        fp2_sub3(&tmp1, &M[i].x, &M[i].z);
            fp2_add3(&tmp0, &M[i].x, &M[i].z);
        fp2_mul2(&prod.x, &tmp1);
        fp2_mul2(&prod.z, &tmp0);
            fp2_mul2(&tmp1, &Psum);
            fp2_mul2(&tmp0, &Pdif);
            fp2_add3(&tmp2, &tmp0, &tmp1);
        fp2_mul2(&Q.x, &tmp2);
            fp2_sub3(&tmp2, &tmp0, &tmp1);
        fp2_mul2(&Q.z, &tmp2);

    }


    // point evaluation
    fp2_sq1(&Q.x);
    fp2_sq1(&Q.z);
    fp2_mul2(&P->x, &Q.x);
    fp2_mul2(&P->z, &Q.z);

    //compute Aed.x^k, Aed.z^k
    exp_by_squaring_(&Aed.x, &Aed.z, k);

    //compute prod.x^8, prod.z^8
    fp2_sq1(&prod.x);
    fp2_sq1(&prod.x);
    fp2_sq1(&prod.x);
    fp2_sq1(&prod.z);
    fp2_sq1(&prod.z);
    fp2_sq1(&prod.z);

    //compute image curve parameters
    fp2_mul2(&Aed.z, &prod.x);
    fp2_mul2(&Aed.x, &prod.z);

    //compute Montgomery params
    fp2_add3(&A->x, &Aed.x, &Aed.z);
    fp2_sub3(&A->z, &Aed.x, &Aed.z);
    fp2_add2(&A->x, &A->x);
}

bool is_on_curve(const proj *P, const proj *A) {
  fp2 xz, tmp1, tmp2;

  // (X^2 + Z^2) C
  fp2_sq2(&tmp1, &P->x);
  fp2_sq2(&tmp2, &P->z);
  fp2_add2(&tmp1, &tmp2);
  fp2_mul2(&tmp1, &A->z);

  fp2_mul3(&xz, &P->x, &P->z);

  // X^2C + AXZ + Z^2C
  fp2_mul3(&tmp2, &xz, &A->x);
  fp2_add2(&tmp1, &tmp2);

  fp2_mul2(&tmp1, &xz);
  fp2_mul2(&tmp1, &A->z);

  return fp2_issquare(&tmp1);
}

bool mont_equal(proj const *P1, proj const *P2) {
    fp2 x1z2, x2z1;
    fp2 x1z2_x2z1;

    fp2_mul3(&x1z2, &P1->x, &P2->z);
    fp2_mul3(&x2z1, &P2->x, &P1->z);
    fp2_sub3(&x1z2_x2z1, &x1z2, &x2z1);

    return fp2_iszero(&x1z2_x2z1);
}

bool mont_iszero(proj const *P) {
    return fp2_iszero(&P->z);
}

void mont_add(proj *Q, proj const *A, proj const *P1, proj const *P2) {
  // y^2 = x^3 + (a/c)*x^2 + x
  fp2 y1, y2, x3, z3;
  bool twist1 = xLIFT(&y1, A, P1);
  bool twist2 = xLIFT(&y2, A, P2);
  assert(twist1==twist2);
  assert(!twist1);

  fp2_mul2(&y1,&P2->x);
  fp2_mul2(&y2,&P1->x);

  fp2_sub3(&x3, &y1, &y2);
  fp2_sq1(&x3);
  fp2_mul2(&x3,&P1->z);
  fp2_mul2(&x3,&P2->z);

  fp2_mul3(&y1,&P2->x,&P1->z);
  fp2_mul3(&y2,&P1->x,&P2->z);
  fp2_sub3(&z3, &y1, &y2);
  fp2_sq1(&z3);
  fp2_mul2(&z3,&P1->x);
  fp2_mul2(&z3,&P2->x);

  Q->x = x3;
  Q->z = z3;
}


void xyADD(proj2 *Q, proj const *A, proj2 const *P1, proj2 const *P2) {
  if(xy_is_zero(P1)) { *Q = *P2; return; }
  if(xy_is_zero(P2)) { *Q = *P1; return; }
  if(xy_equal(P1,P2)) { xyDBL(Q, A, P1); return; }

  fp2 a,b;
  fp2 U,V,x3,y3,z3,X21,Y21;


  assert(xy_is_on_curve(A, P1));
  assert(xy_is_on_curve(A, P2));

  a = A->z;
  assert(!fp2_iszero(&a));
  fp2_inv(&a);
  fp2_mul2(&a,&A->x);

  b = fp2_1;

  fp2_mul3(&U,&P2->x,&P1->z);
  fp2_mul3(&V,&P1->x,&P2->z);
  fp2_sub3(&X21, &U, &V);

  fp2_mul3(&U,&P2->y,&P1->z);
  fp2_mul3(&V,&P1->y,&P2->z);
  fp2_sub3(&Y21, &U, &V);

  // z3 = z1z2 * x1x2 * (x2z1 - x1z2)^3
  fp2_sq2(&z3,&X21);
  fp2_mul2(&z3,&X21);
  fp2_mul2(&z3,&P1->x);
  fp2_mul2(&z3,&P2->x);
  fp2_mul2(&z3,&P1->z);
  fp2_mul2(&z3,&P2->z);

  // x3 = b * (x2y1 - x1y2)^2 * (z1z2)^2 * (x2z1 - x1z2)
  fp2_mul3(&U,&P2->x,&P1->y);
  fp2_mul3(&V,&P1->x,&P2->y);
  fp2_sub3(&x3, &U, &V);
  fp2_sq1(&x3);
  fp2_mul2(&x3,&b);
  fp2_mul2(&x3,&P1->z);
  fp2_mul2(&x3,&P2->z);
  fp2_mul2(&x3,&P1->z);
  fp2_mul2(&x3,&P2->z);
  fp2_mul2(&x3,&X21);

  fp2 y3_A, y3_B, y3_C;

  fp2_mul3(&U,&P1->x,&P2->z);
  fp2_add2(&U,&U);
  fp2_mul3(&V,&P2->x,&P1->z);
  fp2_add2(&U,&V);
  fp2_mul3(&V,&P2->z,&P1->z);
  fp2_mul2(&V,&a);
  fp2_add2(&U,&V);

  fp2_sq2(&y3_A,&X21);
  fp2_mul2(&y3_A,&U);
  fp2_mul2(&y3_A,&Y21);

  fp2_sq2(&y3_B,&Y21);
  fp2_mul2(&y3_B,&Y21);
  fp2_mul2(&y3_B,&b);
  fp2_mul2(&y3_B,&P1->z);
  fp2_mul2(&y3_B,&P2->z);

  fp2_sq2(&y3_C,&X21);
  fp2_mul2(&y3_C,&X21);
  fp2_mul2(&y3_C,&P1->y);
  fp2_mul2(&y3_C,&P2->z);

  fp2_sub3(&y3, &y3_A, &y3_B);
  fp2_sub2(&y3, &y3_C);
  fp2_mul2(&y3,&P1->x);
  fp2_mul2(&y3,&P2->x);

  Q->x = x3;
  Q->y = y3;
  Q->z = z3;


  assert(xy_is_on_curve(A, Q));
}

void xyDBL(proj2 *Q, proj const *A, proj2 const *P1){

  if(xy_is_zero(P1)) { *Q = *P1; return; }

  fp2 a,b;
  fp2 U,x3,y3,z3, Dx, Dy, xz,x2, z2, By2;
  const fp2 x = P1->x, y = P1->y, z = P1->z;

  assert(xy_is_on_curve(A, P1));

  a = A->z;
  assert(!fp2_iszero(&a));
  fp2_inv(&a);
  fp2_mul2(&a,&A->x);

  b = fp2_1;

  fp2_sq2(&x2,&x);
  fp2_sq2(&z2,&z);
  fp2_mul3(&xz,&x,&z);

  // Dx = 4*x*(x^2+A*x*z+z^2)
  fp2_mul3(&Dx,&a,&xz);
  fp2_add2(&Dx,&x2);
  fp2_add2(&Dx,&z2);
  fp2_mul2(&Dx,&x);
  fp2_add2(&Dx,&Dx);
  fp2_add2(&Dx,&Dx);

  // Dy = (2*B*y)^3
  fp2_mul3(&By2,&b,&y);
  fp2_add2(&By2,&By2);
  fp2_sq2(&Dy,&By2);
  fp2_mul2(&Dy,&By2);


  // z3 = Dx*Dy*z^3
  fp2_mul3(&z3,&Dx,&Dy);
  fp2_mul2(&z3,&z);
  fp2_mul2(&z3,&z);
  fp2_mul2(&z3,&z);


  // x3 = Dy*(x^2-z^2)^2*z^2
  fp2_sub3(&x3,&x2,&z2);
  fp2_sq1(&x3);
  fp2_mul2(&x3,&z);
  fp2_mul2(&x3,&z);
  fp2_mul2(&x3,&Dy);


  // y3

  fp2 y3_A, y3_B, y3_C;

  fp2_mul3(&U,&a,&xz);
  fp2_add2(&U,&U);
  fp2_add2(&U,&z2);
  fp2_add2(&U,&x2);
  fp2_add2(&U,&x2);
  fp2_add2(&U,&x2);

  fp2_mul3(&y3_A, &a, &z);
  fp2_add2(&y3_A, &x);
  fp2_add2(&y3_A, &x);
  fp2_add2(&y3_A, &x);

  fp2_mul2(&y3_A,&U);
  fp2_mul2(&y3_A,&By2);
  fp2_mul2(&y3_A,&By2);
  fp2_mul2(&y3_A,&z);

  fp2_sq2(&y3_B, &U);
  fp2_mul2(&y3_B, &U);
  fp2_mul2(&y3_B, &b);


  fp2_mul3(&y3_C, &Dy, &y);
  fp2_mul2(&y3_C, &z2);


  fp2_sub3(&y3, &y3_A, &y3_B);
  fp2_sub2(&y3, &y3_C);
  fp2_mul2(&y3,&Dx);

  Q->x = x3;
  Q->y = y3;
  Q->z = z3;


  assert(xy_is_on_curve(A, Q));
}


void xyNEG(proj2 *Q, proj2 const *P1) {
  *Q = *P1;
  fp2_neg1(&Q->y);
}
//void xyMUL(proj2 *Q, proj const *A, proj2 const *P, uintbig const *k);


bool xy_is_zero(const proj2 *P) {
  return (fp2_iszero(&P->z) && fp2_iszero(&P->x) && !fp2_iszero(&P->y));
}

bool xy_equal(const proj2 *P1, const proj2 *P2) {
    fp2 x1z2,x2z1;
    fp2 x1z2_x2z1;
    fp2 y1z2,y2z1;
    fp2 y1z2_y2z1;
    fp2_mul3(&x1z2, &P1->x, &P2->z);
    fp2_mul3(&x2z1, &P2->x, &P1->z);
    fp2_sub3(&x1z2_x2z1, &x1z2, &x2z1);

    fp2_mul3(&y1z2, &P1->y, &P2->z);
    fp2_mul3(&y2z1, &P2->y, &P1->z);
    fp2_sub3(&y1z2_y2z1, &y1z2, &y2z1);

    return fp2_iszero(&x1z2_x2z1) && fp2_iszero(&y1z2_y2z1);
}

// check if it is on the curve (with B = 1, i.e., not the twist)
bool xy_is_on_curve(const proj *A, const proj2 *P) {

    if (fp2_iszero(&P->z) && fp2_iszero(&P->x) && !fp2_iszero(&P->y)) return true;

    fp2 z2, tmp1, tmp2, tmp3;

    // (X^2 + Z^2) C
    fp2_sq2(&tmp1, &P->x);
    fp2_sq2(&z2, &P->z);
    fp2_add2(&tmp1, &z2);
    fp2_mul2(&tmp1, &A->z);

    // X^2C + AXZ + Z^2C
    fp2_mul3(&tmp2, &P->x, &P->z);
    fp2_mul2(&tmp2, &A->x);
    fp2_add2(&tmp1, &tmp2);

    // X^3C + AX^2Z + XZ^2C = Z^3(Cx^3 + Ax^2 + Cx) = Z^3 C (B*y^2) = Z C (B*Y^2) // x = X/Z, B = 1
    fp2_mul2(&tmp1, &P->x);

    // (ZC)
    fp2_mul3(&tmp2, &A->z, &P->z);

    fp2_sq2(&tmp3,&P->y);
    fp2_mul2(&tmp3, &tmp2);

    return fp2_equal(&tmp1,&tmp3);

}


void xyMUL(proj2 *Q, proj const *A, proj2 const *P, uintbig const *k)
{
    proj2 R = *P;

    Q->x = fp2_0;
    Q->y = fp2_1;
    Q->z = fp2_0;

    unsigned long i = BITS;
    while (--i && !uintbig_bit(k, i));

    do {

        bool bit = uintbig_bit(k, i);

        if (bit) { proj2 T = *Q; *Q = R; R = T; }

        xyADD(&R, A, Q, &R);
        xyDBL(Q, A, Q);

        if (bit) { proj2 T = *Q; *Q = R; R = T; }

    } while (i--);
}


void xtoxy(proj2 *Q, const proj *A, const proj *P) {
  if (fp2_iszero(&P->z)) { Q->x = fp2_0; Q->y = fp2_1; Q->z = fp2_0; }
  else {
    Q->x = P->x; Q->z = P->z;
    assert(!xLIFT(&Q->y, A, P)); // doesn't handle the twist
  }
}

void xytox(proj *Q, const proj2 *P) {
  if (fp2_iszero(&P->z)) { Q->x = fp2_1; Q->z = fp2_0; }
  else { Q->x = P->x; Q->z = P->z; }
}

void normalize_proj(proj *A){
  fp2 tmpfp2;
  tmpfp2= A->z;
  fp2_inv(&tmpfp2);
  fp2_mul3(&(A->x),&(A->x),&tmpfp2);
  fp2_mul3(&(A->z),&(A->z),&tmpfp2);
}

bool naive_two_DLP(uint64_t *a, const proj*Q, const proj *P){
  if (mont_iszero(Q)) {*a = 0; return true;}
  if (mont_equal(Q,P)) {*a = 1; return true;}
  return false;
}

bool mont_two_DLP_slow(uint64_t *a,const proj *A, const proj *Q, const proj *P,const proj *PQ,long e){


  uint64_t log=0,ell_pow=1;
  long elle=pow(2,e);
  uintbig x_big,uintbig1;
  uintbig_set(&uintbig1,1);
  uint64_t x;
  proj Ps,R, Ri;

  Ps = *P;
  for (int i = 1; i < e; ++i) {
      xDBL(&Ps, A, &Ps);
      assert(!mont_iszero(&Ps));
  }

  R = *Q;
  for (int i = 0; i < e; ++i) {
      Ri = R;
      for (int j = 0; j < e-1-i; ++j) {
        xDBL(&Ri,A,&Ri);
      }
      if(!naive_two_DLP(&x, &Ri, &Ps))
          {return false; }
      log += x*ell_pow;

      if (x ==1 ) {
        uintbig_set(&x_big,elle - log);
        xBIDIM(&R,A,Q,&uintbig1,P,&x_big,PQ);

      }
      ell_pow*=2;

  }
  *a = log;
  return true;
}


// translation by the two-torsion point T = (0:0:1)
void xTRANST(proj *R, const proj *P) {
  fp2 tmp = P->x;
  R->x = P->z;
  R->z = tmp;
}

void xADD_safe(proj *R, const proj *A, const proj *P, const proj *Q, const proj *PQ) {
  proj res;
  if (mont_iszero(PQ)) { 
    assert(mont_equal(P,Q));
    xDBL(&res, A, P);
  }
  else if (fp2_iszero(&PQ->x)) { // P-Q = T, P+Q = 2Q + T
    xDBL(&res, A, P);
    xTRANST(&res, &res);
  }
  else {
    xADD(&res, P, Q, PQ);
  }
  *R = res;
}

// PQ = P-Q
bool mont_two_DLP_rec(uint64_t *a, const proj *A, long len, proj *Q, proj *P, proj *PQ, long stacklen){
  if (len == 0) {
    *a = 0;
    return true;
  }
  else if (len == 1) {
    if (mont_iszero(&Q[stacklen-1])) {
      *a = 0;
      for (int i = 0; i < stacklen-1; ++i) {
        xADD_safe(&PQ[i], A, &P[i], &PQ[i], &Q[i]); // newPQ = P + (P-Q) = 2P - Q
        xDBL(&P[i], A, &P[i]); // newP = 2P (we get  newP - newQ = 2P - Q = newPQ)
      }
      return true;
    }
    else if (mont_equal(&Q[stacklen-1],&P[stacklen-1])) {
      *a = 1;
      proj tmp;
      for (int i = 0; i < stacklen-1; ++i) {
        xADD_safe(&tmp, A, &P[i], &Q[i], &PQ[i]); // tmp = P+Q
        xDBL(&P[i], A, &P[i]); // newP = 2P
        Q[i] = PQ[i]; // new Q = Q-P
        if (fp2_iszero(&PQ[i].x)){
          PQ[i] = tmp;
        }
        else {
          xADD_safe(&PQ[i], A, &P[i], &Q[i], &tmp); // newPQ = newP-newQ = 2P + (P-Q) (diff P+Q = tmp)
        }
      }
      return true;
    }
    else { return false; }
  }
  else {
    long right = (double)len * 0.5;
    long left = len - right;
    Q[stacklen] = Q[stacklen-1];
    P[stacklen] = P[stacklen-1];
    PQ[stacklen] = PQ[stacklen-1];
    for (int i = 0; i < left; i++) {
      xDBL(&Q[stacklen], A, &Q[stacklen]);
      xDBL(&P[stacklen], A, &P[stacklen]);
      xDBL(&PQ[stacklen], A, &PQ[stacklen]);
    }
    uint64_t dlp1 = 0, dlp2 = 0;
    bool ok;
    ok = mont_two_DLP_rec(&dlp1, A, right, Q, P, PQ, stacklen+1);
    if (!ok) return false;
    ok = mont_two_DLP_rec(&dlp2, A, left, Q, P, PQ, stacklen);
    if (!ok) return false;
    *a = (dlp2 << right) + dlp1;
    return true;
  }
}

// PQ = P+Q
bool mont_two_DLP(uint64_t *a,const proj *A, const proj *Q, const proj *P,const proj *PQ,long e){
  long log, len = e;
  for (log = 0; len > 1; len >>= 1) log++;
  log += 1;

  proj stack_Q[log], stack_P[log], stack_PQ[log];

  stack_Q[0] = *Q;
  stack_P[0] = *P;
  xADD_safe(&stack_PQ[0], A, P, Q, PQ);

  *a = 0;

  bool ok = mont_two_DLP_rec(a, A, e, stack_Q, stack_P,stack_PQ, 1);



  return ok;
}

