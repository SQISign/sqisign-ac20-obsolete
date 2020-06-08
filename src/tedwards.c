#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <gmp.h>
#define FP_LIMBS (4 * 64 / GMP_LIMB_BITS)

#include "steps.h"
#include "tedwards.h"
#include "uintbig.h"
#include "poly.h"

// a*x^2+y^2=1+d*x^2*y^2
// a = A.x, d = A.z

bool ted_is_on_curve(point const *P, proj const *E) {
    fp2 x2, y2, z2, tmp1, tmp2;
    fp2_sq2(&x2,&P->x);
    fp2_sq2(&y2,&P->y);
    fp2_sq2(&z2,&P->z);

    fp2_mul3(&tmp1, &x2, &E->x);
    fp2_add2(&tmp1, &y2);
    fp2_mul2(&tmp1, &z2);

    fp2_mul3(&tmp2, &x2, &y2);
    fp2_mul2(&tmp2, &E->z);
    fp2_sq1(&z2);
    fp2_add2(&tmp2, &z2);

    fp2_sub2(&tmp1, &tmp2);

    return fp2_iszero(&tmp1);
}




void ted_double(point *Q, proj const *E, point const *P) {
    // A = X1^2
    // B = Y1^2
    // C = 2*Z1^2
    // D = a*A
    // K = (X1+Y1)^2-A-B
    // G = D+B
    // F = G-C
    // H = D-B
    // X3 = K*F
    // Y3 = G*H
    // T3 = K*H
    // Z3 = F*G

    // TODO: neutral element
    fp2 A,B,C,D,K,G,F,H;

    fp2_sq2(&A,&P->x);
    fp2_sq2(&B,&P->y);
    fp2_sq2(&C,&P->z);
    fp2_add2(&C,&C);
    fp2_mul3(&D,&A,&E->x);
    fp2_add3(&K,&P->x,&P->y);
    fp2_sq1(&K);
    fp2_sub2(&K,&A);
    fp2_sub2(&K,&B);
    fp2_add3(&G,&D,&B);
    fp2_sub3(&F,&G,&C);
    fp2_sub3(&H,&D,&B);


    fp2_mul3(&Q->x,&K,&F);
    fp2_mul3(&Q->y,&G,&H);
    fp2_mul3(&Q->t,&K,&H);
    fp2_mul3(&Q->z,&F,&G);

}

void ted_add(point *S, proj const *E, point const *P, point const *Q) {
    // A = X1*X2
    // B = Y1*Y2
    // C = Z1*T2
    // D = T1*Z2
    // K = D+C
    // F = (X1-Y1)*(X2+Y2)+B-A
    // G = B+a*A
    // H = D-C
    // X3 = K*F
    // Y3 = G*H
    // T3 = K*H
    // Z3 = F*G

    // TODO: neutral element
    point res;

    fp2 A,B,C,D,K,F,G,H,tmp;

    fp2_mul3(&A,&P->x, &Q->x);
    fp2_mul3(&B,&P->y, &Q->y);
    fp2_mul3(&C,&P->z, &Q->t);
    fp2_mul3(&D,&P->t, &Q->z);
    fp2_add3(&K,&D,&C);
    fp2_add3(&F,&Q->x, &Q->y);
    fp2_sub3(&tmp,&P->x, &P->y);
    fp2_mul2(&F,&tmp);
    fp2_add2(&F,&B);
    fp2_sub2(&F,&A);
    fp2_mul3(&G,&A,&E->x);
    fp2_add2(&G,&B);
    fp2_sub3(&H,&D,&C);


    fp2_mul3(&res.x,&K,&F);
    fp2_mul3(&res.y,&G,&H);
    fp2_mul3(&res.t,&K,&H);
    fp2_mul3(&res.z,&F,&G);

    if (fp2_iszero(&res.x) && fp2_iszero(&res.y) && fp2_iszero(&res.z)) {
        ted_double(S, E, P);
    }
    else *S = res;
}

void ted_neg(point *Q, point const *P) {
    fp2_neg2(&Q->x, &P->x);
    Q->y = P->y;
    Q->z = P->z;
    fp2_neg2(&Q->t, &P->t);
}

void ted_mul(point *res, point const *P, proj const *E, uintbig const *k)
{
    const point Pcopy = *P;
    res->x = fp2_0;
    res->y = fp2_1;
    res->z = fp2_1;

    unsigned long i = BITS;
    while (--i && !uintbig_bit(k, i));

    do {
        ted_double(res,E,res);
        if (uintbig_bit(k, i)) {
            ted_add(res, E, res, &Pcopy);
        }
    } while (i--);
}



bool ted_iszero(point const *P) {
    if (fp2_iszero(&P->x)) {
        fp2 a;
        fp2_sub3(&a, &P->y, &P->z);
        return fp2_iszero(&a);
    }
    else return false;
}

void mont_to_ted(proj *E, proj const *A, bool twist) {
    fp2 tmp, two;
    tmp = A->z;
    fp2_inv(&tmp);
    fp2_mul2(&tmp,&A->x);
    fp2_set(&two,2);
    fp2_add3(&E->x, &tmp, &two);
    fp2_sub3(&E->z, &tmp, &two);
    if (twist) {
        // B = Fp2_inv(fp2_non_residue)
        tmp = fp2_non_residue();
        fp2_mul2(&E->x,&tmp);
        fp2_mul2(&E->z,&tmp);
    }
}

void mont_to_ted_point(point *Q, proj const *A, proj const *P) {
    if (fp2_iszero(&P->z)) {
        fp2_set(&Q->x, 0);
        fp2_set(&Q->y, 1);
        fp2_set(&Q->z, 1);
        fp2_set(&Q->t, 0);
    }
    else {
        fp2 tmp, y;
        xLIFT(&y, A, P);

        fp2_add3(&tmp,&P->x,&P->z);
        fp2_mul3(&Q->x,&P->x,&tmp);

        fp2_sub3(&Q->y,&P->x,&P->z);
        fp2_mul2(&Q->y,&y);

        fp2_mul3(&Q->z,&tmp,&y);

        Q->t = Q->z;
        fp2_inv(&Q->t);
        fp2_mul2(&Q->t,&Q->x);
        fp2_mul2(&Q->t,&Q->y);
    }
}

void ted_to_mont_point(proj *Q, point const *P) {
    fp2_add3(&Q->x, &P->z, &P->y);
    fp2_sub3(&Q->z, &P->z, &P->y);
}

bool ted_equal(point const *P1, point const *P2) {
    fp2 x1z2, y1z2;
    fp2 y2z1, x2z1;
    fp2 x1z2_x2z1, y1z2_y2z1;

    fp2_mul3(&x1z2, &P1->x, &P2->z);
    fp2_mul3(&y1z2, &P1->y, &P2->z);
    fp2_mul3(&y2z1, &P2->y, &P1->z);
    fp2_mul3(&x2z1, &P2->x, &P1->z);
    fp2_sub3(&x1z2_x2z1, &x1z2, &x2z1);
    fp2_sub3(&y1z2_y2z1, &y1z2, &y2z1);

    return fp2_iszero(&x1z2_x2z1) && fp2_iszero(&y1z2_y2z1);
}


void ted_miller_dou(fp2 *cz2,fp2 *cxy,fp2 *cxz, point *P3, proj const *E, point const *P1) {
    fp2 A,B,C,D,E_,F,G,H,I,J,K;

    fp2_sq2(&A,&P1->x);
    fp2_sq2(&B,&P1->y);
    fp2_sq2(&C,&P1->z);
    fp2_add3(&D,&P1->x,&P1->y);
    fp2_sq2(&D,&D);
    fp2_add3(&E_,&P1->y,&P1->z);
    fp2_sq2(&E_,&E_);
    fp2_add3(&F,&A,&B);
    fp2_sub3(&F,&D,&F);
    fp2_add3(&G,&B,&C);
    fp2_sub3(&G,&E_,&G);
    fp2_mul3(&H,&A,&E->x);
    fp2_add3(&I,&H,&B);
    fp2_sub3(&J,&C,&I);
    fp2_add3(&K,&J,&C);

    // coefficients of the conic
    fp2_sub3(cz2,&P1->t,&P1->x);
    fp2_mul2(cz2,&P1->y);
    fp2_add2(cz2,cz2);
    fp2_add3(cxy,&J,&J);
    fp2_add2(cxy,&G);
    fp2_mul3(cxz,&P1->x,&P1->t);
    fp2_mul2(cxz,&E->x);
    fp2_sub2(cxz,&B);
    fp2_add2(cxz,cxz);

    // compute P3 = 2*P1
    fp2_mul3(&P3->x,&F,&K);
    fp2_sub3(&P3->y,&B,&H);
    fp2_mul2(&P3->y,&I);
    fp2_mul3(&P3->z,&I,&K);
    fp2_sub3(&P3->t,&B,&H);
    fp2_mul2(&P3->t,&F);
}

void ted_miller_add(fp2 *cz2,fp2 *cxy,fp2 *cxz, point *P3, proj const *E, point const *P1, point const *P2) {
    fp2 A,B,C,D,E_,F,G,H,I,tmp;

    fp2_mul3(&A,&P1->x,&P2->x);
    fp2_mul3(&B,&P1->y,&P2->y);
    fp2_mul3(&C,&P1->z,&P2->t);
    fp2_mul3(&D,&P1->t,&P2->z);
    fp2_add3(&E_,&D,&C);
    fp2_add3(&F,&P2->x,&P2->y);
    fp2_sub3(&tmp,&P1->x,&P1->y);
    fp2_mul2(&F,&tmp);
    fp2_add2(&F,&B);
    fp2_sub2(&F,&A);
    fp2_mul3(&G,&A,&E->x);
    fp2_add2(&G,&B);
    fp2_sub3(&H,&D,&C);
    fp2_mul3(&I,&P1->t,&P2->t);

    // coefficients of the conic
    fp2_sub3(cz2,&P1->t,&P1->x);
    fp2_add3(&tmp,&P2->t,&P2->x);
    fp2_mul2(cz2,&tmp);
    fp2_sub2(cz2,&I);
    fp2_add2(cz2,&A);
    fp2_mul3(cxy,&P1->x,&P2->z);
    fp2_mul3(&tmp,&P2->x,&P1->z);
    fp2_sub2(cxy,&tmp);
    fp2_add2(cxy,&F);
    fp2_sub3(cxz,&P1->y,&P1->t);
    fp2_add3(&tmp,&P2->y,&P2->t);
    fp2_mul2(cxz,&tmp);
    fp2_sub2(cxz,&B);
    fp2_add2(cxz,&I);
    fp2_sub2(cxz,&H);

    // compute P3 = 2*P1
    fp2_mul3(&P3->x,&E_,&F);
    fp2_mul3(&P3->y,&G,&H);
    fp2_mul3(&P3->t,&E_,&H);
    fp2_mul3(&P3->z,&F,&G);
}

// eta_Q = (Q.z+Q.y)/Q.x and Y_Q = Q.y/Q.z
void ted_phi_l1l2(fp2 *f, fp2 *g, const fp2 *cz2, const fp2 *cxy, const fp2 *cxz, const fp2 *eta_Q, const fp2 *Y_Q, const point *P3){
    fp2 tmp, f0, g0;
    fp2_mul3(&f0, cz2, eta_Q);
    fp2_mul3(&tmp, cxy, Y_Q);
    fp2_add2(&f0, &tmp);
    fp2_add2(&f0, cxz);
    fp2_mul3(&g0, Y_Q, &P3->z);
    fp2_sub2(&g0, &P3->y);

    *f = f0;
    *g = g0;
}

void ted_miller(fp2 *res, fp2 *res2, proj const *E, point const *P, point const *Q, point const *Q2, uintbig const *k) {
    point R = *P;
    fp2 f, g,f2,g2, cz2, cxy, cxz, f0,g0, eta_Q, Y_Q, eta_Q2, Y_Q2, tmp;
    fp2_set(&f, 1);
    fp2_set(&g, 1);

    Y_Q = Q->z;
    fp2_inv(&Y_Q);
    fp2_mul2(&Y_Q,&Q->y);

    fp2_add3(&eta_Q,&Q->z,&Q->y);
    tmp = Q->x;
    fp2_inv(&tmp);
    fp2_mul2(&eta_Q,&tmp);

    if (Q2) {
        fp2_set(&f2, 1);
        fp2_set(&g2, 1);

        Y_Q2 = Q2->z;
        fp2_inv(&Y_Q2);
        fp2_mul2(&Y_Q2,&Q2->y);

        fp2_add3(&eta_Q2,&Q2->z,&Q2->y);
        tmp = Q2->x;
        fp2_inv(&tmp);
        fp2_mul2(&eta_Q2,&tmp);
    }

    unsigned long i = BITS;
    while (--i && !uintbig_bit(k, i));
    i--;
    do {
        ted_miller_dou(&cz2,&cxy,&cxz, &R, E, &R);


        ted_phi_l1l2(&f0,&g0, &cz2,&cxy,&cxz, &eta_Q, &Y_Q, &R);
        fp2_sq1(&f);
        fp2_sq1(&g);
        fp2_mul2(&f,&f0);
        fp2_mul2(&g,&g0);

        if (Q2) {
            ted_phi_l1l2(&f0,&g0, &cz2,&cxy,&cxz, &eta_Q2, &Y_Q2, &R);
            fp2_sq1(&f2);
            fp2_sq1(&g2);
            fp2_mul2(&f2,&f0);
            fp2_mul2(&g2,&g0);
        }

        if (uintbig_bit(k, i)) {
            ted_miller_add(&cz2,&cxy,&cxz, &R, E, &R, P);

            ted_phi_l1l2(&f0,&g0, &cz2,&cxy,&cxz, &eta_Q, &Y_Q, &R);
            fp2_mul2(&f,&f0);
            fp2_mul2(&g,&g0);

            if (Q2) {
                ted_phi_l1l2(&f0,&g0, &cz2,&cxy,&cxz, &eta_Q2, &Y_Q2, &R);
                fp2_mul2(&f2,&f0);
                fp2_mul2(&g2,&g0);
            }
        }

    } while (i--);

    // // testing
    // point test;
    // ted_mul(&test, P, E, k);
    // ted_neg(&test,&test);
    // ted_add(&test, E, &test, &R);
    // assert(ted_iszero(&test));

    if (fp2_iszero(&f) || fp2_iszero(&g)) {
        *res = fp2_0;
    }
    else {
        *res = g;
        fp2_inv(res);
        fp2_mul2(res,&f);
    }

    if (Q2) {
        if (fp2_iszero(&f2) || fp2_iszero(&g2)) {
            *res2 = fp2_0;
        }
        else {
            *res2 = g2;
            fp2_inv(res2);
            fp2_mul2(res2,&f2);
        }
    }
}

void ted_weil(fp2 *res, const proj *E, const point *P, const point *Q, const uintbig *k) {
    fp2 fQT, fPQT, fPT, fQPT;
    point S,R,T;

    ted_neg(&S,Q);
    ted_add(&T, E, P, &S);


    ted_add(&S, E, P, &T);
    ted_miller(&fQT, &fQPT, E, Q, &T, &S, k);
    if (fp2_iszero(&fQT) || fp2_iszero(&fQPT)) { *res = fp2_1; return; }

    //ted_miller(&fQPT, NULL, E, Q, &S, NULL, k);
    //if (fp2_iszero(&fQPT)) { *res = fp2_1; return; }

    ted_neg(&S,&T);
    ted_add(&R, E, Q, &S);
    ted_miller(&fPT, &fPQT, E, P, &S, &R, k);
    if (fp2_iszero(&fPT) || fp2_iszero(&fPQT)) { *res = fp2_1; return; }


    fp2_mul3(res, &fPT, &fQPT);
    fp2_inv(res);
    fp2_mul2(res, &fQT);
    fp2_mul2(res, &fPQT);
}





bool ted_bidim_log_weil(long *a, long *b, const proj *E, const point *Q, const point *P1, const point *P2, long ell) {
    uintbig ell_big;
    uintbig_set(&ell_big, ell);

    fp2 weil_12, weil_Q1, weil_Q2;

    ted_weil(&weil_12, E, P1, P2, &ell_big);
    ted_weil(&weil_Q1, E, Q,  P1, &ell_big);
    ted_weil(&weil_Q2, E, Q,  P2, &ell_big);

    if (!fp2_dlp_naive(a, &weil_Q2, &weil_12, ell)) { return false; }
    if (!fp2_dlp_naive(b, &weil_Q1, &weil_12, ell)) { return false; }
    *b = (ell -  *b) % ell;
    return true;
}

bool ted_bidim_log(GEN *a, GEN *b, const proj *E, const point *Q, const point *P1, const point *P2, long ell, long e) {
    pari_sp ltop = avma;

    uintbig ell_big, x_big, y_big;
    uintbig_set(&ell_big, ell);

    GEN log_1 = gen_0, log_2 = gen_0, ell_pow = gen_1, gtmp;
    long x,y;
    point Ps1[e], Ps2[e], R, Ri, tmp;

    Ps1[0] = *P1;
    Ps2[0] = *P2;

    for (int i = 1; i < e; ++i) {
        ted_mul(&Ps1[i], &Ps1[i-1], E, &ell_big);
        ted_mul(&Ps2[i], &Ps2[i-1], E, &ell_big);
    }

    R = *Q;

    for (int i = 0; i < e; ++i) {
        Ri = R;
        for (int j = 0; j < e-1-i; ++j) {
            ted_mul(&Ri, &Ri, E, &ell_big);
        }
        if(!ted_bidim_log_weil(&x,&y, E, &Ri, &Ps1[e-1], &Ps2[e-1], ell))
            { return false; }

        uintbig_set(&x_big, x);
        ted_mul(&tmp, &Ps1[i], E, &x_big);
        ted_neg(&tmp, &tmp);
        ted_add(&R, E, &R, &tmp);

        uintbig_set(&y_big, y);
        ted_mul(&tmp, &Ps2[i], E, &y_big);
        ted_neg(&tmp, &tmp);
        ted_add(&R, E, &R, &tmp);

        gtmp = mului(x, ell_pow);
        log_1 = gadd(log_1, gtmp);

        gtmp = mului(y, ell_pow);
        log_2 = gadd(log_2, gtmp);

        ell_pow = muliu(ell_pow,ell);

    }

    *a = log_1;
    *b = log_2;

    gerepileall(ltop, 2, a, b);
    return true;
}
