#define _XOPEN_SOURCE

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <pari/pari.h>
#include <math.h>
#include <assert.h>
#include <gmp.h>

#include "ideal.h"
#include "toolbox.h"
#include "klpt.h"

#include "mont.h"
#include "tedwards.h"
#include "constants.h"
#include "precomputed.h"

#define FP_LIMBS (4 * 64 / GMP_LIMB_BITS)




const uintbig cof_p_plus1 = {517434778561ULL,0,0,0};
const uintbig cof_p_plus2 = {26602537156291ULL,0,0,0};

// ted0_twister_frob = B^((p-1)/2) where By^2 = x^3 + x is the quadratic twist of y^2 = x^3 + x
// B is set to be Fp2_inv(fp2_non_residue) (see implementation of twisted edwards, notably xLIFT)
const fp2 ted0_twister_frob = { { 7694077626149142205ULL, 2250429401854817681ULL, 12195039634258678503ULL, 4190125920647857ULL },
  { 14451885770804936219ULL, 10930435212066601406ULL, 15495326356676002808ULL, 5880257503665917138ULL } };



proj torsion_basis[12][3];
proj torsion_basis_sum[3];
point torsion_basis_ted_sum[3];
proj torsion_basis_twist[19][3];
proj torsion_basis_twist_sum[3];
point torsion_basis_twist_ted_sum[3];
proj torsion_basis_two[3];

struct precomp_struct global_setup;



uintbig stobig(long long x) {
    uintbig x_big;
    uintbig_set(&x_big, x);
    return x_big;
}

char* pari_int_code(GEN i) {
    if (is_bigint(i)) return pari_sprintf("strtoi(\"%Ps\")", i);
    else return pari_sprintf("stoi(%PsULL)", i);
}

char* pari_2x2_matrix_code(GEN M) {
    return  pari_sprintf("mkmat2(mkcol2(%s,%s),mkcol2(%s,%s))", pari_int_code(gcoeff(M,1,1)),
            pari_int_code(gcoeff(M,2,1)),
            pari_int_code(gcoeff(M,1,2)),
            pari_int_code(gcoeff(M,2,2)));
}


char* fp_code(const fp *x) {
    return pari_sprintf("{ %luULL, %luULL, %luULL, %luULL }", x->x.c[0], x->x.c[1], x->x.c[2], x->x.c[3]);
}

char* fp2_code(const fp2 *x) {
    return pari_sprintf("{ %s,\n %s }", fp_code(&x->re), fp_code(&x->im));
}

char* proj_code(const proj *P) {
    return pari_sprintf("{ %s,\n %s }", fp2_code(&P->x), fp2_code(&P->z));
}

char* ted_code(const point *P) {
    return pari_sprintf("{ %s,\n %s,\n %s,\n %s }", fp2_code(&P->x), fp2_code(&P->y), fp2_code(&P->z), fp2_code(&P->t));
}

GEN norm0(GEN x) {
    return algnorm(global_setup.B, x,0);
}

void fp2_random_replayable(fp2 *x) {
  uint64_t thrash;
  x->re.x.c[0] = random_Fl(0xffffffffffffffff);
  x->re.x.c[1] = random_Fl(0xffffffffffffffff);
  x->re.x.c[2] = random_Fl(0xffffffffffffffff);
  x->re.x.c[3] = random_Fl(0xffffffffffffffff);
  mpn_tdiv_qr((mp_ptr)&thrash, (mp_ptr)x->re.x.c, 0, (mp_srcptr)x->re.x.c, FP_LIMBS, (mp_srcptr)p.c, FP_LIMBS);
  x->im.x.c[0] = random_Fl(0xffffffffffffffff);
  x->im.x.c[1] = random_Fl(0xffffffffffffffff);
  x->im.x.c[2] = random_Fl(0xffffffffffffffff);
  x->im.x.c[3] = random_Fl(0xffffffffffffffff);
  mpn_tdiv_qr((mp_ptr)&thrash, (mp_ptr)x->im.x.c, 0, (mp_srcptr)x->im.x.c, FP_LIMBS, (mp_srcptr)p.c, FP_LIMBS);
}


void random_point(proj *P, proj const *A, long ell, long e, bool twist) {
    uintbig cofactor;
    uintbig_add3(&cofactor, &p, &uintbig_1);
    uintbig ell_big;
    uintbig_set(&ell_big, ell);
    for (int i = 0; i < e; ++i) {
        uintbig_div3_64(&cofactor, &cofactor, ell); 
    }
    proj Z;

    while (1) {
        fp2_random_replayable(&P->x); fp2_random_replayable(&P->z);
        if (twist == is_on_curve(P, A)) continue;
        xMUL(P, A, P, &cofactor);
        Z = *P;
        for (int i = 0; i < e-1; ++i) {
            xMUL(&Z, A, &Z, &ell_big);
        }
        if (!fp2_iszero(&Z.z)) { 
            // a final test
            xMUL(&Z, A, &Z, &ell_big);
            assert(fp2_iszero(&Z.z));
            return;
        }
    }
}

void random_basis(proj *P1, proj *P2, point *P1_ted, point *P2_ted, proj const *A, long ell, long e, bool twist) {
    point P1_mul, P2_mul, tmp;
    uintbig ell_big;
    uintbig_set(&ell_big, ell);
    fp2 weil;

    proj E;
    mont_to_ted(&E, A, twist);

    random_point(P1, A, ell, e, twist);

    mont_to_ted_point(P1_ted, A, P1);
    P1_mul = *P1_ted;
    for (int i = 0; i < e-1; ++i) {
        ted_mul(&P1_mul, &P1_mul, &E, &ell_big);
    }

    assert(ted_is_on_curve(&P1_mul,&E));
    assert(!ted_iszero(&P1_mul));
    ted_mul(&tmp, &P1_mul, &E, &ell_big);
    assert(ted_iszero(&tmp));

    do {
        random_point(P2, A, ell, e, twist);
        mont_to_ted_point(P2_ted, A, P2);
        P2_mul = *P2_ted;
        for (int i = 0; i < e-1; ++i) {
            ted_mul(&P2_mul, &P2_mul, &E, &ell_big);
        }
        assert(ted_is_on_curve(&P2_mul,&E));
        assert(!ted_iszero(&P2_mul));
        ted_mul(&tmp, &P2_mul, &E, &ell_big);
        assert(ted_iszero(&tmp));

        ted_weil(&weil, &E, &P1_mul, &P2_mul, &ell_big);
        fp2_sub2(&weil, &fp2_1);
    } while (fp2_iszero(&weil));

}

void random_basis_2e(proj *P1, proj *P2, proj const *A, long e, bool twist) {
    proj P1_mul, P2_mul, tmp;
    uintbig ell_big;
    long ell = 2;
    uintbig_set(&ell_big, ell);

    random_point(P1, A, ell, e, twist);

    P1_mul = *P1;
    for (int i = 0; i < e-1; ++i) {
        xMUL(&P1_mul, A, &P1_mul, &ell_big);
    }

    assert(is_on_curve(&P1_mul,A));
    assert(!mont_iszero(&P1_mul));
    xMUL(&tmp, A, &P1_mul, &ell_big);
    assert(mont_iszero(&tmp));

    do {
        random_point(P2, A, ell, e, twist);
        P2_mul = *P2;
        for (int i = 0; i < e-1; ++i) {
            xMUL(&P2_mul, A, &P2_mul, &ell_big);
        }
        assert(is_on_curve(&P2_mul,A));
        assert(!mont_iszero(&P2_mul));
        xMUL(&tmp, A, &P2_mul, &ell_big);
        assert(mont_iszero(&tmp));

    } while (mont_equal(&P1_mul,&P2_mul));

}

bool fp2_ord_trivial(long *res, const fp2 *g, long bound) {
    long order = 1;
    fp2 x = *g;

    for (int i = 0; i < bound; ++i) {
        if (fp2_equal(&fp2_1,&x)) { *res = order; return true; }
        order++;
        fp2_mul2(&x,g);
    }

    return false;
}

// log of Q is base P1,P2, a basis of the 2-torsion
// ASSUMES Q is in the 2-torsion
void bidim_log_2(long long *a, long long *b, const proj *Q, const proj *P1, const proj *P2) {
    if (mont_iszero(Q))         { *a = 0; *b = 0; }
    else if (mont_equal(Q,P1))  { *a = 1; *b = 0; }
    else if (mont_equal(Q,P2))  { *a = 0; *b = 1; }
    else                        { *a = 1; *b = 1; }
}

// Far from optimal
// P12 = P1 + P2
bool bidim_log_2e(long long *a, long long *b, const proj *A, const proj *Q, const proj *P1, const proj *P2, const proj *P12, long e) {
    uintbig big_2, x, y, tmp_big;
    uintbig_set(&big_2, 2);

    long long log_1 = 0, log_2 = 0, tmp;
    proj Ps1[e], Ps2[e], Ps12[e], Qs[e], R;

    Ps1[0] = *P1;
    Ps2[0] = *P2;
    Ps12[0] = *P12;
    Qs[0] = *Q;

    for (int i = 1; i < e; ++i) {
        xMUL(&Ps1[i], A, &Ps1[i-1], &big_2);
        xMUL(&Ps2[i], A, &Ps2[i-1], &big_2);
        xMUL(&Ps12[i], A, &Ps12[i-1], &big_2);
        xMUL(&Qs[i], A, &Qs[i-1], &big_2);
    }

    // first bit
    bidim_log_2(&log_1, &log_2, &Qs[e-1], &Ps1[e-1], &Ps2[e-1]);

    // next bits
    for (int i = 1; i < e; ++i) {

        uintbig_set(&x, log_1);
        uintbig_set(&y, log_2);

        tmp = (1ULL << i);
        uintbig_set(&tmp_big, tmp);

        xBIDIM(&R, A, &Ps1[e-i-1], &x, &Ps2[e-i-1], &y, &Ps12[e-i-1]);
        if (mont_equal(&Qs[e-i-1], &R)) {
            // do nothing
        }
        else {
            uintbig_add3(&x, &x, &tmp_big);
            xBIDIM(&R, A, &Ps1[e-i-1], &x, &Ps2[e-i-1], &y, &Ps12[e-i-1]);
            if (mont_equal(&Qs[e-i-1], &R)) {
                log_1 += tmp;
            }
            else {
                uintbig_set(&x, log_1);
                uintbig_add3(&y, &y, &tmp_big);
                xBIDIM(&R, A, &Ps1[e-i-1], &x, &Ps2[e-i-1], &y, &Ps12[e-i-1]);
                if (mont_equal(&Qs[e-i-1], &R)) {
                    log_2 += tmp;
                }
                else {
                    log_1 += tmp;
                    log_2 += tmp;
                }
            }
                
        }
    }
    *a = log_1;
    *b = log_2;

    // test
    uintbig_set(&x, *a);
    uintbig_set(&y, *b);
    xBIDIM(&R, A, P1, &x, P2, &y, P12);
    assert(mont_equal(Q, &R));
    return true;
}


// distorsion map on E0, twisted Edwards form
void ted0_dist(point *Q, const point *P) {
    point Pcopy = *P;
    fp2_mul3(&Q->x, &Pcopy.t, &fp2_i);
    Q->y = Pcopy.z;
    Q->z = Pcopy.y;
    fp2_mul3(&Q->t, &Pcopy.x, &fp2_i);
}

// distorsion map on E0, montgomery form
void mont0_dist(proj *Q, const proj *P) {
    proj Pcopy = *P;
    fp2_neg2(&Q->x, &Pcopy.x);
    Q->z = Pcopy.z;
}

// distorsion map on E0, montgomery form
void montxy0_dist(proj2 *Q, const proj2 *P) {
    proj2 Pcopy = *P;
    fp2_neg2(&Q->x, &Pcopy.x);
    fp2_mul3(&Q->y, &Pcopy.y, &fp2_i);
    Q->z = Pcopy.z;
}

// frobenius map on E0, montgomery form
void ted0_frob(point *Q, const point *P) {
    point Pcopy = *P;
    fp2_frob2(&Q->x, &Pcopy.x);
    fp2_frob2(&Q->y, &Pcopy.y);
    fp2_frob2(&Q->z, &Pcopy.z);
    fp2_frob2(&Q->t, &Pcopy.t);
}

// frobenius map on the twist of E0, montgomery form
void ted0_frob_twist(point *Q, const point *P) {
    point Pcopy = *P;
    fp2_frob2(&Q->x, &Pcopy.x);
    fp2_frob2(&Q->y, &Pcopy.y);
    fp2_frob2(&Q->z, &Pcopy.z);
    fp2_frob2(&Q->t, &Pcopy.t);
    fp2_mul2(&Q->x, &ted0_twister_frob);
    fp2_mul2(&Q->t, &ted0_twister_frob);
}

// frobenius map on E0, montgomery form
void mont0_frob(proj *Q, const proj *P) {
    proj Pcopy = *P;
    fp2_frob2(&Q->x, &Pcopy.x);
    fp2_frob2(&Q->z, &Pcopy.z);
}

// distorsion map on E0, montgomery form
void montxy0_frob(proj2 *Q, const proj2 *P) {
    proj2 Pcopy = *P;
    fp2_frob2(&Q->x, &Pcopy.x);
    fp2_frob2(&Q->y, &Pcopy.y);
    fp2_frob2(&Q->z, &Pcopy.z);
}


void compute_action(GEN *M, void (*endo)(point*, const point*), const point *P1, const point *P2, const proj *E, long ell, long e) {
    point Q;
    GEN a,b,c,d;
    endo(&Q,P1); 
    assert(ted_is_on_curve(&Q,E));

    assert(ted_bidim_log(&a, &c, E, &Q, P1, P2, ell, e));

    endo(&Q,P2);
    assert(ted_bidim_log(&b, &d, E, &Q, P1, P2, ell, e));

    *M = mkmat2(mkcol2(a,c),mkcol2(b,d));
}

// P12 = P1 + P2
void compute_action_2e(GEN *M, void (*endo)(proj*, const proj*), const proj *P1, const proj *P2, const proj *P12, const proj *A, long e) {
    proj Q, R;
    long long a,b,c,d,x,y;
    uintbig biga, bigc;

    endo(&Q,P1); 
    assert(is_on_curve(&Q,A));

    assert(bidim_log_2e(&a, &c, A, &Q, P1, P2, P12, e));

    uintbig_set(&biga, a);
    uintbig_set(&bigc, c);
    xBIDIM(&R, A, P1, &biga, P2, &bigc, P12);
    assert(mont_equal(&Q, &R));

    endo(&Q,P2);
    assert(bidim_log_2e(&b, &d, A, &Q, P1, P2, P12, e));

    uintbig_set(&biga, b);
    uintbig_set(&bigc, d);
    xBIDIM(&R, A, P1, &biga, P2, &bigc, P12);
    assert(mont_equal(&Q, &R));


    endo(&Q,P12);
    assert(bidim_log_2e(&x, &y, A, &Q, P1, P2, P12, e));

    uintbig_set(&biga, x);
    uintbig_set(&bigc, y);
    xBIDIM(&R, A, P1, &biga, P2, &bigc, P12);
    assert(mont_equal(&Q, &R));

    if (((x != (a+b)%(1LL<<e)) && ((1LL<<e)-x != (a+b)%(1LL<<e))) || ((y != (c+d)%(1LL<<e)) && ((1LL<<e)-y != (c+d)%(1LL<<e)))) {
        b = (1LL<<e)-b;
        d = (1LL<<e)-d;
    }

    assert((x == (a+b)%(1LL<<e)) || ((1LL<<e)-x == (a+b)%(1LL<<e)));
    assert((y == (c+d)%(1LL<<e)) || ((1LL<<e)-y == (c+d)%(1LL<<e)));

    *M = mkmat2(mkcol2(stoi(a),stoi(c)),mkcol2(stoi(b),stoi(d)));
}

GEN action_two_3_4(GEN m_i, GEN m_j, long e) {
    GEN m_ij2,m_1ji2;
    GEN gelle = stoi(1LL<<e);
    GEN gelle1 = stoi(1LL<<(e-1));
    GEN iMi, MM, test1, test2;
    GEN pplus14 = gdiv(gadd(global_setup.p,gen_1),stoi(4));
    GEN id = mkmat2(mkcol2s(1,0),mkcol2s(0,1));


    for (int i11 = 0; i11 < 2; ++i11){
        for (int i12 = 0; i12 < 2; ++i12){
            for (int i21 = 0; i21 < 2; ++i21){
                for (int i22 = 0; i22 < 2; ++i22){
        
                    m_ij2 = gdiv(gadd(m_i, m_j),gen_2);
                    if (i11) gcoeff(m_ij2,1,1) = gadd(gcoeff(m_ij2,1,1),gelle1);
                    if (i12) gcoeff(m_ij2,1,2) = gadd(gcoeff(m_ij2,1,2),gelle1);
                    if (i21) gcoeff(m_ij2,2,1) = gadd(gcoeff(m_ij2,2,1),gelle1);
                    if (i22) gcoeff(m_ij2,2,2) = gadd(gcoeff(m_ij2,2,2),gelle1);

                    iMi = gmod(gmul(gmul(m_i,m_ij2),m_i),gelle);
                    test1 = gmod(gsub(m_ij2,m_i),gelle);

                    MM = gmod(gmul(m_ij2,m_ij2),gelle);
                    test2 = gmod(gneg(gmul(id,pplus14)),gelle);



                    if (gequal(iMi,test1) && gequal(MM,test2)) { // a candidate
                        m_1ji2 = gmod(gneg(gmul(m_ij2,m_i)),gelle);



                        printf("\t/* candidate %d%d%d%d */\n", i11,i12,i21,i22);
                        printf("\taction_two_3 = %s;\n", pari_2x2_matrix_code(m_1ji2));
                        printf("\taction_two_4 = %s;\n", pari_2x2_matrix_code(m_ij2));
                    }
                }
            }
        }
    }


    return NULL;
}

void check_action(GEN M, void (*endo)(point*, const point*), void (*mont_endo)(proj*, const proj*), const point *P1, const point *P2, const proj *E, const proj *E_mont, long ell, long e) {
    long x,y,u,v;
    uintbig X,Y,U,V;
    point Q, endoQ, endoQ_combination, tmp;

    point P12;
    proj M1,M2,M12,MQ,MendoQ,MendoQ_combination;

    ted_to_mont_point(&M1, P1);
    ted_to_mont_point(&M2, P2);
    ted_add(&P12, E, P1, P2);
    ted_to_mont_point(&M12, &P12);

    if (ell == 3) return; // TODO: ell^e does not fit in a long for ell = 3

    for (int i = 0; i < 10; ++i) {    
        x = random_Fl(ell);
        y = random_Fl(ell);
        uintbig_set(&X, x);
        uintbig_set(&Y, y);

        ted_mul(&Q, P1, E, &X);
        ted_mul(&tmp, P2, E, &Y);
        ted_add(&Q, E, &Q, &tmp);

        endo(&endoQ,&Q);

        long A,B,C,D;
        A = itos_or_0(gcoeff(M,1,1));
        B = itos_or_0(gcoeff(M,1,2));
        C = itos_or_0(gcoeff(M,2,1));
        D = itos_or_0(gcoeff(M,2,2));

        u = (A*x+B*y);
        v = (C*x+D*y);
        uintbig_set(&U, u);
        uintbig_set(&V, v);

        ted_mul(&endoQ_combination, P1, E, &U);
        ted_mul(&tmp, P2, E, &V);
        ted_add(&endoQ_combination, E, &endoQ_combination, &tmp);

        assert(ted_equal(&endoQ,&endoQ_combination));


        // same check in montgomery form
        ted_to_mont_point(&MQ, &Q);
        mont_endo(&MendoQ,&MQ);

        xBIDIM(&MendoQ_combination, E_mont, &M1, &U, &M2, &V, &M12);

        assert(mont_equal(&MendoQ,&MendoQ_combination));
    }
}

bool check_action_2e(GEN M, void (*endo)(proj*, const proj*), void (*endoxy)(proj2*, const proj2*), const proj2 *P1xy, const proj2 *P2xy, const proj *E, long e) {
    long long x,y;
    uintbig X,Y,U,V;
    proj2 Qxy,Rxy,exy1,exy2,exy12,eQxy_combination,eQxy;

    proj2 P12xy;

    GEN t;

    for (int i = 0; i < 20; ++i) {  
        x = random_Fl(1ULL<<e);
        y = random_Fl(1ULL<<e);
        uintbig_set(&X, x);
        uintbig_set(&Y, y);


        endoxy(&exy1,P1xy);
        endoxy(&exy2,P2xy);
        endoxy(&exy12,&P12xy);


        xyMUL(&Qxy, E, P1xy, &X);
        xyMUL(&Rxy, E, P2xy, &Y);
        xyADD(&Qxy, E, &Qxy, &Rxy);

        endoxy(&eQxy,&Qxy);

        t = gmod(RgM_RgC_mul(M, mkcol2s(x,y)), stoi(1ULL<<e));
        uintbig_set(&U, itos_or_0(gel(t,1)));
        uintbig_set(&V, itos_or_0(gel(t,2)));

        xyMUL(&eQxy_combination, E, P1xy, &U);
        xyMUL(&Rxy, E, P2xy, &V);
        xyADD(&eQxy_combination, E, &eQxy_combination, &Rxy);

        proj2 neg;
        xyNEG(&neg, &eQxy);

        assert(xy_equal(&eQxy,&eQxy_combination) || xy_equal(&neg,&eQxy_combination));

        if (xy_equal(&neg,&eQxy_combination)) return false; // flip sign!
    }
    return true;
}

static void gentobig(uintbig *res, GEN a) {
    pari_sp ltop = avma;
    GEN b;
    res->c[0] = umodi2n(a,32);
    b = shifti(a, -32);
    res->c[0] += (umodi2n(b,32) << 32);
    b = shifti(b, -32);
    res->c[1] = umodi2n(b,32);
    b = shifti(b, -32);
    res->c[1] += (umodi2n(b,32) << 32);
    b = shifti(b, -32);
    res->c[2] = umodi2n(b,32);
    b = shifti(b, -32);
    res->c[2] += (umodi2n(b,32) << 32);
    b = shifti(b, -32);
    res->c[3] = umodi2n(b,32);
    b = shifti(b, -32);
    res->c[3] += (umodi2n(b,32) << 32);
    b = shifti(b, -32);
    avma = ltop;
}

// argv[1] is the random seed; default = 1
int main(int argc, char *argv[]){
    pari_init(80000000, 1<<18);

    setrand(stoi(1));
    srand48(1);
    if( argc > 1 ) {
      setrand(strtoi(argv[1]));
      srand48(atoi(argv[1]));
    }

    long var = fetch_var();
    GEN nf = nfinit(pol_x(fetch_var()),LOWDEFAULTPREC);
    
    GEN a = stoi(-1),
        p_pari = strtoi("73743043621499797449074820543863456997944695372324032511999999999999999999999"),
        b = negi(p_pari);

    GEN B = alg_hilbert(nf, a, b, var, 0);

    GEN B_1 = mkcol4s(1,0,0,0);
    GEN B_i = mkcol4s(0,1,0,0);
    GEN B_j = mkcol4s(0,0,1,0);
    GEN B_ji = mkcol4s(0,0,0,1);
    //GEN B_ij = mkcol4s(0,0,0,-1);

    GEN B_1k_2 = mkcol4(ghalf,gen_0,gen_0,gneg(ghalf)); // (1-ji)/2
    GEN B_ij_2 = mkcol4(gen_0,ghalf,ghalf,gen_0); // (i+j)/2

    GEN B_O0 = alglathnf(B,mkmat4(B_1, B_i, B_1k_2, B_ij_2), gen_0);

    global_setup.p = p_pari;
    global_setup.B = B; // the quaternion algebra
    global_setup.qf = mkmat4(mkcol4s(1,0,0,0),
                             mkcol4s(0,1,0,0),
                             mkcol4(gen_0,gen_0,p_pari,gen_0),
                             mkcol4(gen_0,gen_0,gen_0,p_pari)); // quadratic form defined by the reduced norm

    global_setup.torsion_fm = Z_factor_limit(strtoi(
        "197530174297949459837634878151545563369632855190375548677707409417459236752253845947265965991865263091519488000000000000000000000"
        ), 30000);

    global_setup.O0 = B_O0; // the cannonical maximal order
    global_setup.one = B_1;
    global_setup.i = B_i;
    global_setup.j = B_j;
    global_setup.ji = B_ji;

    global_setup.O0_b1 = B_1;
    global_setup.O0_b2 = B_i;
    global_setup.O0_b3 = B_1k_2;
    global_setup.O0_b4 = B_ij_2;
    global_setup.O0_to_standard = mkmat4(B_1, B_i, B_1k_2, B_ij_2);
    global_setup.standard_to_O0 = RgM_inv(global_setup.O0_to_standard);

    global_setup.E0.x = fp2_0;
    global_setup.E0.z = fp2_1;






    // ted0_twister_frob = B^((p-1)/2) where By^2 = x^3 + x is the quadratic twist of y^2 = x^3 + x
    // B is set to be Fp2_inv(fp2_non_residue) (see implementation of twisted edwards, notably xLIFT)
    fp2 ted0_twister_frob_ = fp2_non_residue();
    uintbig p12;
    uintbig_set(&p12,1);
    uintbig_sub3(&p12, &p, &p12);
    uintbig_div3_64(&p12, &p12, 2);
    fp2_exp(&ted0_twister_frob_, &ted0_twister_frob_, &p12);

    assert(fp2_equal(&ted0_twister_frob,&ted0_twister_frob_));










    proj E, E_twist;
    mont_to_ted(&E, &global_setup.E0, false);
    mont_to_ted(&E_twist, &global_setup.E0, true);

    GEN m_frob, m_dist;
    GEN gelle, inv2;

    proj basis[p_plus_len][3];
    point basis_ted[p_plus_len][3];

    GEN m1 = mkmat2(mkcol2s(1,0),mkcol2s(0,1));
    GEN action_2[p_plus_len],action_3[p_plus_len],action_4[p_plus_len];

    for (int i = 0; i < p_plus_len; i++) {
        long ell = p_plus_fact[i], e = p_plus_mult[i];

        gelle = powuu(ell,e);
        inv2 = Fp_inv(gen_2, gelle);

        random_basis(&basis[i][0], &basis[i][1], &basis_ted[i][0], &basis_ted[i][1], &global_setup.E0, ell, e, false);

        ted_add(&basis_ted[i][2], &E, &basis_ted[i][0], &basis_ted[i][1]);
        ted_to_mont_point(&basis[i][2], &basis_ted[i][2]);


        point test;
        proj test2, test3;

        // TEST that mont_add works properly
        mont_add(&test2, &global_setup.E0, &basis[i][0], &basis[i][1]);
        ted_neg(&test, &basis_ted[i][1]);
        ted_add(&test, &E, &basis_ted[i][0], &test);
        ted_to_mont_point(&test3, &test);
        assert(mont_equal(&test2,&basis[i][2]) || mont_equal(&test2,&test3));
        // END TEST



        compute_action(&m_dist, ted0_dist, &basis_ted[i][0], &basis_ted[i][1], &E, ell, e);
        check_action(m_dist, ted0_dist, mont0_dist, &basis_ted[i][0], &basis_ted[i][1], &E, &global_setup.E0, ell, e);

        compute_action(&m_frob, ted0_frob, &basis_ted[i][0], &basis_ted[i][1], &E, ell, e);
        check_action(m_frob, ted0_frob, mont0_frob, &basis_ted[i][0], &basis_ted[i][1], &E, &global_setup.E0, ell, e);

        action_2[i] = m_dist;
        action_3[i] = gmod(gmul(gsub(m1,gmul(m_frob,m_dist)),inv2),gelle); // (1-ji)/2
        action_4[i] = gmod(gmul(gadd(m_dist,m_frob),inv2),gelle); //(i+j)/2



    }


    proj basis_twist[p_minus_len][3];
    point basis_twist_ted[p_minus_len][3];

    GEN action_twist_2[p_minus_len],action_twist_3[p_minus_len],action_twist_4[p_minus_len];

    for (int i = 0; i < p_minus_len; i++) {
        long ell = p_minus_fact[i], e = p_minus_mult[i];

        gelle = powuu(ell,e);
        inv2 = Fp_inv(gen_2, gelle);

        random_basis(&basis_twist[i][0], &basis_twist[i][1], &basis_twist_ted[i][0], &basis_twist_ted[i][1], &global_setup.E0, ell, e, true);


        point test;
        proj test2;
        uintbig ell_big;
        gentobig(&ell_big, powuu(ell,e));
        ted_mul(&test, &basis_twist_ted[i][0], &E_twist, &ell_big);
        assert(ted_iszero(&test));
        ted_mul(&test, &basis_twist_ted[i][1], &E_twist, &ell_big);
        assert(ted_iszero(&test));

        ted_to_mont_point(&test2, &basis_twist_ted[i][0]);
        assert(!is_on_curve(&test2, &global_setup.E0));
        xMUL(&test2, &global_setup.E0, &test2, &ell_big);
        assert(mont_iszero(&test2));

        ted_to_mont_point(&test2, &basis_twist_ted[i][1]);
        assert(!is_on_curve(&test2, &global_setup.E0));
        xMUL(&test2, &global_setup.E0, &test2, &ell_big);
        assert(mont_iszero(&test2));


        ted_add(&basis_twist_ted[i][2], &E_twist, &basis_twist_ted[i][0], &basis_twist_ted[i][1]);
        ted_to_mont_point(&basis_twist[i][2], &basis_twist_ted[i][2]);

        ted_mul(&test, &basis_twist_ted[i][2], &E_twist, &ell_big);
        assert(ted_iszero(&test));
        ted_to_mont_point(&test2, &basis_twist_ted[i][2]);
        assert(!is_on_curve(&test2, &global_setup.E0));
        xMUL(&test2, &global_setup.E0, &test2, &ell_big);
        assert(mont_iszero(&test2));
        
        compute_action(&m_dist, ted0_dist, &basis_twist_ted[i][0], &basis_twist_ted[i][1], &E_twist, ell, e);
        check_action(m_dist, ted0_dist, mont0_dist, &basis_twist_ted[i][0], &basis_twist_ted[i][1], &E_twist, &global_setup.E0, ell, e);

        compute_action(&m_frob, ted0_frob_twist, &basis_twist_ted[i][0], &basis_twist_ted[i][1], &E_twist, ell, e);
        check_action(m_frob,       ted0_frob_twist, mont0_frob, &basis_twist_ted[i][0], &basis_twist_ted[i][1], &E_twist, &global_setup.E0, ell, e);


        action_twist_2[i] = m_dist;
        action_twist_3[i] = gmod(gmul(gsub(m1,gmul(m_frob,m_dist)),inv2),gelle); // (1-ji)/2
        action_twist_4[i] = gmod(gmul(gadd(m_dist,m_frob),inv2),gelle); //(i+j)/2
    }


    point basis_ted_sum[3], basis_twist_ted_sum[3];
    proj basis_sum[3], basis_twist_sum[3];

    ted_add(&basis_ted_sum[0], &E, &basis_ted[0][0], &basis_ted[1][0]);
    ted_add(&basis_ted_sum[1], &E, &basis_ted[0][1], &basis_ted[1][1]);
    for (int i = 2; i < p_plus_len; i++) {
        ted_add(&basis_ted_sum[0], &E, &basis_ted_sum[0], &basis_ted[i][0]);
        ted_add(&basis_ted_sum[1], &E, &basis_ted_sum[1], &basis_ted[i][1]);
    }
    ted_add(&basis_ted_sum[2], &E, &basis_ted_sum[0], &basis_ted_sum[1]);
    ted_to_mont_point(&basis_sum[0], &basis_ted_sum[0]);
    ted_to_mont_point(&basis_sum[1], &basis_ted_sum[1]);
    ted_to_mont_point(&basis_sum[2], &basis_ted_sum[2]);

    ted_add(&basis_twist_ted_sum[0], &E_twist, &basis_twist_ted[0][0], &basis_twist_ted[1][0]);
    ted_add(&basis_twist_ted_sum[1], &E_twist, &basis_twist_ted[0][1], &basis_twist_ted[1][1]);
    for (int i = 2; i < p_minus_len; i++) {
        ted_add(&basis_twist_ted_sum[0], &E_twist, &basis_twist_ted_sum[0], &basis_twist_ted[i][0]);
        ted_add(&basis_twist_ted_sum[1], &E_twist, &basis_twist_ted_sum[1], &basis_twist_ted[i][1]);
    }
    ted_add(&basis_twist_ted_sum[2], &E_twist, &basis_twist_ted_sum[0], &basis_twist_ted_sum[1]);
    ted_to_mont_point(&basis_twist_sum[0], &basis_twist_ted_sum[0]);
    ted_to_mont_point(&basis_twist_sum[1], &basis_twist_ted_sum[1]);
    ted_to_mont_point(&basis_twist_sum[2], &basis_twist_ted_sum[2]);


    printf("const proj torsion_basis_sum[3] = \n");
    printf("{%s,\n %s,\n %s };\n", proj_code(&basis_sum[0]), proj_code(&basis_sum[1]), proj_code(&basis_sum[2]));

    printf("const proj torsion_basis_twist_sum[3] = \n");
    printf("{%s,\n %s,\n %s };\n", proj_code(&basis_twist_sum[0]), proj_code(&basis_twist_sum[1]), proj_code(&basis_twist_sum[2]));

    printf("\n");




    printf("const point torsion_basis_ted_sum[3] = \n");
    printf("{%s,\n %s,\n %s };\n", ted_code(&basis_ted_sum[0]), ted_code(&basis_ted_sum[1]), ted_code(&basis_ted_sum[2]));

    printf("const point torsion_basis_twist_ted_sum[3] = \n");
    printf("{%s,\n %s,\n %s };\n", ted_code(&basis_twist_ted_sum[0]), ted_code(&basis_twist_ted_sum[1]), ted_code(&basis_twist_ted_sum[2]));

    printf("\n");


    printf("void init_action() {\n");
    for (int i = 0; i < p_plus_len; i++) {
        printf("\tglobal_setup.action_2[%d] = %s;\n", i, pari_2x2_matrix_code(action_2[i]));
        printf("\tglobal_setup.action_3[%d] = %s;\n", i, pari_2x2_matrix_code(action_3[i]));
        printf("\tglobal_setup.action_4[%d] = %s;\n", i, pari_2x2_matrix_code(action_4[i]));
    }
    for (int i = 0; i < p_minus_len; i++) {
        printf("\tglobal_setup.action_twist_2[%d] = %s;\n", i, pari_2x2_matrix_code(action_twist_2[i]));
        printf("\tglobal_setup.action_twist_3[%d] = %s;\n", i, pari_2x2_matrix_code(action_twist_3[i]));
        printf("\tglobal_setup.action_twist_4[%d] = %s;\n", i, pari_2x2_matrix_code(action_twist_4[i]));
    }



    // two-torsion
    
    proj basis_two[3];
    random_basis_2e(&basis_two[0], &basis_two[1], &global_setup.E0, two_tors_height, false);

    mont_add(&basis_two[2], &global_setup.E0, &basis_two[0], &basis_two[1]);
    assert(is_on_curve(&basis_two[2], &global_setup.E0));




    proj2 basisxy_two[3];

    xtoxy(&basisxy_two[0],&global_setup.E0,&basis_two[0]);
    xtoxy(&basisxy_two[1],&global_setup.E0,&basis_two[1]);

    xyADD(&basisxy_two[2], &global_setup.E0, &basisxy_two[0], &basisxy_two[1]);

    proj pt;
    xytox(&pt, &basisxy_two[2]);

    assert(mont_equal(&pt, &basis_two[2]));



    printf("const proj torsion_basis_two[3] = \n");
    printf("{%s,\n %s,\n %s };\n", proj_code(&basis_two[0]), proj_code(&basis_two[1]), proj_code(&basis_two[2]));

    bool correct_sign;

    compute_action_2e(&m_dist, mont0_dist, &basis_two[0], &basis_two[1], &basis_two[2], &global_setup.E0, two_tors_height);
    do {
        correct_sign = check_action_2e(m_dist, mont0_dist, montxy0_dist, &basisxy_two[0], &basisxy_two[1], &global_setup.E0, two_tors_height);
        if (!correct_sign) { m_dist = gmod(gneg(m_dist), powuu(2,two_tors_height));}
    } while (!correct_sign);

    compute_action_2e(&m_frob, mont0_frob, &basis_two[0], &basis_two[1], &basis_two[2], &global_setup.E0, two_tors_height);
    do {
        correct_sign = check_action_2e(m_frob, mont0_frob, montxy0_frob, &basisxy_two[0], &basisxy_two[1], &global_setup.E0, two_tors_height);
        if (!correct_sign) { m_frob = gmod(gneg(m_frob), powuu(2,two_tors_height));}
    } while (!correct_sign);

    printf("\tglobal_setup.action_two_2 = %s;\n", pari_2x2_matrix_code(m_dist));
    action_two_3_4(m_dist, m_frob, two_tors_height);

    printf("}\n");



    return 0;
}



