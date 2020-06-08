#define _XOPEN_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <pari/pari.h>
#include <math.h>
#include <assert.h>
#include <gmp.h>

#define FP_LIMBS (4 * 64 / GMP_LIMB_BITS)


#include "ideal.h"
#include "toolbox.h"
#include "klpt.h"

#include "mont.h"
#include "idiso.h"
#include "constants.h"
#include "precomputed.h"

/*
static char* fp2_hash(fp2 x) {
    return pari_sprintf("h%lu", (x.re.x.c[0]+3*x.re.x.c[1]+5*x.re.x.c[2]+7*x.re.x.c[3]
                    +11*x.im.x.c[0]+13*x.im.x.c[1]+17*x.im.x.c[2]+23*x.im.x.c[3]) % 100003);
}
static fp2 fp2_ratio(fp2 *x, fp2 *y) {
    fp2 tmp;
    tmp = *y;
    assert(!fp2_iszero(&tmp));
    fp2_inv(&tmp);
    fp2_mul2(&tmp, x);
    //printf("%lld %lld %lld %lld\n", tmp.re.x.c[0], tmp.re.x.c[1], tmp.im.x.c[0], tmp.im.x.c[1]);
    return tmp;
}
static char* proj_hash(proj x) {
    return fp2_hash(fp2_ratio(&x.x,&x.z));
}



static GEN norm0(GEN x) {
    return algnorm(global_setup.B, x,0);
}
*/

int test_kertoid_action() {

    float accumulated_time_ms = 0.;
    int repetitions = 1000;
    clock_t t;
    long ell, e;
    GEN gelle, m_1, m_i, m_j, m_ji, m_2, m_3, m_4, v, v1, v2, w;
    GEN ideal;

    ell = 6983;
    e = 1;
    gelle = gpowgs(stoi(ell),e);

    m_1 = mkmat2(mkcol2s(1,0), mkcol2s(0,1));
    m_i = mkmat2(mkcol2s(0,1), mkcol2s(-1,0));
    m_j = mkmat2(mkcol2s(1296,525), mkcol2s(525,-1296));
    m_ji = gmul(m_j, m_i);

 // B_1k_2, B_ij_2
    m_2 = m_i;
    m_3 = gmod(gmul(gsub(m_1, m_ji), Fp_inv(gen_2,gelle)),gelle); // (1 - ji) / 2
    m_4 = gmod(gmul(gadd(m_i, m_j), Fp_inv(gen_2,gelle)),gelle); // (i + j) / 2


    ell = 4283;
    e = 1;
    gelle = gpowgs(stoi(ell),e);
    m_2 = mkmat2(mkcol2s(-1298618089,2597427807), mkcol2s(-2442345774,1298618089));
    m_3 = mkmat2(mkcol2s(3594923693,3500346076), mkcol2s(-3142991081,-3594923692));
    m_4 = mkmat2(mkcol2s(-3020381676,2095475696), mkcol2s(-348475943,3020381676));


    t = tic();
    for (int i = 0; i < repetitions; i++) {
        do { 
            v1 = randomi(gelle);
            v2 = randomi(gelle);
        } while (!(umodiu(v1,ell) || umodiu(v2,ell)));

        v = mkcol2(v1,v2);

        t = tic();
        ideal = kernel_to_ideal_action_O0(v, m_1, m_2, m_3, m_4, ell, e);
        w = ideal_to_kernel_action_O0(ideal, m_1, m_2, m_3, m_4, ell, e);
        accumulated_time_ms += toc(t);

        // output(gmod(v,stoi(ell)));
        // output(gmod(w,stoi(ell)));
        assert(gcmp(gmod(QM_det(mkmat2(v,w)),gelle),gen_0) == 0);

    }
    printf("average kertoid\t [%f ms]\n",  (accumulated_time_ms / repetitions));


    // printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));

    return 0;
}

int test_kertoid_O0_ell() {
    float accumulated_time_ms = 0.;
    int repetitions = 1000;
    clock_t t;
    long ell, e;
    GEN I, v = mkcol2s(0,0), w, gelle;

    for (int i = 0; i < repetitions; i++) {
        if (i % 2) {
            ell = p_plus_fact[(i>>1)%p_plus_len];
            e = p_plus_mult[(i>>1)%p_plus_len];
        }
        else {
            ell = p_minus_fact[(i>>1)%p_minus_len];
            e = p_minus_mult[(i>>1)%p_minus_len];
        }

        gelle = powuu(ell,e);

        gel(v,1) = randomi(gelle);
        gel(v,2) = randomi(gelle);

        t = tic();
        I = kernel_to_ideal_O0_ell(v, ell);
        w = ideal_to_kernel_O0_ell(I, ell);
        accumulated_time_ms += toc(t);

        GEN gcd1 = ggcd(gelle,ggcd(gel(v,1),gel(v,2)));
        GEN gcd2 = ggcd(gelle,ggcd(gel(w,1),gel(w,2)));
        assert(gcmp(gcd1,gcd2) == 0);

        assert(gcmp(gmod(QM_det(mkmat2(v,w)),gelle),gen_0) == 0);
    }

    printf("average kertoid_ell\t [%f ms]\n",  (accumulated_time_ms / repetitions));
    return 1;
}

int test_kertoid_O0_two() {
    float accumulated_time_ms = 0.;
    int repetitions = 1000;
    clock_t t;
    long ell = 2, e = two_tors_height;
    GEN I, v = mkcol2s(0,0), w, gelle = powuu(ell,e);

    for (int i = 0; i < repetitions; i++) {


        gel(v,1) = randomi(gelle);
        gel(v,2) = randomi(gelle);

        t = tic();
        I = kernel_to_ideal_O0_ell(v, ell);
        w = ideal_to_kernel_O0_ell(I, ell);
        accumulated_time_ms += toc(t);


        GEN gcd1 = ggcd(gelle,ggcd(gel(v,1),gel(v,2)));
        GEN gcd2 = ggcd(gelle,ggcd(gel(w,1),gel(w,2)));
        assert(gcmp(gcd1,gcd2) == 0);


        assert(gcmp(gmod(QM_det(mkmat2(v,w)),gelle),gen_0) == 0);
    }

    printf("average kertoid_two\t [%f ms]\n",  (accumulated_time_ms / repetitions));
    return 1;
}

int test_kertoid_O0_T() {
    float accumulated_time_kertoid = 0., accumulated_time_idtoker = 0.;
    int repetitions = 2;
    clock_t t;
    GEN I, w, gelle, fact_norm;

    GEN coeff_1 = cgetg(p_plus_len+1, t_VEC), coeff_2 = cgetg(p_plus_len+1, t_VEC);
    GEN coeff_twist_1 = cgetg(p_minus_len+1, t_VEC), coeff_twist_2 = cgetg(p_minus_len+1, t_VEC);

    odd_isogeny phi;
    proj P;

    pari_sp ltop = avma;
    for (int i = 0; i < repetitions; i++) {

        for (int i = 0; i < p_plus_len; ++i) {
            gelle = powuu(p_plus_fact[i],p_plus_mult[i]);
            gel(coeff_1,i+1) = randomi(gelle);
            gel(coeff_2,i+1) = randomi(gelle);
        }
        for (int i = 0; i < p_minus_len; ++i) {
            gelle = powuu(p_minus_fact[i],p_minus_mult[i]);
            gel(coeff_twist_1,i+1) = randomi(gelle);
            gel(coeff_twist_2,i+1) = randomi(gelle);
        }
        t = tic();
        I = kernel_to_ideal_O0_T(mkvec2(mkvec2(coeff_1,coeff_2),mkvec2(coeff_twist_1,coeff_twist_2)));
        accumulated_time_kertoid += toc(t);

        fact_norm = Z_factor(lideal_norm(I));

        w = ideal_to_kernel_O0_T(I, fact_norm);
        GEN comp1 = torsion_crt_compose (gel(w,1), false);
        gel(w,1) = torsion_crt_decompose (comp1, false);
        GEN comp2 = torsion_crt_compose (gel(w,2), true);
        gel(w,2) = torsion_crt_decompose (comp2, true);


        t = tic();

        phi = ideal_to_isogeny_O0_T(I, fact_norm);

        fp2_random(&P.x); fp2_random(&P.z);
        proj A = global_setup.E0;
        eval(&A, &phi, &P);
        A = global_setup.E0;
        dual(&A, &phi);

        accumulated_time_idtoker += toc(t);

        for (int i = 0; i < p_plus_len; ++i) {
            gelle = powuu(p_plus_fact[i],p_plus_mult[i]);


            GEN c1 = mkcol2(gel(coeff_1,i+1),gel(coeff_2,i+1));
            GEN c2 = mkcol2(gel(gel(gel(w,1),1),i+1),gel(gel(gel(w,1),2),i+1));
            GEN gcd1 = ggcd(gelle,ggcd(gel(c1,1),gel(c1,2)));
            GEN gcd2 = ggcd(gelle,ggcd(gel(c2,1),gel(c2,2)));
            GEN M = mkmat2(c1,c2);
            assert(gcmp(gcd1,gcd2) == 0);
            assert(gcmp(gmod(QM_det(M),gelle),gen_0) == 0);
        }
        for (int i = 0; i < p_minus_len; ++i) {
            gelle = powuu(p_minus_fact[i],p_minus_mult[i]);
            GEN c1 = mkcol2(gel(coeff_twist_1,i+1),gel(coeff_twist_2,i+1));
            GEN c2 = mkcol2(gel(gel(gel(w,2),1),i+1),gel(gel(gel(w,2),2),i+1));
            GEN gcd1 = ggcd(gelle,ggcd(gel(c1,1),gel(c1,2)));
            GEN gcd2 = ggcd(gelle,ggcd(gel(c2,1),gel(c2,2)));
            GEN M = mkmat2(c1,c2);
            assert(gcmp(gcd1,gcd2) == 0);
            assert(gcmp(gmod(QM_det(M),gelle),gen_0) == 0);
        }

        avma = ltop;
    }

    //TOC(t, "kertoid_T");
    printf("average kertoid\t [%f ms]\n",  (accumulated_time_kertoid / repetitions));
    printf("average idtoker\t [%f ms]\n",  (accumulated_time_idtoker / repetitions));

    return 1;
}


void random_odd_isogeny(odd_isogeny *phi){
    pari_sp ltop = avma;

    GEN coeff_plus_1 = cgetg(p_plus_len+1, t_VEC);
    GEN coeff_plus_2 = cgetg(p_plus_len+1, t_VEC);
    GEN coeff_minus_1 = cgetg(p_minus_len+1, t_VEC);
    GEN coeff_minus_2 = cgetg(p_minus_len+1, t_VEC);

    for (int i = 0; i < p_plus_len; ++i) {
        gel(coeff_plus_1,i+1) = gen_1;
        gel(coeff_plus_2,i+1) = randomi(powuu(p_plus_fact[i],p_plus_mult_com[i]));
    }

    for (int i = 0; i < p_minus_len; ++i) {
        gel(coeff_minus_1,i+1) = gen_1;
        gel(coeff_minus_2,i+1) = randomi(powuu(p_minus_fact[i],p_minus_mult_com[i]));
    }

    GEN coeff_plus = mkvec2(coeff_plus_1, coeff_plus_2);
    GEN coeff_minus = mkvec2(coeff_minus_1, coeff_minus_2);
    GEN coeff = mkvec2(coeff_plus, coeff_minus);

    GEN I = kernel_to_ideal_O0_T(coeff);

    proj ker = coeff_to_E0(gel(coeff,1), false);
    proj ker_twist = coeff_to_E0(gel(coeff,2), true);


    isog_degree deg, deg_twist;

    //deg = mont_order(&ker, &global_setup.E0, false);
    //deg_twist = mont_order(&ker_twist, &global_setup.E0, true);


    famat_to_degree(&deg, &deg_twist, Z_factor(lideal_norm(I)));

    phi->kernel_plus = ker;
    phi->kernel_minus = ker_twist;
    phi->deg_plus = deg;
    phi->deg_minus = deg_twist;
    avma = ltop;
}

void random_two_walk(two_walk *phi){
    pari_sp ltop = avma;
    
    phi->A = global_setup.E0;

    phi->len = two_tors_height;

    const proj *basis = torsion_basis_two;
    uintbig x, y;
    proj P;

    do {
        gentobig(&x, gen_1);
        gentobig(&y, randomi(powuu(2,two_tors_height)));

        xBIDIM(&phi->ker, &global_setup.E0, &basis[0], &x, &basis[1], &y, &basis[2]);

        gentobig(&x, powuu(2,two_tors_height-1));
        xMUL(&P, &global_setup.E0, &phi->ker, &x);
        assert(!mont_iszero(&P));
    } while (!fp2_iszero(&P.x));

    gentobig(&x, gen_2);
    xMUL(&P, &global_setup.E0, &P, &x);
    assert(mont_iszero(&P));

    avma = ltop;
}


void random_point(proj *P, proj const *A, bool twist) {
    do {
        fp2_random(&P->x); P->z = fp2_1;
    } while (twist == is_on_curve(P, A));
}



isog_degree mont_order(const proj *P, const proj *A, bool twist) {
  
    const long *fact, *mult;
    long len;
      

    if(!twist) {
        fact = p_plus_fact; mult = p_plus_mult;
        len = p_plus_len;
    }
    else {
        fact = p_minus_fact; mult = p_minus_mult;
        len = p_minus_len;
    }

  proj tmp;
  uintbig cof;
  isog_degree deg = degree_co((isog_degree){ 0 }, mult, len);
  for (int j = 0; j < len; j++) {
    degree_unset(&deg, j);
    degree_to_uint(&cof, deg, fact, len);
    xMUL(&tmp, A, P, &cof);
    uintbig_set(&cof, fact[j]);
    uint8_t v = 0;
    for ( ; !mont_iszero(&tmp); v++) {
      xMUL(&tmp, A, &tmp, &cof);
    }
    degree_set(&deg, j, v);
  }
  return deg;
}

int test_push_isogenies() {
    float accumulated_time = 0.;
    int repetitions = 10;
    clock_t t;

    odd_isogeny phi_odd1, phi_odd2;
    proj phi_odd_source1, phi_odd_source2, P, A1, A2, P1, P2;
    two_walk phi_two1, phi_two2;

    for (int i = 0; i < repetitions; i++) {
        random_odd_isogeny(&phi_odd1);
        phi_odd2 = phi_odd1;
        phi_odd_source1 = global_setup.E0;
        phi_odd_source2 = global_setup.E0;
        random_two_walk(&phi_two1);
        phi_two2 = phi_two1;

        t = tic();

        phi_odd1 = push_odd_isogeny_through_two_walk(&phi_odd1, &phi_odd_source1, &phi_two1);
        
        phi_two2 = push_two_walk_through_odd_isogeny(&phi_two2, &phi_odd2, &phi_odd_source2);

        accumulated_time += toc(t);

        random_point(&P, &global_setup.E0, false);

        P1 = P;
        A1 = global_setup.E0;

        isomorphism isom;

        eval_walk_isom(&isom, &phi_two1, &A1, &P1, &phi_two1, &P1);

        assert(mont_equal(&A1,&phi_odd_source1));
        eval(&A1, &phi_odd1, &P1);

        P2 = P;
        A2 = global_setup.E0;
        eval(&A2, &phi_odd2, &P2);
        assert(mont_equal(&A2,&phi_two2.A));

        eval_walk_isom(&isom, &phi_two2, &A2, &P2, &phi_two2, &P2);

        assert(mont_equal(&A1,&A2));
        assert(mont_equal(&P1,&P2));

    }

    printf("average push\t [%f ms]\n",  (accumulated_time / repetitions));
    return 1;
}


// static GEN alg_O0_to_standard(GEN elt) {
//     return RgM_RgC_mul(global_setup.O0_to_standard, elt);
// }

// static GEN alg_standard_to_O0(GEN elt) {
//     return RgM_RgC_mul(global_setup.standard_to_O0, elt);
// }


// distorsion map on E0, twisted Edwards form
static void ted0_dist(point *Q, const point *P) {
    point Pcopy = *P;
    fp2_mul3(&Q->x, &Pcopy.t, &fp2_i);
    Q->y = Pcopy.z;
    Q->z = Pcopy.y;
    fp2_mul3(&Q->t, &Pcopy.x, &fp2_i);
}

// distorsion map on E0, montgomery form
static void mont0_dist(proj *Q, const proj *P) {
    proj Pcopy = *P;
    fp2_neg2(&Q->x, &Pcopy.x);
    Q->z = Pcopy.z;
}

// frobenius map on E0, montgomery form
static void ted0_frob(point *Q, const point *P) {
    point Pcopy = *P;
    fp2_frob2(&Q->x, &Pcopy.x);
    fp2_frob2(&Q->y, &Pcopy.y);
    fp2_frob2(&Q->z, &Pcopy.z);
    fp2_frob2(&Q->t, &Pcopy.t);
}

// frobenius map on the twist of E0, montgomery form
static void ted0_frob_twist(point *Q, const point *P) {
    fp2 ted0_twister_frob;
    fp_enc(&ted0_twister_frob.re, &(uintbig){ 7694077626149142205ULL, 2250429401854817681ULL, 12195039634258678503ULL, 4190125920647857ULL });
    fp_enc(&ted0_twister_frob.im, &(uintbig){ 14451885770804936219ULL, 10930435212066601406ULL, 15495326356676002808ULL, 5880257503665917138ULL });
  
    point Pcopy = *P;
    fp2_frob2(&Q->x, &Pcopy.x);
    fp2_frob2(&Q->y, &Pcopy.y);
    fp2_frob2(&Q->z, &Pcopy.z);
    fp2_frob2(&Q->t, &Pcopy.t);
    fp2_mul2(&Q->x, &ted0_twister_frob);
    fp2_mul2(&Q->t, &ted0_twister_frob);
}

// frobenius map on E0, montgomery form
static void mont0_frob(proj *Q, const proj *P) {
    proj Pcopy = *P;
    fp2_frob2(&Q->x, &Pcopy.x);
    fp2_frob2(&Q->z, &Pcopy.z);
}


static void mont0_id_plus_dist(proj *Q, const proj *P, const proj *A) {
    proj Pcopy = *P;
    fp2 x,x2,x3;
    x = Pcopy.x;
    fp2_sq2(&x2, &Pcopy.x);
    fp2_mul3(&x3, &x2, &Pcopy.x);

    fp2_mul3(&x, &x, &Pcopy.z);
    fp2_mul3(&x, &x, &Pcopy.z); // x*z^2
    fp2_mul3(&x2, &x2, &Pcopy.z); // x^2*z

    fp2_mul3(&Q->z, &x2, &A->z);
    fp2_add2(&Q->z, &Q->z); // multiplication by 2

    fp2_mul2(&x, &A->z);
    fp2_mul2(&x2, &A->x);
    fp2_mul2(&x3, &A->z);
    fp2_add3(&Q->x, &x, &x2);
    fp2_add2(&Q->x, &x3);
    fp2_mul2(&Q->x, &fp2_i);
    fp2_neg1(&Q->x);
}




static point ted_action(GEN coeff, const point *P, const proj *E, bool twist) {
    point A, B, C, D;

    A = B = C = D = *P;

    ted0_dist(&B,&B);

    if (!twist) ted0_frob(&C,&C);
    else ted0_frob_twist(&C,&C);

    D = B;
    if (!twist) ted0_frob(&D,&D);
    else ted0_frob_twist(&D,&D);

    uintbig x;

    gentobig(&x, gel(coeff,1));
    ted_mul(&A, &A, E, &x);
    gentobig(&x, gel(coeff,2));
    ted_mul(&B, &B, E, &x);
    gentobig(&x, gel(coeff,3));
    ted_mul(&C, &C, E, &x);
    gentobig(&x, gel(coeff,4));
    ted_mul(&D, &D, E, &x);

    point Q = A;
    ted_add(&Q, E, &Q, &B);
    ted_add(&Q, E, &Q, &C);
    ted_add(&Q, E, &Q, &D);

    return Q;
}















static bool is_in_kernel_of_endo_two(GEN endo, GEN v, proj pt) {

    proj P0 = torsion_basis_two[0], P1 = torsion_basis_two[1], P2 = torsion_basis_two[2];
    uintbig x,y;

    proj K;

    gentobig(&x, gel(v,1));
    gentobig(&y, gel(v,2));
    xBIDIM(&K, &global_setup.E0, &P0, &x, &P1, &y, &P2);

    if (!(mont_equal(&K,&pt))) return false;


    proj iP0, iP1, iP2;



    mont0_dist(&iP0, &P0);
    mont0_dist(&iP1, &P1);
    mont0_dist(&iP2, &P2);


    proj P0_iP0, P1_iP1, P2_iP2;

    mont0_id_plus_dist(&P0_iP0, &P0, &global_setup.E0);
    mont0_id_plus_dist(&P1_iP1, &P1, &global_setup.E0);
    mont0_id_plus_dist(&P2_iP2, &P2, &global_setup.E0);


    proj a_ib_0, a_ib_1, a_ib_2, a_ib_K;

    gentobig(&x, gel(endo,1));
    gentobig(&y, gel(endo,2));
    xBIDIM(&a_ib_0, &global_setup.E0, &P0, &x, &iP0, &y, &P0_iP0);
    xBIDIM(&a_ib_1, &global_setup.E0, &P1, &x, &iP1, &y, &P1_iP1);
    xBIDIM(&a_ib_2, &global_setup.E0, &P2, &x, &iP2, &y, &P2_iP2);

    gentobig(&x, gel(v,1));
    gentobig(&y, gel(v,2));
    xBIDIM(&a_ib_K, &global_setup.E0, &a_ib_0, &x, &a_ib_1, &y, &a_ib_2);


    proj c_id_0, c_id_1, c_id_2, c_id_K;

    gentobig(&x, gel(endo,3));
    gentobig(&y, gel(endo,4));
    xBIDIM(&c_id_0, &global_setup.E0, &P0, &x, &iP0, &y, &P0_iP0);
    xBIDIM(&c_id_1, &global_setup.E0, &P1, &x, &iP1, &y, &P1_iP1);
    xBIDIM(&c_id_2, &global_setup.E0, &P2, &x, &iP2, &y, &P2_iP2);

    gentobig(&x, gel(v,1));
    gentobig(&y, gel(v,2));
    xBIDIM(&c_id_K, &global_setup.E0, &c_id_0, &x, &c_id_1, &y, &c_id_2);

    proj jc_jid_K;
    mont0_frob(&jc_jid_K, &c_id_K);


    return mont_equal(&jc_jid_K, &a_ib_K);


}




static bool is_in_kernel_of_endo_ell(GEN endo, GEN v, proj *pt, long ell) {
    long e = ell_to_e(ell);
    bool twist;
    //unsigned long index = ell_to_index(ell, &twist);
    ell_to_index(ell, &twist);
    uintbig x,y;

    GEN gelle = powuu(ell,e);

    GEN cof = (twist) ? famat_prod(global_setup.gen_p_minus_fact) : famat_prod(global_setup.gen_p_plus_fact);
    cof = gdiv(cof, gelle);

    uintbig cof_big;

    gentobig(&cof_big,cof);

    const proj *basis = (twist) ? torsion_basis_twist_sum : torsion_basis_sum;

    proj P0 = basis[0], P1 = basis[1], P2 = basis[2];
    xMUL(&P0, &global_setup.E0, &P0, &cof_big);
    xMUL(&P1, &global_setup.E0, &P1, &cof_big);
    xMUL(&P2, &global_setup.E0, &P2, &cof_big);


    proj ker_alt;
    gentobig(&x, gel(v,1));
    gentobig(&y, gel(v,2));
    xBIDIM(&ker_alt, &global_setup.E0, &P0, &x, &P1, &y, &P2);

    assert(mont_equal(&ker_alt, pt));


    // GEN m_id, m_i, m_j, m_ji;
    // action_from_elle(&m_id, &m_i, &m_j, &m_ji, ell, e);


    proj iP0, iP1, iP2; 

    mont0_dist(&iP0, &P0);
    mont0_dist(&iP1, &P1);
    mont0_dist(&iP2, &P2);


    proj P0_iP0, P1_iP1, P2_iP2;


    mont0_id_plus_dist(&P0_iP0, &P0 , &global_setup.E0);
    mont0_id_plus_dist(&P1_iP1, &P1 , &global_setup.E0);
    mont0_id_plus_dist(&P2_iP2, &P2 , &global_setup.E0);


    proj a_ib_0, a_ib_1, a_ib_2, a_ib_K;


    //endo = gmod(algmul(global_setup.B, endo, mkcol4s(0,1,0,0)), gelle);


    gentobig(&x, gel(endo,1));
    gentobig(&y, gel(endo,2));
    xBIDIM(&a_ib_0, &global_setup.E0, &P0, &x, &iP0, &y, &P0_iP0);
    xBIDIM(&a_ib_1, &global_setup.E0, &P1, &x, &iP1, &y, &P1_iP1);
    xBIDIM(&a_ib_2, &global_setup.E0, &P2, &x, &iP2, &y, &P2_iP2);

    gentobig(&x, gel(v,1));
    gentobig(&y, gel(v,2));
    xBIDIM(&a_ib_K, &global_setup.E0, &a_ib_0, &x, &a_ib_1, &y, &a_ib_2);

    proj c_id_0, c_id_1, c_id_2, c_id_K;

    gentobig(&x, gel(endo,3));
    gentobig(&y, gel(endo,4));
    xBIDIM(&c_id_0, &global_setup.E0, &P0, &x, &iP0, &y, &P0_iP0);
    xBIDIM(&c_id_1, &global_setup.E0, &P1, &x, &iP1, &y, &P1_iP1);
    xBIDIM(&c_id_2, &global_setup.E0, &P2, &x, &iP2, &y, &P2_iP2);

    gentobig(&x, gel(v,1));
    gentobig(&y, gel(v,2));
    xBIDIM(&c_id_K, &global_setup.E0, &c_id_0, &x, &c_id_1, &y, &c_id_2);

    proj jc_jid_K;
    mont0_frob(&jc_jid_K, &c_id_K);

    //printf("compare: %s %s\n", proj_hash(a_ib_K), proj_hash(jc_jid_K));

    return (mont_equal(&jc_jid_K, &a_ib_K));

}















static bool is_in_kernel_of_endo_odd(GEN endo, GEN v, proj *pt, bool twist) {

    long len = (twist) ? p_minus_len : p_plus_len;
    const long *fact = (twist) ? p_minus_fact : p_plus_fact;
    const long *mult = (twist) ? p_minus_mult : p_plus_mult;


    for (int i = 0; i < len; ++i) {
        long ell = fact[i];
        long e = mult[i];

        GEN gelle = powuu(ell,e);
        GEN cof = (twist) ? famat_prod(global_setup.gen_p_minus_fact) : famat_prod(global_setup.gen_p_plus_fact);
        cof = gdiv(cof, gelle);
        uintbig cof_big;

        gentobig(&cof_big,cof);

        proj P;
        xMUL(&P, &global_setup.E0, pt, &cof_big);

        if (!(is_in_kernel_of_endo_ell(gmod(endo,gelle), gmod(v,gelle), &P, ell))) {
            printf("is_in_kernel_of_endo_odd failed at ell = %ld\n", ell);
            return false;
        }
    }

    // the above method does not detect distorsions. The one below should be finer

    proj E;
    mont_to_ted(&E, &global_setup.E0, twist);
    point K, K1, K2, K_endo;

    uintbig x,y;
    proj K_mont, K_mont_bis;


    const point *basis_ted = (twist) ? torsion_basis_twist_ted_sum : torsion_basis_ted_sum;
    const proj *basis = (twist) ? torsion_basis_twist_sum : torsion_basis_sum;
    GEN order = (twist) ? famat_prod(global_setup.gen_p_minus_fact) : famat_prod(global_setup.gen_p_plus_fact);

    K1 = basis_ted[0];
    K2 = basis_ted[1];
    gentobig(&x, gel(v,1));
    ted_mul(&K1, &K1, &E, &x);
    gentobig(&y, gel(v,2));
    ted_mul(&K2, &K2, &E, &y);

    ted_add(&K, &E, &K1, &K2);
    endo = gmod(endo,order);
    K_endo = ted_action(endo, &K, &E, twist);
    ted_to_mont_point(&K_mont, &K);
    xBIDIM(&K_mont_bis, &global_setup.E0, &basis[0], &x, &basis[1], &y, &basis[2]);

    assert(mont_equal(&K_mont, &K_mont_bis));

    return ted_iszero(&K_endo);


}




int test_odd_isog() {
    float accumulated_time = 0.;
    int repetitions = 20;
    clock_t t;
    GEN A = global_setup.B;
    GEN order = global_setup.O0;
    GEN fm = famat_mul(global_setup.gen_p_plus_fact, global_setup.gen_p_minus_fact);

    for (int i = 0; i < repetitions; i++) {
        t = tic();

        GEN Nodd = gsqr(global_setup.gen_odd_torsion);
        GEN Ntwo = powiu(gen_2, 20);

        GEN gamma = NULL;
        while (!gamma) {
            // parity option is 1 so the 2-walk is not backtracking
            gamma = norm_equation_special(global_setup.p, gmul(Nodd,Ntwo), 1, true);
        }
        gamma = gtrans(gamma);
        GEN n;

        gamma = alg_primitive(&n, A, order, gamma);

        GEN N = algnorm(A,gamma,0);
        Nodd = ggcd(Nodd, N);

        assert(gcmp(N, gmul(Nodd,Ntwo)) == 0);

        GEN N_H1 = ggcd(global_setup.gen_odd_torsion, Nodd);

        GEN H1_odd = lideal_create(A, order, gamma, N_H1);


        assert(gcmp(N_H1, lideal_norm(H1_odd)) == 0);

        GEN coeff = ideal_to_kernel_O0_T(H1_odd, famat_Z_gcd(fm,lideal_norm(H1_odd)));

        proj E, E_twist;
        mont_to_ted(&E, &global_setup.E0, false);
        mont_to_ted(&E_twist, &global_setup.E0, true);




        point K, K1, K2, K_endo;


        // CHECK KERNEL PLUS
        GEN coeff_plus = gel(coeff,1);

        GEN endo;
        uintbig x,y;
        GEN v_plus, v_minus;
        proj K_mont, K_mont_bis;

        K1 = torsion_basis_ted_sum[0];
        K2 = torsion_basis_ted_sum[1];

        v_plus = torsion_crt_compose(coeff_plus, false);



        gentobig(&x, gel(v_plus,1));
        ted_mul(&K1, &K1, &E, &x);
        gentobig(&y, gel(v_plus,2));
        ted_mul(&K2, &K2, &E, &y);

        ted_add(&K, &E, &K1, &K2);
        endo = gmod(gmul(gamma, gen_2), gadd(global_setup.p, gen_1));
        K_endo = ted_action(endo, &K, &E, false);
        assert(ted_iszero(&K_endo));

        ted_to_mont_point(&K_mont, &K);
        xBIDIM(&K_mont_bis, &global_setup.E0, &torsion_basis_sum[0], &x, &torsion_basis_sum[1], &y, &torsion_basis_sum[2]);

        assert(mont_equal(&K_mont, &K_mont_bis));


        // CHECK KERNEL MINUS
        GEN coeff_minus = gel(coeff,2);


        K1 = torsion_basis_twist_ted_sum[0];
        K2 = torsion_basis_twist_ted_sum[1];

        v_minus = torsion_crt_compose(coeff_minus, true);

        gentobig(&x, gel(v_minus,1));
        ted_mul(&K1, &K1, &E_twist, &x);
        gentobig(&y, gel(v_minus,2));
        ted_mul(&K2, &K2, &E_twist, &y);

        ted_add(&K, &E_twist, &K1, &K2);
        endo = gmod(gmul(gamma, gen_2), gsub(global_setup.p, gen_1));
        K_endo = ted_action(endo, &K, &E_twist, true);

        ted_to_mont_point(&K_mont, &K);
        xBIDIM(&K_mont_bis, &global_setup.E0, &torsion_basis_twist_sum[0], &x, &torsion_basis_twist_sum[1], &y, &torsion_basis_twist_sum[2]);

        assert(mont_equal(&K_mont, &K_mont_bis));
        assert(ted_iszero(&K_endo));



        odd_isogeny psi_1 = ideal_to_isogeny_O0_T(H1_odd, famat_Z_gcd(fm,lideal_norm(H1_odd)));
        endo = gmod(gamma, gadd(global_setup.p, gen_1));
        assert(is_in_kernel_of_endo_odd(endo, v_plus, &psi_1.kernel_plus, false));
        endo = gmod(gamma, gsub(global_setup.p, gen_1));
        assert(is_in_kernel_of_endo_odd(endo, v_minus, &psi_1.kernel_minus, true));




        printf("passed\n");

        accumulated_time += toc(t);
    }

    printf("average time\t [%f ms]\n",  (accumulated_time / repetitions));
    return 1;
}









int test_two_isog() {
    float accumulated_time = 0.;
    int repetitions = 20;
    clock_t t;
    GEN A = global_setup.B;
    GEN order = global_setup.O0;

    for (int i = 0; i < repetitions; i++) {
        t = tic();

        GEN Nodd = gsqr(global_setup.gen_odd_torsion);
        long e1 = 33;
        GEN Ntwo = powiu(gen_2, e1 + 20);

        GEN gamma = NULL;
        while (!gamma) {
            // parity option is 1 so the 2-walk is not backtracking
            gamma = norm_equation_special(global_setup.p, gmul(Nodd,Ntwo), 1, true);
        }
        gamma = gtrans(gamma);
        GEN n;
        gamma = alg_primitive(&n, A, order, gamma);



        GEN H1_two = lideal_create(A, order, gamma, powuu(2, e1));
        



        // GEN v = ideal_to_kernel_O0_ell(H1_two, 2);


        // proj P0 = torsion_basis_two[0], P1 = torsion_basis_two[1], P2 = torsion_basis_two[2];
        // proj iP0, iP1, iP2;//, iP0_alt, iP1_alt, iP2_alt;


        // //proj jP0, jP1, jP2, jP0_alt, jP1_alt, jP2_alt;
        // //proj jiP0, jiP1, jiP2, jiP0_alt, jiP1_alt, jiP2_alt;

        // //GEN gelle = powuu(2,33);

        // GEN m_id = mkmat2(mkcol2s(1,0),mkcol2s(0,1));
        // GEN m_i = global_setup.action_two_2;
        // // GEN m_1ji2 = global_setup.action_two_3;
        // // GEN m_ij2 = global_setup.action_two_4;

        // //GEN m_j = gmod(gsub(gmul(m_ij2,gen_2), m_i),gelle);
        // //GEN m_ji = gmod(gsub(m_id, gmul(m_1ji2,gen_2)),gelle);

        // //assert(isexactzero(gmod(gsub(m_ji, gmul(m_j,m_i)), gelle)));


        // mont0_dist(&iP0, &P0);
        // mont0_dist(&iP1, &P1);
        // mont0_dist(&iP2, &P2);

        // // mont0_frob(&jP0, &P0);
        // // mont0_frob(&jP1, &P1);
        // // mont0_frob(&jP2, &P2);

        // // mont0_frob(&jiP0, &iP0);
        // // mont0_frob(&jiP1, &iP1);
        // // mont0_frob(&jiP2, &iP2);

        // // gentobig(&x, gcoeff(m_i,1,1));
        // // gentobig(&y, gcoeff(m_i,2,1));
        // // xBIDIM(&iP0_alt, &global_setup.E0, &P0, &x, &P1, &y, &P2);
        // // assert(mont_equal(&iP0, &iP0_alt));

        // // gentobig(&x, gcoeff(m_i,1,2));
        // // gentobig(&y, gcoeff(m_i,2,2));
        // // xBIDIM(&iP1_alt, &global_setup.E0, &P0, &x, &P1, &y, &P2);
        // // assert(mont_equal(&iP1, &iP1_alt));
        
        // // gentobig(&x, gadd(gcoeff(m_i,1,1),gcoeff(m_i,1,2)));
        // // gentobig(&y, gadd(gcoeff(m_i,2,1),gcoeff(m_i,2,2)));
        // // xBIDIM(&iP2_alt, &global_setup.E0, &P0, &x, &P1, &y, &P2);
        // // assert(mont_equal(&iP2, &iP2_alt));


        // // gentobig(&x, gcoeff(m_j,1,1));
        // // gentobig(&y, gcoeff(m_j,2,1));
        // // xBIDIM(&jP0_alt, &global_setup.E0, &P0, &x, &P1, &y, &P2);
        // // assert(mont_equal(&jP0, &jP0_alt));

        // // gentobig(&x, gcoeff(m_j,1,2));
        // // gentobig(&y, gcoeff(m_j,2,2));
        // // xBIDIM(&jP1_alt, &global_setup.E0, &P0, &x, &P1, &y, &P2);
        // // assert(mont_equal(&jP1, &jP1_alt));

        // // gentobig(&x, gadd(gcoeff(m_j,1,1),gcoeff(m_j,1,2)));
        // // gentobig(&y, gadd(gcoeff(m_j,2,1),gcoeff(m_j,2,2)));
        // // xBIDIM(&jP2_alt, &global_setup.E0, &P0, &x, &P1, &y, &P2);
        // // assert(mont_equal(&jP2, &jP2_alt));

        // // gentobig(&x, gcoeff(m_ji,1,1));
        // // gentobig(&y, gcoeff(m_ji,2,1));
        // // xBIDIM(&jiP0_alt, &global_setup.E0, &P0, &x, &P1, &y, &P2);
        // // assert(mont_equal(&jiP0, &jiP0_alt));

        // // gentobig(&x, gcoeff(m_ji,1,2));
        // // gentobig(&y, gcoeff(m_ji,2,2));
        // // xBIDIM(&jiP1_alt, &global_setup.E0, &P0, &x, &P1, &y, &P2);
        // // assert(mont_equal(&jiP1, &jiP1_alt));

        // // gentobig(&x, gadd(gcoeff(m_ji,1,1),gcoeff(m_ji,1,2)));
        // // gentobig(&y, gadd(gcoeff(m_ji,2,1),gcoeff(m_ji,2,2)));
        // // xBIDIM(&jiP2_alt, &global_setup.E0, &P0, &x, &P1, &y, &P2);
        // // assert(mont_equal(&jiP2, &jiP2_alt));




        // GEN endo = gmod(gmul(gamma, gen_2), powuu(2,33));


        // // GEN m_endo, m_endo_ab, m_endo_cd;
        // // m_endo_ab = gmul(gel(endo,1), m_id);
        // // m_endo_ab = gadd(m_endo_ab, gmul(gel(endo,2), m_i));
        // // m_endo_cd = gmul(gel(endo,3), m_j);
        // // m_endo_cd = gadd(m_endo_cd, gmul(gel(endo,4), m_ji));

        // // m_endo = gadd(m_endo_ab, m_endo_cd);

        // // proj R;
        // // proj a_ib_K_alt, jc_jid_K_alt;


        // // GEN w, wab,wcd;

        // // w = gmul(m_endo, v);
        // // gentobig(&x, gel(w,1));
        // // gentobig(&y, gel(w,2));
        // // xBIDIM(&R, &global_setup.E0, &P0, &x, &P1, &y, &P2);

        // // assert(mont_iszero(&R));

        // // wab = gmul(m_endo_ab, v);
        // // gentobig(&x, gel(wab,1));
        // // gentobig(&y, gel(wab,2));
        // // xBIDIM(&a_ib_K_alt, &global_setup.E0, &P0, &x, &P1, &y, &P2);

        // // wcd = gmul(m_endo_cd, v);
        // // gentobig(&x, gel(wcd,1));
        // // gentobig(&y, gel(wcd,2));
        // // xBIDIM(&jc_jid_K_alt, &global_setup.E0, &P0, &x, &P1, &y, &P2);

        // // assert(mont_equal(&jc_jid_K_alt, &a_ib_K_alt));




        // proj P0_iP0_alt, P1_iP1_alt, P2_iP2_alt;

        // proj P0_iP0, P1_iP1, P2_iP2;


        // mont0_id_plus_dist(&P0_iP0, &P0 , &global_setup.E0);
        // mont0_id_plus_dist(&P1_iP1, &P1 , &global_setup.E0);
        // mont0_id_plus_dist(&P2_iP2, &P2 , &global_setup.E0);


        // // GEN m_id_plus_i = gadd(m_id, m_i);

        // // gentobig(&x, gcoeff(m_id_plus_i,1,1));
        // // gentobig(&y, gcoeff(m_id_plus_i,2,1));
        // // xBIDIM(&P0_iP0_alt, &global_setup.E0, &P0, &x, &P1, &y, &P2);


        // // gentobig(&x, gcoeff(m_id_plus_i,1,2));
        // // gentobig(&y, gcoeff(m_id_plus_i,2,2));
        // // xBIDIM(&P1_iP1_alt, &global_setup.E0, &P0, &x, &P1, &y, &P2);

        // // gentobig(&x, gadd(gcoeff(m_id_plus_i,1,1),gcoeff(m_id_plus_i,1,2)));
        // // gentobig(&y, gadd(gcoeff(m_id_plus_i,2,1),gcoeff(m_id_plus_i,2,2)));
        // // xBIDIM(&P2_iP2_alt, &global_setup.E0, &P0, &x, &P1, &y, &P2);

        // // assert(mont_equal(&P0_iP0, &P0_iP0_alt));
        // // assert(mont_equal(&P1_iP1, &P1_iP1_alt));
        // // assert(mont_equal(&P2_iP2, &P2_iP2_alt));



        // proj a_ib_0, a_ib_1, a_ib_2, a_ib_K;
        // // proj a_ib_0_alt, a_ib_1_alt, a_ib_2_alt;

        // gentobig(&x, gel(endo,1));
        // gentobig(&y, gel(endo,2));
        // xBIDIM(&a_ib_0, &global_setup.E0, &P0, &x, &iP0, &y, &P0_iP0);
        // xBIDIM(&a_ib_1, &global_setup.E0, &P1, &x, &iP1, &y, &P1_iP1);
        // xBIDIM(&a_ib_2, &global_setup.E0, &P2, &x, &iP2, &y, &P2_iP2);


        // // gentobig(&x, gcoeff(m_endo_ab,1,1));
        // // gentobig(&y, gcoeff(m_endo_ab,2,1));
        // // xBIDIM(&a_ib_0_alt, &global_setup.E0, &P0, &x, &P1, &y, &P2);
        // // assert(mont_equal(&a_ib_0, &a_ib_0_alt));
        // // gentobig(&x, gcoeff(m_endo_ab,1,2));
        // // gentobig(&y, gcoeff(m_endo_ab,2,2));
        // // xBIDIM(&a_ib_1_alt, &global_setup.E0, &P0, &x, &P1, &y, &P2);
        // // assert(mont_equal(&a_ib_1, &a_ib_1_alt));
        // // gentobig(&x, gadd(gcoeff(m_endo_ab,1,1),gcoeff(m_endo_ab,1,2)));
        // // gentobig(&y, gadd(gcoeff(m_endo_ab,2,1),gcoeff(m_endo_ab,2,2)));
        // // xBIDIM(&a_ib_2_alt, &global_setup.E0, &P0, &x, &P1, &y, &P2);
        // // assert(mont_equal(&a_ib_2, &a_ib_2_alt));



        // gentobig(&x, gel(v,1));
        // gentobig(&y, gel(v,2));
        // xBIDIM(&a_ib_K, &global_setup.E0, &a_ib_0, &x, &a_ib_1, &y, &a_ib_2);


        // //assert(mont_equal(&a_ib_K, &a_ib_K_alt));

        // proj c_id_0, c_id_1, c_id_2, c_id_K;

        // gentobig(&x, gel(endo,3));
        // gentobig(&y, gel(endo,4));
        // xBIDIM(&c_id_0, &global_setup.E0, &P0, &x, &iP0, &y, &P0_iP0);
        // xBIDIM(&c_id_1, &global_setup.E0, &P1, &x, &iP1, &y, &P1_iP1);
        // xBIDIM(&c_id_2, &global_setup.E0, &P2, &x, &iP2, &y, &P2_iP2);

        // gentobig(&x, gel(v,1));
        // gentobig(&y, gel(v,2));
        // xBIDIM(&c_id_K, &global_setup.E0, &c_id_0, &x, &c_id_1, &y, &c_id_2);

        // proj jc_jid_K;
        // mont0_frob(&jc_jid_K, &c_id_K);


        // assert(mont_equal(&jc_jid_K, &a_ib_K));

        // two_walk phi = ideal_to_isogeny_O0_two(H1_two);
        // proj ker_alt;
        // gentobig(&x, gel(v,1));
        // gentobig(&y, gel(v,2));
        // xBIDIM(&ker_alt, &global_setup.E0, &P0, &x, &P1, &y, &P2);

        // assert(mont_equal(&ker_alt, &phi.ker));

        two_walk phi = ideal_to_isogeny_O0_two(H1_two);
        GEN v = ideal_to_kernel_O0_ell(H1_two, 2);
        GEN endo = gmod(gmul(gamma, gen_1), powuu(2,33));

        assert(is_in_kernel_of_endo_two(endo, v, phi.ker));

        printf("passed\n");

        accumulated_time += toc(t);
    }

    printf("average time\t [%f ms]\n",  (accumulated_time / repetitions));
    return 1;
}









int test_diamond() {
    float accumulated_time = 0.;
    int repetitions = 5;
    clock_t t;
    GEN A = global_setup.B;
    GEN order = global_setup.O0;
    GEN fm = famat_mul(global_setup.gen_p_plus_fact, global_setup.gen_p_minus_fact);

    for (int i = 0; i < repetitions; i++) {
        t = tic();

        long e1 = 33;
        long e2 = 33;
        GEN Nodd = gsqr(global_setup.gen_odd_torsion);
        GEN Ntwo = powiu(gen_2, e1+e2);


        // STEP 1: Find an endomorphism gamma of norm Nodd*Ntwo = T^2 * 2^e1 * 2^e2,
        //         where T = global_setup.gen_odd_torsion

        GEN gamma = NULL;
        while (!gamma) {
            // parity option is 1 so the 2-walk is not backtracking
            gamma = norm_equation_special(global_setup.p, gmul(Nodd,Ntwo), 1, true);
        }
        gamma = gtrans(gamma);
        GEN n;
        gamma = alg_primitive(&n, A, order, gamma);

        GEN N = algnorm(A,gamma,0);
        Nodd = ggcd(Nodd, N);

        assert(gcmp(N, gmul(Nodd,Ntwo)) == 0);
        assert(Z_lval(N, 2) == e1+e2);

        // long m = remove_1_i(A,order,&gamma);
        // assert(m < e1);
        // e1 -= m;

        GEN gamma_conj = alg_conj(A, gamma);
        
        // m = remove_1_i(A,order,&gamma_conj);
        // assert(m < e2);
        // e2 -= m;
        // gamma = alg_conj(A, gamma_conj);
        // printf("e1 = %ld; e2 = %ld\n", e1, e2);

        // GEN gamma_i = algmul(A, gamma, mkcol4s(0,1,0,0));
        // GEN gamma_conj_i = algmul(A, gamma_conj, mkcol4s(0,1,0,0));


        // STEP 2: Compute the kernel of the degree T and 2^e1 initial isogenies of 
        //         gamma, and the degree T and 2^e2 initial isogenies of gamma_conj
        //         The corresponding isogenies are psi_1 (T), phi_1 (2^e1), 
        //         psi_2 (T), and phi_2 (2^e2), 

        GEN H1_odd = lideal_create(A, order, gamma, ggcd(global_setup.gen_odd_torsion, Nodd));
        GEN H1_two = lideal_create(A, order, gamma, powuu(2, e1));

        GEN H2_odd = lideal_create(A, order, gamma_conj, gdiv(Nodd, lideal_norm(H1_odd)));
        GEN H2_two = lideal_create(A, order, gamma_conj, powuu(2, e2));

        odd_isogeny psi_1 = ideal_to_isogeny_O0_T(H1_odd, famat_Z_gcd(fm,lideal_norm(H1_odd)));
        odd_isogeny psi_2 = ideal_to_isogeny_O0_T(H2_odd, famat_Z_gcd(fm,lideal_norm(H2_odd)));

        proj psi_1_source = global_setup.E0;
        proj psi_2_source = global_setup.E0;
        two_walk phi_1 = ideal_to_isogeny_O0_two(H1_two);
        two_walk phi_2 = ideal_to_isogeny_O0_two(H2_two);



        // STEP 3: Check the isogeny correctness of the kernels from previous step

        GEN coeff_1 = ideal_to_kernel_O0_T(H1_odd, famat_Z_gcd(fm,lideal_norm(H1_odd)));
        GEN v_plus_1 = torsion_crt_compose(gel(coeff_1,1), false);
        GEN v_minus_1 = torsion_crt_compose(gel(coeff_1,2), true);

        assert(is_in_kernel_of_endo_odd(gamma, v_plus_1, &psi_1.kernel_plus, false));
        assert(is_in_kernel_of_endo_odd(gamma, v_minus_1, &psi_1.kernel_minus, true));

        GEN coeff_2 = ideal_to_kernel_O0_T(H2_odd, famat_Z_gcd(fm,lideal_norm(H2_odd)));
        GEN v_plus_2 = torsion_crt_compose(gel(coeff_2,1), false);
        GEN v_minus_2 = torsion_crt_compose(gel(coeff_2,2), true);

        assert(is_in_kernel_of_endo_odd(gamma_conj, v_plus_2, &psi_2.kernel_plus, false));
        assert(is_in_kernel_of_endo_odd(gamma_conj, v_minus_2, &psi_2.kernel_minus, true));

        GEN v1 = ideal_to_kernel_O0_ell(H1_two, 2);
        GEN endo1 = gmod(gamma, powuu(2,33));
        assert(is_in_kernel_of_endo_two(endo1, v1, phi_1.ker));

        GEN v2 = ideal_to_kernel_O0_ell(H2_two, 2);
        GEN endo2 = gmod(gamma_conj, powuu(2,33));
        assert(is_in_kernel_of_endo_two(endo2, v2, phi_2.ker));


        // STEP 3: Compute psi_1*phi_1, phi_1*psi_1, psi_2*phi_2, phi_2*psi_2

        two_walk phi_1_pushed = push_two_walk_through_odd_isogeny(&phi_1, &psi_1, &psi_1_source);
        two_walk phi_2_pushed = push_two_walk_through_odd_isogeny(&phi_2, &psi_2, &psi_2_source);
        
        proj psi_1_pushed_source = psi_1_source;
        proj psi_2_pushed_source = psi_2_source;
        odd_isogeny psi_1_pushed = push_odd_isogeny_through_two_walk(&psi_1, &psi_1_pushed_source, &phi_1);
        odd_isogeny psi_2_pushed = push_odd_isogeny_through_two_walk(&psi_2, &psi_2_pushed_source, &phi_2);

        proj A1, A2, pt1 = phi_1_pushed.ker, pt2 = phi_2_pushed.ker;
        isomorphism isom1, isom2;
        two_walk phi_1_pushed_adjusted;
        two_walk phi_2_pushed_adjusted;
        eval_walk_isom(&isom1, &phi_1_pushed_adjusted, &A1, &pt1, &phi_1_pushed, &pt1);
        eval_walk_isom(&isom2, &phi_2_pushed_adjusted, &A2, &pt2, &phi_2_pushed, &pt2);

        assert(fp2_iszero(&pt1.z) && !fp2_iszero(&pt1.x));
        assert(fp2_iszero(&pt2.z) && !fp2_iszero(&pt2.x));

        pt1 = psi_1_pushed.kernel_plus;
        pt2 = psi_2_pushed.kernel_plus;

        proj B1 = psi_1_pushed_source;
        eval(&B1, & psi_1_pushed, &pt1);
        proj B2 = psi_2_pushed_source;
        eval(&B2, & psi_2_pushed, &pt2);

        assert(fp2_iszero(&pt1.z) && !fp2_iszero(&pt1.x));
        assert(fp2_iszero(&pt2.z) && !fp2_iszero(&pt2.x));

        // STEP 3: Check that psi_1*phi_1, phi_1*psi_1, psi_2*phi_2, phi_2*psi_2
        //         all have the same target curve

        proj jA1,jA2;
        jinv256(&jA1, &A1);
        jinv256(&jA2, &A2);

        proj jB1,jB2;
        jinv256(&jB1, &B1);
        jinv256(&jB2, &B2);

        assert(mont_equal(&jA1,&jB1));
        assert(mont_equal(&jA2,&jB2));
        assert(mont_equal(&jA1,&jA2));

        accumulated_time += toc(t);
    }

    printf("average time\t [%f ms]\n",  (accumulated_time / repetitions));
    return 1;
}


// argv[1] is the random seed; default = 1
int main(int argc, char *argv[]){
    pari_init(80000000, 1<<18);
    init_precomputations();

    setrand(stoi(1));
    srand48(1);
    if( argc > 1 ) {
      setrand(strtoi(argv[1]));
      srand48(atoi(argv[1]));
    }

    //test_kertoid_action();
    //test_kertoid_O0_ell();
    //test_kertoid_O0_two();
    //test_kertoid_O0_T();
    //test_push_isogenies();
    //test_odd_isog();
    //test_two_isog();
    test_diamond();
    //test_push_isogenies_long();

    printf("    \033[1;32mAll tests passed\033[0m\n");
    exit(0);
}



