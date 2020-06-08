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
#include "tedwards.h"
#include "constants.h"

struct quaternion_setup_t {
    GEN p; // the prime
    GEN B; // the quaternion algebra
    GEN qf; // the quaternion algebra
    GEN O0; // the cannonical maximal order
    GEN one;
    GEN i;
    GEN j;
    GEN ji;
    GEN torsion_fm; // factorisation matrix of the available torsion

    GEN O0_b1;
    GEN O0_b2;
    GEN O0_b3;
    GEN O0_b4;
    GEN O0_to_standard;
    GEN standard_to_O0;
};

struct quaternion_setup_t global_setup;


GEN norm0(GEN x) {
    return algnorm(global_setup.B, x,0);
}

/*
static GEN alg_O0_to_standard(GEN elt) {
    return RgM_RgC_mul(global_setup.O0_to_standard, elt);
}
*/

static GEN alg_standard_to_O0(GEN elt) {
    return RgM_RgC_mul(global_setup.standard_to_O0, elt);
}



proj random_point(proj const *A, long ell, long e) {
    proj P;
    uintbig cofactor;
    uintbig_add3(&cofactor, &p, &uintbig_1);
    uintbig ell_big;
    uintbig_set(&ell_big, ell);
    for (int i = 0; i < e; ++i) {
        uintbig_div3_64(&cofactor, &cofactor, ell); 
    }
    proj Z;

    while (1) {
        fp2_random(&P.x); fp2_random(&P.z);
        if (!is_on_curve(&P, A)) continue;
        xMUL(&P, A, &P, &cofactor);
        Z = P;
        for (int i = 0; i < e-1; ++i) {
            xMUL(&Z, A, &Z, &ell_big);
        }
        if (!fp2_iszero(&Z.z)) { 
            //xMUL(&Z, A, &Z, &ell_big);
            //assert(fp2_iszero(&Z.z));
            return P;
        }
    }

}

int test_random() {
    // float accumulated_time_ms = 0.;
    int repetitions = 1000;
    clock_t t;

    proj A = { fp2_0, fp2_1 };
    long len = p_plus_len;
    proj P;

    t = tic();
    for (int i = 0; i < repetitions; i++) {
        long ell = p_plus_fact[i % len], e = p_plus_mult[i % len];
        P = random_point(&A, ell, e);
    }
    TOC(t, "rand");

    // printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));

    return 0;
}


int test_ted() {
    float accumulated_time_ms = 0.;
    clock_t t;

    uintbig k;

    proj A = { fp2_0, fp2_1 };
    long len = p_plus_len;
    int repetitions = p_plus_len*10;
    proj P,Q;

    proj E;
    mont_to_ted(&E, &A, false);

    point tP,tQ,tPpQ,tPmQ,t2P, tsum;

    for (int i = 0; i < repetitions; i++) {
        t = tic();
        long ell = p_plus_fact[i % len], e = p_plus_mult[i % len];

        P = random_point(&A, ell, e);
        Q = random_point(&A, ell, e);

        mont_to_ted_point(&tP, &A, &P);
        mont_to_ted_point(&tQ, &A, &Q);

        ted_add(&tPpQ, &E, &tP, &tQ);
        ted_neg(&tQ, &tQ);
        ted_add(&tPmQ, &E, &tP, &tQ);
        ted_neg(&tQ, &tQ);

        ted_double(&t2P, &E, &tP);

        assert(ted_is_on_curve(&tPpQ,&E));
        assert(ted_is_on_curve(&tPmQ,&E));
        assert(ted_is_on_curve(&t2P,&E));

        ted_neg(&tsum, &t2P);
        ted_add(&tsum, &E, &tsum, &tPpQ);
        ted_add(&tsum, &E, &tsum, &tPmQ);

        assert(ted_iszero(&tsum));


        ted_double(&tsum, &E, &tP);
        ted_add(&tsum, &E, &tsum, &tP);
        ted_add(&tsum, &E, &tsum, &tP);
        ted_add(&tsum, &E, &tsum, &tP);
        ted_add(&tsum, &E, &tsum, &tP);
        ted_add(&tsum, &E, &tsum, &tP);
        ted_add(&tsum, &E, &tsum, &tP);
        ted_add(&tsum, &E, &tsum, &tP);
        ted_add(&tsum, &E, &tsum, &tP);
        ted_add(&tsum, &E, &tsum, &tP);
        ted_add(&tsum, &E, &tsum, &tP);
        ted_add(&tsum, &E, &tsum, &tP); // 13*P


        uintbig_set(&k,13);
        ted_mul(&tQ, &tP, &E, &k);

        ted_neg(&tsum, &tsum);
        ted_add(&tsum, &E, &tsum, &tQ);

        assert(ted_iszero(&tsum));


        accumulated_time_ms += toc(t);

    }

    //printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));

    return 0;
}


int test_weil() {
    float accumulated_time_ms = 0.;
    clock_t t;

    proj A = { fp2_0, fp2_1 };
    long len = p_plus_len;
    int repetitions = p_plus_len*10;
    proj montP,montQ, montR, montS;
    uintbig ellbig;

    proj E;
    mont_to_ted(&E, &A, false);

    point P,Q,R,S,T;
    fp2 a, b, tmp;


    for (int i = 0; i < repetitions; i++) {

        long ell = p_plus_fact[i % len], e = p_plus_mult[i % len];
        //ell = 5;
        e = 1;

        uintbig_set(&ellbig, ell);


        montP = random_point(&A, ell, e);
        montQ = random_point(&A, ell, e);
        montR = random_point(&A, ell, e);
        montS = random_point(&A, ell, e);

        mont_to_ted_point(&P, &A, &montP);
        mont_to_ted_point(&Q, &A, &montQ);
        mont_to_ted_point(&R, &A, &montR);
        mont_to_ted_point(&S, &A, &montS);

        assert(ted_is_on_curve(&P,&E));
        assert(ted_is_on_curve(&Q,&E));
        assert(ted_is_on_curve(&R,&E));
        assert(ted_is_on_curve(&S,&E));

        ted_mul(&T, &P, &E, &ellbig);
        assert(ted_iszero(&T));
        ted_mul(&T, &Q, &E, &ellbig);
        assert(ted_iszero(&T));
        ted_mul(&T, &R, &E, &ellbig);
        assert(ted_iszero(&T));
        ted_mul(&T, &S, &E, &ellbig);
        assert(ted_iszero(&T));

        t = tic();
        ted_weil(&a, &E, &P, &Q, &ellbig);
        fp2_exp(&tmp, &a, &ellbig);
        fp2_sub2(&tmp, &fp2_1);

        assert(fp2_iszero(&tmp));

        ted_add(&R, &E, &Q, &P);
        ted_weil(&b, &E, &P, &R, &ellbig);

        fp2_sub3(&tmp,&a,&b);

        assert(fp2_iszero(&tmp));
        accumulated_time_ms += toc(t);

    }

    //printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));

    return 0;
}


int test_bidim() {
    // float accumulated_time_ms = 0.;
    clock_t t;

    proj A = { fp2_0, fp2_1 };
    long len = p_plus_len;
    int repetitions = p_plus_len*10;
    proj P,Q,PQ,R1,R2;
    uintbig k_big, a_big, b_big, k_1, a_kb;
    long k, a, b;
    fp2 x1z2,x2z1,diff;

    t = tic();
    for (int i = 0; i < repetitions; i++) {
        long ell = p_plus_fact[i % len], e = p_plus_mult[i % len];
        k = random_Fl(ell);
        a = random_Fl(ell);
        b = random_Fl(ell);
        uintbig_set(&k_big, k);
        uintbig_set(&a_big, a);
        uintbig_set(&b_big, b);

        uintbig_set(&k_1, k+1);

        P = random_point(&A, ell, e);
        //printf("%ld^%ld\n",ell,e);
        xMUL(&Q, &A, &P, &k_big);
        xMUL(&PQ, &A, &P, &k_1);

        xBIDIM(&R1, &A, &P, &a_big, &Q, &b_big, &PQ);


        uintbig_set(&a_kb, a+k*b);
        xMUL(&R2, &A, &P, &a_kb);

        fp2_mul3(&x1z2, &R1.x, &R2.z);
        fp2_mul3(&x2z1, &R2.x, &R1.z);
        fp2_sub3(&diff, &x1z2, &x2z1);
        
        assert(fp2_iszero(&diff));
    }
    TOC(t, "bidim");

    // printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));

    return 0;
}



GEN kerner_to_ideal(GEN v, GEN m1, GEN m2, GEN m3, GEN m4, long ell, long e) {
    pari_sp ltop = avma;

    GEN gelle = gpowgs(stoi(ell),e);

    GEN v1 = gmul(m1, v);
    GEN v2 = gmul(m2, v);
    GEN v3 = gmul(m3, v);
    GEN v4 = gmul(m4, v);

    GEN matsys = mkmat4(v1,v2,v3,v4);
    GEN ker = matkermod(matsys, gelle, NULL); // flag = 1 because integral entries


    GEN sol, sol_reduced;

    int dim_ker = lg(ker)-1;
    for (int i = 1; i <= dim_ker; ++i){
        sol = gel(ker,i);
        sol_reduced = gmodgs(sol,ell);
        if (!isexactzero(sol_reduced)) break;
    }

    GEN generator = gmul(gel(sol,1),global_setup.O0_b1);
    generator = gadd(generator, gmul(gel(sol,2),global_setup.O0_b2));
    generator = gadd(generator, gmul(gel(sol,3),global_setup.O0_b3));
    generator = gadd(generator, gmul(gel(sol,4),global_setup.O0_b4));

    GEN ideal = lideal_create(global_setup.B, global_setup.O0, generator, gelle);

    return gerepilecopy(ltop, ideal);
}

GEN ideal_to_kernel(GEN I, GEN m1, GEN m2, GEN m3, GEN m4, long ell, long e) {
    pari_sp ltop = avma;

    GEN gelle = gpowgs(stoi(ell),e);

    GEN generator = lideal_generator(I);
    GEN generator_O0 = alg_standard_to_O0(generator);
    GEN endo = gmul(gel(generator_O0,1),m1);
    endo = gadd(endo, gmul(gel(generator_O0,2),m2));
    endo = gadd(endo, gmul(gel(generator_O0,3),m3));
    endo = gadd(endo, gmul(gel(generator_O0,4),m4));

    GEN ker = matkermod(endo, gelle, NULL);

    GEN sol, sol_reduced;

    int dim_ker = lg(ker)-1;
    for (int i = 1; i <= dim_ker; ++i){
        sol = gel(ker,i);
        sol_reduced = gmodgs(sol,ell);
        if (!isexactzero(sol_reduced)) break;
    }

    // remains to compute sol[1]*P1 + sol[2]*P2 where P1,P2 is a basis of the torsion

    return gerepilecopy(ltop, sol);
}

int test_kertoid() {

    // float accumulated_time_ms = 0.;
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


    ell = 2;
    e = 33;
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

        ideal = kerner_to_ideal(v, m_1, m_2, m_3, m_4, ell, e);

        // n1 = lideal_norm(ideal);
        // gel(ideal,2) = gen_0;
        // n2 = lideal_norm(ideal);
        // printf("%d",gcmp(n1,n2));

        w = ideal_to_kernel(ideal, m_1, m_2, m_3, m_4, ell, e);

        assert(gcmp(gmod(QM_det(mkmat2(v,w)),gelle),gen_0) == 0);

    }
    TOC(t, "kertoid");

    // printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));

    return 0;
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
        p = strtoi("73743043621499797449074820543863456997944695372324032511999999999999999999999"),
        b = negi(p);

    GEN B = alg_hilbert(nf, a, b, var, 0);

    GEN B_1 = mkcol4s(1,0,0,0);
    GEN B_i = mkcol4s(0,1,0,0);
    GEN B_j = mkcol4s(0,0,1,0);
    GEN B_ji = mkcol4s(0,0,0,1);
    //GEN B_ij = mkcol4s(0,0,0,-1);

    GEN B_1k_2 = mkcol4(ghalf,gen_0,gen_0,gneg(ghalf)); // (1-ji)/2
    GEN B_ij_2 = mkcol4(gen_0,ghalf,ghalf,gen_0); // (i+j)/2

    GEN B_O0 = alglathnf(B,mkmat4(B_1, B_i, B_1k_2, B_ij_2), gen_0);

    global_setup.p = p;
    global_setup.B = B; // the quaternion algebra
    global_setup.qf = mkmat4(mkcol4s(1,0,0,0),
                             mkcol4s(0,1,0,0),
                             mkcol4(gen_0,gen_0,p,gen_0),
                             mkcol4(gen_0,gen_0,gen_0,p)); // quadratic form defined by the reduced norm

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
    
    test_weil();

    printf("    \033[1;32mAll tests passed\033[0m\n");
    exit(0);
}



