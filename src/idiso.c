
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pari/pari.h>
#include <assert.h>


#include "ideal.h"
#include "idiso.h"
#include "constants.h"
#include "precomputed.h"
#include "isogenies.h"
#include "klpt.h"
#include "toolbox.h"




#ifndef NDEBUG
// static isog_degree mont_order(const proj *P, const proj *A, const long *fact, const long *mult, long len) {
//   proj tmp;
//   uintbig cof;
//   isog_degree deg = degree_co((isog_degree){ 0 }, mult, len);
//   for (int j = 0; j < len; j++) {
//     degree_unset(&deg, j);
//     degree_to_uint(&cof, deg, fact, len);
//     xMUL(&tmp, A, P, &cof);
//     uintbig_set(&cof, fact[j]);
//     uint8_t v = 0;
//     for ( ; !mont_iszero(&tmp); v++) {
//       xMUL(&tmp, A, &tmp, &cof);
//     }
//     degree_set(&deg, j, v);
//   }
//   return deg;
// }
// static char* fp2_hash(fp2 x) {
//     return pari_sprintf("h%lu", (x.re.x.c[0]+3*x.re.x.c[1]+5*x.re.x.c[2]+7*x.re.x.c[3]
//                     +11*x.im.x.c[0]+13*x.im.x.c[1]+17*x.im.x.c[2]+23*x.im.x.c[3]) % 100003);
// }
// static fp2 fp2_ratio(fp2 *x, fp2 *y) {
//     fp2 tmp;
//     tmp = *y;
//     assert(!fp2_iszero(&tmp));
//     fp2_inv(&tmp);
//     fp2_mul2(&tmp, x);
//     //printf("%lld %lld %lld %lld\n", tmp.re.x.c[0], tmp.re.x.c[1], tmp.im.x.c[0], tmp.im.x.c[1]);
//     return tmp;
// }
// static char* proj_hash(proj x) {
//     if (fp2_iszero(&x.z)) return "infty";
//     return fp2_hash(fp2_ratio(&x.x,&x.z));
// }
// // distorsion map on E0, montgomery form
// static void mont0_dist(proj *Q, const proj *P) {
//     proj Pcopy = *P;
//     fp2_neg2(&Q->x, &Pcopy.x);
//     Q->z = Pcopy.z;
// }
#endif


odd_isogeny trivial_odd_isogeny() {
    odd_isogeny phi;

    degree_one(&phi.deg_plus);
    degree_one(&phi.deg_minus);
    phi.kernel_plus.x = fp2_1;
    phi.kernel_plus.z = fp2_0;
    phi.kernel_minus.x = fp2_1;
    phi.kernel_minus.z = fp2_0;

    return phi;
}

special_isogeny trivial_special_isogeny() {
    special_isogeny phi;

    phi.source = global_setup.E0;
    phi.target = global_setup.E0;
    phi.middle = global_setup.E0;

    phi.phi2_set = true;
    phi.phi2_dual_set = true;

    phi.phi1 = trivial_odd_isogeny();
    phi.phi2 = trivial_odd_isogeny();
    phi.phi2_dual = trivial_odd_isogeny();

    return phi;
}


void free_two_walk_long(two_walk_long *phi) {
    if (phi->len) { free(phi->phi); phi->len = 0; }
}

void init_trivial_two_walk_long(two_walk_long *phi) {
    phi->phi = NULL;
    phi->len = 0;
}

static GEN alg_standard_to_O0(GEN elt) {
    return RgM_RgC_mul(global_setup.standard_to_O0, elt);
}


GEN kernel_to_ideal_gen_action(GEN v, GEN m1, GEN m2, GEN m3, GEN m4, long ell, long e) {

    if (e < 1) {
        return mkcol4s(1,0,0,0);
    }

    pari_sp ltop = avma;

    GEN gelle = powuu(ell,e);

    GEN v1 = gmul(m1, v);
    GEN v2 = gmul(m2, v);
    GEN v3 = gmul(m3, v);
    GEN v4 = gmul(m4, v);

    GEN matsys = mkmat4(v1,v2,v3,v4);
    GEN ker = matkermod(matsys, gelle, NULL);

    GEN sol, sol_reduced;

    int dim_ker = lg(ker)-1;


    if (dim_ker == 0) {
        avma = ltop;
        return mkcol4s(1,0,0,0);
    }
    for (int i = 1; i <= dim_ker; ++i){
        sol = gel(ker,i);
        sol_reduced = gmodgs(sol,ell);
        if (!isexactzero(sol_reduced)) break;
    }

    return gerepilecopy(ltop, sol);
}


GEN kernel_to_ideal_action_O0(GEN v, GEN m1, GEN m2, GEN m3, GEN m4, long ell, long e) {

    pari_sp ltop = avma;

    GEN gelle = powuu(ell,e);
    GEN gen_v = kernel_to_ideal_gen_action(v, m1, m2, m3, m4, ell, e);



    GEN generator = gmul(gel(gen_v,1),global_setup.O0_b1);
    generator = gadd(generator, gmul(gel(gen_v,2),global_setup.O0_b2));
    generator = gadd(generator, gmul(gel(gen_v,3),global_setup.O0_b3));
    generator = gadd(generator, gmul(gel(gen_v,4),global_setup.O0_b4));

    GEN ideal = lideal_create(global_setup.B, global_setup.O0, generator, gelle);

    return gerepilecopy(ltop, ideal);
}

// endo is an endomorphism expressed in the basis of the order whose action on the torsion correspond to m1 m2 m3 m4
// i.e. endo acts on the torsion as endo[0]*m1 + endo[1]*m2 + endo[2]*m3 + endo[3]*m4
GEN endo_to_kernel_action(GEN endo, GEN m1, GEN m2, GEN m3, GEN m4, long ell, long e) {
    pari_sp ltop = avma;



    GEN gelle = gpowgs(stoi(ell),e);

    GEN endo_m = gmul(gel(endo,1),m1);
    endo_m = gadd(endo_m, gmul(gel(endo,2),m2));
    endo_m = gadd(endo_m, gmul(gel(endo,3),m3));
    endo_m = gadd(endo_m, gmul(gel(endo,4),m4));



    GEN ker = matkermod(endo_m, gelle, NULL);

    GEN sol, sol_reduced;

    int dim_ker = lg(ker)-1;

    if (dim_ker == 0) {
        return mkcol2(gen_0,gen_0);
    }

    for (int i = 1; i <= dim_ker; ++i){
        sol = gel(ker,i);
        sol_reduced = gmodgs(sol,ell);
        if (!isexactzero(sol_reduced)) break;
        assert(i < dim_ker);
    }


    // remains to compute sol[1]*P1 + sol[2]*P2 where P1,P2 is a basis of the torsion

    return gerepilecopy(ltop, sol);
}

GEN ideal_to_kernel_action_O0(GEN I, GEN m1, GEN m2, GEN m3, GEN m4, long ell, long e) {
    pari_sp ltop = avma;

    GEN generator = lideal_generator(I);
    GEN endo = alg_standard_to_O0(generator);
    GEN sol = endo_to_kernel_action(endo, m1, m2, m3, m4, ell, e);
    // remains to compute sol[1]*P1 + sol[2]*P2 where P1,P2 is a basis of the torsion

    return gerepilecopy(ltop, sol);
}

void action_from_elle(GEN *m1, GEN *m2, GEN *m3, GEN *m4, long ell, long e) {
    pari_sp ltop = avma;

    GEN gelle = powuu(ell,e);

    *m1 = mkmat2(mkcol2s(1,0),mkcol2s(0,1));
    if (ell == 2) {
        *m2 = gmod(global_setup.action_two_2, gelle);
        *m3 = gmod(global_setup.action_two_3, gelle);
        *m4 = gmod(global_setup.action_two_4, gelle);

    }
    else {
        bool twist;
        unsigned long index = ell_to_index(ell, &twist);
        if (!twist) {
            *m2 = gmod(global_setup.action_2[index], gelle);
            *m3 = gmod(global_setup.action_3[index], gelle);
            *m4 = gmod(global_setup.action_4[index], gelle);
        }
        else {
            *m2 = gmod(global_setup.action_twist_2[index], gelle);
            *m3 = gmod(global_setup.action_twist_3[index], gelle);
            *m4 = gmod(global_setup.action_twist_4[index], gelle);
        }
    }

    gerepileall(ltop, 4, m1,m2,m3,m4);
}

GEN kernel_to_ideal_gen_O0_ell(GEN v, long ell, long *e) {
    pari_sp ltop = avma;
    long e_max = ell_to_e(ell);
    GEN gcd = ggcd(gel(v,1),gel(v,2));
    long e_diff = (isexactzero(gcd)) ? e_max : Z_lval(gcd, ell);

    if (e_diff >= e_max) {
        avma = ltop;
        *e = 0;
        return mkcol4s(1,0,0,0);
    }

    GEN generator;
    GEN m1,m2,m3,m4;
    *e = (e_diff < e_max) ? e_max - e_diff : 0;
    action_from_elle(&m1, &m2, &m3, &m4, ell, *e);
    v = gdiv(v, powuu(ell,e_diff)); // from ell^e_max torsion to ell^e torsion
    generator = kernel_to_ideal_gen_action(v, m1, m2, m3, m4, ell, *e);
    return gerepilecopy(ltop, generator);
}

GEN kernel_to_ideal_O0_ell(GEN v, long ell) {
    pari_sp ltop = avma;
    long e_max = ell_to_e(ell);
    GEN gcd = ggcd(gel(v,1),gel(v,2));
    long e_diff = (isexactzero(gcd)) ? e_max : Z_lval(gcd, ell);

    if (e_diff >= e_max) {
        avma = ltop;
        return lideal_create(global_setup.B, global_setup.O0, alg_scalar(global_setup.B,gen_1), gen_1);
    }

    GEN I;
    GEN m1,m2,m3,m4;
    long e = (e_diff < e_max) ? e_max - e_diff : 0;
    action_from_elle(&m1, &m2, &m3, &m4, ell, e);
    v = gdiv(v, powuu(ell,e_diff)); // from ell^e_max torsion to ell^e torsion
    I = kernel_to_ideal_action_O0(v, m1, m2, m3, m4, ell, e);
    return gerepilecopy(ltop, I);
}

// assume I is a cyclic ideal
GEN ideal_to_kernel_O0_ell(GEN I, long ell) {
    pari_sp ltop = avma;
    long e = Z_lval(lideal_norm(I), ell);
    long e_max = ell_to_e(ell);
    long e_diff = e_max - e;
    GEN v;
    GEN m1,m2,m3,m4;
    action_from_elle(&m1, &m2, &m3, &m4, ell, e);
    v = ideal_to_kernel_action_O0(I, m1, m2, m3, m4, ell, e);
    v = gmul(v, powuu(ell,e_diff)); // from ell^e torsion to ell^e_max torsion
    return gerepilecopy(ltop, v);
}

// assume I is a cyclic ideal
GEN endo_to_kernel_O0_ell(GEN alpha, long ell, long e) {
    pari_sp ltop = avma;
    long e_max = ell_to_e(ell);
    long e_diff = e_max - e;
    GEN v;
    GEN m1,m2,m3,m4;
    action_from_elle(&m1, &m2, &m3, &m4, ell, e);

    GEN endo = alg_standard_to_O0(alpha);

    v = endo_to_kernel_action(endo, m1, m2, m3, m4, ell, e);
    v = gmul(v, powuu(ell,e_diff)); // from ell^e torsion to ell^e_max torsion
    return gerepilecopy(ltop, v);
}



GEN ideal_to_kernel_O0_T(GEN I, GEN fact_norm) {
    pari_sp ltop = avma;
    long len = lg(gel(fact_norm,1)), ell, e_max, e;
    GEN coeff_1 = zerovec(p_plus_len), coeff_2 = zerovec(p_plus_len);
    GEN coeff_twist_1 = zerovec(p_minus_len), coeff_twist_2 = zerovec(p_minus_len);
    GEN alpha = lideal_generator(I);
    long index;
    bool twist;
    GEN v;

    for (int i = 1; i < len; ++i) {
        ell = itos(gel(gel(fact_norm,1),i));
        e_max = itos(gel(gel(fact_norm,2),i));
        index = ell_to_index(ell, &twist);
        e = Z_lval(lideal_norm(I), ell);
        v = endo_to_kernel_O0_ell(alpha, ell, e);
        if (!twist) {
            gel(coeff_1,index+1) = gel(v,1);
            gel(coeff_2,index+1) = gel(v,2);
        }
        else {
            gel(coeff_twist_1,index+1) = gel(v,1);
            gel(coeff_twist_2,index+1) = gel(v,2);
        }
    }

    return gerepilecopy(ltop, mkvec2(mkvec2(coeff_1,coeff_2),mkvec2(coeff_twist_1,coeff_twist_2)));
}

GEN torsion_crt_compose (GEN coeff, bool twist) {
    pari_sp ltop = avma;
    GEN p_primary = (twist) ? global_setup.gen_p_minus_primary : global_setup.gen_p_plus_primary;

    GEN a = ZV_chinese(gel(coeff,1), p_primary, NULL);
    GEN b = ZV_chinese(gel(coeff,2), p_primary, NULL);

    return gerepilecopy(ltop, mkcol2(a,b));
}

GEN torsion_crt_decompose (GEN v, bool twist) {
    pari_sp ltop = avma;
    long len = (twist) ? p_minus_len : p_plus_len;
    GEN p_primary = (twist) ? global_setup.gen_p_minus_primary : global_setup.gen_p_plus_primary;

    GEN coeff_1 = cgetg(len+1, t_VEC);
    GEN coeff_2 = cgetg(len+1, t_VEC);

    for (int i = 0; i < len; ++i) {
        gel(coeff_1,i+1) = gmod(gel(v,1), gel(p_primary,i+1));
        gel(coeff_2,i+1) = gmod(gel(v,2), gel(p_primary,i+1));
    }

    return gerepilecopy(ltop, mkvec2(coeff_1,coeff_2));
}

void gentobig(uintbig *res, GEN a) {
    assert(gsigne(a) >= 0);
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
    assert(isexactzero(b));
    avma = ltop;
}

proj vec_to_E0(GEN v, bool twist) {
    proj pt;
    const proj *basis = (twist) ? torsion_basis_twist_sum : torsion_basis_sum;
    uintbig x, y;

    gentobig(&x, gel(v,1));
    gentobig(&y, gel(v,2));

    // TODO: replace the following with a small multiplication and a combination of the form P + yQ (i.e., x = 1)
    xBIDIM(&pt, &global_setup.E0, &basis[0], &x, &basis[1], &y, &basis[2]);

    return pt;
}


void famat_to_degree(isog_degree *deg, isog_degree *deg_twist, GEN f) {
    degree_one(deg);
    degree_one(deg_twist);

    if (lg(f) == 1) return;

    long m = lg(gel(f,1)) - 1;
    long ell, e, index;
    bool twist;

    for (long i = 1; i <= m; ++i) {
        ell = itos_or_0(gel(gel(f,1),i));
        e = itos_or_0(gel(gel(f,2),i));
        index = ell_to_index(ell, &twist);
        if (!twist) {
            degree_set(deg, index, e);
        }
        else {
            degree_set(deg_twist, index, e);
        }
    }
}

proj coeff_to_E0(GEN coeff, bool twist) {
    GEN v = torsion_crt_compose (coeff, twist);
    return vec_to_E0(v, twist);
}

odd_isogeny ideal_to_isogeny_O0_T(GEN I, GEN factorisation_norm) {
    odd_isogeny phi;
    GEN coeff = ideal_to_kernel_O0_T(I,factorisation_norm);
    proj ker = coeff_to_E0(gel(coeff,1), false);
    proj ker_twist = coeff_to_E0(gel(coeff,2), true);
    isog_degree deg, deg_twist;
    famat_to_degree(&deg, &deg_twist, factorisation_norm);

    phi.kernel_plus = ker;

    phi.kernel_minus = ker_twist;
    phi.deg_plus = deg;
    phi.deg_minus = deg_twist;

    return phi;
}



two_walk ideal_to_isogeny_O0_two(GEN I) {
    pari_sp ltop = avma;
    two_walk phi;
    phi.A = global_setup.E0;

    GEN v = ideal_to_kernel_O0_ell(I, 2);

    const proj *basis = torsion_basis_two;
    uintbig x, y;

    gentobig(&x, gel(v,1));
    gentobig(&y, gel(v,2));

    // TODO: replace the following with a small multiplication and a combination of the form P + yQ (i.e., x = 1)
    xBIDIM(&phi.ker, &global_setup.E0, &basis[0], &x, &basis[1], &y, &basis[2]);



    phi.len = Z_lval(lideal_norm(I), 2);

    avma = ltop;
    return phi;
}




GEN kernel_to_ideal_O0_T(GEN coeff) {
    pari_sp ltop = avma;

    GEN I, v, N = gen_1;
    long ell, e, e_max;

    GEN x, ideal_generator_decomposed = mkvec4(
        cgetg(p_plus_len + p_minus_len + 1, t_VEC),
        cgetg(p_plus_len + p_minus_len + 1, t_VEC),
        cgetg(p_plus_len + p_minus_len + 1, t_VEC),
        cgetg(p_plus_len + p_minus_len + 1, t_VEC));

    for (int i = 0; i < p_plus_len; ++i) {

        ell = p_plus_fact[i];
        e_max = p_plus_mult[i];
        v = mkcol2(gel(gel(gel(coeff,1),1),i+1), gel(gel(gel(coeff,1),2),i+1));
        x = kernel_to_ideal_gen_O0_ell(v, ell, &e);
        N = gmul(N, powuu(ell,e));
        gel(gel(ideal_generator_decomposed,1),i+1) = gel(x,1);
        gel(gel(ideal_generator_decomposed,2),i+1) = gel(x,2);
        gel(gel(ideal_generator_decomposed,3),i+1) = gel(x,3);
        gel(gel(ideal_generator_decomposed,4),i+1) = gel(x,4);
    }

    for (int i = 0; i < p_minus_len; ++i) {
        ell = p_minus_fact[i];
        e_max = p_minus_mult[i];
        v = mkcol2(gel(gel(gel(coeff,2),1),i+1), gel(gel(gel(coeff,2),2),i+1));
        x = kernel_to_ideal_gen_O0_ell(v, ell, &e);
        N = gmul(N, powuu(ell,e));
        gel(gel(ideal_generator_decomposed,1),p_plus_len+i+1) = gel(x,1);
        gel(gel(ideal_generator_decomposed,2),p_plus_len+i+1) = gel(x,2);
        gel(gel(ideal_generator_decomposed,3),p_plus_len+i+1) = gel(x,3);
        gel(gel(ideal_generator_decomposed,4),p_plus_len+i+1) = gel(x,4);
    }

    GEN p_primary = shallowconcat(global_setup.gen_p_plus_primary, global_setup.gen_p_minus_primary);

    GEN a = ZV_chinese(gel(ideal_generator_decomposed,1), p_primary, NULL);
    GEN b = ZV_chinese(gel(ideal_generator_decomposed,2), p_primary, NULL);
    GEN c = ZV_chinese(gel(ideal_generator_decomposed,3), p_primary, NULL);
    GEN d = ZV_chinese(gel(ideal_generator_decomposed,4), p_primary, NULL);

    GEN generator = gmul(a,global_setup.O0_b1);
    generator = gadd(generator, gmul(b,global_setup.O0_b2));
    generator = gadd(generator, gmul(c,global_setup.O0_b3));
    generator = gadd(generator, gmul(d,global_setup.O0_b4));

    I = lideal_create(global_setup.B, global_setup.O0, generator, N);

    return gerepilecopy(ltop, I);
}



// Evaluate special isogeny phi : source -> target at point P.
// sets P to the image point
proj eval_special(proj *range, special_isogeny *phi, const proj *P) {
    proj A = phi->source;
    proj Q = *P;
    isomorphism isom;

    eval(&A, &phi->phi1, &Q);

    if (!phi->phi2_set) {
        assert(phi->phi2_dual_set);
        phi->phi2 = phi->phi2_dual;
        phi->middle = phi->target;
        dual(&phi->middle, &phi->phi2);
        phi->phi2_set = true;
    }

    mont_isom(&isom, &A, &phi->middle);
    mont_isom_apply(&isom, &Q);
    A = phi->middle;
    eval(&A, &phi->phi2, &Q);

    *range = A;

    return Q;
}



void two_walk_stol(two_walk_long *phil, const two_walk *phi) {
    if (phi->len  == 0) { free_two_walk_long(phil); init_trivial_two_walk_long(phil); return; }
    two_walk_long res;
    res.len = 1;
    res.phi = malloc(sizeof(two_walk)*(1));
    res.phi[0] = *phi;
    free_two_walk_long(phil);
    *phil = res;
}

void two_walk_composition_ll(two_walk_long *phi, const two_walk_long *phi2, const two_walk_long *phi1){
    two_walk_long res;
    res.len = phi1->len + phi2->len;
    res.phi = malloc(sizeof(two_walk)*(res.len));
    for (int i = 0; i < phi1->len; ++i) {
        res.phi[i] = phi1->phi[i];
    }
    for (int i = 0; i < phi2->len; ++i) {
        res.phi[i+phi1->len] = phi2->phi[i];
    }
    free_two_walk_long(phi);
    *phi = res;
}

void copy_two_walk_long(two_walk_long *copy, const two_walk_long *phi) {
    two_walk_long triv;
    init_trivial_two_walk_long(&triv);
    two_walk_composition_ll(copy, phi, &triv);
}

void two_walk_composition_ls(two_walk_long *phi, const two_walk_long *phi2, const two_walk *phi1) {
    if (phi1->len == 0) { copy_two_walk_long(phi, phi2); return; }
    two_walk_long phi1l;
    init_trivial_two_walk_long(&phi1l);
    two_walk_stol(&phi1l, phi1);
    two_walk_composition_ll(phi, phi2, &phi1l);
    free_two_walk_long(&phi1l);
}

void two_walk_composition_sl(two_walk_long *phi, const two_walk *phi2, const two_walk_long *phi1) {
    if (phi2->len == 0) { copy_two_walk_long(phi, phi1); return; }
    two_walk_long phi2l;
    init_trivial_two_walk_long(&phi2l);
    two_walk_stol(&phi2l, phi2);
    two_walk_composition_ll(phi, &phi2l, phi1);
    free_two_walk_long(&phi2l);
}

void two_walk_composition_ss(two_walk_long *phi, const two_walk *phi2, const two_walk *phi1) {
    if (phi2->len == 0) {
        if (phi1->len == 0) { free_two_walk_long(phi); init_trivial_two_walk_long(phi); }
        else { two_walk_stol(phi, phi1); return; }
    }
    if (phi1->len == 0) { two_walk_stol(phi, phi2); return; }
    two_walk_long phi1l;
    init_trivial_two_walk_long(&phi1l);
    two_walk_stol(&phi1l, phi1);
    two_walk_long phi2l;
    init_trivial_two_walk_long(&phi2l);
    two_walk_stol(&phi2l, phi2);
    two_walk_composition_ll(phi, &phi2l, &phi1l);
    free_two_walk_long(&phi1l);
    free_two_walk_long(&phi2l);
}

// tested
odd_isogeny push_odd_isogeny_through_two_walk(const odd_isogeny *phi_odd, proj *phi_odd_source, const two_walk *phi_two) {
    odd_isogeny phi_odd_isom = *phi_odd;

    if (phi_two->len == 0) { return phi_odd_isom;}
    two_walk phi_two_isom = *phi_two;
    isomorphism isom1, isom2;

    #ifndef NDEBUG
    proj j1,j2;
    jinv256(&j1, phi_odd_source);
    jinv256(&j2, &phi_two->A);
    assert(mont_equal(&j1,&j2));
    #endif

    if (!mont_equal(phi_odd_source,&phi_two_isom.A)) {

        #ifndef NDEBUG
        // Check the j-invariant is not 0 or 1728 (otherwise the isomorphism is ambiguous)
        proj j1728;
        fp2_set(&j1728.x,1728);
        fp2_set(&j1728.z,256);
        assert(!mont_equal(&j1,&j1728) && !mont_iszero(&j1));
        #endif

        mont_isom(&isom2, phi_odd_source, &phi_two_isom.A);
        mont_isom_apply(&isom2, &phi_odd_isom.kernel_plus);
        mont_isom_apply(&isom2, &phi_odd_isom.kernel_minus);

    }


    eval_walk_isom(&isom1, &phi_two_isom, NULL, NULL, &phi_two_isom, NULL);

    mont_isom_apply(&isom1, &phi_odd_isom.kernel_plus);
    mont_isom_apply(&isom1, &phi_odd_isom.kernel_minus);


    #ifndef NDEBUG
    uintbig k;
    proj P = phi_two_isom.ker;
    gentobig(&k, powuu(2,phi_two_isom.len-1));
    xMUL(&P, &phi_two_isom.A, &P, &k);
    assert(!mont_iszero(&P));
    assert(!fp2_iszero(&P.x));
    gentobig(&k, powuu(2,1));
    xMUL(&P, &phi_two_isom.A, &P, &k);
    assert(mont_iszero(&P));

    jinv256(&j1, phi_odd_source);
    jinv256(&j2, &phi_two_isom.A);
    assert(mont_equal(&j1,&j2));
    #endif


    odd_isogeny res = phi_odd_isom;
    eval_walk(&phi_two_isom, phi_odd_source, &res.kernel_plus);
    eval_walk(&phi_two_isom, phi_odd_source, &res.kernel_minus);
    return res;
}

// tested
odd_isogeny push_odd_isogeny_through_two_walk_long(const odd_isogeny *phi_odd, proj *phi_odd_source, const two_walk_long *phi_two) {
    odd_isogeny res = *phi_odd;
    for (int i = 0; i < phi_two->len; ++i) {
        res = push_odd_isogeny_through_two_walk(&res, phi_odd_source, &phi_two->phi[i]);
    }
    return res;
}

// tested
two_walk push_two_walk_through_odd_isogeny(const two_walk *phi_two, const odd_isogeny *phi_odd, const proj *phi_odd_source) {
    two_walk res = *phi_two;
    isomorphism isom;

    mont_isom(&isom, &phi_two->A, phi_odd_source);
    mont_isom_apply(&isom, &res.ker);
    //printf("isom %s %s %s\n", fp2_hash(isom.Nx), fp2_hash(isom.Nz), fp2_hash(isom.D));

    res.A = *phi_odd_source;

    eval(&res.A, phi_odd, &res.ker);

    return res;
}

void push_two_walk_long_through_odd_isogeny(two_walk_long *phi, const two_walk_long *phi_two, const odd_isogeny *phi_odd, const proj *phi_odd_source) {
    two_walk_long res;

    init_trivial_two_walk_long(&res);
    copy_two_walk_long(&res, phi_two);
    odd_isogeny phi_odd_i = *phi_odd;
    proj phi_odd_source_i = *phi_odd_source;


    for (int i = 0; i < phi_two->len; ++i) {
        res.phi[i] = push_two_walk_through_odd_isogeny(&res.phi[i], &phi_odd_i, &phi_odd_source_i);
        if (i+1 < phi_two->len) {
            #ifndef NDEBUG
            proj j1,j2;
            jinv256(&j1, &phi_odd_source_i);
            jinv256(&j2, &phi_two->phi[i].A);
            assert(mont_equal(&j1,&j2));
            #endif
            phi_odd_i = push_odd_isogeny_through_two_walk(&phi_odd_i, &phi_odd_source_i, &phi_two->phi[i]);
        }
    }

    free_two_walk_long(phi);
    copy_two_walk_long(phi, &res);
}

two_walk push_two_walk_through_special_isogeny(const two_walk *phi_two, special_isogeny *phi_special) {
    two_walk res = *phi_two;
    proj E0 = phi_special->source;
    isomorphism isom;

    mont_isom(&isom, &phi_two->A, &E0);
    mont_isom_apply(&isom, &res.ker);

    res.ker = eval_special(&res.A, phi_special, &res.ker);

    return res;
}


special_isogeny special_ideal_to_isogeny(GEN J, GEN I, const two_walk_long *phi_I) {
    pari_sp ltop = avma;
    special_isogeny phi;
    GEN A = lideal_algebra(J);
    GEN order = lideal_order(J);
    GEN H1 = lideal_create(A, order, lideal_generator(J), ggcd(global_setup.gen_odd_torsion,lideal_norm(J)));
    GEN beta = lideal_isom(I, J); // I*beta = J
    GEN alpha = gmul(alg_conj(A,beta), lideal_norm(I)); // J = chi_I(alpha) = I conj(alpha)/n(I)
    GEN H2 = lideal_create(A, order, alpha, gdiv(lideal_norm(J),lideal_norm(H1)));

    GEN fm = famat_mul(global_setup.gen_p_plus_fact,global_setup.gen_p_minus_fact);
    GEN fm_norm1 = famat_Z_gcd(fm, lideal_norm(H1));
    phi.phi1 = ideal_to_isogeny_O0_T(H1, fm_norm1);
    phi.source = global_setup.E0;

    GEN fm_norm2 = famat_Z_gcd(fm, lideal_norm(H2));
    phi.phi2_dual = ideal_to_isogeny_O0_T(H2, fm_norm2);
    phi.target = global_setup.E0;

    // push phi->phi2 though phi_I
    phi.phi2_dual = push_odd_isogeny_through_two_walk_long(&phi.phi2_dual, &phi.target, phi_I);

    #ifndef NDEBUG
    odd_isogeny X = phi.phi1,Y = phi.phi2_dual;
    proj E1 = phi.source, E2 = phi.target;
    dual(&E1, &X);
    dual(&E2, &Y);
    assert(mont_equal(&E1,&E2));
    #endif


    phi.phi2_dual_set = true;
    phi.phi2_set = false;
    avma = ltop;
    return phi;
}

void eval_walk_long_mult(const two_walk_long *phi, proj *B, proj *P, long cardinality) {
    isomorphism isom;

    two_walk dummy;
    for (int i = 0; i < phi->len; ++i) {
        mont_isom(&isom, B, &phi->phi[i].A);
        for (int i = 0; i < cardinality; ++i) {
            mont_isom_apply(&isom, &P[i]);
        }
        eval_walk_isom_mult(&isom, &dummy, B, &phi->phi[i], P, cardinality);
    }
}


GEN dual_coeff(GEN coeff, isog_degree deg_plus, isog_degree deg_minus) {
    pari_sp ltop = avma;
    GEN coeff_dual = coeff;
    for (int i = 0; i < p_plus_len; ++i) {
        long ell = p_plus_fact[i];
        long e = p_plus_mult[i];
        GEN c1 = gel(gel(gel(coeff,1), 1), i+1);
        GEN c2 = gel(gel(gel(coeff,1), 2), i+1);
        // printf("coeff %ld^%ld\n", ell, e);
        // output(mkvec2(c1,c2));
        long d = degree_get(deg_plus, i);
        if (Z_lval(c1, ell) == e-d) {
            gel(gel(gel(coeff_dual,1), 1), i+1) = gen_0;
            gel(gel(gel(coeff_dual,1), 2), i+1) = (d == 0) ? gen_1 : powuu(ell,e-d);
        }
        else if (Z_lval(c2, ell) == e-d) {
            gel(gel(gel(coeff_dual,1), 1), i+1) = (d == 0) ? gen_1 : powuu(ell,e-d);
            gel(gel(gel(coeff_dual,1), 2), i+1) = gen_0;
        }
        else {
            assert(isexactzero(c1) && isexactzero(c2));
            gel(gel(gel(coeff_dual,1), 1), i+1) = gen_0;
            gel(gel(gel(coeff_dual,1), 2), i+1) = gen_0;
        }

    }


    for (int i = 0; i < p_minus_len; ++i) {
        long ell = p_minus_fact[i];
        long e = p_minus_mult[i];
        GEN c1 = gel(gel(gel(coeff,2), 1), i+1);
        GEN c2 = gel(gel(gel(coeff,2), 2), i+1);
        long d = degree_get(deg_minus, i);
        if (Z_lval(c1, ell) == e-d) {
            gel(gel(gel(coeff_dual,2), 1), i+1) = gen_0;
            gel(gel(gel(coeff_dual,2), 2), i+1) = (d == 0) ? gen_1 : powuu(ell,e-d);
        }
        else if (Z_lval(c2, ell) == e-d) {
            gel(gel(gel(coeff_dual,2), 1), i+1) = (d == 0) ? gen_1 : powuu(ell,e-d);
            gel(gel(gel(coeff_dual,2), 2), i+1) = gen_0;
        }
        else {
            assert(isexactzero(c1) && isexactzero(c2));
            gel(gel(gel(coeff_dual,2), 1), i+1) = gen_0;
            gel(gel(gel(coeff_dual,2), 2), i+1) = gen_0;
        }
    }
    return gerepilecopy(ltop, coeff_dual);
}

// T = global_setup.gen_odd_torsion
// f = two_tors_height
// I is a left O_0-ideal of norm dividing T^2 \ell^{2f+delta}
// J is a left O_0-ideal containing I of norm gcd(T^2,n(I))
// K is a left O_0-ideal equivalent to J of norm a power of 2
// Finds phi such that phi_I = phi * phi_J
// Finds L equivalent to I of norm dividing T^2
// phi_K_basis is the image through phi_K of a basis of the odd torsion, then
// phi is applied to it (a basis = 6 points: A, B and A+B for the curve and its twist)
void ideal_to_isogeny_two_2f_delta(two_walk_long *phi, GEN *L,
    special_isogeny *phi_L, GEN I, GEN J, GEN K,
    proj *phi_K_basis, proj *phi_K_target, int delta, GEN I_long) {
    pari_sp ltop = avma;
    GEN A = lideal_algebra(J);
    GEN order = lideal_order(J);
    isomorphism isom;

    #ifndef NDEBUG
    proj j1,j2;
    #endif


    long len_step = Z_lval(lideal_norm(I), 2);
    long e1 = (len_step < two_tors_height) ? len_step : two_tors_height;
    if (e1 > len_step) e1 = len_step;
    long e2 = len_step - delta - e1;
    if (e2 < 0) e2 = 0;

    long dist = len_step - e1 - e2;
    if (dist < 0) dist = 0;

    if (dist > 22) fprintf(stderr,"Warning: MITM distance is %ld\n", dist);



    GEN L_;

    if (len_step + Z_lval(lideal_norm(K), 2) < dbllog2r(itor(global_setup.p,10))/2. + 10 ) {
        GEN M;
        GEN alpha = lideal_isom(J, K); // J*alpha = K

        M = lideal_mul(I, alpha); // I*alpha, equivalent to I, but norm a power of 2

        L_ = klpt_special_smooth_small_2e_input(M, famat_sqr(famat_mul(global_setup.gen_p_plus_fact, global_setup.gen_p_minus_fact)));
    }
    else {
        L_ = klpt_special_smooth(I, famat_sqr(famat_mul(global_setup.gen_p_plus_fact, global_setup.gen_p_minus_fact)));
    }

    assert(lideal_isom(I, L_));

    GEN a = lideal_isom(J, K); // J*a = K
    if (gcmp(lideal_norm(K), gen_1) == 0) { a = alg_scalar(A,gen_1); /* make sure we don't apply a distorsion */ }

    GEN M = lideal_mul(I, a);
    assert(lideal_isom(L_,M));
    GEN b = lideal_isom(L_,M); // L_*gamma = M
    GEN gamma = gmul(b, lideal_norm(L_));

    assert(alglatcontains(A, lideal_lattice(K), gamma, NULL));
    assert(alglatcontains(A, lideal_lattice(L_), alg_conj(A, gamma), NULL));
    assert(gcmp(algnorm(A,gamma,0), gmul(powuu(2,len_step),gmul(lideal_norm(K),lideal_norm(L_)))) == 0);

    GEN n;
    alg_primitive(&n, A, order, gamma);
    assert(gcmp(n,gen_1) == 0);

    // GEN H1_two = lideal_create(A, order, gamma, gmul(lideal_norm(K),powuu(2, e1)));


    GEN H1_odd = lideal_create(A, order, gamma, ggcd(global_setup.gen_odd_torsion, lideal_norm(L_)));



    odd_isogeny psi_1; //= ideal_to_isogeny_O0_T(H1_odd, famat_Z_gcd(famat_mul(global_setup.gen_p_plus_fact, global_setup.gen_p_minus_fact),lideal_norm(H1_odd)));


    GEN factorisation_norm = famat_Z_gcd(famat_mul(global_setup.gen_p_plus_fact, global_setup.gen_p_minus_fact),lideal_norm(H1_odd));
    GEN coeff = ideal_to_kernel_O0_T(H1_odd,factorisation_norm);

    GEN v_plus = torsion_crt_compose(gel(coeff,1), false);
    GEN v_minus = torsion_crt_compose(gel(coeff,2), true);

    proj ker_plus, ker_minus;

    uintbig x, y;
    gentobig(&x, gel(v_plus,1));
    gentobig(&y, gel(v_plus,2));
    xBIDIM(&ker_plus, phi_K_target, &phi_K_basis[0], &x, &phi_K_basis[1], &y, &phi_K_basis[2]);
    gentobig(&x, gel(v_minus,1));
    gentobig(&y, gel(v_minus,2));
    xBIDIM(&ker_minus, phi_K_target, &phi_K_basis[3], &x, &phi_K_basis[4], &y, &phi_K_basis[5]);

    isog_degree deg_plus, deg_minus;
    famat_to_degree(&deg_plus, &deg_minus, factorisation_norm);

    psi_1.kernel_plus = ker_plus;
    psi_1.kernel_minus = ker_minus;
    psi_1.deg_plus = deg_plus;
    psi_1.deg_minus = deg_minus;
    proj psi_1_source = *phi_K_target;


    odd_isogeny psi_1_dual = psi_1;
    proj psi_1_dual_source = psi_1_source;


    if ((psi_1_dual.deg_plus.val != 0) || (psi_1_dual.deg_minus.val != 0)) {

        proj ker_dual[2];

        GEN coeff_dual = dual_coeff(coeff, psi_1.deg_plus, psi_1.deg_minus);

        GEN v_plus_dual = torsion_crt_compose(gel(coeff_dual,1), false);
        gentobig(&x, gel(v_plus_dual,1));
        gentobig(&y, gel(v_plus_dual,2));
        xBIDIM(&ker_dual[0], &psi_1_source, &phi_K_basis[0], &x, &phi_K_basis[1], &y, &phi_K_basis[2]);


        GEN v_minus_dual = torsion_crt_compose(gel(coeff_dual,2), true);
        gentobig(&x, gel(v_minus_dual,1));
        gentobig(&y, gel(v_minus_dual,2));
        xBIDIM(&ker_dual[1], &psi_1_source, &phi_K_basis[3], &x, &phi_K_basis[4], &y, &phi_K_basis[5]);


        eval_mult(&psi_1_dual_source, &psi_1, ker_dual, 2);

        psi_1_dual.kernel_plus = ker_dual[0];
        psi_1_dual.kernel_minus = ker_dual[1];

    }


    GEN gamma_conj = alg_conj(A, gamma);
    GEN H2_odd = lideal_create(A, order, gamma_conj, gdiv(lideal_norm(L_), lideal_norm(H1_odd)));
    GEN H2_two = lideal_create(A, order, gamma_conj, powuu(2, e2));

    odd_isogeny psi_2 = ideal_to_isogeny_O0_T(H2_odd, famat_Z_gcd(famat_mul(global_setup.gen_p_plus_fact, global_setup.gen_p_minus_fact),lideal_norm(H2_odd)));




    // TODO: redundant computation
    GEN beta = lideal_isom(I, L_); // I*beta = L_
    GEN I_long_next = lideal_mul(I_long, beta); // I_i_long*beta
    GEN I1_next = lideal_create(A, order, lideal_generator(I_long_next), powuu(2, two_tors_height));

    two_walk phi_1_next = ideal_to_isogeny_O0_two(I1_next);



    two_walk phi_2 = ideal_to_isogeny_O0_two(H2_two);

    proj pt;


    proj pt2[2];
    pt2[0] = phi_1_next.ker;
    pt2[1] = phi_2.ker;

    eval_mult(&phi_1_next.A, &psi_2, pt2, 2);

    phi_1_next.ker = pt2[0];
    phi_2.ker = pt2[1];
    phi_2.A = phi_1_next.A;


    two_walk phi_2_adjusted;// = phi_2; // push kernel away from (0:1)
    proj phi_2_adjusted_target; // = phi_2.A;

    pt = phi_2.ker;
    eval_walk_isom(&isom, &phi_2_adjusted, &phi_2_adjusted_target, &pt, &phi_2, &pt);


    two_walk eta;
    two_walk phi_2_dual = phi_2_adjusted;

    dual_walk(&phi_2_dual);


    // Check dual
    #ifndef NDEBUG
    proj pt_save;
    proj B = phi_2_adjusted.A;
    pt.z = fp2_1;
    do {
        fp2_random(&pt.x);
    } while(!is_on_curve(&pt, &phi_2_adjusted.A));

    pt_save = pt;


    eval_walk(&phi_2_adjusted, &B, &pt);

    mont_isom(&isom, &B, &phi_2_dual.A);
    mont_isom_apply(&isom, &pt);

    B = phi_2_dual.A;
    eval_walk(&phi_2_dual, &B, &pt);

    mont_isom(&isom, &B, &phi_2_adjusted.A);
    mont_isom_apply(&isom, &pt);

    jinv256(&j1, &phi_2_adjusted_target);
    jinv256(&j2, &phi_2_dual.A);
    assert(mont_equal(&j1,&j2));
    #endif


    proj from = psi_1_dual_source, from0 = psi_1_dual_source;
    proj to = phi_2_adjusted_target; // phi_2 source
    two_walk_long phi_2_dual_eta;
    init_trivial_two_walk_long(&phi_2_dual_eta);

    // float accumulated_time = 0.;
    // clock_t t;
    // t = tic();
    if (dist > 0) {
        bool done;
        //rand_isom(&isom, &to);
        assert(dist != 1); // for now, the case dist == 1 crashes
        done = MITM2(&eta, &from, &to, dist);
        assert(done);
        two_walk_composition_ss(&phi_2_dual_eta, &phi_2_dual, &eta);
    }
    else {
        two_walk_stol(&phi_2_dual_eta, &phi_2_dual);
    }
    //TOC(t, "MITM");





    phi_L->source = global_setup.E0;
    phi_L->phi1 = psi_2;

    phi_L->phi2_dual_set = false;

    // since psi_1_dual has already been computed...
    phi_L->middle = psi_1_dual_source; // = psi_1_dual_source
    phi_L->phi2 = push_odd_isogeny_through_two_walk_long(&psi_1_dual, &phi_L->middle, &phi_2_dual_eta);

    phi_L->phi2_set = true;



    // push phi_2 and phi_1_next simultaneously through phi_L->phi2
    two_walk phi2_pushed = phi_2;


    isomorphism isom3;
    mont_isom(&isom3, &phi_1_next.A, &phi_L->middle);
    mont_isom_apply(&isom3, &phi_1_next.ker);
    phi_1_next.A = phi_L->middle;



    mont_isom_apply(&isom3, &phi2_pushed.ker);
    phi2_pushed.A = phi_L->middle;

    pt2[0] = phi2_pushed.ker;
    pt2[1] = phi_1_next.ker;
    eval_mult(&phi2_pushed.A, &phi_L->phi2, pt2, 2);

    phi2_pushed.ker = pt2[0];
    phi_1_next.ker = pt2[1];
    phi_1_next.A = phi2_pushed.A;


    rand_isom(&isom, &phi2_pushed.A);
    mont_isom_apply(&isom, &phi2_pushed.ker);

    proj phi2_pushed_target, proj_tmp;
    eval_walk(&phi2_pushed, &phi2_pushed_target, &proj_tmp);

    two_walk phi2_pushed_dual = phi2_pushed;
    dual_walk(&phi2_pushed_dual);



    // float accumulated_time = 0.;
    // clock_t t;
    // t = tic();

    from = psi_1_source, from0 = psi_1_source;
    to = phi2_pushed_target; // phi_2 source
    two_walk mitm_up;
    if (dist > 0) {
        bool done;
        done = MITM2(&mitm_up, &from, &to, dist);
        assert(done);
        two_walk_composition_ss(phi, &phi2_pushed_dual, &mitm_up);
    }
    else {
        two_walk_stol(phi, &phi2_pushed_dual);
    }
    //TOC(t, "MITM");


    two_walk_composition_sl(phi, &phi_1_next, phi);

    eval_walk_long_mult(phi, phi_K_target, phi_K_basis, 6);

    *L = gerepilecopy(ltop, L_);
    assert(lideal_isom(I, *L));

    free_two_walk_long(&phi_2_dual_eta);
}








// T = global_setup.gen_odd_torsion
// I is a left O0-ideal of norm dividing T^2 2^e for some positive integer e
// J = I + O0*T^2
// K is a left O0-ideal equivalent to J of norm a power of 2
// Finds phi such that phi_I = phi * phi_J
// Finds L equivalent to I of norm dividing T^2
void ideal_to_isogeny_two(two_walk_long *phi_res, GEN *L, special_isogeny *phi_L,
    GEN I, GEN J, GEN K, const special_isogeny *phi_J, const two_walk_long *phi_K,
    bool endpoint_close_to_E0){
    pari_sp ltop = avma;
    GEN A = lideal_algebra(J);
    GEN order = lideal_order(J);

    assert(lideal_isom(J, K));


    #ifndef NDEBUG
    GEN X = lideal_create(A, order, lideal_generator(I), lideal_norm(J));
    assert(gcmp(algnorm(A,lideal_isom(J, X),0),gen_1) == 0);
    assert(lideal_isom(J, K));

    if (phi_K->len > 0) {
        proj j1,j2;
        proj P = phi_K->phi[phi_K->len-1].ker;
        proj E = phi_K->phi[phi_K->len-1].A;
        jinv256(&j1, &phi_J->target);
        eval_walk(&phi_K->phi[phi_K->len-1], &E, &P);
        jinv256(&j2, &E);
        assert(mont_equal(&j1,&j2));
        assert(phi_J->phi2_set);
        assert(phi_J->phi2_dual_set);
    }

    long len = 0;
    for (int i = 0; i < phi_K->len; ++i) {
        len += phi_K->phi[i].len;
    }
    #endif


    long delta = 14;
    long len_phi = Z_lval(lideal_norm(I), 2);
    long len_step = 2*two_tors_height + delta;

    if (len_phi == 0) {
        avma = ltop; *L = J; *phi_L = *phi_J;
        free_two_walk_long(phi_res);
        init_trivial_two_walk_long(phi_res);
        return;
    }

    long steps = len_phi / len_step;
    if (steps*len_step < len_phi) ++steps;


    proj phi_K_basis[6], phi_K_target = global_setup.E0;
    phi_K_basis[0] = torsion_basis_sum[0];
    phi_K_basis[1] = torsion_basis_sum[1];
    phi_K_basis[2] = torsion_basis_sum[2];
    phi_K_basis[3] = torsion_basis_twist_sum[0];
    phi_K_basis[4] = torsion_basis_twist_sum[1];
    phi_K_basis[5] = torsion_basis_twist_sum[2];

    special_isogeny phi_J_i;
    two_walk_long phi_K_i, phi[steps];
    init_trivial_two_walk_long(&phi_K_i);

    long len_first_step = (endpoint_close_to_E0) ? (len_phi % len_step) : len_step;

    GEN I_i_long = I;
    GEN I_i_short = lideal_create(A, order, lideal_generator(I), gmul(lideal_norm(J),powuu(2, len_first_step)));
    GEN J_i = J;
    GEN K_i = K;
    phi_J_i = *phi_J;
    copy_two_walk_long(&phi_K_i, phi_K);


    GEN alpha, beta;


    int len_phi1 = (len_first_step < two_tors_height) ? len_first_step : two_tors_height;
    GEN I1 = lideal_create(A, order, lideal_generator(I_i_short), powuu(2, len_phi1));
    two_walk phi_1_0 = ideal_to_isogeny_O0_two(I1);
    two_walk phi_1 = push_two_walk_through_special_isogeny(&phi_1_0, &phi_J_i);

    two_walk_composition_sl(&phi_K_i, &phi_1, &phi_K_i);


    eval_walk_long_mult(&phi_K_i, &phi_K_target, phi_K_basis, 6);



    for (int i = 0; i < steps; ++i) {
        alpha = lideal_isom(J_i, K_i); // J_i*alpha = K_i
        if (gcmp(lideal_norm(K_i), gen_1) == 0) {
            alpha = alg_scalar(A, gen_1); /* make sure we don't apply a distorsion */

        }

        init_trivial_two_walk_long(&phi[i]);


        ideal_to_isogeny_two_2f_delta(&phi[i], &J_i, &phi_J_i,
                                      I_i_short, J_i, K_i,
                                      phi_K_basis, &phi_K_target, delta,
                                      I_i_long);

        assert(lideal_isom(J_i, I_i_short));

        two_walk_composition_ll(&phi_K_i, &phi[i], &phi_K_i);
        // update K, phi_K and I_i_*

        K_i = lideal_mul(I_i_short, alpha); // I_i_short*alpha
        beta = lideal_isom(I_i_short, J_i); // I_i*beta = J_i
        I_i_long = lideal_mul(I_i_long, beta); // I_i_long*beta
        I_i_short = lideal_create(A, order, lideal_generator(I_i_long), gmul(lideal_norm(J_i),powuu(2, len_step)));

    }

    *L = gerepilecopy(ltop, J_i);
    assert(lideal_isom(I, *L));
    *phi_L = phi_J_i;



    // reconstruct phi
    two_walk_long phi_full;
    init_trivial_two_walk_long(&phi_full);
    two_walk_stol(&phi_full, &phi_1);
    for (int i = 0; i < steps; ++i) {
        two_walk_composition_ll(&phi_full, &phi[i], &phi_full);
        free_two_walk_long(&phi[i]);
    }


    copy_two_walk_long(phi_res, &phi_full);
    free_two_walk_long(&phi_full);
    free_two_walk_long(&phi_K_i);

}






void ideal_to_isogeny_O0_two_long(two_walk_long *phi, GEN *L, special_isogeny *phi_L, GEN I,bool endpoint_close_to_E0) {
    pari_sp ltop = avma;

    GEN trivial_ideal = lideal_create(global_setup.B, global_setup.O0, mkcol4s(1,0,0,0), gen_1);
    special_isogeny triv_special = trivial_special_isogeny();
    two_walk_long triv_two;
    init_trivial_two_walk_long(&triv_two);
    //output(lideal_norm(I));

    ideal_to_isogeny_two(phi, L, phi_L, I, trivial_ideal, trivial_ideal, &triv_special, &triv_two,endpoint_close_to_E0);
    *L = gerepilecopy(ltop, *L);
}
