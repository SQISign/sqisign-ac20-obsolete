#define _XOPEN_SOURCE

#include <stdint.h>
//#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pari/pari.h>
#include <assert.h>


#include "ideal.h"
#include "idiso.h"
#include "constants.h"
#include "precomputed.h"
#include "tedwards.h"
#include "isogenies.h"
#include "klpt.h"
#include "toolbox.h"
#include "sqisign.h"
// #include "isomorphism.h"
// #include "uintbig.h"


void init_compressed_sig(compressed_signature *comp_sigma) {
  comp_sigma->zip=malloc(sizeof(uint64_t)*signing_length_two_tors_height_step);
}

void free_compressed_sig(compressed_signature *comp_sigma) {
  free(comp_sigma->zip);
}

void keygen(public_key *pk, secret_key *sk) {

    unsigned int length_NI = security_level/2;
    unsigned int e_tau = (double)security_level*2;
    GEN NI = NULL;
    GEN NJ = powiu(gen_2, e_tau);

    do {
        NI = randomprime(powiu(gen_2, length_NI));
    } while (Fp_issquare(gen_2,NI));

    GEN gamma = NULL;
    while (!gamma) {
        // parity option is 1 so the 2-walk is not backtracking
        gamma = norm_equation_special(global_setup.p, gmul(NI,NJ), 1, true);
    }
    gamma = gtrans(gamma);

    sk->I_large_prime = lideal_create(global_setup.B, global_setup.O0, gamma, NI);
    sk->I_two = lideal_create(global_setup.B, global_setup.O0, alg_conj(global_setup.B, gamma), NJ);

    init_trivial_two_walk_long(&sk->phi_two);
    ideal_to_isogeny_O0_two_long(&sk->phi_two, &sk->I_T, &sk->phi_T, sk->I_two, true);


    if (!sk->phi_T.phi2_set) {
        assert(sk->phi_T.phi2_dual_set);
        sk->phi_T.phi2 = sk->phi_T.phi2_dual;
        sk->phi_T.middle = sk->phi_T.target;
        dual(&sk->phi_T.middle, &sk->phi_T.phi2);
        sk->phi_T.phi2_set = true;
    }

    if (!sk->phi_T.phi2_dual_set) {
        assert(sk->phi_T.phi2_set);
        sk->phi_T.phi2_dual = sk->phi_T.phi2;
        sk->phi_T.target = sk->phi_T.middle;
        dual(&sk->phi_T.target, &sk->phi_T.phi2_dual);
        sk->phi_T.phi2_dual_set = true;
    }

    pk->E = sk->phi_T.target;

}

void commitment(GEN *coeff, GEN *I, odd_isogeny *phi_com){
    pari_sp ltop = avma;

    GEN coeff_plus_1 = cgetg(p_plus_len+1, t_VEC);
    GEN coeff_plus_2 = cgetg(p_plus_len+1, t_VEC);
    GEN coeff_minus_1 = cgetg(p_minus_len+1, t_VEC);
    GEN coeff_minus_2 = cgetg(p_minus_len+1, t_VEC);

    for (int i = 0; i < p_plus_len; ++i) {
        gel(coeff_plus_1,i+1) = powuu(p_plus_fact[i], p_plus_mult[i] - p_plus_mult_com[i]);
        gel(coeff_plus_2,i+1) = randomi(powuu(p_plus_fact[i],p_plus_mult_com[i]));
        gel(coeff_plus_2,i+1) = gmul(gel(coeff_plus_2,i+1), gel(coeff_plus_1,i+1));

        gel(coeff_plus_1,i+1) = gmod(gel(coeff_plus_1,i+1), powuu(p_plus_fact[i],p_plus_mult[i]));
        gel(coeff_plus_2,i+1) = gmod(gel(coeff_plus_2,i+1), powuu(p_plus_fact[i],p_plus_mult[i]));
    }

    for (int i = 0; i < p_minus_len; ++i) {
        gel(coeff_minus_1,i+1) = powuu(p_minus_fact[i], p_minus_mult[i] - p_minus_mult_com[i]);;
        gel(coeff_minus_2,i+1) = randomi(powuu(p_minus_fact[i],p_minus_mult_com[i]));
        gel(coeff_minus_2,i+1) = gmul(gel(coeff_minus_2,i+1), gel(coeff_minus_1,i+1));

        gel(coeff_minus_1,i+1) = gmod(gel(coeff_minus_1,i+1), powuu(p_minus_fact[i],p_minus_mult[i]));
        gel(coeff_minus_2,i+1) = gmod(gel(coeff_minus_2,i+1), powuu(p_minus_fact[i],p_minus_mult[i]));
    }

    GEN coeff_plus = mkvec2(coeff_plus_1, coeff_plus_2);
    GEN coeff_minus = mkvec2(coeff_minus_1, coeff_minus_2);
    *coeff = mkvec2(coeff_plus, coeff_minus);

    *I = kernel_to_ideal_O0_T(*coeff);


    proj ker = coeff_to_E0(gel(*coeff,1), false);
    proj ker_twist = coeff_to_E0(gel(*coeff,2), true);
    isog_degree deg, deg_twist;
    famat_to_degree(&deg, &deg_twist, Z_factor(lideal_norm(*I)));

    phi_com->kernel_plus = ker;
    phi_com->kernel_minus = ker_twist;
    phi_com->deg_plus = deg;
    phi_com->deg_minus = deg_twist;

    gerepileall(ltop, 2, I, coeff);
}



void deterministic_point(proj *P, proj const *A, long ell, long e, bool twist, GEN seed) {
    pari_sp ltop = avma;
    GEN rand_state = getrand();
    setrand(seed);

    unsigned short newseed[3] = {1,2,3};
    unsigned short *oldptr = seed48(newseed);

    uintbig cofactor;
    uintbig_add3(&cofactor, &p, &uintbig_1);
    uintbig ell_big;
    uintbig_set(&ell_big, ell);
    for (int i = 0; i < e; ++i) {
        uintbig_div3_64(&cofactor, &cofactor, ell);
    }
    proj Z;

    while (1) {
        fp2_random(&P->x); fp2_random(&P->z);
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

    setrand(rand_state);
    seed48(oldptr);
    avma = ltop;
}


void deterministic_basis(point *P1_ted, point *P2_ted, proj const *A, long ell, long e, bool twist) {
    proj P1, P2;
    point P1_mul, P2_mul, tmp;
    uintbig ell_big;
    uintbig_set(&ell_big, ell);
    fp2 weil;

    proj E;
    mont_to_ted(&E, A, twist);

    GEN seed = gen_1;

    deterministic_point(&P1, A, ell, e, twist, seed);

    seed = gadd(seed, gen_1);

    mont_to_ted_point(P1_ted, A, &P1);
    P1_mul = *P1_ted;
    for (int i = 0; i < e-1; ++i) {
        ted_mul(&P1_mul, &P1_mul, &E, &ell_big);
    }

    assert(ted_is_on_curve(&P1_mul,&E));
    assert(!ted_iszero(&P1_mul));
    ted_mul(&tmp, &P1_mul, &E, &ell_big);
    assert(ted_iszero(&tmp));

    do {
        deterministic_point(&P2, A, ell, e, twist, seed);
        mont_to_ted_point(P2_ted, A, &P2);
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
        seed = gadd(seed, gen_1);
    } while (fp2_iszero(&weil));
}




bool mont_dlp(GEN *a, GEN *b, const proj *A, const proj *P, const proj *P1,
    const proj *P2, const proj *P1plusP2, long ell, long e, bool twist) {
    proj Q;
    proj E;
    point basis_ted[3], P_ted;
    mont_to_ted(&E, A, twist);

    mont_to_ted_point(&P_ted, A, P);
    mont_to_ted_point(&basis_ted[0], A, P1);
    mont_to_ted_point(&basis_ted[1], A, P2);

    ted_add(&basis_ted[2], &E, &basis_ted[0], &basis_ted[1]);
    ted_to_mont_point(&Q, &basis_ted[2]);
    if (!mont_equal(&Q, P1plusP2)) { ted_neg(&basis_ted[1], &basis_ted[1]); }
    ted_add(&basis_ted[2], &E, &basis_ted[0], &basis_ted[1]);
    ted_to_mont_point(&Q, &basis_ted[2]);
    assert(mont_equal(&Q, P1plusP2));

    return ted_bidim_log(a, b, &E, &P_ted, &basis_ted[0], &basis_ted[1], ell, e);
}


//in fact only work for ell=2
bool mont_power_dlp(uint64_t *a,const proj *A, const proj *Q, const proj *P,const proj *PQ,long e) {
  return mont_two_DLP(a,A,Q,P,PQ,e);
}


// PRF to generate points
static void hash(proj *P, int i) {
  uintbig_set(&P->x.re.x, 3 * i + 13);
  uintbig_set(&P->z.re.x, 5 * i * i + 17);
  uintbig_set(&P->x.im.x, 7 * i * i * i + 19);
  uintbig_set(&P->z.im.x, 11 * i * i * i + 23);
}
// Find a basis of the 2-torsion of A, deterministically
//
// Outputs x(P), x(Q) and x(P-Q) of a basis (P,Q) such that [2^(n-1)]P
// = (0,0).
//
// Assumes the curve A has order p+1
static void find_basis(proj *P, proj *Q, proj *PQ, proj *A) {
  bool oncurve = class_mod_4 == 3;
  proj P2, Q2, tmp;
  long cnt = 1;

    //normalize for deterministic computation
    normalize_proj(A);

  // Get first point
  while (true) {
    hash(P, cnt++);
    if (is_on_curve(P, A) != oncurve)
      continue;
    // multiply by cofactor
    xMUL(P, A, P, &p_even_cofactor);
    // check it has maximal order
    P2 = *P;
    for (int i = 1; i < two_tors_height; i++)
      xDBL(&P2, A, &P2);
    if (!mont_iszero(&P2))
      break;
  }

  // Get linearly independent point
  while (true) {
    hash(Q, cnt++);
    if (is_on_curve(Q, A) != oncurve)
      continue;
    // multiply by cofactor
    xMUL(Q, A, Q, &p_even_cofactor);
    // check it has maximal order
    Q2 = *Q;
    for (int i = 1; i < two_tors_height; i++)
      xDBL(&Q2, A, &Q2);
    if (!mont_iszero(&Q2) && !mont_equal(&Q2, &P2))
      break;
  }

  // Compute P-Q
  xBILIFT(PQ, &tmp, P, Q, A);

  // Shuffle to satisfy constraint
  if (fp2_iszero(&P2.x)) {
  } else if (fp2_iszero(&Q2.x)) {
    fp2_cswap(&P->x, &Q->x, true);
    fp2_cswap(&P->z, &Q->z, true);
  } else {
    fp2_cswap(&P->x, &PQ->x, true);
    fp2_cswap(&P->z, &PQ->z, true);
  }
}


//renormalize the two walk into a walk of step having length two_tors_height
long normalized_walk(two_walk *w,uint64_t *zip, long *n) {
  long tot = 0;
  for (int i=0;i<*n;i++) {
    tot+=w[i].len;
  }
  long step= tot / two_tors_height;
  long add_step=1;
  if (tot == step*two_tors_height) {
    add_step=0;
  }
  long index_w = 0;
  two_walk norm_walk[step+add_step];
  // proj Pc[step+add_step],Qc[step+add_step],PQc[step+add_step],Ac[step+add_step];
  while ( w[index_w].len == two_tors_height ) {
    norm_walk[index_w]=w[index_w];
    //this computes a basis <P,Q> and the value log such that the kernel is generated by
    // P + log*Q
    // currently the algorithm is pushing the basis through and computing a DLP
    // for this one we already have the kernel but doing a bidimensionnal DLP appears difficult
    // as we would need some extra point (the difference between the kernel and the points of the basis ?)
    proj P1,P2,P3,dummy_A;
    find_basis(&P1,&P2,&P3,&w[index_w].A);
    // Pc[0]=P1;
    // Qc[0]=P2;
    // PQc[0]=P3;
    // Ac[0]=w[index_w].A;
    proj push_points[3];
    push_points[0]=P2;
    push_points[1]=P1;
    push_points[2]=P3;
    eval_walk_mult(&w[index_w],&dummy_A,push_points,3);
    // what we do when P is the kernel generator ?
    bool dlp=mont_power_dlp(&(zip[index_w]),&dummy_A,&(push_points[0]),&(push_points[1]),&(push_points[2]),two_tors_height);
    zip[index_w]=pow(2,two_tors_height) - zip[index_w];
    assert(dlp);

    index_w++;
    isomorphism isom;
    mont_isom(&isom,&w[index_w].A,&dummy_A);
    mont_isom_apply(&isom,&w[index_w].ker);
    // w[index_w].A= dummy_A;
  }

  proj A = w[index_w].A;
  proj P = w[index_w].ker;
  long order_P = w[index_w].len;
  assert(mont_equal(&A,&w[index_w].A));


  for (int index=index_w; index < step; index++ ){
    norm_walk[index].A=A;
    norm_walk[index].len=two_tors_height;
    proj P1,P2,P3;
    proj dummy_A;
    find_basis(&P1,&P2,&P3,&A);
    // Pc[index]=P1;
    // Qc[index]=P2;
    // PQc[index]=P3;
    // Ac[index]=A;
    #ifndef NDEBUG
    proj test_order;
    test_order=P1;
    for (int i=1;i<two_tors_height; i++) {
      xDBL(&test_order,&A,&test_order);
    }
    assert(!mont_iszero(&test_order));
    xDBL(&test_order,&A,&test_order);
    assert(mont_iszero(&test_order));
    uint64_t a;
    assert(!mont_power_dlp(&a,&A,&P2,&P1,&P3,two_tors_height));
    uintbig log_test;
    long mult_fac=55;
    uint64_t x_test;
    uintbig_set(&log_test,mult_fac);
    proj test_point,test_pt2;
    xMUL(&test_point,&A,&P1,&log_test);
    uintbig_set(&log_test,mult_fac+1);
    xMUL(&test_pt2,&A,&P1,&log_test);

    assert(mont_power_dlp(&x_test,&A,&test_point,&P1,&test_pt2,two_tors_height));
    #endif

    long w_step=1;
    long remain_step = order_P;
    while (remain_step+w[index_w + w_step].len <two_tors_height) {
      remain_step+=w[index_w+w_step].len;
      w_step++;
    }
    two_walk loc_phi[w_step+1];
    loc_phi[0].A = A;
    loc_phi[0].ker=P;
    loc_phi[0].len =order_P;


    if (w_step > 1) {
      for (int i=1; i < w_step ;i++){
        loc_phi[i] = w[index_w+i];
      }
      // remain_step+=w[index_w+1].len;
    }
    remain_step = two_tors_height - remain_step;

    if (remain_step == w[index_w+w_step].len) {
        loc_phi[w_step] = w[index_w+w_step];
        index_w += w_step+1;
        order_P = w[index_w].len;
        P = w[index_w].ker;
    }
    else {
      loc_phi[w_step].A = w[index_w+w_step].A;
      loc_phi[w_step].len=remain_step;
      xDBL(&loc_phi[w_step].ker,&w[index_w+w_step].A,&w[index_w+w_step].ker);
      for (int i=0; i < w[index_w+w_step].len - remain_step -1 ; i++ ) {
        xDBL(&loc_phi[w_step].ker,&w[index_w+w_step].A,&loc_phi[w_step].ker);
      }
      index_w += w_step;
      order_P = w[index_w].len - remain_step;
      P=w[index_w].ker;
      isomorphism isom;
      eval_walk_isom_mult(&isom,&(loc_phi[w_step]), &dummy_A,&(loc_phi[w_step]),&P,1);
    }
    proj push_points[3];
    push_points[0]=P2;push_points[1]=P1;push_points[2]=P3;

    for (int i=0; i < w_step +1 ;i ++) {
      long log;
      isomorphism isom;
      eval_walk_isom_mult(&isom,&(loc_phi[i]),&A,&(loc_phi[i]),push_points,3);
      if (i != w_step) {
        isomorphism isom;
        mont_isom(&isom,&loc_phi[i+1].A,&A);
        mont_isom_apply(&isom,&loc_phi[i+1].ker);
        loc_phi[i+1].A=A;
      }
      if (i == w_step) {
        isomorphism isom;
        mont_isom(&isom,&dummy_A,&A);
        mont_isom_apply(&isom,&P);
      }
    }
    uintbig x,y;

    if (mont_iszero(&(push_points[0])) ) {
      uintbig_set(&x,1);
      uintbig_set(&y,0);
    }
    if (mont_iszero(&(push_points[1])) ){
      uintbig_set(&x,0);
      uintbig_set(&y,1);
    }
    else {
      uintbig_set(&x,1);
      bool dlp =mont_power_dlp(&(zip[index]),&A,&(push_points[0]),&(push_points[1]),&(push_points[2]),two_tors_height );
      zip[index]=pow(2,two_tors_height) - zip[index];
      uintbig_set(&y,zip[index]);
      assert(dlp);
      #ifndef NDEBUG
      proj verif_dlp;
      verif_dlp=push_points[1];
      xMUL(&verif_dlp,&A,&verif_dlp,&y);
      assert(mont_equal(&verif_dlp,&(push_points[0])));
      #endif
    }

    xBIDIM(&norm_walk[index].ker,&(loc_phi[0].A),&P2,&x,&P1,&y,&P3);
    isomorphism isom;
    proj A_target;
    eval_walk_isom(&isom,&norm_walk[index],&A_target,NULL,&norm_walk[index],NULL);
    #ifndef NDEBUG
    proj verif=P;
    for (int i=1; i < order_P; i++) {
      xDBL(&verif,&A,&verif);
    }
    assert(!mont_iszero(&verif));
    xDBL(&verif,&A,&verif);
    assert(mont_iszero(&verif));
    proj proj_tmp = {fp2_1, fp2_0};
    proj j1,j2;
    jinv256(&j1,&A);
    // eval_walk(&norm_walk[index],&A_target,&proj_tmp);;
    jinv256(&j2,&A_target);
    assert(mont_equal(&j1, &j2));
    #endif

  }
  if (add_step == 1){
    norm_walk[step].A=A;
    norm_walk[step].ker=P;
    norm_walk[step].len=order_P;

    //compute the compression coeff for the last one
    isomorphism isom;

    proj P1,P2,P3,A_target;
    // eval_walk_isom(&isom,&norm_walk[step],&A_target,NULL,&norm_walk[step],NULL);
    find_basis(&P1,&P2,&P3,&A);
    two_walk phi_test;
    phi_test.A=A;
    phi_test.len=order_P;
    proj push_points[3];
    push_points[0]=P2;
    push_points[1]=P1;
    push_points[2]=P3;
    for (int i=0;i < two_tors_height - order_P; i++) {
      for (int j=0; j<3; j++){
        xDBL(&(push_points[j]),&A,&(push_points[j]));
      }
    }
    P1=push_points[1];
    P2=push_points[0];
    P3=push_points[2];
    eval_walk_isom_mult(&isom,&norm_walk[step],&A_target,&norm_walk[step],push_points,3);
    //what we do when P1 is the kernel generator ?
    bool dlp=mont_power_dlp(&zip[step],&A_target,&(push_points[0]),&(push_points[1]),&(push_points[2]),order_P);
    zip[step]= pow(2,order_P) - zip[step];
    assert(dlp);
    #ifndef NDEBUG
    proj A_test;
    uintbig a;
    uintbig_set(&a,zip[step]);
    xBIDIM(&phi_test.ker,&phi_test.A,&P1,&a,&P2,&uintbig_1,&P3);
    isomorphism isom2;
    eval_walk_isom(&isom2,&phi_test,&A_test,NULL,&phi_test,NULL);
    proj j1,j2;
    jinv256(&j1,&A_target);
    jinv256(&j2,&A_test);
    assert(mont_equal(&j1,&j2));
    #endif

  }
  for ( int i=0; i<step+add_step; i++){
    w[i]=norm_walk[i];
  }
  *n=step+add_step;
  #ifndef NDEBUG
  uintbig a;
  proj Q, PQ;
  A=w[0].A;
  two_walk walk[step+add_step];
  for (int i = 0; i < step; i++) {
    uintbig_set(&a, zip[i] );
    // long hint = (zip[i] & hint_mask) >> two_tors_height;
    // get the next kernel
    find_basis(&P, &Q, &PQ, &A);  // TODO: use point Q from previous step + hint

    // assert(mont_equal(&A,&Ac[i]));
    // assert(mont_equal(&P,&Pc[i]));
    // assert(mont_equal(&Q,&Qc[i]));
    // assert(mont_equal(&PQ,&PQc[i]));
    xBIDIM(&(walk[i].ker), &A, &P, &a, &Q, &uintbig_1, &PQ);
    // assert(mont_equal(&walk[i].ker,&w[i].ker));
    walk[i].A = A;
    walk[i].len = two_tors_height;
    // take the next step
    isomorphism isom;
    eval_walk_isom(&isom,&walk[i], &A, NULL,&walk[i],NULL);
    // A=walk[i].A;
    proj j1,j2;
    jinv256(&j1,&A);
    jinv256(&j2,&w[i+1].A);
    assert(mont_equal(&j1,&j2));
  }
  #endif

  return order_P;

}



void challenge(proj *E_cha, const uintbig *m, const proj *E_com, const proj *basis_plus, const proj *basis_minus, GEN *dlog, proj *basis_two){
    pari_sp ltop = avma;
    //unsigned short newseed[3] = {1,2,3};
    //unsigned short *oldptr = seed48(newseed);

    odd_isogeny phi;
    proj A = *E_com;
    long index;
    bool twist;
    uintbig k;
    long ell;
    proj H, P, Z;

    uintbig_set(&H.x.re.x, m->c[0]);
    uintbig_set(&H.x.im.x, m->c[1]);
    uintbig_set(&H.z.re.x, m->c[2]);
    uintbig_set(&H.z.im.x, m->c[3]);


    uintbig cofactor_plus = {1,0,0,0}, cofactor_minus = {1,0,0,0};
    uintbig order_plus = {1,0,0,0}, order_minus = {1,0,0,0};
    isog_degree deg;
    degree_one(&deg);


    // compute cofactor and degree of the 'minus' part
    for (int i = 0; i < p_plus_len; ++i) {
        ell = p_plus_fact[i];
        for (int j = 0; j < p_plus_mult[i] - p_plus_mult_cha[i]; ++j){
            uintbig_mul3_64(&cofactor_plus, &cofactor_plus, ell);
        }
        for (int j = 0; j < p_plus_mult_cha[i]; ++j){
            uintbig_mul3_64(&order_plus, &order_plus, ell);
            index = ell_to_index(p_plus_fact[i], &twist);
            degree_set(&deg, index, p_plus_mult_cha[i]);
        }

    }


    bool bad;


    // find the 'plus' part of the kernel
    while (1) {
        //fp2_random(&P->x); fp2_random(&P->z);
        fp_add2(&H.x.re, &fp_1);
        if (!is_on_curve(&H, E_com)) continue;
        xMUL(&P, E_com, &H, &cofactor_plus);
        xMUL(&P, E_com, &P, &p_plus_odd_cofactor);

        bad = false;
        for (int i = 0; i < p_plus_len; ++i) {
            if (p_plus_mult_cha[i] > 0) {
                ell = p_plus_fact[i];
                Z = P;
                uintbig_div3_64(&k, &order_plus, ell);

                xMUL(&Z, E_com, &Z, &k);
                if (mont_iszero(&Z)) { bad = true; break; }

                #ifndef NDEBUG
		uintbig ell_big;
                uintbig_set(&ell_big, ell);
                xMUL(&Z, E_com, &Z, &ell_big);
                assert(mont_iszero(&Z));
                #endif
            }
        }
        if (bad) continue;
        else break;

    }

    GEN dlog_plus = NULL;
    if (basis_plus) {
        long len = p_plus_len;
        GEN coeff_1 = cgetg(len+1, t_VEC);
        GEN coeff_2 = cgetg(len+1, t_VEC);
        proj Q;
        for (int i = 0; i < len; ++i) {

            gel(coeff_1,i+1) = gen_0;
            gel(coeff_2,i+1) = gen_0;

            if (0 < p_plus_mult_cha[i]) {
                GEN a,b;
                bool dlp = mont_dlp(&a, &b, &A, &P, &basis_plus[0], &basis_plus[1], &basis_plus[2],
                    p_plus_fact[i], p_plus_mult_cha[i], false);
                assert(dlp);
                gel(coeff_1,i+1) = a;
                gel(coeff_2,i+1) = b;
            }
        }

        //dlog_plus = torsion_crt_compose(mkvec2(coeff_1, coeff_2), false);
        dlog_plus = mkvec2(coeff_1, coeff_2);

        GEN dlog_plus_composed = torsion_crt_compose(dlog_plus, false);

        uintbig a_big, b_big;
        gentobig(&a_big, gel(dlog_plus_composed, 1));
        gentobig(&b_big, gel(dlog_plus_composed, 2));
        xBIDIM(&Q, &A, &basis_plus[0], &a_big, &basis_plus[1], &b_big, &basis_plus[2]);
        assert(mont_equal(&Q,&P));

    }




    phi.kernel_plus = P;
    phi.kernel_minus.x = fp2_1;
    phi.kernel_minus.z = fp2_0;
    phi.deg_plus = deg;
    phi.deg_minus.val = 0;


    long len_points = (basis_plus) ? 3 : 1;
    if (basis_two) len_points += 3;

    proj points[len_points];
    points[0] = P;

    if (basis_minus) {
        points[0] = basis_minus[0];
        points[1] = basis_minus[1];
        points[2] = basis_minus[2];
    }
    if (basis_two) {
        points[3] = basis_two[0];
        points[4] = basis_two[1];
        points[5] = basis_two[2];
    }


    eval_mult(&A, &phi, points, len_points);


    // compute cofactor and degree of the 'minus' part
    for (int i = 0; i < p_minus_len; ++i) {
        ell = p_minus_fact[i];
        for (int j = 0; j < p_minus_mult[i] - p_minus_mult_cha[i]; ++j){
            uintbig_mul3_64(&cofactor_minus, &cofactor_minus, ell);
        }
        for (int j = 0; j < p_minus_mult_cha[i]; ++j){
            uintbig_mul3_64(&order_minus, &order_minus, ell);
            index = ell_to_index(p_minus_fact[i], &twist);
            degree_set(&deg, index, p_minus_mult_cha[i]);
        }

    }

    // find the 'minus' part of the kernel
    while (1) {
        fp_add2(&H.x.re, &fp_1);
        //fp2_random(&P->x); fp2_random(&P->z);
        if (is_on_curve(&H, &A)) continue;
        xMUL(&P, &A, &H, &cofactor_minus);
        xMUL(&P, &A, &P, &p_minus_odd_cofactor);


        bad = false;
        for (int i = 0; i < p_minus_len; ++i) {
            if (p_minus_mult_cha[i] > 0) {
                ell = p_minus_fact[i];
                Z = P;
                uintbig_div3_64(&k, &order_minus, ell);

                xMUL(&Z, &A, &Z, &k);
                if (mont_iszero(&Z)) { bad = true; break; }

                #ifndef NDEBUG
		uintbig ell_big;
                uintbig_set(&ell_big, ell);
                xMUL(&Z, &A, &Z, &ell_big);
                assert(mont_iszero(&Z));
                #endif
            }
        }
        if (bad) continue;
        else break;
    }


    GEN dlog_minus = NULL;
    if (basis_minus) {
        long len = p_minus_len;
        GEN coeff_1 = cgetg(len+1, t_VEC);
        GEN coeff_2 = cgetg(len+1, t_VEC);
        proj Q;
        for (int i = 0; i < len; ++i) {

            gel(coeff_1,i+1) = gen_0;
            gel(coeff_2,i+1) = gen_0;

            if (0 < p_minus_mult_cha[i]) {
                GEN a,b;


                bool dlp = mont_dlp(&a, &b, &A, &P, &points[0], &points[1], &points[2],
                    p_minus_fact[i], p_minus_mult_cha[i], true);
                assert(dlp);
                gel(coeff_1,i+1) = a;
                gel(coeff_2,i+1) = b;
            }
        }


        //dlog_minus = torsion_crt_compose(mkvec2(coeff_1, coeff_2), true);

        dlog_minus = mkvec2(coeff_1, coeff_2);

        GEN dlog_minus_composed = torsion_crt_compose(dlog_minus, true);

        uintbig a_big, b_big;
        gentobig(&a_big, gel(dlog_minus_composed, 1));
        gentobig(&b_big, gel(dlog_minus_composed, 2));
        xBIDIM(&Q, &A, &points[0], &a_big, &points[1], &b_big, &points[2]);
        assert(mont_equal(&Q,&P));

    }


    len_points = (basis_two) ? 3 : 1;

    if (basis_two) {
        points[0] = points[3];
        points[1] = points[4];
        points[2] = points[5];
    }

    phi.kernel_minus = P;
    phi.kernel_plus.x = fp2_1;
    phi.kernel_plus.z = fp2_0;
    phi.deg_minus = deg;
    phi.deg_plus.val = 0;
    eval_mult(&A, &phi, points, len_points);

    if (basis_two) {
        basis_two[0] = points[0];
        basis_two[1] = points[1];
        basis_two[2] = points[2];
    }

    *E_cha = A;

    if (dlog) { *dlog = gerepilecopy(ltop, mkvec2(dlog_plus, dlog_minus)); }

    //seed48(oldptr);
    //avma = ltop;
}

//function to test the decompression
void decompress(two_walk *walk, proj *A, const uint64_t *zip, long len,long last_step) {
  long mask = (1 << two_tors_height) - 1;
  long hint_mask = (0xf << two_tors_height);
  uintbig a;
  proj P, Q, PQ;
  for (int i = 0; i < len-1; i++) {
    // printf("zip %ld ",zip[i]);

    uintbig_set(&a, zip[i]);
    long hint = (zip[i] & hint_mask) >> two_tors_height;
    // get the next kernel
    find_basis(&P, &Q, &PQ, A);  // TODO: use point Q from previous step + hint
    xBIDIM(&(walk[i].ker), A, &P, &a, &Q, &uintbig_1, &PQ);
    // printf(" k %ld ",walk[i].ker.x.re.x.c[0]);
    walk[i].A = *A;
    walk[i].len = two_tors_height;
    // take the next step
    isomorphism isom;
    eval_walk_isom(&isom,&walk[i], A, NULL,&walk[i],NULL);
    // printf("\n");
  }
  //last step of smaller size
  uintbig_set(&a, zip[len-1]);
  long hint = (zip[len-1] & hint_mask) >> two_tors_height;
  // get the next kernel
  find_basis(&P, &Q, &PQ, A);
  // printf(" A %ld \n",A->x.re.x.c[0]);
  for (int i=0; i < two_tors_height - last_step ; i++){
    xDBL(&P,A,&P);
    xDBL(&Q,A,&Q);
    xDBL(&PQ,A,&PQ);
  }
  xBIDIM(&walk[len-1].ker, A, &P, &a, &Q, &uintbig_1, &PQ);
  walk[len-1].A = *A;
  walk[len-1].len = last_step;
  isomorphism isom;
  eval_walk_isom(&isom,&walk[len-1],A,NULL,&walk[len-1], NULL);
  // printf(" A %ld \n",A->x.re.x.c[0]);
}

// basis_two is the image of the basis of the two torsion through phi_com and phi_cha
void response(two_walk_long *sigma,uint64_t *zip,  GEN coeff_ker_challenge_commitment, const secret_key *sk, const proj *basis_two, const proj *E_cha){
    pari_sp ltop = avma;
    GEN A = global_setup.B;
    GEN order = global_setup.O0;

    GEN I = sk->I_large_prime;
    GEN I_two = sk->I_two;

    // STEP 1: compute the ideal of phi_challenge_commitment
    GEN I_phipsi = kernel_to_ideal_O0_T(coeff_ker_challenge_commitment);
    // STEP 2: find J of norm a power of two such that J inter I is equivalent to I_phipsi
    GEN beta = lideal_isom(I_two, I); // I_two*beta = I
    GEN alpha = gmul(beta, lideal_norm(I_two)); // n(alpha) = n(I)
    GEN gamma = lideal_generator_coprime(I_phipsi, gmul(gen_2, lideal_norm(I)));
    GEN generator = algmul(A, gamma, alpha);
    GEN norm = gmul(lideal_norm(I_two), lideal_norm(I_phipsi));
    GEN K = lideal_create(A, order, generator, norm);

    assert(lideal_isom(I_phipsi, lideal_inter(K,I))); // K inter I is equivalent to I_phipsi


    GEN J;
    GEN n;
    do{
        J = klpt_general_power(I, K, gen_2); // J inter I is equivalent to I_phipsi
        // printf("path length %ld\n", Z_lval(lideal_norm(J),2));
        alg_primitive(&n, A, order, algmul(A, lideal_generator(J), lideal_generator(I_two)));

        // backtracking?
    } while(gcmp(n,gen_1) != 0);
    // printf("path length %ld\n", Z_lval(lideal_norm(J),2));

    // STEP 3: compute L of norm n(J) such that L inter sk->I_T is equivalent to I_phipsi

    GEN L = lideal_inter(J,I);
    assert(lideal_isom(I_phipsi, L));

    beta = lideal_isom(I, sk->I_T); // I*beta = I_T;

    L = lideal_mul(L,beta);

    assert(lideal_isom(I_phipsi, L));
    assert(gcmp(lideal_norm(L), gmul(lideal_norm(sk->I_T), lideal_norm(J))) == 0);


    GEN dummy_ideal;
    special_isogeny dummy_isogeny;

    // STEP 4: convert to an isogeny


    #ifndef NDEBUG
    alg_primitive(&n, A, order, algmul(A, lideal_generator(L), lideal_generator(sk->I_two)));
    assert(gcmp(n,gen_1) == 0);
    #endif

    long delta = 14;
    long len_tail = two_tors_height + delta;
    GEN L_cut = lideal_create(A, order, lideal_generator(L), shifti(lideal_norm(L), -len_tail));

    ideal_to_isogeny_two(sigma, &dummy_ideal, &dummy_isogeny, L_cut, sk->I_T, sk->I_two, &sk->phi_T, &sk->phi_two, false);

    GEN gen_tail = lideal_isom(L, I_phipsi); // L*gen_tail = I_phipsi;
    gen_tail = gmul(gen_tail, lideal_norm(L));
    GEN L_tail = lideal_create(A, order, gen_tail, powuu(2,two_tors_height));

    two_walk phi_tail_dual;
    phi_tail_dual.A = *E_cha;
    phi_tail_dual.len = two_tors_height;

    GEN v_tail = ideal_to_kernel_O0_ell(L_tail, 2);
    uintbig x, y;
    gentobig(&x, gel(v_tail,1));
    gentobig(&y, gel(v_tail,2));

    xBIDIM(&phi_tail_dual.ker, &phi_tail_dual.A, &basis_two[0], &x, &basis_two[1], &y, &basis_two[2]);


    isomorphism isom;
    proj L_cut_target, phi_tail_dual_target, proj_tmp = {fp2_1, fp2_0};
    eval_walk(&sigma->phi[sigma->len-1], &L_cut_target, &proj_tmp);

    eval_walk_isom(&isom, &phi_tail_dual, &phi_tail_dual_target, NULL, &phi_tail_dual, NULL);

    two_walk eta;

    bool done;

    proj from = L_cut_target;
    proj to = phi_tail_dual_target; // phi_2 source

    done = MITM(&eta, &from, &to, delta);
    assert(done);

    two_walk phi_tail = phi_tail_dual;
    dual_walk(&phi_tail);

    two_walk_composition_sl(sigma, &eta, sigma);
    two_walk_composition_sl(sigma, &phi_tail, sigma);

    normalized_walk(sigma->phi,zip,&(sigma->len));
    // sigma->len = normalized_n;
    // assert(mont_equal(&sigma->phi[0].A)    //

    avma = ltop;
}


bool simple_check_signature(const two_walk_long *sigma, const uint64_t *zip, const public_key *pk, proj *E_cha) {
    proj j1,j2,j3;
    jinv256(&j1, &pk->E);
    jinv256(&j2, &sigma->phi[0].A);
    assert(mont_equal(&j1, &j2));

    proj sigma_target = sigma->phi[0].A, proj_tmp = {fp2_1, fp2_0};
    for (int i = 0; i < sigma->len; ++i) {
        jinv256(&j1, &sigma_target);
        jinv256(&j2, &sigma->phi[i].A);
        if (!mont_equal(&j1, &j2)) return false;
        sigma_target = sigma->phi[i].A;
        eval_walk(&sigma->phi[i], &sigma_target, &proj_tmp);
    }

    jinv256(&j1, &sigma_target);
    jinv256(&j2, E_cha);
    normalize_proj(&j2);
    proj A_target=(sigma->phi)[0].A;
    two_walk check[sigma->len];
    long last_step= sigma->phi[sigma->len-1].len;
    decompress(check,&A_target,zip,(sigma->len),last_step);
    jinv256(&j3, &A_target);
    // printf("%ld %ld \n",j2.x.re.x.c[0],j2.x.re.x.c[1]);
    return (mont_equal(&j1, &j2) && mont_equal(&j1,&j3));
}

// void zip_copy(compressed_signature *comp_sigma,uint64_t *zip, long len) {
//   compressed_signature res;
//   res.E_com=comp_sigma->E_com;
//   res.zip=malloc(sizeof(uint64_t)*len);
//   for (int i=0;i <len;i++){
//     res.zip[i]=zip[i];
//   }
//   free(comp_sigma->zip);
//   // comp_sigma->E_com = { fp2_0, fp2_1 };
//   *comp_sigma= res;


void sign(compressed_signature *comp_sigma, const secret_key *sk, const public_key *pk, const uintbig *m) {
    pari_sp ltop = avma;

    GEN coeff_com, I_com;
    odd_isogeny phi_com;
    proj E_cha;

    //printf("Commitment\n");

    commitment(&coeff_com, &I_com, &phi_com);

    comp_sigma->E_com = global_setup.E0;


    // compute the image of a basis of the torsion used for the challenge

    proj points[9];
    points[0] = torsion_basis_sum[0];
    points[1] = torsion_basis_sum[1];
    points[2] = torsion_basis_sum[2];
    points[3] = torsion_basis_twist_sum[0];
    points[4] = torsion_basis_twist_sum[1];
    points[5] = torsion_basis_twist_sum[2];
    points[6] = torsion_basis_two[0];
    points[7] = torsion_basis_two[1];
    points[8] = torsion_basis_two[2];

    eval_mult(&comp_sigma->E_com, &phi_com, points, 9);

    proj basis_plus[3], basis_minus[3], basis_two[3];
    basis_plus[0] = points[0];
    basis_plus[1] = points[1];
    basis_plus[2] = points[2];
    basis_minus[0] = points[3];
    basis_minus[1] = points[4];
    basis_minus[2] = points[5];
    basis_two[0] = points[6];
    basis_two[1] = points[7];
    basis_two[2] = points[8];

    uintbig ell_big;
    for (int i = 0; i < p_plus_len; ++i) {
        uintbig_set(&ell_big, p_plus_fact[i]);
        for (int j = 0; j < p_plus_mult[i] - p_plus_mult_cha[i]; ++j){
            for (int l = 0; l < 3; ++l) {
                xMUL(&basis_plus[l], &comp_sigma->E_com, &basis_plus[l], &ell_big);
            }
        }
    }

    for (int i = 0; i < p_minus_len; ++i) {
        uintbig_set(&ell_big, p_minus_fact[i]);
        for (int j = 0; j < p_minus_mult[i] - p_minus_mult_cha[i]; ++j){
            for (int l = 0; l < 3; ++l) {
                xMUL(&basis_minus[l], &comp_sigma->E_com, &basis_minus[l], &ell_big);
            }
        }
    }

    //printf("Challenge\n");

    // comp_sigma->E_com= Sigma->E_com;

    GEN dlog;

    challenge(&E_cha, m, &comp_sigma->E_com, basis_plus, basis_minus, &dlog, basis_two);

    GEN coeff_ker_challenge_commitment = gadd(coeff_com,dlog);


    #ifndef NDEBUG

    odd_isogeny phi_test;
    GEN I = kernel_to_ideal_O0_T(coeff_ker_challenge_commitment);
    proj ker = coeff_to_E0(gel(coeff_ker_challenge_commitment,1), false);
    proj ker_twist = coeff_to_E0(gel(coeff_ker_challenge_commitment,2), true);
    isog_degree deg, deg_twist;
    famat_to_degree(&deg, &deg_twist, Z_factor(lideal_norm(I)));

    phi_test.kernel_plus = ker;
    phi_test.kernel_minus = ker_twist;
    phi_test.deg_plus = deg;
    phi_test.deg_minus = deg_twist;

    proj A = global_setup.E0, H = {fp2_1, fp2_0};
    eval(&A, &phi_test, &H);

    assert(mont_equal(&A, &E_cha));

    #endif


    uint64_t zip[signing_length_two_tors_height_step];
    two_walk_long sigma;
    init_trivial_two_walk_long(&sigma);

    response(&sigma, comp_sigma->zip, coeff_ker_challenge_commitment, sk, basis_two, &E_cha);

    // Sigma->sigma = sigma;
    // zip_copy(comp_sigma,zip,signing_length_two_tors_height_step);


    // for (int i=0;i<signing_length_two_tors_height_step ;i ++ ){
    //     comp_sigma->zip[i]=zip[i];
    // }


    assert(simple_check_signature(&sigma,comp_sigma->zip, pk, &E_cha));

    free_two_walk_long(&sigma);

    avma = ltop;
}





bool verif(compressed_signature *comp_sigma, const public_key *pk,const uintbig *m){
    proj A_chall = comp_sigma->E_com;
    challenge(&A_chall, m, &A_chall, NULL, NULL, NULL, NULL);
    proj A_check=pk->E;
    two_walk walk_check[signing_length_two_tors_height_step];
    // for (int i=0;i< 30;i++){
    //   printf("%ld ",comp_sigma->zip[i]);
    // }
    // printf("\n");

    decompress(walk_check,&A_check,comp_sigma->zip,signing_length_two_tors_height_step,last_step_length);

    proj j1,j2;
    jinv256(&j1,&A_chall);
    normalize_proj(&j1);
    // printf("%ld %ld \n",j1.x.re.x.c[0],j1.x.re.x.c[1]);
    jinv256(&j2,&A_check);
    normalize_proj(&j2);
    // printf("%ld %ld \n",j2.x.re.x.c[0],j2.x.re.x.c[1]);
    return mont_equal(&j1,&j2);
}
