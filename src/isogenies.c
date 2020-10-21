#include "isogenies.h"
#include "isomorphism.h"
#include <assert.h>

// Internal function: evaluates an isogeny proceeding from the largest
// to the smallest degrees.
//
// To be improved: use SIDH-like strategies for the smallest degree
void eval_rtl_old(proj *A, proj *Q, int n, isog_degree deg, const long* fact, long len) {
  uintbig cof;
  degree_to_uint(&cof, deg, fact, len);
  for (int i = len - 1; i >= 0; i--) {
    int v = degree_get(deg, i);
    for (int j = 0; j < v; j++) {
      proj K;
      uintbig_div3_64(&cof, &cof, fact[i]);
      xMUL(&K, A, Q+1, &cof);
      assert(!mont_iszero(&K));
      xISOG_many(A, Q, n, &K, fact[i]);
    }
  }
  assert(mont_iszero(Q+1));
}

//Internal function, recursively building a balanced strategy to evaluate points
void eval_rtl_rec(proj *A, proj *K, int len, bool advance,proj *P, int stacklen, long fact) {

  if (len==0) return;

  if (len==1) {
    xISOG_many(A, P, stacklen, K, fact);

  } else {
    long right = len / 2;
    long left = len - right;

    P[stacklen]=*K;

    uintbig fac_mul;
    uintbig_set(&fac_mul,fact);

    for (int i=0; i<left;i++)
      xMUL(K,A,K,&fac_mul);


    eval_rtl_rec(A, K, right,advance, P, stacklen+1,fact);
    K[right*advance]=P[stacklen];

    eval_rtl_rec(A,K+right*advance,left,advance,P,stacklen,fact);
  }


}

//Internal function
//supersedes eval_rtl_old using the recursive function for bettre isogeny computation
void eval_rtl(proj *A, proj *Q, int n, isog_degree deg, const long* fact, long len) {
  uintbig cof;
  degree_to_uint(&cof, deg, fact, len);

  //starting by the 3^53 or 5^21 isogeny reduces the number of scalar multiplication
  int start=0;
  int v = degree_get(deg, 0);
  if (fact[0]<6 && v>8) {
    start=1;
    proj K;
    for (int j = 0; j < v; j++) {
      uintbig_div3_64(&cof, &cof, fact[0]);
    }
    xMUL(&K, A, Q+1, &cof);


    int log,leni = v;
    for (log =0 ; leni > 1; leni >>=1 ) log++;

    assert(!mont_iszero(&K));

    proj stack[n+log];
    for (int i=0; i<n;i++)
      stack[i]=*(Q+i);
    eval_rtl_rec(A, &K, v, false, stack, n, fact[0]);
    for (int i=0;i<n;i++)
        *(Q+i)=(stack[i]);

  }
  for (int i = len - 1; i >= start; i--) {
    int v = degree_get(deg, i);
    if (fact[i]<6 && v!=0) {
      proj K;
      for (int j = 0; j < v; j++) {
        uintbig_div3_64(&cof, &cof, fact[i]);
      }
      xMUL(&K, A, Q+1, &cof);


      int log,leni = v;
      for (log =0 ; leni > 1; leni >>=1 ) log++;

      assert(!mont_iszero(&K));

      proj stack[n+log];
      for (int i=0; i<n;i++)
        stack[i]=*(Q+i);
      eval_rtl_rec(A, &K, v, false, stack, n, fact[i]);
      assert(mont_iszero(&stack[1]));
      for (int i=0;i<n;i++)
          *(Q+i)=(stack[i]);
    }
    else {
      for (int j = 0; j < v; j++) {
        proj K;
        assert(!mont_iszero(Q+1));
        uintbig_div3_64(&cof, &cof, fact[i]);
        xMUL(&K, A, Q+1, &cof);

        assert(!mont_iszero(&K));
        xISOG_many(A, Q, n, &K, fact[i]);
      }
    }

  }
  assert(mont_iszero(Q+1));
}



void eval(proj *A, const odd_isogeny *phi, proj *P) {
  proj Q[3];
  Q[0] = *P;

  // The p+1 part
  Q[1] = phi->kernel_plus;
  Q[2] = phi->kernel_minus;
  eval_rtl(A, Q, 3, phi->deg_plus, p_plus_fact, p_plus_len);

  // The p-1 part
  Q[1] = Q[2];
  eval_rtl(A, Q, 2, phi->deg_minus, p_minus_fact, p_minus_len);

  *P = Q[0];
}

void eval_mult(proj *A, const odd_isogeny *phi, proj *P, int n) {
  proj Q[n+2];
  Q[0] = P[0];

  // The p+1 part
  Q[1] = phi->kernel_plus;

  for (int i = 1; i < n; ++i) {
    Q[i+1] = P[i];
  }

  Q[n+1] = phi->kernel_minus;

  eval_rtl(A, Q, n+2, phi->deg_plus, p_plus_fact, p_plus_len);

  // The p-1 part
  Q[1] = Q[n+1];
  eval_rtl(A, Q, n+1, phi->deg_minus, p_minus_fact, p_minus_len);

  P[0] = Q[0];
  for (int i = 1; i < n; ++i) {
    P[i] = Q[i+1];
  }
}

// Quite naive way to compute a dual isogeny of degree dividing p±1,
// but it seems difficult to do better with XZ-coordinates
//
// res must point to an array of at least 2 proj, of length given by reslen.
// If reslen > 2, any point in res[2:] is pushed through the isogeny
void dual_one_side(proj *A, const proj *K, bool K_is_on_curve, isog_degree deg,
		   proj *res, long reslen,
		   const uintbig* cof, const long* fact, const long *mult, long len) {
  isog_degree deg_co;
  uintbig ord, cof2;
  proj Abak, tmp;

  while (true) {
    fp2_random(&res[0].x);
    res[0].z = fp2_1;
    if (is_on_curve(res, A) != K_is_on_curve)
      continue;
    xMUL(res, A, res, cof);
    deg_co = degree_co(deg, mult, len);
    degree_to_uint(&cof2, deg_co, fact, len);
    xMUL(res, A, res, &cof2);
    // Now we have a point of order dividing deg
    res[1] = *K;
    Abak = *A;
    eval_rtl(A, res, reslen, deg, fact, len);
    reslen = 2; // stop pushing the input points
    // Now we need to test the order of res[0]
    degree_to_uint(&ord, deg, fact, len);
    long i;
    for (i = 0; i < len; i++) {
      if (degree_get(deg, i)) {
	uintbig_div3_64(&cof2, &ord, fact[i]);
	xMUL(&tmp, A, res, &cof2);
	if (mont_iszero(&tmp)) break;
      }
    }
    // res[0] generates the dual isogeny!
    if (i == len) break;
    else *A = Abak;
  }
}

void dual(proj *A, odd_isogeny *phi) {


  proj Q[3];

  // The p+1 part
  bool K_is_on_curve = true; // WARNING: this is assuming that all
			     // curves are in the isogeny class of
			     // order (p+1)
  Q[2] = phi->kernel_minus;
  long two_tors = class_mod_4 == 3 ? two_tors_height : 1;
  dual_one_side(A, &phi->kernel_plus, K_is_on_curve, phi->deg_plus,
		Q, 3, &p_plus_odd_cofactor, p_plus_fact, p_plus_mult, p_plus_len);

  // The p-1 part
  phi->kernel_minus = Q[2];
  Q[2] = Q[0];
  two_tors = class_mod_4 == 3 ? 1 : two_tors_height;
  dual_one_side(A, &phi->kernel_minus, !K_is_on_curve, phi->deg_minus,
		Q, 3, &p_minus_odd_cofactor, p_minus_fact, p_minus_mult, p_minus_len);

  phi->kernel_plus = Q[2];
  phi->kernel_minus = Q[0];

}
