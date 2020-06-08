#include "two_walks.h"
#include <assert.h>

void two_isog(const proj *K, proj *P) {
  fp2 tmp1, tmp2, tmp3;
  assert(!fp2_iszero(&K->x));
  fp2_mul3(&tmp1, &P->x, &K->x);
  fp2_mul3(&tmp2, &P->z, &K->z);
  fp2_sub2(&tmp1, &tmp2);
  fp2_mul3(&tmp2, &P->x, &K->z);
  fp2_mul3(&tmp3, &P->z, &K->x);
  fp2_sub2(&tmp2, &tmp3);
  fp2_mul2(&P->x, &tmp1);
  fp2_mul2(&P->z, &tmp2);
}

void two_isog_dual(const proj *K, proj *P) {
  fp2 tmp;
  fp2_add3(&tmp, &P->x, &P->z);
  fp2_mul2(&P->z, &P->x);
  fp2_mul2(&P->z, &K->x);
  fp2_add2(&P->z, &P->z);
  fp2_add2(&P->z, &P->z);
  fp2_sq2(&P->x, &tmp);
  fp2_mul2(&P->x, &K->z);
}

// Internal function.
//
// Careful! P is an array of points, filled up to stacklen, and it
// must contain space for up to ⎡log₂(len)⎤ points after that.
//
// advance indicates whether K is a single element or an array of
// length len. In the former case, K will contain the kernels of the
// isogeny walk. In the latter, it will only contain the kernel of the
// last step.
void eval_walk_rec(proj *A, proj *K, long len, bool advance, proj *P, long stacklen) {
  if (len == 0)
    return;
  if (len == 1) {
    // push points
    for (int i = 0; i < stacklen; i++)
      two_isog(K, P+i);
    // push curve
    fp2_sq2(&A->z, &K->z);
    fp2_sq2(&A->x, &K->x);
    fp2_add2(&A->x, &A->x);
    fp2_sub3(&A->x, &A->z, &A->x);
    fp2_add2(&A->x, &A->x);
  } else {
    long right = len / 2;
    long left = len - right;
    P[stacklen] = *K;
    for (int i = 0; i < left; i++)
      xDBL(K, A, K);
    eval_walk_rec(A, K, right, advance, P, stacklen+1);
    K[right*advance] = P[stacklen];
    eval_walk_rec(A, K+right*advance, left, advance, P, stacklen);
  }
}

void eval_walk(const two_walk *phi, proj *B, proj *P) {
  *B = phi->A;
  proj K = phi->ker;
  long log, len = phi->len;
  for (log = 0; len > 1; len >>= 1) log++;
  proj stack[1+log];
  stack[0] = *P;
  eval_walk_rec(B, &K, phi->len, false, stack, 1);
  *P = stack[0];
}

void eval_walk_mult(const two_walk *phi, proj *B, proj *P, long cardinality) {
  *B = phi->A;
  proj K = phi->ker;
  long log, len = phi->len;
  for (log = 0; len > 1; len >>= 1) log++;
  proj stack[cardinality+log];
  for (int i = 0; i < cardinality; ++i) {
    stack[i] = P[i];
  }
  eval_walk_rec(B, &K, phi->len, false, stack, cardinality);
  for (int i = 0; i < cardinality; ++i) {
    P[i] = stack[i];
  }
}


void eval_dual(const two_walk *phi, const proj *B, proj *P) {
  proj A = phi->A;
  proj K[phi->len];
  K[0] = phi->ker;
  long log, len = phi->len;
  for (log = 0; len > 1; len >>= 1) log++;
  proj stack[log];
  eval_walk_rec(&A, K, phi->len, true, stack, 0);
  // Move P to A if needed
  if (!mont_equal(&A, B)) {
    isomorphism isom;
    mont_isom(&isom, B, &A);
    mont_isom_apply(&isom, P);
  }
  // Eval the chain backwards
  for (int i = phi->len - 1; i >= 0; i--) {
    two_isog_dual(K+i, P);
  }
}

void dual_walk(two_walk* phi) {
  proj Kd, K2, tmp, tmp2, A;
  bool oncurve = class_mod_4 == 3;

  K2 = phi->ker;
  for (int i = 1; i < phi->len; i++)
    xDBL(&K2, &phi->A, &K2);

  assert(!fp2_iszero(&K2.x));
  
  while (true) {
    fp2_random(&Kd.x); Kd.z = fp2_1;
    if (is_on_curve(&Kd, &phi->A) != oncurve)
      continue;
    // multiply by cofactor
    xMUL(&Kd, &phi->A, &Kd, &p_even_cofactor);
    // check order
    tmp = Kd;
    for (int i = 1; i < phi->len; i++)
      xDBL(&tmp, &phi->A, &tmp);
    long cof;
    xDBL(&tmp2, &phi->A, &tmp);
    for (cof = 0; !mont_iszero(&tmp2); cof++) {
      tmp = tmp2;
      xDBL(&tmp2, &phi->A, &tmp2);
    }
    // check orthogonaliy to phi->ker
    if (mont_iszero(&tmp) || mont_equal(&tmp, &K2))
      continue;
    // adjust order
    for (; cof > 0; cof--)
      xDBL(&Kd, &phi->A, &Kd);
    break;
  }

  // Push dual generator
  eval_walk(phi, &A, &Kd);
  // Adjust direction
  isomorphism isom;
  rand_isom(&isom, &A);
  mont_isom_apply(&isom, &Kd);

  phi->ker = Kd;
  phi->A = A;
}


// finds isom such that phi*(isom inverse) can be evaluated
// set phi_new phi*(isom inverse)
// then set P to phi_new(isom(P)), and A to the target curve
void eval_walk_isom(isomorphism *isom, two_walk *phi_new, proj *B, proj *R, const two_walk *phi, const proj *P) {
    proj A0 = phi->A;
    *phi_new = *phi;
    if (R) { assert(P); *R = *P; }
    trivial_isom(isom);

    if (phi_new->len) {
        proj Q = phi_new->ker;

        uintbig k;
        uintbig_set(&k, 1ULL<<(phi_new->len-1));
        xMUL(&Q, &phi_new->A, &Q, &k);
        assert(!mont_iszero(&Q));
        if(fp2_iszero(&Q.x)) {
            // need to push (0,0) away
            rand_isom(isom, &phi_new->A);
            mont_isom_apply(isom, &phi_new->ker);
            if (R) mont_isom_apply(isom, R);
        }
        xDBL(&Q, &A0, &Q); // Q is still on the first curve
        assert(mont_iszero(&Q));

        if (R || B) { 
            proj B0;
            proj* B1 = (B) ? B : &B0;

            proj R0 = {fp2_1, fp2_0};
            proj* R1 = (R) ? R : &R0;

            eval_walk(phi_new, B1, R1);
        } 
    }
    else { if (B) *B = phi->A; }
}

void eval_walk_isom_mult(isomorphism *isom, two_walk *phi_new, proj *B, const two_walk *phi, proj *P, long cardinality) {
    proj A0 = phi->A;
    *phi_new = *phi;
    trivial_isom(isom);

    if (phi_new->len) {
        proj Q = phi_new->ker;

        uintbig k;
        uintbig_set(&k, 1ULL<<(phi_new->len-1));
        xMUL(&Q, &phi_new->A, &Q, &k);
        assert(!mont_iszero(&Q));
        if(fp2_iszero(&Q.x)) {
            // need to push (0,0) away
            rand_isom(isom, &phi_new->A);
            mont_isom_apply(isom, &phi_new->ker);
            for (int i = 0; i < cardinality; ++i) {
              mont_isom_apply(isom, &P[i]);
            }
        }
        xDBL(&Q, &A0, &Q); // Q is still on the first curve
        assert(mont_iszero(&Q));

        eval_walk_mult(phi_new, B, P, cardinality);
        
    }
    else { *B = phi->A; }
}

