#include "two_walks.h"
#include <stdlib.h>

// An entry of the hash table
typedef struct entry {
  uint16_t hash[3];
  uint16_t a;
} entry;

static uint64_t hash(proj *j, const proj *A) {
  jinv256(j, A);
  fp2_inv(&j->z); // ouch!
  fp2_mul2(&j->x, &j->z);
  j->z = fp2_1;
  // Quite arbitrary hash function mixing some words of the j-invariant
  return (j->x.re.x.c[0] + j->x.re.x.c[3] + j->x.im.x.c[1] + j->x.im.x.c[2]) | (1l << 47);
}

static void insert(uint64_t hash, uint16_t a, entry *table, long tab_size) {
  long i;
  for (i = hash % tab_size; table[i].hash[2]; i = (i+1) % tab_size);
  table[i].hash[0] = hash;
  table[i].hash[1] = hash >> 16;
  table[i].hash[2] = hash >> 32;
  table[i].a = a;
}

static long lookup(uint64_t hash, entry *table, long tab_size, long start) {
  start = (start < 0 ? hash : (start + 1)) % tab_size;
  while (table[start].hash[2]) {
    if (table[start].hash[0] == (hash & 0xffff)
	&& table[start].hash[1] == (hash >> 16 & 0xffff)
	&& table[start].hash[2] == (hash >> 32 & 0xffff))
      return start;
    start = (start + 1) % tab_size;
  }
  return -1;
}

// Find a basis of the 2-torsion of A
//
// Outputs x(P), x(Q) and x(P-Q) of a basis (P,Q) such that [2^(n-1)]P
// = (0,0).
//
// Assumes the curve A has order p+1
static void find_basis(proj *P, proj *Q, proj *PQ, const proj *A) {
  bool oncurve = class_mod_4 == 3;
  proj P2, Q2, tmp;
  // Get first point
  while (true) {
    fp2_random(&P->x); P->z = fp2_1;
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
    fp2_random(&Q->x); Q->z = fp2_1;
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


bool MITM_cutoff(two_walk *phi, const proj *from, const proj *to, long len, long tab_size) {
  long dual_len = len - tab_size;
  proj Pf, Qf, PQf, Pt, Qt, PQt, B, C, j;
  proj dual_kers[dual_len], phi_kers[tab_size];

  // Prepare stack for isogeny evaluations
  long stacksize, ts = len - tab_size;
  ts = ts < tab_size ? tab_size : ts;
  for (stacksize = 0; ts > 1; ts >>= 1) stacksize++;
  proj stack[stacksize];

  // Compute bases  
  find_basis(&Pf, &Qf, &PQf, from);
  find_basis(&Pt, &Qt, &PQt, to);

  // This may be too large for the stack, better malloc
  long h_tab_size = (1 << tab_size) + (1 << (tab_size - 1));
  entry *h_tab = malloc(sizeof(entry) * h_tab_size);
  for (long i = 0; i < h_tab_size; i++)
    h_tab[i].hash[2] = 0;

  // Fill the table
  uintbig a;
  long cof = two_tors_height - tab_size;
  for (long i = 0; i < cof; i++) {
    xDBL(&Pf, from, &Pf);
    xDBL(&Qf, from, &Qf);
    xDBL(&PQf, from, &PQf);
  }

  for (uint64_t i = 0; i < 1 << tab_size; i++) {
    uintbig_set(&a, i);
    xBIDIM(phi_kers, from, &Pf, &a, &Qf, &uintbig_1, &PQf);
    B = *from;
    eval_walk_rec(&B, phi_kers, tab_size, false, stack, 0);
    long h = hash(&j, &B);
    insert(h, i, h_tab, h_tab_size);
  }

  // Prepare phi
  cof = two_tors_height - len;
  for (long i = 0; i < cof; i++) {
    xDBL(&Pt, to, &Pt);
    xDBL(&Qt, to, &Qt);
    xDBL(&PQt, to, &PQt);
  }
  phi->A = *from;
  phi->len = len;
  phi->ker = Pt;
  
  // Search for collisions
  for (long i = 0; i < tab_size; i++) {
    xDBL(&Pt, to, &Pt);
    xDBL(&Qt, to, &Qt);
    xDBL(&PQt, to, &PQt);
  }
  for (uint64_t i = 0; i < 1 << dual_len; i++) {
    uintbig_set(&a, i);
    xBIDIM(dual_kers, to, &Pt, &a, &Qt, &uintbig_1, &PQt);
    B = *to;
    eval_walk_rec(&B, dual_kers, dual_len, true, stack, 0);
    long h = hash(&j, &B);
    long start = -1;
    do {
      start = lookup(h, h_tab, h_tab_size, start);
      if (start >= 0) {
      	// Potential collision, recompute the walk and check
      	proj j1;
      	uintbig_set(&a, h_tab[start].a);
      	xBIDIM(phi_kers, from, &Pf, &a, &Qf, &uintbig_1, &PQf);
      	C = *from;
      	eval_walk_rec(&C, phi_kers, tab_size, true, stack, 0);
      	jinv256(&j1, &C);
      	if (mont_equal(&j, &j1)) {
      	  // Collision found, reconstruct kernel
      	  // Push phi->ker through dual_kers
      	  for (long j = 0; j < dual_len; j++) {
      	    two_isog(dual_kers+j, &phi->ker);
      	  }
      	  // Compute and evaluate isomorphism at meeting point
      	  isomorphism isom;
      	  mont_isom(&isom, &B, &C);
      	  mont_isom_apply(&isom, &phi->ker);
      	  // Push phi->ker through the dual of phi_kers
      	  for (long j = tab_size - 1; j >= 0; j--) {
      	    two_isog_dual(phi_kers+j, &phi->ker);
      	  }
      	  
      	  free(h_tab);
      	  return true;
      	}
      }
    } while (start >= 0);
  }
  
  free(h_tab);
  return false;
}

























static uint64_t hash_fp2(const fp2 *p) {
    return (p->re.x.c[0] + p->re.x.c[3] + p->im.x.c[1] + p->im.x.c[2]) | (1l << 47);
}

static void push_points(proj *curve, proj *new_stack_1, proj *new_stack_2, proj *new_stack_3,
  const proj *stack_1, const proj *stack_2, const proj *stack_3, long stacklen) {

  fp2_sq2(&curve->z, &stack_1[stacklen-1].z);
  fp2_sq2(&curve->x, &stack_1[stacklen-1].x);

  fp2_add2(&curve->x, &curve->x);
  fp2_sub3(&curve->x, &curve->z, &curve->x);
  fp2_add2(&curve->x, &curve->x);

  proj A24;
  if (1 < stacklen) {
    fp2_add3(&A24.x, &curve->z, &curve->z);    //precomputation of A24=(A+2C:4C)
    fp2_add3(&A24.z, &A24.x, &A24.x);
    fp2_add2(&A24.x, &curve->x);
  }

  for (int i = 0; i < stacklen-1; i++) {
    new_stack_1[i] = stack_1[i];
    new_stack_2[i] = stack_2[i];
    new_stack_3[i] = stack_3[i];
    two_isog(&stack_1[stacklen-1], &new_stack_1[i]);
    two_isog(&stack_1[stacklen-1], &new_stack_2[i]);
    two_isog(&stack_1[stacklen-1], &new_stack_3[i]);
    xDBLADD(&new_stack_3[i], &new_stack_2[i], &new_stack_3[i], &new_stack_2[i], &new_stack_1[i], &A24);
  }


}

static long build_list_rec(proj *A, long lenA, long len, proj **stack_1, proj **stack_2, proj **stack_3, long stacklen) {

  if (len == 0) {
    return lenA;
  }
  if (len == 1) {
    for (int j = 0; j < lenA; ++j) {

      // push points though second isogeny
      push_points(&A[lenA+j], stack_2[lenA+j], stack_1[lenA+j], stack_3[lenA+j],
        stack_2[j], stack_1[j], stack_3[j], stacklen);

      // push points though first isogeny
      push_points(&A[j], stack_1[j], stack_2[j], stack_3[j],
        stack_1[j], stack_2[j], stack_3[j], stacklen);

    }

    return lenA*2;

  } else {
    long right = (double)len * 0.5;
    long left = len - right;
    for (int j = 0; j < lenA; ++j) {
      stack_1[j][stacklen] = stack_1[j][stacklen-1];
      stack_2[j][stacklen] = stack_2[j][stacklen-1];
      stack_3[j][stacklen] = stack_3[j][stacklen-1];
      for (int i = 0; i < left; i++) {
        xDBL(&stack_1[j][stacklen], &A[j], &stack_1[j][stacklen]);
        xDBL(&stack_2[j][stacklen], &A[j], &stack_2[j][stacklen]);
        xDBL(&stack_3[j][stacklen], &A[j], &stack_3[j][stacklen]);
      }
    }
    lenA = build_list_rec(A, lenA, right, stack_1, stack_2, stack_3, stacklen+1);
    return build_list_rec(A, lenA, left, stack_1, stack_2, stack_3, stacklen);
  }
}

// compute p.x[i]/p.z[i] for i = 0..(len-1)
static void simultaneous_ratio(fp2 *res, proj* p, long len) {
  fp2 prod = fp2_1;
  for (int i = 0; i < len; ++i) {
    fp2_mul3(&res[i], &prod, &p[i].x);
    fp2_mul2(&prod, &p[i].z);
  }
  fp2_inv(&prod);
  for (int i = len-1; i >= 0; --i) {
    fp2_mul2(&res[i], &prod);
    fp2_mul2(&prod, &p[i].z);
  }
}

// P3 = P1-P2 corresponds to the (0,0) direction
// returns a list A[2^length] such that A[i] is the target of
// the isogeny of kernel P1 + a*P3
static void build_list(long length, const proj *from, proj *A, const proj *P1, const proj *P2, const proj *P3) {
  // B is an array of length 2^length
  A[0] = *from;

  long lenA = (1ull<<length);
  proj *stack_1[lenA], *stack_2[lenA], *stack_3[lenA];

  long log, len = length;
  for (log = 0; len > 1; len >>= 1) log++;
  log += 1;

  for (int i = 0; i < lenA; ++i) {
    // TODO: this is too much
    stack_1[i] = malloc(sizeof(proj)*log);
    stack_2[i] = malloc(sizeof(proj)*log);
    stack_3[i] = malloc(sizeof(proj)*log);
  }

  stack_1[0][0] = *P1;
  stack_2[0][0] = *P2;
  stack_3[0][0] = *P3;

  long lenB = build_list_rec(A, 1, length, stack_1, stack_2, stack_3, 1);

  for (int i = 0; i < lenB; ++i) {
    free(stack_1[i]);
    free(stack_2[i]);
    free(stack_3[i]);
  }
}


// Qi is a basis of the 2^length torsion of B
// K is a point in B, then pushed through the isogeny found
static bool search_collision(long *i1, proj *K, proj *C, entry *h_tab, long h_tab_size,
  const fp2 *j1_fp, long length, const proj *B,
  const proj *Q1, const proj *Q2, const proj *Q3) {

  long lenA = 1ul << length;

  proj A[lenA];

  build_list(length, B, A, Q1, Q2, Q3);

  proj j2_proj[lenA];
  for (int i = 0; i < lenA; ++i) {
    jinv256(j2_proj+i, A+i);
  }
  fp2 j2_fp[lenA];
  simultaneous_ratio(j2_fp, j2_proj, lenA);

  uintbig a;
  uint64_t h;
  two_walk phi2;
  phi2.len = length;
  phi2.A = *B;

  proj Q13;
  xADD(&Q13, Q1, Q3, Q2);

  for (int i2 = 0; i2 < lenA; ++i2) {
    h = hash_fp2(j2_fp + i2);      
    long start = -1;
    do {
      start = lookup(h, h_tab, h_tab_size, start);
      if (start >= 0) {
        // Potential collision, check
        *i1 = h_tab[start].a;


        if (fp2_equal(&j1_fp[*i1], &j2_fp[i2])) {
          // Collision, reconstitute the walk and return
          uintbig_set(&a, i2);
          xBIDIM(&phi2.ker, B, Q1, &uintbig_1, Q3, &a, &Q13);
          
          eval_walk(&phi2, C, K);

          return true;
        }
      }
    } while (start >= 0);
  }
  return false;

}


bool MITM2(two_walk *eta, const proj *from, const proj *to, long length) {
  long len2 = length/2;
  long len1 = length-len2;

  if (length < 6) {
    len1 = length;
    len2 = 0;
  }
  if (length == 0) {
    eta->len = 0;
    eta->A = *from;
    eta->ker = (proj){fp2_1,fp2_0};
    proj j1, j2;

    jinv256(&j1, from);
    jinv256(&j2, to);

    return mont_equal(&j1,&j2);
  }


  long lenA1 = 1ul << len1;

  proj A1[lenA1], P1, P2, P3;
  find_basis(&P3, &P1, &P2, from); // P3 has the (0,0) direction

  for (long i = 0; i < two_tors_height-len1; i++) {
    xDBL(&P1, from, &P1);
    xDBL(&P2, from, &P2);
    xDBL(&P3, from, &P3);
  }

  build_list(len1, from, A1, &P1, &P2, &P3);

  // Fill table

  proj j1_proj[lenA1];
  for (long i = 0; i < lenA1; ++i) {
    jinv256(j1_proj+i, A1+i);
  }
  fp2 j1_fp[lenA1];
  simultaneous_ratio(j1_fp, j1_proj, lenA1);

  long h_tab_size = (1 << len1) + (1 << (len1 - 1));
  entry *h_tab = malloc(sizeof(entry) * h_tab_size);
  for (long i = 0; i < h_tab_size; i++)
    h_tab[i].hash[2] = 0;

  uint64_t h;
  for (long i = 0; i < lenA1; ++i) {
    h = hash_fp2(j1_fp + i);
    insert(h, i, h_tab, h_tab_size);
  }

  proj K;
  bool found;
  long i1;

  uintbig a;
  proj C;
  two_walk phi1;
  phi1.len = len1;
  phi1.A = *from;
  proj P13;
  xADD(&P13, &P1, &P3, &P2);

  int head_len = (len2 < 4) ? 0 : 3;
  two_walk head_walk;
  head_walk.len = head_len;
  head_walk.A = *to;
  proj head_curve;

  proj Q1_long, Q2_long, Q3_long, Q1, Q2, Q3;
  find_basis(&Q3_long, &Q1_long, &Q2_long, to); // P3 has the (0,0) direction

  for (long i = 0; i < two_tors_height-length; i++) {
    xDBL(&Q1_long, to, &Q1_long);
    xDBL(&Q2_long, to, &Q2_long);
    xDBL(&Q3_long, to, &Q3_long);
  }
  Q1 = Q1_long;
  Q2 = Q2_long;
  Q3 = Q3_long;
  for (long i = 0; i < length-head_len; i++) {
    xDBL(&Q1, to, &Q1);
    xDBL(&Q2, to, &Q2);
    xDBL(&Q3, to, &Q3);
  }
  proj Q13;
  xADD(&Q13, &Q1, &Q3, &Q2);


  for (int head = 0; head < (1 << head_len); ++head) {
    uintbig_set(&a, head);
    xBIDIM(&head_walk.ker, to, &Q1, &uintbig_1, &Q3, &a, &Q13);
    K = Q3_long;
    eval_walk(&head_walk, &head_curve, &K);

    // compute a basis of 2^(len2-head_len) torsion in head_curve with R3 the (0,0) direction

    proj R1 = Q1_long, R2 = Q2_long, R3 = K;
    eval_walk(&head_walk, &head_curve, &R1);
    eval_walk(&head_walk, &head_curve, &R2);

    for (int i = 0; i < head_len; ++i) {
      if (head & (1 << i))
        xADD(&R1, &R1, &R3, &R2);
      else
        xADD(&R2, &R2, &R3, &R1);
      xDBL(&R3, &head_curve, &R3);
    }

    for (long i = 0; i < length - len2; i++) {
      xDBL(&R1, &head_curve, &R1);
      xDBL(&R2, &head_curve, &R2);
      xDBL(&R3, &head_curve, &R3);
    }


    found = search_collision(&i1, &K, &C, h_tab, h_tab_size, j1_fp,
      len2-head_len, &head_curve, &R1, &R2, &R3);


    if (found) {
      eta->ker = K;
      uintbig_set(&a, i1);
      xBIDIM(&phi1.ker, from, &P1, &uintbig_1, &P3, &a, &P13);

      eval_dual(&phi1, &C, &eta->ker);

      eta->len = length;
      eta->A = *from;

      free(h_tab);
      return true;
    }
  }
  free(h_tab);
  return false;
}









