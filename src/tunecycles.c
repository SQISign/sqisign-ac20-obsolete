#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "isogenies.h"
#include "steps.h"
#include "cycle.h"
#include "constants.h"

int comparelonglong(const void *uptr,const void *vptr)
{
  long long u = *(const long long *) uptr;
  long long v = *(const long long *) vptr;
  if (u < v) return -1;
  if (u > v) return 1;
  return 0;
}

long long median(long long *x,long long xlen)
{
  if (xlen <= 0) return 0;
  qsort(x,xlen,sizeof(long long),comparelonglong);
  if (xlen&1) return x[xlen/2];
  return (x[xlen/2-1]+x[xlen/2])/2;
}

void isog_setup(proj *A, proj *P, proj *K,
		isog_degree deg,
		const uintbig *cof, const long *fact, const long *mult, long len,
		bool twist)
{
  A->x = fp2_0;
  A->z = fp2_1;

  for (;;) {
    do {
      fp2_random(&K->x); fp2_random(&K->z);
    } while (is_on_curve(K, A) == twist);
    xMUL(K,A,K,cof);
    uintbig cof2;
    degree_to_uint(&cof2, degree_co(deg, mult, len), fact, len);
    xMUL(K,A,K,&cof2);

    if (!mont_iszero(K)) break;
  }

  uintbig cof2;
  degree_to_uint(&cof2, deg, fact, len);
  xMUL(P,A,K,&cof2);
  if (!mont_iszero(P)) abort();

  fp2_random(&P->x);
  fp2_random(&P->z);
}

#define TIMINGS 31

int main()
{
  proj A[TIMINGS];
  proj P[TIMINGS];
  proj K[TIMINGS];
  long long baseline[TIMINGS];
  long long bench[TIMINGS];

  long len = p_minus_len;
  const long *fact = p_minus_fact,
    *mult = p_minus_mult;
  const uintbig *cofactor = &p_minus_odd_cofactor;
  int twist = 1;

  do {
    for (long long lpos = 0;lpos < len;++lpos) {
      long long l = fact[lpos];
      isog_degree deg = { 0 };
      degree_set(&deg, lpos, 1);

      for (long long t = 0;t < TIMINGS;++t)
	isog_setup(&A[t],&P[t],&K[t],deg,cofactor,fact,mult,len,twist);
      for (long long t = 0;t < TIMINGS;++t) {
	baseline[t] = getticks();
	xISOG_old(&A[t],&P[t],&K[t],l);
	baseline[t] = getticks() - baseline[t];
      }
      long long baselinemedian = median(baseline,TIMINGS);

      long long bestbs = 0;
      long long bestgs = 0;
      long long bestbenchmedian = -1;

      for (long long bs = 0;bs <= 100;bs += 2) {
	for (long long gs = 0;;++gs) {
	  if (!gs) if (bs) continue;
	  if (!bs) if (gs) break;
	  if (2*bs*gs > (l-1)/2) break;
	  if (gs > bs*2) continue;
	  if (bs > gs*3) continue;

	  steps_override(bs,gs);

	  for (long long t = 0;t < TIMINGS;++t)
	    isog_setup(&A[t],&P[t],&K[t],deg,cofactor,fact,mult,len,twist);

	  for (long long t = 0;t < TIMINGS;++t) {
	    bench[t] = getticks();
	    xISOG(&A[t],&P[t],&K[t],l);
	    bench[t] = getticks() - bench[t];
	  }
	  /* XXX: check for stability, re-run if necessary */

	  long long benchmedian = median(bench,TIMINGS);

	  if (benchmedian > 0)
	    if (bestbenchmedian < 0 || benchmedian < bestbenchmedian) {
	      bestbs = bs;
	      bestgs = gs;
	      bestbenchmedian = benchmedian;
	    }
	}
      }

      printf("%lld %lld %lld %lld %lld\n",l,bestbs,bestgs,bestbenchmedian,baselinemedian);
      fflush(stdout);
    }

    len = p_plus_len;
    fact = p_plus_fact;
    mult = p_plus_mult;
    cofactor = &p_plus_odd_cofactor;
    twist--;
  } while (twist >= 0);
  
  return 0;
}
