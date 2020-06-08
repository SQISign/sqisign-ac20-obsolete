#include <stdio.h>
#include <assert.h>
#include <inttypes.h>
#include "mont.h"

static inline void print_big(const uintbig *x) {
  printf("%" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 "\n",
	 x->c[0], x->c[1], x->c[2], x->c[3]);
}


static inline void print_fp2(const fp2 *x) {
  printf("re ");
  print_big(&x->re.x);
  printf("im ");
  print_big(&x->im.x);
}

static inline unsigned long fp2_hash(fp2 x) {
  return (x.re.x.c[0]+3*x.re.x.c[1]+5*x.re.x.c[2]+7*x.re.x.c[3]
	  +11*x.im.x.c[0]+13*x.im.x.c[1]+17*x.im.x.c[2]+23*x.im.x.c[3]) % 100003;
}
static inline fp2 fp2_ratio(fp2 *x, fp2 *y) {
  fp2 tmp;
  tmp = *y;
  assert(!fp2_iszero(&tmp));
  fp2_inv(&tmp);
  fp2_mul2(&tmp, x);
  return tmp;
}
static inline void proj2_print(proj2 x) {
  if (fp2_iszero(&x.z)) { printf("(infty)"); }
  else { printf("(%lu,%lu)", fp2_hash(fp2_ratio(&x.x,&x.z)), fp2_hash(fp2_ratio(&x.y,&x.z))); }
}
