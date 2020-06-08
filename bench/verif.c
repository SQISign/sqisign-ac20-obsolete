#define _XOPEN_SOURCE

#include <getopt.h>
#include <inttypes.h>
#include <stdio.h>
#include <time.h>
#include <pari/pari.h>

#include "precomputed.h"
#include "sqisign.h"

static __inline__ uint64_t rdtsc(void)
{
    uint32_t hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return lo | (uint64_t) hi << 32;
}

int main(int argc, char **argv) {
  int keys = 5, sigs = 5, samples = 10, seed = 1;

  int opt;
  while ((opt = getopt(argc, argv, "k:m:s:r:h")) != -1) {
    switch (opt) {
    case 'k':
      keys = atoi(optarg);
      break;
    case 'm':
      sigs = atoi(optarg);
      break;
    case 's':
      samples = atoi(optarg);
      break;
    case 'r':
      seed = atoi(optarg);
      break;
    default:
      fprintf(stderr,
	      "Usage: %s [-k keys] [-m signatures] [-s samples] [-r seed]\n",
	      argv[0]);
      exit(-1);
    }
  }

  pari_init(800000000, 1<<18);
  init_precomputations();

  setrand(mkintn(1, seed));
  srand48(seed);

  printf("### Verify\n");
  printf("# key\tmessage\tcycles\t\tms\t\tlength\n");
  for (int k = 0; k < keys; k++) {
    public_key pk;
    secret_key sk;
    keygen(&pk, &sk);

    for (int s = 0; s < sigs; s++) {
      uintbig m;
      randombytes(m.c, 32);
      signature Sigma;
      compressed_signature comp_sigma;
      
      sign(&Sigma, &comp_sigma, &sk, &pk, &m);      
    
      for (int i = 0; i < samples; i++) {
	clock_t t = -clock();
	uint64_t c = -rdtsc();
	verif(&comp_sigma, &pk, &m, 10, 31);
	c += rdtsc();
	t += clock();

	printf("%d\t%d\t%" PRIu64 "\t%.3lf\n", k, s, c, 1000. * t / CLOCKS_PER_SEC);
      }
    }
  }

  return 0;
}
