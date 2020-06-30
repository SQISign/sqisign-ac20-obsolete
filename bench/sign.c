#define _XOPEN_SOURCE

#include <getopt.h>
#include <inttypes.h>
#include <stdio.h>
#include <time.h>
#include <pari/pari.h>

#include "precomputed.h"
#include "sqisign.h"
#include "constants.h"

static __inline__ uint64_t rdtsc(void)
{
    uint32_t hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return lo | (uint64_t) hi << 32;
}

int main(int argc, char **argv) {
  int keys = 10, samples = 10, seed = 1;

  int opt;
  while ((opt = getopt(argc, argv, "k:s:r:h")) != -1) {
    switch (opt) {
    case 'k':
      keys = atoi(optarg);
      break;
    case 's':
      samples = atoi(optarg);
      break;
    case 'r':
      seed = atoi(optarg);
      break;
    default:
      fprintf(stderr,
	      "Usage: %s [-k keys] [-s samples] [-r seed]\n",
	      argv[0]);
      exit(-1);
    }
  }

  pari_init(800000000, 1<<18);
  init_precomputations();

  setrand(mkintn(1, seed));
  srand48(seed);

  printf("### Sign\n");
  printf("# key\tcycles\t\tms\t\tlength\n");
  for (int k = 0; k < keys; k++) {
    public_key pk;
    secret_key sk;
    keygen(&pk, &sk);

    for (int i = 0; i < samples; i++) {
      // signature Sigma;
      compressed_signature comp_sigma;
      init_compressed_sig(&comp_sigma);
      uintbig m;
      randombytes(m.c, 32);

      clock_t t = -clock();
      uint64_t c = -rdtsc();
      sign(&comp_sigma ,&sk, &pk, &m);
      c += rdtsc();
      t += clock();
      free_compressed_sig(&comp_sigma);
  //     int len = 0;
  //     for (int j = 0; j < Sigma.sigma.len; j++)
	// len += Sigma.sigma.phi[j].len;

      printf("%d\t%" PRIu64 "\t%.3lf\t%ld\n", k, c, 1000. * t / CLOCKS_PER_SEC, signing_length);
    }
  }

  return 0;
}
