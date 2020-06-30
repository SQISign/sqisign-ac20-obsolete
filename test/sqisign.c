#define _XOPEN_SOURCE

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
#include "tedwards.h"
#include "isogenies.h"
#include "klpt.h"
#include "toolbox.h"
#include "sqisign.h"
#include "mont.h"

// argv[1] is the random seed; default = 1
int main(int argc, char *argv[]){
    pari_init(800000000, 1<<18);
    init_precomputations();

    setrand(stoi(1));
    srand48(1);
    if( argc > 1 ) {
      setrand(strtoi(argv[1]));
      srand48(atoi(argv[1]));
    }

    for (int i = 0; i < 10; i++) {
      uintbig m;
      randombytes(m.c, 32);

      public_key pk;
      secret_key sk;

      // printf("Key generation\n");
      keygen(&pk, &sk);

      // printf("sk->I_T\n");
      // output(sk.I_T);
      // sk.I_T = gcopy(sk.I_T);
      compressed_signature comp_sigma;
      init_compressed_sig(&comp_sigma);
      // signature Sigma;


      sign(&comp_sigma ,&sk, &pk, &m);

      assert(verif(&comp_sigma,&pk,&m));

      randombytes(m.c, 32);
      assert(!verif(&comp_sigma,&pk,&m));
      free_compressed_sig(&comp_sigma);
    }

    printf("    \033[1;32mAll tests passed\033[0m\n");
    exit(0);

    return 0;
}
