#define _XOPEN_SOURCE
#include "rng.h"
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <assert.h>

static void urandom(void *x, size_t l)
{
	// printf("WARNING: irreproducible randomness\n");
    static int fd = -1;
    ssize_t n;
    if (fd < 0 && 0 > (fd = open("/dev/urandom", O_RDONLY)))
        exit(1);
    for (size_t i = 0; i < l; i += n)
        if (0 >= (n = read(fd, (char *) x + i, l - i)))
            exit(2);
}

static void drand(void *x, size_t l)
{
  for (size_t i = 0; i < l; i += 4) {
    long b = mrand48();
    for (int j = 0; j < 4 && i+j < l; j++)
      ((char *)x)[i+j] = (b >> 8*j) & 0xff;
  }
}

void _randombytes(void *x, size_t l) {
  drand(x, l);
}

/* Ridiculous hack for cross-platform assembly */
void randombytes(void *x, size_t l)
{
  _randombytes(x, l);
}

