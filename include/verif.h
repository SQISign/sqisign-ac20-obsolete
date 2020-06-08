#include "two_walks.h"
#include "isogenies.h"

// Compress a sequence of 2-walks of length len to a sequence of
// integers (n bits + 4 hint bits).
//
// zip must have space for len words
void compress(uint64_t *zip, const two_walk *walk, long len);

// Inverse of the above
//
// A is the starting curve. At the end of the routine it is the
// arrival curve.
void decompress_old(two_walk *walk, proj *A, const uint64_t *zip, long len);

// Deterministically apply an isogeny of degree 3^a·5^b from A
// TODO: merge with challenge in sqisign
void challenge_alt(proj *A, const uintbig *m);
