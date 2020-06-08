#include "isomorphism.h"
#include <assert.h>

void jinv256(proj *j, const proj *A) {
  // j(A) / 256 = (A²-3)³/(A²-4)
  fp2 X2, Z2;
  fp2_sq2(&X2, &A->x);
  fp2_sq2(&Z2, &A->z);
  fp2_add3(&j->x, &Z2, &Z2);
  fp2_sub3(&j->z, &X2, &j->x);
  fp2_sub2(&j->z, &j->x);
  fp2_add3(&j->x, &j->z, &Z2);
  fp2_mul2(&j->z, &Z2);
  fp2_mul2(&j->z, &Z2);
  fp2_mul3(&X2, &j->x, &j->x);
  fp2_mul2(&j->x, &X2);
}

void mont_isom(isomorphism *isom, const proj *A, const proj *B) {
  proj A2, B2;
  fp2 tmp, tmp2;
  
  fp2_mul3(&tmp, &A->x, &B->z);
  fp2_mul3(&A2.x, &A->z, &B->x);

  fp2_sub3(&tmp2, &tmp, &A2.x);
  fp2_add2(&tmp, &A2.x);
  // The A = B case, mapping to x -> x
  if (fp2_iszero(&tmp2)) {
    isom->D = isom->Nx = fp2_1;
    isom->Nz = fp2_0;
    return;
  }

  // The A = -B case, mapping to x -> -x
  if (fp2_iszero(&tmp)) {
    isom->D = isom->Nx = fp2_1;
    fp2_neg1(&isom->D);
    isom->Nz = fp2_0;
  }
  // The A = 0, B = ±3/√2 (j = 1728) case, mapping x -> ±(x-i)/√-2
  //
  // Watch out: in this case, the output depends on an arbitrary
  // choice for i.
  else if (fp2_iszero(&A->x)) {
    fp2_neg2(&isom->Nz, &B->x);
    fp2_mul3(&isom->Nx, &fp2_i, &B->x);
    fp2_add3(&isom->D, &B->z, &B->z);
    fp2_add2(&isom->D, &B->z);
    fp2_neg1(&isom->D);
    // Issue a warning, nevertheless
    fprintf(stderr, "WARNING: calling mont_isom on j=1728\n");
  }
  else {
    fp2_sq2(&B2.x, &B->x); fp2_sq2(&B2.z, &B->z);
    fp2_sq2(&A2.x, &A->x); fp2_sq2(&A2.z, &A->z);
    fp2_add3(&isom->Nx, &B2.z, &B2.z);
    fp2_add2(&isom->Nx, &B2.z);
    fp2_sub3(&isom->Nx, &B2.x, &isom->Nx);
    
    // We should never arrive here: if B = √3, then j=0 and A=±√3, so
    // this is either a mistake (A=B), or it has been caught earlier
    // (A=-B).
    //
    // Of course, one should probably never call this function when
    // j=0, as the isomorphism is ambiguous.
    if (fp2_iszero(&isom->Nx))
      assert(false);

    fp2_mul2(&isom->Nx, &A->x);  // Ax(Bx²-3Bz²)
      
    fp2_mul3(&tmp, &A2.z, &B2.z);
    fp2_add3(&isom->Nz, &tmp, &tmp);
    fp2_add2(&isom->Nz, &isom->Nz);
    fp2_add2(&isom->Nz, &isom->Nz);
    fp2_add2(&isom->Nz, &tmp);
    fp2_mul3(&tmp, &A2.x, &B2.z);
    fp2_sub2(&isom->Nz, &tmp);
    fp2_mul3(&tmp, &A2.z, &B2.x);
    fp2_sub2(&isom->Nz, &tmp);
    fp2_sub2(&isom->Nz, &tmp);  // 9Az²Bz² - Ax²Βz² - 2Az²Bx²

    fp2_mul3(&isom->D, &isom->Nx, &A->x);
    fp2_add2(&isom->D, &isom->Nz);
    fp2_add2(&isom->D, &isom->Nz);
    fp2_add2(&isom->D, &isom->Nz);  // 3(9Az²Bz² - Ax²Βz² - 2Az²Bx²) + Ax²(Bx²-3Bz²)

    fp2_mul2(&isom->Nx, &A->z);
    
    fp2_mul2(&isom->Nx, &B->x); // Bx ···
    fp2_mul2(&isom->Nz, &B->x); // Bx ···
    fp2_mul2(&isom->D,  &B->z);  // Bz ···
  }
}

void rand_isom(isomorphism *isom, proj *A) {
  fp2_add3(&isom->Nx, &A->z, &A->z);    // 2 Az
  fp2_mul3(&isom->D, &isom->Nx, &A->z); // 2 Az²
  fp2_mul3(&isom->Nz, &A->x, &A->x);    // Ax²
  fp2_sub2(&isom->Nz, &isom->D);        // Ax² - 2 Az²
  fp2_sub2(&isom->Nz, &isom->D);        // Ax² - 4 Az²
  fp2_sqrt(&isom->Nz);                  // √(Ax² - 4 Az²)
  fp2_sub2(&isom->Nz, &A->x);           // (α:β) = (-Ax + √(Ax² - 4 Az²) : 2 Az)
  
  fp2_mul3(&A->x, &isom->Nz, &isom->Nz);    // α²
  fp2_mul3(&isom->D, &isom->Nx, &isom->Nx); // β²
  fp2_sub3(&isom->D, &A->x, &isom->D);      // α² - β²
  fp2_add2(&A->x, &isom->D);                // 2α² - β²
  fp2_sqrt(&isom->D);                       // √(α²-β²)
  fp2_mul3(&A->z, &isom->Nz, &isom->D);     // α √(α²-β²)
}
