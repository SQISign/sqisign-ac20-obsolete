
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pari/pari.h>


#include "klpt.h"
#include "ideal.h"
#include "toolbox.h"
#include "constants.h"

// return NULL if L*N(beta) is not a quadratic residue mod N
GEN klpt_strong_approximation(GEN A, GEN order, GEN p, GEN N, GEN beta, GEN gamma, GEN L, GEN lambda) {
    pari_sp ltop = avma;


    GEN C0 = gel(beta,3);
    GEN D0 = gel(beta,4);

    // the following should not happen: make sure in klpt that both C0 and D0 are invertible
    int swap = 0;
    if (gcmp(ggcd(D0,N),gen_1) != 0) {
        swap = 1;
        C0 = gel(beta,4);
        D0 = gel(beta,3);
    }

    GEN beta_norm = algnorm(A,beta,0);

    GEN p_lambda_2 = gmul(gmul(p,lambda),gen_2);
    GEN coeff_c = gmul(p_lambda_2,C0);
    GEN coeff_d = gmul(p_lambda_2,D0);
    GEN cst_term = diviiexact(gsub(L, gmul(gsqr(lambda),beta_norm)), N); // (L-lambda^2*beta_norm)/N

    GEN coeff_d_inv = Fp_inv(coeff_d,N);

    GEN cp0 = gen_0;
    GEN dp0 = Fp_mul(cst_term,coeff_d_inv,N);

    GEN lattice_basis = gmul(mkmat2(mkcol2(gen_1, gneg(Fp_mul(coeff_c,coeff_d_inv, N))), mkcol2(gen_0,N)),N);
    GEN target = gadd(gmul(mkcol2(C0,D0),lambda), gmul(mkcol2(cp0,dp0),N));

    // find a vector in lattice_basis close to target

    GEN closest = lattice_nearest_plane(lattice_basis, target, 1);
    GEN gram = gmul(gtrans(lattice_basis), lattice_basis);

    // find all the vectors *close enough* to the target

    GEN diff = gsub(target, closest);
    GEN distance = RgV_dotproduct(diff,diff);

    GEN bound0 = gdiv(L,p);
    GEN bound = gadd(gadd(bound0, distance), gmul(stoi(2),gsqrt(gmul(bound0,distance),20) ));

    GEN small_vectors = qfminim0(gram, gfloor(bound) , stoi(2000), 2, 10); // returns at most 2000 vectors...

    GEN short_v, close_v;
    for (int i = 0; i < lg(gel(small_vectors,3)); ++i) {

        if (i == 0) {
            close_v = closest;
        }
        else {
            short_v = gmul(lattice_basis, gel(gel(small_vectors,3),i));
            close_v = gadd(closest, short_v);
        }

        diff = gsub(target, close_v);
        distance = gmul(RgV_dotproduct(diff,diff), p);
        if (gcmp(distance, L) <= 0) {
            GEN rhs = diviiexact(gsub(L, distance), gsqr(N));

            GEN ap, bp;
            if (cornacchia_extended(rhs, &ap, &bp) ){
                // Solution found! Finalizing...

                GEN a = gmul(N,ap);
                GEN b = gmul(N,bp);

                GEN c = gel(diff,1);
                GEN d = gel(diff,2);
                if (swap) {
                    c = gel(diff,2);
                    d = gel(diff,1);
                }
                GEN betap = mkcol4(a,b,c,d);

                if ( Z_lval(L,2) ==0  ) {
                  return gerepilecopy(ltop, betap);
                }
                else {
                  GEN X;
                  GEN Xg;
                  alglatcontains(A,order,gamma,&Xg);
                  GEN ng = content(Xg);
                  alglatcontains(A,order,algmul(A,gamma,betap),&X);
                  GEN n =content(X);
                  //checking that the constant factor is of the desired length to ensure fixed length
                  if (  Z_lval(n,2) == 2 + Z_lval(ng,2) ) {
                    return gerepilecopy(ltop, betap);
                  }
                }


            }
        }
    }

    avma = ltop;
    return NULL;
}

// as above but N is a power of 2
GEN klpt_strong_approximation_2e(GEN A, GEN p, GEN N, GEN beta, GEN L, GEN lambda) {
    pari_sp ltop = avma;

    GEN C0 = gel(beta,3);
    GEN D0 = gel(beta,4);

    int swap = 0;
    if (gcmp(ggcd(D0,N),gen_1) != 0) {
        swap = 1;
        C0 = gel(beta,4);
        D0 = gel(beta,3);
    }

    GEN beta_norm = algnorm(A,beta,0);

    GEN p_lambda = gmul(p,lambda);
    GEN coeff_c_half = gmul(p_lambda,C0);
    GEN coeff_d_half = gmul(p_lambda,D0);
    GEN cst_term_half = diviiexact(gsub(L, gmul(gsqr(lambda),beta_norm)), gmul(N,gen_2)); // (L-lambda^2*beta_norm)/(2*N)

    GEN coeff_d_half_inv = Fp_inv(coeff_d_half,N);

    GEN cp0 = gen_0;
    GEN dp0 = Fp_mul(cst_term_half,coeff_d_half_inv,N);

    GEN lattice_basis = gmul(mkmat2(mkcol2(gen_1, gneg(Fp_mul(coeff_c_half,coeff_d_half_inv, N))), mkcol2(gen_0,N)),N);
    GEN target = gadd(gmul(mkcol2(C0,D0),lambda), gmul(mkcol2(cp0,dp0),N));

    // find a vector in lattice_basis close to target

    GEN closest = lattice_nearest_plane(lattice_basis, target, 1);
    GEN gram = gmul(gtrans(lattice_basis), lattice_basis);

    // find all the vectors *close enough* to the target

    GEN diff = gsub(target, closest);
    GEN distance = RgV_dotproduct(diff,diff);

    GEN bound0 = gdiv(L,p);
    GEN bound = gadd(gadd(bound0, distance), gmul(stoi(2),gsqrt(gmul(bound0,distance),20) ));

    GEN small_vectors = qfminim0(gram, gfloor(bound) , stoi(2000), 2, 10); // returns at most 2000 vectors...

    GEN short_v, close_v;
    for (int i = 0; i < lg(gel(small_vectors,3)); ++i) {

        if (i == 0) {
            close_v = closest;
        }
        else {
            short_v = gmul(lattice_basis, gel(gel(small_vectors,3),i));
            close_v = gadd(closest, short_v );
        }

        diff = gsub(target, close_v);
        distance = gmul(RgV_dotproduct(diff,diff), p);
        if (gcmp(distance, L) <= 0) {
            GEN rhs = diviiexact(gsub(L, distance), gsqr(N));

            GEN ap, bp;
            if (cornacchia_extended(rhs, &ap, &bp)){
                // Solution found! Finalizing...

                GEN a = gmul(N,ap);
                GEN b = gmul(N,bp);

                GEN c = gel(diff,1);
                GEN d = gel(diff,2);
                if (swap) {
                    c = gel(diff,2);
                    d = gel(diff,1);
                }

                GEN betap = mkcol4(a,b,c,d);

                return gerepilecopy(ltop, betap);
            }
        }
    }

    avma = ltop;
    return NULL;
}

// compute beta in span(j,j*i) such that gamma*beta in J with GCD(N(beta),N) = 1
GEN klpt_solve_beta(GEN A, GEN gamma, GEN J, GEN N) {
    pari_sp ltop = avma;
    GEN A_i = mkcol4s(0,1,0,0);
    GEN A_j = mkcol4s(0,0,1,0);

    GEN J_basis = lideal_basis(J);
    GEN J_scalar = lideal_scalar(J);
    GEN gamma_j = algmul(A,gamma,A_j);

    GEN matsys = gmul(mkmat2(gamma_j, algmul(A,gamma_j, A_i)), Q_denom(J_scalar));
    matsys = gconcat(matsys, gmul(J_basis, Q_remove_denom(J_scalar, NULL)));

    GEN ker;

    if (mpodd(N)) // prime case
        ker = matkermod(matsys, N, NULL); // flag = 1 because integral entries
    else // power of 2 case
        ker = matkermod(matsys, gmul(N,Q_denom(J_scalar)), NULL);

    unsigned long i = lg(ker)-1;
    GEN C0,D0;
    do {
        C0 = gel(gel(ker, i), 1);
        D0 = gel(gel(ker, i), 2);

        i--;
        if (i < 1) break;
    } while ((gcmp(C0,gen_0) == 0) || (gcmp(D0,gen_0) == 0));

    if ((gcmp(C0,gen_0) == 0) || (gcmp(D0,gen_0) == 0)) {
        avma = ltop;
        return NULL;
    }

    GEN beta = mkcol4(gen_0, gen_0, C0, D0); // beta in span(j,j*i) and gamma*beta in J
    return gerepilecopy(ltop, beta);
}


// compute beta in span(j,j*i) such that there is a zeta \in (Z/N(J)Z)^* with
// gamma*beta*delta - zeta in J*N
// N and N(J) are coprime
GEN klpt_solve_beta_special(GEN A, GEN gamma, GEN delta, GEN N, GEN J) {
    pari_sp ltop = avma;
    GEN A_1 = mkcol4s(1,0,0,0);
    GEN A_i = mkcol4s(0,1,0,0);
    GEN A_j = mkcol4s(0,0,1,0);
    GEN NJ = lideal_norm(J);

    GEN J_basis = lideal_basis(J);
    GEN J_scalar = lideal_scalar(J);

    GEN gamma_j = algmul(A,gamma,A_j);

    GEN gamma_j_delta = algmul(A,gamma_j,delta);
    GEN gamma_j_i_delta = algmul(A,algmul(A,gamma_j,A_i),delta);


    GEN matsys = gmul(mkmat3(gamma_j_delta, gamma_j_i_delta, A_1), Q_denom(J_scalar));

    matsys = gconcat(matsys, gmul(J_basis, gmul(Q_remove_denom(J_scalar, NULL),N)));


    GEN ker = matkermod(matsys, NJ, NULL); // flag = 1 because integral entries

    unsigned long i = lg(ker)-1;
    GEN C0,D0, zeta;
    do {
        C0 = gel(gel(ker, i), 1);
        D0 = gel(gel(ker, i), 2);
        zeta = gel(gel(ker, i), 3);
        i--;
        if (i < 1) break;
    } while (gcmp(zeta,gen_0) == 0);

    if (gcmp(zeta,gen_0) == 0) {
        avma = ltop;
        return NULL;
    }

    GEN beta = mkcol4(gen_0, gen_0, C0, D0); // beta in span(j,j*i) and gamma*beta in J
    return gerepilecopy(ltop, beta);
}

// runs KLPT for the left ideal I in the special order of the quaternion algebra A
// the result is an ideal equivalent to I of norm dividing the integer whose factorisation matrix is fm
// Assumes the basis of A is 1, i, j, j*i, where i^2 = -1 and j^2 = -p
GEN klpt_special_smooth_given_nearprime_J(GEN J0, GEN L0, GEN N, GEN fm) {
    pari_sp ltop = avma;

    GEN A = lideal_algebra(J0);
    GEN order = lideal_order(J0);

    GEN A_j = mkcol4s(0,0,1,0);
    GEN p = algnorm(A,A_j,0);

    GEN J = lideal_create(A, order, lideal_generator(J0), N);
    GEN p_div_N = gadd(truedivii(p, N), gen_1);

    // try 100 times then abandon
    for (int n = 1; n <= 100; ++n) {

        // find a quaternion gamma in (1,i,j,ji) of norm N*L1 where L1 divides fm
        GEN fm_1, L1;

        GEN coeff = NULL;
        unsigned long ctr = 1;

        while(!coeff) {
            fm_1 = famat_random(fm, gmulgs(p_div_N,(ctr/10) + 3));
            L1 = famat_prod(fm_1);
            coeff = norm_equation_special(p, gmul(L1,N), 0, false);
            ctr += 1;
        }

        GEN gamma = gtrans(coeff); // norm N*L1
        GEN fm_rem = famat_div(fm,fm_1); // remaining factors

        // compute beta in span(j,j*i) such that gamma*beta in J with GCD(N(beta),N) = 1
        GEN beta = klpt_solve_beta(A, gamma, J, N);
        if (!beta) continue;
        GEN beta_norm = algnorm(A,beta,0);

        if (gcmp(ggcd(beta_norm, N),gen_1) != 0) {
            /// A problem with this beta!
            continue;
        }

        // generate L2 such that L2/beta_norm is a square modulo N
        int beta_norm_is_square = Fp_issquare(beta_norm,N), L2_is_square;
        GEN bound_L2 = gmul(gmul(gmul(p,N),gmul(N,N)),stoi(80)), fm_2, L2; // (p*N^3)*MARGIN

        int safety_ctr = 0;
        GEN list_L2 = const_vec(0, gen_0); // a list of values already tried for L2
        GEN betap = NULL;

        do {
            do {
                fm_2 = famat_random(fm_rem, bound_L2);
                L2 = famat_prod(fm_2);
                L2_is_square = Fp_issquare(L2,N);
                safety_ctr++;
            } while (((L2_is_square != beta_norm_is_square) || RgV_isin(list_L2, L2)) && (safety_ctr < 40));
            if (safety_ctr >= 40) { continue; } // not enough entropy in fm_rem to ensure a solution

            list_L2 = vec_append(list_L2,L2);

            // strong approximation
            GEN lambda = Fp_sqrt(Fp_div(L2,beta_norm,N), N);

            betap = klpt_strong_approximation(A, order, p, N, beta,gamma, L2, lambda);

            if (betap) {

                // TODO: implement and use lideal_mul_elt instead
                GEN norm = gmul(gmul(L1,L2),L0);
                GEN alpha = algmul(A,gamma,betap);
                GEN generator = algmul(A,lideal_generator_coprime(J0,norm),alg_conj(A,alpha));
                GEN newideal = lideal_create(A, order, generator, norm);

                return gerepilecopy(ltop, newideal);
            }
        } while (!betap && (safety_ctr < 40));

    }


    avma = ltop;
    return NULL;
}

GEN klpt_special_smooth(GEN I, GEN fm) {
    pari_sp ltop = avma;


    // find an equivalent nearprime ideal
    GEN J0 = lideal_equiv_nearprime(I,fm,0);

    GEN fm_0 = famat_Z_gcd(fm, lideal_norm(J0));
    GEN L0 = famat_prod(fm_0);
    GEN N = diviiexact(lideal_norm(J0), L0);

    GEN fm_without_L0 = famat_div(fm,fm_0);

    GEN klpt_sol = klpt_special_smooth_given_nearprime_J(J0,L0,N,fm_without_L0);

    long ctr = 0;
    while (!klpt_sol && ctr < 50) { // fallback strategy if first J0 failed (negligible probability)
        do {
            ctr++;
            // TODO: implement and use a randomised lideal_equiv_nearprime instead
            J0 = lideal_equiv_prime_random(I,NULL,stoi(5+ctr));
        } while (!J0);

        N = lideal_norm(J0);
        klpt_sol = klpt_special_smooth_given_nearprime_J(J0,gen_1,N,fm);
    }

    if (!klpt_sol) { avma = ltop; fprintf(stderr, "klpt_special_smooth did not find a solution\n"); return NULL; }
    else return gerepilecopy(ltop, klpt_sol);
}


GEN klpt_special_smooth_given_2e_J(GEN J, GEN fm) {
    pari_sp ltop = avma;

    GEN A = lideal_algebra(J);
    GEN order = lideal_order(J);

    GEN N = lideal_norm(J);
    GEN N_2 = gmul(N,gen_2);

    GEN A_j = mkcol4s(0,0,1,0);
    GEN p = algnorm(A,A_j,0);
    GEN p_div_N = gadd(truedivii(p, N), gen_1);

    // try 100 times then abandon
    for (int n = 1; n <= 100; ++n) {

        // find a quaternion gamma in (1,i,j,ji) of norm N*L1 where L1 divides fm
        GEN fm_1, L1;

        GEN coeff = NULL;
        unsigned long ctr = 1;
        while(!coeff) {
            fm_1 = famat_random(fm, gmulgs(p_div_N,(ctr/20) + 1));
            L1 = famat_prod(fm_1);
            coeff = norm_equation_special(p, gmul(L1,N), 1, false);
            ctr += 1;
        }

        GEN gamma = gtrans(coeff); // norm N*L1
        GEN fm_rem = famat_div(fm,fm_1); // remaining factors

        // compute beta in span(j,j*i) such that gamma*beta in J with GCD(N(beta),N) = 1
        GEN beta = klpt_solve_beta(A, gamma, J, N);
        if (!beta) continue;


        if (!mpodd(gel(beta,3)) && !mpodd(gel(beta,4))) continue; // beta is divisible by 2

        GEN beta_norm = algnorm(A,beta,0);
        if (gcmp(ggcd(beta_norm, N),gen_1) != 0) {
            // printf("Problem with beta!\n");
            continue;
        }
        // beta_norm is 1 mod 4 (exactly one coefficient of beta is odd)

        // generate L2 such that L2/beta_norm is a square modulo N
        GEN lambda = NULL;
        GEN bound_L2 = gmul(gmul(gmul(p,N),gmul(N,N)),stoi(80)), fm_2, L2; // (p*N^3)*MARGIN

        int safety_ctr = 0;
        GEN list_L2 = const_vec(0, gen_0); // a list of values already tried for L2
        GEN betap = NULL;

        do {
            do {
                fm_2 = famat_random(fm_rem, bound_L2);
                L2 = famat_prod(fm_2);
                lambda = Zn_sqrt(Fp_div(L2,beta_norm,N_2), N_2);  // lambda is found mod 2*N instead of N

                safety_ctr++;
            } while (((!lambda) || RgV_isin(list_L2, L2)) && (safety_ctr < 40));
            if (safety_ctr >= 40) { continue; } // not enough entropy in fm_rem to ensure a solution

            list_L2 = vec_append(list_L2,L2);

            // strong approximation
            GEN betap = klpt_strong_approximation_2e(A, p, N, beta, L2, lambda);
            if (betap) {

                // TODO: implement and use lideal_mul_elt instead
                GEN norm = gmul(L1,L2);
                GEN alpha = algmul(A,gamma,betap);
                GEN generator = algmul(A,lideal_generator_coprime(J,norm),alg_conj(A,alpha));
                GEN newideal = lideal_create(A, order, generator, norm);

                return gerepilecopy(ltop, newideal);
            }
        }
        while ((!betap) && (safety_ctr < 40));
    }


    avma = ltop;
    return NULL;
}


GEN klpt_special_smooth_small_2e_input(GEN I, GEN fm){
    pari_sp ltop = avma;

    GEN A = lideal_algebra(I);
    GEN order = lideal_order(I);
    GEN i_plus_1 = mkcol4s(1,1,0,0);



    GEN I_red = lideal_create(A, order, lideal_generator(I), gen_2);
    GEN O_i_1 = lideal_create(A, order, i_plus_1, gen_2);
    GEN J;

    // we make sure the input ideal is a fixed point for the action of (R/2R)^*
    if (lideal_equals(I_red, O_i_1)) {
        J = I;
    }
    else {
        J = lideal_create(A, order, algmul(A,lideal_generator(I),i_plus_1), gmul(lideal_norm(I),gen_2));
    }

    GEN klpt_sol = klpt_special_smooth_given_2e_J(J,fm);

    if (!klpt_sol) { avma = ltop; fprintf(stderr, "klpt_special_smooth_small_2e_input did not find a solution\n"); return NULL; }
    else return gerepilecopy(ltop, klpt_sol);

}



GEN klpt_general_power_given_J(GEN I, GEN l, GEN J, GEN delta) {
    pari_sp ltop = avma;

    GEN A = lideal_algebra(I);
    GEN order = lideal_order(I);

    GEN A_j = mkcol4s(0,0,1,0);
    GEN p = algnorm(A,A_j,0);

    GEN N;
    N = lideal_norm(J);

    int l_is_square_N = Fp_issquare(l,N);
    int l_is_square_NI = Fp_issquare(l,lideal_norm(I));

    // find a quaternion gamma in (1,i,j,ji) of norm N*L1 where L1 divides fm
    GEN p_div_N = gadd(truedivii(p, N), gen_1);

    long L1_exp = floor(dbllog2r(itor(p_div_N,10)) / dbllog2r(itor(l,10)) );
    GEN L1 = powis(l, L1_exp);

    GEN coeff, gamma, beta1, beta1_norm;

    int beta_norm_is_square_N, parity = -1;
    unsigned int safety_ctr = 0;
    do {
        coeff = NULL;
        while((!coeff) && safety_ctr < 16) {
            L1 = gmul(L1, l);
            safety_ctr++;
            coeff = norm_equation_special(p, gmul(L1,N), 0, false);
        }

        if (!coeff) break;

        gamma = gtrans(coeff); // norm N*L1

        // compute beta1 in span(j,j*i) such that gamma*beta in J with GCD(N(beta),N) = 1
        beta1 = klpt_solve_beta(A, gamma, J, N);

        // make sure that there exists an e such that l^e*norm(beta) is a quadratic residue mod N
        if (beta1) {
            beta1_norm = algnorm(A,beta1,0);
            beta_norm_is_square_N = Fp_issquare(beta1_norm,N);

            if (beta_norm_is_square_N) {
                if (l_is_square_N) {}
                else if (!l_is_square_N) parity = 0;
            }
            else if (!beta_norm_is_square_N) {
                if (l_is_square_N) { beta1 = NULL; } // try again...
                else if (!l_is_square_N) parity = 1;
            }
        }
    }
    // sometimes a gamma is found such that gamma*j is in J, in which case beta1 is not found... a few repetitions should avoid that
    while ((!beta1) && (safety_ctr < 16)); // max 16 because we do not want L1 to get too big

    if (!beta1) {
        avma = ltop;
        // beta1 not found
        return NULL;
    }


    if (gcmp(ggcd(beta1_norm, N),gen_1) != 0) {
        // Problem with beta1
        avma = ltop;
        return NULL;
    }

    // compute beta2 in span(j,j*i) such that gamma*beta*delta/N - zeta in I with GCD(zeta,N(I)) = 1

    GEN beta2 = klpt_solve_beta_special(A, gamma, delta, N, I);

    if (!beta2) {
        avma = ltop;
        // beta2 not found
        return NULL;
    }

    GEN beta2_norm = algnorm(A,beta2,0);

    int beta_norm_is_square_NI = Fp_issquare(beta2_norm,lideal_norm(I));

    if (gcmp(ggcd(beta2_norm, N),gen_1) != 0) {
        // Problem with beta2
        avma = ltop;
        return NULL;
    }
    // recover beta as a CRT combination of beta1 and beta2
    GEN C0 = Z_chinese(gel(beta1,3), gel(beta2,3), N, lideal_norm(I));
    GEN D0 = Z_chinese(gel(beta1,4), gel(beta2,4), N, lideal_norm(I));

    GEN beta = mkcol4(gen_0,gen_0,C0,D0);

    GEN NNI = gmul(N,lideal_norm(I));
    GEN bound_L2 = gmul(gmul(gmul(p,NNI),gmul(NNI,NNI)),stoi(80)); // (p*N^3)*MARGIN

    long L2_exp = floor(dbllog2r(itor(bound_L2,10)) / dbllog2r(itor(l,10)));

    GEN beta_norm = algnorm(A,beta,0);

    int fail = 0;
    if (beta_norm_is_square_NI) {
        if (l_is_square_NI) {}
        else if (!l_is_square_NI) {
            if (parity == 1) { fail = 1; } // no way
            else { parity = 0; }
        }
    }
    else if (!beta_norm_is_square_NI) {
        if (l_is_square_NI) { fail = 1; } // no way
        else if (!l_is_square_NI) {
            if (parity == 0) { fail = 1; } // no way
            else { parity = 1; }
        }
    }

    if (fail) {
        avma = ltop;
        return NULL;
    }

    if ((parity != -1) && ((L2_exp % 2) != parity)) L2_exp++;

    GEN L2 = powis(l, L2_exp);

    for (int j = 0; j < 10; j++) {
        GEN lambda1 = Fp_sqrt(Fp_div(L2,beta_norm,N), N);
        GEN lambda2 = Fp_sqrt(Fp_div(L2,beta_norm,lideal_norm(I)), lideal_norm(I));
        GEN lambda = Z_chinese(lambda1, lambda2, N, lideal_norm(I));
        GEN mu = klpt_strong_approximation(A, order, p, NNI, beta, gamma, L2, lambda);
        if (mu) {
            GEN alpha = algmul(A,gamma,mu);
            GEN generatorJ = lideal_generator_coprime(J,l);
            GEN generator = algmul(A,generatorJ,alg_conj(A,alpha));
            GEN norm = gmul(L1,L2);
            // N(alpha) = N*norm
            GEN newideal = lideal_create(A, order, generator, norm);
            return gerepilecopy(ltop, newideal);
        }

        L2 = gmul(L2,l);
        if (parity != -1) L2 = gmul(L2,l);
    }

    return NULL;
}

GEN klpt_general_power_given_J_fixed_norm(GEN I, GEN l, GEN J, GEN delta) {
    pari_sp ltop = avma;
    // clock_t t;

    GEN A = lideal_algebra(I);
    GEN order = lideal_order(I);

    GEN A_j = mkcol4s(0,0,1,0);
    GEN p = algnorm(A,A_j,0);

    GEN N;
    N = lideal_norm(J);

    // find a quaternion gamma in (1,i,j,ji) of norm N*L1 where L1 divides fm
    GEN p_div_N = gadd(truedivii(p, N), gen_1);

    long L1_exp = floor(dbllog2r(itor(p_div_N,10)) / dbllog2r(itor(l,10)) )+13;
    GEN L1 = powis(l, L1_exp);




    GEN coeff, gamma, beta1, beta1_norm,beta2,beta2_norm,beta,beta_norm,mu;

    unsigned int safety_ctr = 0;

    GEN NNI = gmul(N,lideal_norm(I));

    long L2_exp;
    GEN L2;
    do {
      safety_ctr ++;
      coeff = norm_equation_special(p,gmul(L1,N), 0, true);
      if (!coeff) break;
      gamma = gtrans(coeff); // norm N*L1

      GEN X1;
      alglatcontains(A,order,gamma,&X1);
      GEN n1 =content(X1);

      //overestimate the size, so that after removing the constant factor, the lenght is exacltly signing_length
      L2_exp = signing_length - L1_exp + 4 + 2*Z_lval(ggcd(n1,L1),2);
      L2 = powis(l, L2_exp);

      // compute beta1 in span(j,j*i) such that gamma*beta in J with GCD(N(beta),N) = 1
      beta1 = klpt_solve_beta(A, gamma, J, N);
      if (!beta1) {
          avma = ltop;
          // beta1 not found
          return NULL;
      }
      beta1_norm = algnorm(A,beta1,0);

      if (gcmp(ggcd(beta1_norm, N),gen_1) != 0) {
          // Problem with beta1
          avma = ltop;
          return NULL;
      }
      // compute beta2 in span(j,j*i) such that gamma*beta*delta/N - zeta in I with GCD(zeta,N(I)) = 1
      beta2 = klpt_solve_beta_special(A, gamma, delta, N, I);

      if (!beta2) {
          avma = ltop;
          return NULL;
          // beta2 not found
      }

      beta2_norm = algnorm(A,beta2,0);


      if (gcmp(ggcd(beta2_norm, N),gen_1) != 0) {
          // Problem with beta2
          avma = ltop;
          return NULL;
      }
      // recover beta as a CRT combination of beta1 and beta2
      GEN C0 = Z_chinese(gel(beta1,3), gel(beta2,3), N, lideal_norm(I));
      GEN D0 = Z_chinese(gel(beta1,4), gel(beta2,4), N, lideal_norm(I));
      beta = mkcol4(gen_0,gen_0,C0,D0);
      beta_norm = algnorm(A,beta,0);
      if ( !(Fp_issquare( Fp_div(L2,beta_norm,N), N ) &&
       Fp_issquare( Fp_div(L2,beta_norm,lideal_norm(I)), lideal_norm(I)) ) ){
        beta=NULL;      }

      if (beta) {
        GEN lambda1 = Fp_sqrt(Fp_div(L2,beta_norm,N), N);
        GEN lambda2 = Fp_sqrt(Fp_div(L2,beta_norm,lideal_norm(I)), lideal_norm(I));
        GEN lambda = Z_chinese(lambda1, lambda2, N, lideal_norm(I));
        mu = klpt_strong_approximation(A, order, p, NNI, beta,gamma, L2, lambda);
      }
    } while ( !beta && !mu && (safety_ctr < 100));
    if (!beta || !mu) {
        avma = ltop;
        return NULL;
        // beta2 not found
    }

    if (mu) {
        GEN alpha = algmul(A,gamma,mu);

        GEN generatorJ = lideal_generator_coprime(J,l);
        GEN generator = algmul(A,generatorJ,alg_conj(A,alpha));
        GEN norm = gmul(L1,L2);
        GEN newideal = lideal_create(A, order, generator, norm);
        return gerepilecopy(ltop, newideal);
    }

    return NULL;
}


GEN klpt_general_power_small_J(GEN I, GEN K, GEN l, GEN list_previous_NJ, GEN* NJ) {
    pari_sp ltop = avma;

    // find prime ideal J equivalent to K
    // delta in K and K*conj(delta)/N(K) = J
    GEN delta, J;

    J = lideal_equiv_prime_except(K,&delta,list_previous_NJ);

    GEN klpt_sol = klpt_general_power_given_J_fixed_norm(I, l, J, delta);
    if (NJ && klpt_sol) {
        *NJ = lideal_norm(J);
        gerepileall(ltop, 2, &klpt_sol, NJ);
    }
    else if (NJ) {
        *NJ = gerepilecopy(ltop, lideal_norm(J));
    }
    else if (klpt_sol) {
        klpt_sol = gerepilecopy(ltop, klpt_sol);
    }
    return klpt_sol;
}


GEN klpt_general_power_random_J(GEN I, GEN K, GEN l, GEN bound_coeff_J) {
    pari_sp ltop = avma;

    // find prime ideal equivalent to K
    // delta in K and K*conj(delta)/N(K) = J
    GEN delta, J;

    J = lideal_equiv_prime_random(K,&delta,bound_coeff_J);
    if (!J) { avma = ltop; return NULL; }

    return klpt_general_power_given_J_fixed_norm(I, l, J, delta);
}


// I and K are both ideals is O0, of coprime norm
// Finds an l-path between the right order of I and the right order of I \cap K
// More precisely, finds an O0-ideal M of norm a power of l such that O_R(I \cap K) = O_R(I \cap M)
// Also assume that I has prime norm

GEN klpt_general_power(GEN I, GEN K, GEN l) {
    GEN list_NJ = const_vec(0, gen_0), NJ = gen_0;
    GEN klpt_sol = NULL;
    float bound_coeff_J = 5;
    int ctr = 0;

    while (!klpt_sol) {
        // start looking for the smallest equivalent primes, then try random ones
        if ((ctr < 3) && false) klpt_sol = klpt_general_power_small_J(I, K, l, list_NJ, &NJ);
        else {
            klpt_sol = klpt_general_power_random_J(I, K, l, stoi(floor(bound_coeff_J)));
            bound_coeff_J += 1.;
        }
        list_NJ = vec_append(list_NJ,NJ);
        ctr++;
    }
    return klpt_sol;
}
