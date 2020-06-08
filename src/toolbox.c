
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <pari/pari.h>
#include <math.h>
#include <limits.h>

#include "ideal.h"
#include "toolbox.h"

static clock_t global_timer;

clock_t
tic()
{
    global_timer = clock();
    return global_timer;
}

float
tac()
{
    float ms = (1000. * (float) (clock() - global_timer) / CLOCKS_PER_SEC);
    return ms;
}

float
TAC(const char *str)
{
    float ms = (1000. * (float) (clock() - global_timer) / CLOCKS_PER_SEC);
    printf("%s [%d ms]\n", str, (int) ms);
    return ms;
}

float
toc(const clock_t t)
{
    float ms = (1000. * (float) (clock() - t) / CLOCKS_PER_SEC);
    return ms;
}

float
TOC(const clock_t t, const char *str)
{
    float ms = (1000. * (float) (clock() - t) / CLOCKS_PER_SEC);
    printf("%s [%d ms]\n", str, (int) ms);
    return ms;
}

// computes the factorisation matrix of x1 by x2 given their factorisation matrices f1 and f2
GEN famat_div(GEN f1, GEN f2) {
    pari_sp ltop = avma;
    GEN f = famat_reduce(famat_div_shallow(f1,f2));
    return gerepilecopy(ltop,f);
}

// returns the first divisor in f1 (first prime in the list with non-zero exponent)
// if f2 is given, it is set to the updated factorisation, where the returned factor has been removed
GEN famat_pop(GEN f1, GEN* f2) {
    if (lg(f1) == 1) return NULL;

    pari_sp ltop = avma;
    GEN new_f1 = gcopy(f1);
    long n = lg(gel(new_f1,1)) - 1;

    for (int i = 1; i <= n; ++i) {
        if (gcmpgs(gel(gel(new_f1,2),i), 0) > 0) {
            if (f2) {
                gel(gel(new_f1,2),i) = gsubgs(gel(gel(new_f1,2),i), 1);
                *f2 = gerepilecopy(ltop, new_f1);
                return gel(gel(*f2,1),i);
            }
            return gerepilecopy(ltop,gel(gel(new_f1,1),i));
        }
    }

    avma = ltop;
    return NULL;
}

GEN famat_degree(GEN f) {
    if (lg(f) == 1) return gen_0;

    pari_sp ltop = avma;

    unsigned long n = lg(gel(f,1)) - 1;
    GEN degree = gen_0;

    for (unsigned long i = 1; i <= n; ++i) {
        degree = gadd(degree, gel(gel(f,2),i));
    }

    return gerepilecopy(ltop,degree);
}

// returns the n-th prime divisor in f, where primes are counted with multiplicity (the 3rd prime of 2*3^2*5 is 3 because 2,3,3,5)
GEN famat_get_nth(GEN f, GEN n, unsigned long *index) {
    if (lg(f) == 1) return NULL;

    pari_sp ltop = avma;

    unsigned long m = lg(gel(f,1)) - 1;

    GEN degree = gen_0;

    for (unsigned long i = 1; i <= m; ++i) {
        degree = gadd(degree, gel(gel(f,2),i));
        if (gcmp(n, degree) <= 0) {
            if (index) *index = i;
            return gerepilecopy(ltop,gel(gel(f,1),i));
        }
    }

    avma = ltop;
    return NULL;
}


GEN famat_random(GEN f1, GEN B) {
    pari_sp ltop = avma;

    int lg_f1 = lg(gel(f1,1)) - 1;
    int m;
    long index, acc, n;

    double logs[lg_f1], logB = dbllog2r(itor(B,10));
    double sum = 0;
    long available_terms[lg_f1], chosen_terms[lg_f1], total_available = 0;

    index = 0;
    for (int i = 1; i <= lg_f1; ++i) {
        if (gcmp(gel(gel(f1,2),i),gen_0) > 0) {
            logs[index] = dbllog2r(itor(gel(gel(f1,1),i),10));
            available_terms[index] = itos_or_0(gel(gel(f1,2),i));
            total_available += available_terms[index];
            chosen_terms[index] = 0;
            index++;
        }
    }
    m = index;

    while (sum < logB && total_available > 0) {
        n = random_Fl(total_available);
        index = -1;
        acc = 0;
        while (acc <= n) {
            index++;
            acc += available_terms[index];
        }
        sum += logs[index];
        available_terms[index]--;
        chosen_terms[index]++;

        total_available--;
    }



    // remove superfluous factors (largest first)
    for (int i = m-1; i >= 0; --i) {
        for (int j = 0; j < chosen_terms[i]; ++j) {
            if (sum - logs[i] > logB ) {
                sum = sum - logs[i];
                chosen_terms[i]--;
            }
            else break;
        }
    }

    // generate result
    GEN result = cgetg(3, t_MAT);
    gel(result,1) = cgetg(m+1, t_COL);
    gel(result,2) = cgetg(m+1, t_COL);
    index = 0;
    for (int i = 1; i <= lg_f1; ++i) {
        if (gcmp(gel(gel(f1,2),i),gen_0) > 0) {
            gel(gel(result,1),index+1) = gcopy(gel(gel(f1,1),i));
            gel(gel(result,2),index+1) = stoi(chosen_terms[index]);
            index++;
        }
    }
    // m == index;

    return gerepileupto(ltop,result);
}

// returns the product
GEN famat_prod(GEN f) {
    if (lg(f) == 1) return gen_1;

    pari_sp ltop = avma;

    long m = lg(gel(f,1)) - 1;

    GEN list = cgetg(m+1, t_VEC);
    for (long i = 1; i <= m; ++i) {
        gel(list,i) = powii(gel(gel(f,1),i), gel(gel(f,2),i));
    }

    return gerepilecopy(ltop,ZV_prod(list));
}

// TODO: currently assumes LONG_IS_64BIT
int cornacchia_extended(GEN N, GEN *x, GEN *y) {

    // test for forbiden factors
    // we could also allow these factors if they appear with an even power, but does not seem to improve performance
    // 11638895555051853627 = 3*7*11*19*23*31*43*47*59*67*71*79*83 = primes that are 3 mod 4 up to 101

    if (ugcd(11638895555051853627ULL, umodiu(N,11638895555051853627ULL)) == 1) {
        // no bad small factor
        pari_sp ltop = avma;
        long valuation_2 = vali(N);
        GEN N_odd = shifti(N, -valuation_2); // remove the even part

        // 10003628061488344205 = 5*13*17*29*37*41*53*61*73*89*97*101 = primes that are 1 mod 4 up to 101
        unsigned long small_factors_1_mod_4 = ugcd(10003628061488344205ULL, umodiu(N_odd,10003628061488344205ULL));
        unsigned long gcd = small_factors_1_mod_4;

        while (gcd != 1) {
            N_odd = diviiexact(N_odd,stoi(gcd));
            gcd = ugcd(gcd, umodiu(N_odd,gcd));
            small_factors_1_mod_4 *= gcd;
        }

        if ((umodiu(N_odd, 4) == 1)) { // we hope the 'unfactored' part is a prime 1 mod 4
            if (ispseudoprime(N_odd,0)) { // the 'unfactored' part is prime, can use Cornacchia

                GEN x0,y0;
                cornacchia(gen_1, N_odd, &x0, &y0);

                GEN small_factors = factoru(small_factors_1_mod_4); // TODO: can improve that...
                GEN sol_2 = gpowgs(mkcomplex(gen_1,gen_1), valuation_2);
                GEN sol_odd = gen_1;
                GEN cx;

                for (int i = 1; i < lg(gel(small_factors,1)); ++i) {
                    switch (gel(small_factors,1)[i]) {
                        case 5: cx = mkcomplex(stoi(2), stoi(1)); break;
                        case 13: cx = mkcomplex(stoi(3), stoi(2)); break;
                        case 17: cx = mkcomplex(stoi(4), stoi(1)); break;
                        case 29: cx = mkcomplex(stoi(5), stoi(2)); break;
                        case 37: cx = mkcomplex(stoi(6), stoi(1)); break;
                        case 41: cx = mkcomplex(stoi(5), stoi(4)); break;
                        case 53: cx = mkcomplex(stoi(7), stoi(2)); break;
                        case 61: cx = mkcomplex(stoi(6), stoi(5)); break;
                        case 73: cx = mkcomplex(stoi(8), stoi(3)); break;
                        case 89: cx = mkcomplex(stoi(8), stoi(5)); break;
                        case 97: cx = mkcomplex(stoi(9), stoi(4)); break;
                        case 101: cx = mkcomplex(stoi(10), stoi(1)); break;
                    }
                    sol_odd = gmul(sol_odd, gpowgs(cx, gel(small_factors,2)[i]));
                }

                GEN sol = gmul(sol_odd, sol_2);
                sol = gmul(sol, mkcomplex(x0,y0));

                *x = gel(sol,1);
                *y = gel(sol,2);
                gerepileall(ltop,2,x,y);
                return 1;
            }
        }
        avma = ltop;
    }
    return 0;
}

// solve x^2 + y^2 + p(u^2 + v^2) = M, with (u,v) != (0,0)
// when parity != 0, ensures that (x+v) and (y+u) are not both even
// (this means that x + y*i + u*j + v*ji is not in 2*Order(1,i,(1-ji)/2, (i+j)/2)
GEN norm_equation_special(GEN p, GEN M, long parity, bool randomized) {
    pari_sp ltop = avma;
    GEN upper_bound_pari = gdivent(M,p); // gdivent is the eauclidean division
    long upper_bound = itos_or_0(upper_bound_pari);
    long u = 1, v = 0, n, delta, A; // sum_of_2_squares;
    GEN N, x, y; // fac, q

    long bound_randomized;
    if (upper_bound == 0) bound_randomized = 1<<15;
    else bound_randomized = (long)(sqrt(upper_bound/2));

    while (1) {
        if (randomized) {
            u = random_Fl(bound_randomized);
            v = random_Fl(bound_randomized);
        }

        n = u*u+v*v;

        // if we are outside the bound, find the next point within the bound if it exists, and quit if it doesn't
        if ((!is_bigint(upper_bound_pari)) && n > upper_bound) {
            if (randomized) continue;

            // we are in the bottom half, outside the bound... find the next point within

            A = u+v;
            delta = 2*upper_bound - A*A;

            if (delta < 0)  {
                avma = ltop;
                return NULL; // no solution!
            }

            u = (long)(((double)A + sqrt(delta))/2.);
            v = A-u;

            n = u*u+v*v;
        }

        N = gsub(M,gmulgs(p,n)); // N = M - p(u^2 + v^2)

        if (cornacchia_extended(N, &x, &y)) { // no bad small factor
            if ((parity == 0) || (smodis(gaddgs(x,v),2) != 0) || (smodis(gaddgs(y,u),2) != 0) ) {
                GEN res = mkvec4(x,y,stoi(u),stoi(v));
                return gerepilecopy(ltop,res);
            }
        }


        // update (u,v)
        if (!randomized) {
            if (v+1 < u) { u--; v++; }
            else { u = v+u+1; v = 0; }
        }
    }

    avma = ltop;
    return NULL;
}

GEN lattice_nearest_plane(GEN lat, GEN target, long flag) {
    pari_sp ltop = avma;

    GEN latlll = gmul(lat, qflll0(lat, flag));
    GEN latgs_sqr_len;
    GEN latgs = RgM_gram_schmidt(latlll, &latgs_sqr_len);
    GEN b = target;
    unsigned long n = lg(latlll) - 1;
    GEN c;
    for (int i = n; i > 0; --i)
    {
        c = ground(gdiv(RgV_dotproduct(b,gel(latgs,i)),gel(latgs_sqr_len,i)));
        b = gsub(b,gmul(c,gel(latlll,i)));
    }
    return gerepilecopy(ltop, gsub(target,b));
}





