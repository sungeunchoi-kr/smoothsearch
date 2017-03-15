#include <iostream>
#include <iomanip>
#include <sstream>
//#include <NTL/ZZ_pX.h>
//#include <NTL/ZZ_pXFactoring.h>
#include "global.h"
#include "get_timestamp.h"
#include "integer_factor.h"
#include "tests.h"

/*----------------------------------------------------------------*/
/*using namespace NTL;
void PowerXModOptimized(NTL::ZZ_pX& hh, const NTL::ZZ& e, const NTL::ZZ_pXModulus& F) {
    long n = NumBits(e);
    long i;
    NTL::ZZ_pX h, h1;
    h.SetMaxLength(F.n);
    set(h);
    for (i=n-1; i>=0; i--) {
        if (bit(e, i)) {
            SqrMod(h1, h, F);
            MulByXMod(h, h1, F);
        } else {
            SqrMod(h, h, F);
        }
    }
    hh = h;
}*/

/*
    a: the polynomial we are trying to factor.
*/
int ctr = 0;
void form_Q(const NTL::ZZ_pX& a, const int deg_a) {
    //NTL::ZZ_pXModulus A;
    //build(A, a);

    //NTL::ZZ_pX h[deg_a];
    NTL::ZZ_pX h;

    ctr++;
    NTL::PowerXMod(h, g_p, a);

    //for (int i=1; i<deg_a; ++i) {
    //    NTL::PowerXMod(h[i], g_p*i, A);
    //}

    //for (int i=1; i<deg_a; ++i) {
    //    std::cout << h[i] << std::endl;
    //}
}
/*----------------------------------------------------------------*/


void conv_monic_poly_time(NTL::ZZ_pX& out_poly, const NTL::ZZ_pX& poly) {
    const ts_t t0 = get_timestamp();

    const NTL::ZZ_p q = NTL::LeadCoeff(poly);
    const long deg = NTL::deg(poly);
    NTL::ZZ_p c; // register
    for (long d=0; d<=deg; ++d) {
        GetCoeff(c, poly, d);
        SetCoeff(out_poly, d, c/q);
    }

    const ts_t t1 = get_timestamp();
    const double time = (t1 - t0);// / 1000.0L;
    std::cerr << "[conv_monic_poly] time: " << time << "ns.\n";
}

void conv_monic_poly(NTL::ZZ_pX& out_poly, const NTL::ZZ_pX& poly) {
    const NTL::ZZ_p q = NTL::LeadCoeff(poly);
    const long deg = NTL::deg(poly);
    NTL::ZZ_p c; // register
    for (long d=0; d<=deg; ++d) {
        NTL::GetCoeff(c, poly, d);
        NTL::SetCoeff(out_poly, d, c/q);
    }
}

struct search_result_t {
    long min_A_offset;
    NTL::ZZ minvalue;
    std::string debug_message;
};

search_result_t
     search_block_doing_nothing(const NTL::ZZ& offset_start,
                                const NTL::ZZ& offset_end,
                                const int log) {
    long block_sz;
    NTL::conv(block_sz, (offset_end - offset_start));
    long A_offset = 0;
    do {
        A_offset++;
    } while(A_offset <= block_sz);
    search_result_t r; {
        r.min_A_offset = 0;
        r.minvalue = 0;
    }
    return r;
}

/* usage:
        search_block(ZZ(0), ZZ(999999), 0);
*/
search_result_t
     search_block_method_1(const NTL::ZZ& offset_start,
                           const NTL::ZZ& offset_end,
                           const int log) {
    /* calculate the size of the block. */
    long block_sz;
    NTL::conv(block_sz, (offset_end - offset_start));

    /* set current modular context to `mod p`. */
    NTL::ZZ_p::init(g_p);    

    /* make a = x^5 + 2 */
    NTL::ZZ_pX a(5, 1);
    NTL::SetCoeff(a, 0, 2);

    /* make F = x^6 + x - 44; */
    NTL::ZZ_pX F(6, 1);
    NTL::SetCoeff(F, 1, 1);
    NTL::SetCoeff(F, 0, -44);

    /* calculate (x^5 + 2)^offset_start `mod` F. */
    NTL::ZZ A_seed = g_A_seed + offset_start;
    NTL::ZZ_pX px_A = NTL::PowerMod(a, A_seed, F); // a^(A_seed) `mod` f

    if (log >= 2) std::cout << px_A << std::endl;
    if (log >= 2) std::cout << "delta is " << block_sz << "." << std::endl;

    ts_t t0 = get_timestamp();

    NTL::vec_pair_ZZ_pX_long factors;
    NTL::ZZ bs[8];
    int bs_len;
    NTL::ZZ_pX px_A_monic;
    NTL::ZZ running_min(NTL::INIT_VAL, "10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    long running_smallest_coeff_A_offset;
    NTL::ZZ register_a(0);

    long A_offset = 0;
    do {
        if (log >= 1) std::cerr << std::hex << "A_offset=" << A_offset << "\n";
        conv_monic_poly(px_A_monic, px_A);

        const bool is_splitting =
            NTL::CanZassShortCircuit(bs, &bs_len, px_A_monic);
        if (is_splitting) {
            //std::cout << "A_offset = " << A_offset << std::endl;
            //for (int k=0; k<bs_len; ++k) {
            //    std::cout << "    bs["<<k<<"] = " << bs[k] << std::endl;
            //}

            /* get the max */
            bool a_is_max = true;
            NTL::ZZ local_max;
            NTL::conv(local_max, NTL::LeadCoeff(px_A));
            //std::cout << "    a     = " << local_max << std::endl;
            for (int k=0; k<bs_len; ++k) {
                if (bs[k] > local_max) {
                    local_max = bs[k];
                    a_is_max = false;
                }
            }

            /* if 'a' is the max, then we factorize it to see if
               we can do better. */
            if (a_is_max) {
                std::stringstream ss;
                ss << NTL::LeadCoeff(px_A);
                std::string st = ss.str();
                char* factor = largest_prime_factor(st.c_str());
                NTL::ZZ new_local_max = NTL::ZZ(NTL::INIT_VAL, factor);
                free(factor);

                //std::cout << "new_local_max = " << new_local_max << std::endl;

                /* now, we try to find max again using new_local_max 
                   as running-max initializer. */
                //std::cout << "    a     = " << local_max << std::endl;
                for (int k=0; k<bs_len; ++k) {
                    if (bs[k] > new_local_max) {
                        new_local_max = bs[k];
                    }
                }
                //std::cout << "new_local_max (after max search) = " << new_local_max << std::endl;
                if (new_local_max < running_min) {
                    running_min = new_local_max;
                    running_smallest_coeff_A_offset = A_offset;
                }
            } else {
                //std::cout << "    lcmax = " << local_max << std::endl;
                if (local_max < running_min) {
                    running_min = local_max;
                    running_smallest_coeff_A_offset = A_offset;
                }
            }
            std::cout << "    glmin = " << running_min << ';' << running_smallest_coeff_A_offset << std::endl;
        }

        px_A = NTL::MulMod(px_A, a, F); // (px_A * a) `mod` F
        A_offset++;
    } while(A_offset <= block_sz);

    if (log >= 1) std::cout << (get_timestamp() - t0) / 1000.0L << "ms." << std::endl;
    if (log >= 1) std::cout << px_A << std::endl;

    search_result_t r;
    r.min_A_offset = running_smallest_coeff_A_offset;
    r.minvalue = running_min;
    return r;
}

search_result_t
     search_block_method_2(const NTL::ZZ& offset_start,
                           const NTL::ZZ& offset_end,
                           const int log) {
    /* calculate the size of the block. */
    unsigned long block_sz;
    NTL::conv(block_sz, (offset_end - offset_start));

    /* set current modular context to `mod p`. */
    NTL::ZZ_p::init(g_p);    

    /* make a = x^5 + 2 */
    NTL::ZZ_pX a(5, 1);
    NTL::SetCoeff(a, 0, 2);

    /* make F = x^6 + x - 44; */
    NTL::ZZ_pX F(6, 1);
    NTL::SetCoeff(F, 1, 1);
    NTL::SetCoeff(F, 0, -44);

    /* calculate (x^5 + 2)^offset_start `mod` F. */
    NTL::ZZ A_seed = g_A_seed + offset_start;
    NTL::ZZ_pX px_A = NTL::PowerMod(a, A_seed, F); // a^(A_seed) `mod` f

    if (log >= 2) std::cout << px_A << std::endl;
    if (log >= 2) std::cout << "delta is " << block_sz << "." << std::endl;

    ts_t t0 = get_timestamp();

    NTL::vec_pair_ZZ_pX_long factors;
    NTL::ZZ_pX px_A_monic;
    NTL::ZZ running_smallest_coeff(0);
    long running_smallest_coeff_A_offset;
    NTL::ZZ register_a(0);

    long A_offset = 0;
    do {
        //std::cout << px_A << std::endl;
        conv_monic_poly(px_A_monic, px_A);
        form_Q(px_A_monic, deg(px_A_monic));

        px_A = NTL::MulMod(px_A, a, F); // (px_A * a) `mod` F
        A_offset++;
    } while(A_offset <= block_sz);

    if (log >= 1) std::cout << (get_timestamp() - t0) / 1000.0L << "ms." << std::endl;
    if (log >= 1) std::cout << px_A << std::endl;

    search_result_t r;
    r.min_A_offset = running_smallest_coeff_A_offset;
    r.minvalue = running_smallest_coeff;
    return r;
}

void test_method_1(void) {
    //const unsigned long BLOCK_SZ = 65536L*65536L;
    const unsigned long BLOCK_SZ = 65536L;
    //const int BLOCK_SZ = 1024;
    //const int BLOCK_SZ = 1;
    const ts_t t0 = get_timestamp();
        NTL::ZZ alpha(NTL::INIT_VAL, "44720719489781720990584648073320583148022340924048099619983541600598964895479283398724586905988448450908371");
        NTL::ZZ omega(alpha + BLOCK_SZ);
        const search_result_t result =
            search_block_method_1(alpha, omega, 0);
    const ts_t t1 = get_timestamp();
    const double ms = (t1 - t0) / 1000.0L;
    std::cerr << "test method_1: searched " << BLOCK_SZ << " in " << ms << "ms."
              << " (" << (65536/BLOCK_SZ)*ms << "ms per standard 2^16 block size.)" << std::endl;
    std::cerr << "\tfound minval: " << result.minvalue << ":" << result.min_A_offset << std::endl
              << "timing message:\n"
              << result.debug_message << std::endl;
    return;
}

void test_method_2(void) {
    //const int BLOCK_SZ = 65536;
    //const int BLOCK_SZ = 1024;
    const int BLOCK_SZ = 1;
    const ts_t t0 = get_timestamp();
        NTL::ZZ alpha(NTL::INIT_VAL, "44720719489781720990584648073320583148022340924048099619983541600598964895479283398724586905988448450908371");
        NTL::ZZ omega(alpha + BLOCK_SZ);
        const search_result_t result =
            search_block_method_2(alpha, omega, 0);
    const ts_t t1 = get_timestamp();
    const double ms = (t1 - t0) / 1000.0L;
    std::cerr << "test method_2: searched " << BLOCK_SZ << " in " << ms << "ms."
              << " (" << (65536/BLOCK_SZ)*ms << "ms per standard 2^16 block size.)" << std::endl;
    std::cerr << "\tfound minval: " << result.minvalue << ":" << result.min_A_offset << std::endl
              << "timing message:\n"
              << result.debug_message << std::endl;
    return;
}

int main(int argc, char** argv) {
    /* initialize all the polynomials, modular context, pari, etc that
       we will be using. see `global.h`. */
    initialize_global_state();
    std::cerr << "program initialized successfully.\n";

    const int BLOCK_SZ = 1000;
    //const int BLOCK_SZ = 65536;

    if (argc < 2) {
        fprintf(stderr, "Provide block address.\n");
        exit(1);
    }

    test_method_1();
    //test_method_2();

    std::cout << "ctr = " << ctr << std::endl;
    return 0;
}

