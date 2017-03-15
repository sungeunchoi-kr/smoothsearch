#include <iostream>
#include <iomanip>
#include <sstream>
#include <ctime>
#include "global.h"
#include "get_timestamp.h"
#include "integer_factor.h"
#include "tests.h"

char iso8601_timestamp_buf__[sizeof "2011-10-08T07:07:09Z"];
const char* get_iso8601_timestamp(void) {
    time_t now;
    time(&now);
    strftime(iso8601_timestamp_buf__, sizeof iso8601_timestamp_buf__, "%FT%TZ", gmtime(&now));
    // this will work too, if your compiler doesn't support %F or %T:
    //strftime(buf, sizeof buf, "%Y-%m-%dT%H:%M:%SZ", gmtime(&now));
    return iso8601_timestamp_buf__;
}

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
};

search_result_t
    search_block_method_1(const NTL::ZZ& page_address,
                          const long page_size) {
    /* calculate the size of the block. */
    const long block_sz = page_size;

    /* NOTE: g_p5 = x^5 + 2 */
    /* NOTE: g_F = x^6 + x - 44; */

    /* calculate (x^5 + 2)^offset_start `mod` g_F. */
    NTL::ZZ A_seed = g_A_seed + (page_address * page_size);
    NTL::ZZ_pX px_A = NTL::PowerMod(g_p5, A_seed, g_F); // g_p5^A_seed) `mod` f

    // std::cerr << px_A << std::endl;
    // std::cerr << "delta is " << block_sz << "." << std::endl;

    ts_t t0 = get_timestamp();

    NTL::vec_pair_ZZ_pX_long factors;
    NTL::ZZ bs[8];
    int bs_len;
    NTL::ZZ_pX px_A_monic;
    NTL::ZZ running_min(NTL::INIT_VAL, "10000000000000000000000000000000000000000000000000");
    long running_smallest_coeff_A_offset;
    NTL::ZZ register_a(0);

    long A_offset = 0;
    do {
        // std::cerr << std::hex << "A_offset=" << A_offset << "\n";
        conv_monic_poly(px_A_monic, px_A);

        const bool is_splitting =
            NTL::CanZassShortCircuit(bs, &bs_len, px_A_monic);
        if (is_splitting) {
            //std::cerr << "A_offset = " << A_offset << std::endl;
            //for (int k=0; k<bs_len; ++k) {
            //    std::cerr << "    bs["<<k<<"] = " << bs[k] << std::endl;
            //}

            /* get the max */
            bool a_is_max = true;
            NTL::ZZ local_max;
            NTL::conv(local_max, NTL::LeadCoeff(px_A));
            //std::cerr << "    g_p5     = " << local_max << std::endl;
            for (int k=0; k<bs_len; ++k) {
                if (bs[k] > local_max) {
                    local_max = bs[k];
                    a_is_max = false;
                }
            }

            /* if 'g_p5' is the max, then we factorize it to see if
               we can do better. */
            if (a_is_max) {
                std::stringstream ss;
                ss << NTL::LeadCoeff(px_A);
                std::string st = ss.str();
                char* factor = largest_prime_factor(st.c_str());
                NTL::ZZ new_local_max = NTL::ZZ(NTL::INIT_VAL, factor);
                free(factor);

                //std::cerr << "new_local_max = " << new_local_max << std::endl;

                /* now, we try to find max again using new_local_max 
                   as running-max initializer. */
                //std::cerr << "    g_p5     = " << local_max << std::endl;
                for (int k=0; k<bs_len; ++k) {
                    if (bs[k] > new_local_max) {
                        new_local_max = bs[k];
                    }
                }
                //std::cerr << "new_local_max (after max search) = " << new_local_max << std::endl;
                if (new_local_max < running_min) {
                    running_min = new_local_max;
                    running_smallest_coeff_A_offset = A_offset;
                }
            } else {
                //std::cerr << "    lcmax = " << local_max << std::endl;
                if (local_max < running_min) {
                    running_min = local_max;
                    running_smallest_coeff_A_offset = A_offset;
                }
            }
            //std::cerr << "    glmin = " << running_min << ';' << running_smallest_coeff_A_offset << std::endl;
        }

        px_A = NTL::MulMod(px_A, g_p5, g_F); // (px_A * g_p5) `mod` g_F
        A_offset++;
    } while(A_offset < block_sz);

    // std::cerr << (get_timestamp() - t0) / 1000.0L << "ms." << std::endl;
    // std::cerr << px_A << std::endl;

    search_result_t r;
    r.min_A_offset = running_smallest_coeff_A_offset;
    r.minvalue = running_min;
    return r;
}

const long PAGE_SZ = 65536L;
void run_search_min(const char* page_address_cstr, const char* proc_name) {
    const ts_t t0 = get_timestamp();

    const NTL::ZZ page_address(NTL::INIT_VAL, page_address_cstr);
    const search_result_t result =
        search_block_method_1(page_address, PAGE_SZ);

    const ts_t t1 = get_timestamp();
    const double ms = (t1 - t0) / 1000.0L;
    std::cerr << "test method_1: searched " << PAGE_SZ << " in " << ms << "ms."
              << " (" << (65536/PAGE_SZ)*ms << "ms per standard 2^16 block size.)" << std::endl;

    /* print result to stdout */
    std::cout << result.minvalue         << ","
              << result.min_A_offset     << ","
              << page_address_cstr       << ","
              << get_iso8601_timestamp() << ","
              << proc_name               << std::endl;
    return;
}

int main(int argc, char** argv) {
    /* initialize all the polynomials, modular context, pari, etc that
       we will be using. see `global.h`. */
    initialize_global_state();
    std::cerr << "program initialized successfully.\n";
    if (argc < 3) {
        fprintf(stderr, "Usage: program (page-address) (process_name)\n");
        exit(1);
    }

    run_search_min(argv[1], argv[2]);
    return 0;
}

