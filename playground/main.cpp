#include <iostream>
#include <vector>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>

#include <sys/time.h>
typedef unsigned long long ts_t;

static ts_t get_timestamp() {
    struct timeval now;
    gettimeofday (&now, NULL);
    return  now.tv_usec + (ts_t)now.tv_sec * 1000000;
}

/* usage */
//ts_t t0 = get_timestamp();
//ts_t t1 = get_timestamp();
//double secs = (t1 - t0) / 1000000.0L;
//double ms = (t1 - t0) / 1000.0L;

NTL_CLIENT
using namespace NTL;

/* factor (x^15 + x^3 + 1) over GF2[x]. */
int changs_poly_factor_test() {
    std::string p_str("2");
    ZZ p(NTL::INIT_VAL, p_str.c_str());
    ZZ_p::init(ZZ(p));    

    // make x^15 + x^3 + 1
    ZZ_pX f(15, 1);
    NTL::SetCoeff(f, 3, 1);
    NTL::SetCoeff(f, 0, 1);

    std::cout << "factoring " << f << std::endl;

    vec_pair_ZZ_pX_long factors;
    NTL::CanZass(factors, f, 0);

    std::cout << "factored: " << factors << std::endl;
    return 0;
}

const char* g_p_str = "837583943092107483758343358937591";
const char* g_A_seed_str = "11302398300000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
NTL::ZZ g_A_seed(NTL::INIT_VAL, g_A_seed_str);
NTL::ZZ g_p = ZZ(NTL::INIT_VAL, g_p_str);

void conv_monic_poly(ZZ_pX& out_poly, const ZZ_pX& poly) {
    //std::cout << "before: " << poly << std::endl;
    NTL::ZZ_p q = NTL::LeadCoeff(poly);
    const long deg = NTL::deg(poly);
    for (long d=0; d<=deg; ++d) {
        ZZ_p c;
        GetCoeff(c, poly, d);
        SetCoeff(out_poly, d, c/q);
    }
    //std::cout << "after: " << out_poly << std::endl;
}
/* usage:
    for (int i=0; i<1; ++i)
        generate_polys(ZZ(0), ZZ(999999), 0);
*/
int work_block(const NTL::ZZ& offset_start,
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

    vec_pair_ZZ_pX_long factors;
    NTL::ZZ_pX px_A_monic;
    NTL::ZZ running_smallest_coeff(0);
    NTL::ZZ register_a(0);

    long A_offset = 0;
    do {
        conv_monic_poly(px_A_monic, px_A);
        NTL::CanZass(factors, px_A_monic, 0); // px_A_monic is only missing the constant
                                              // factor `LeadCoeff(px_A)`. `LeadCoeff(px_A)`
                                              // is in fact the a_1*a_2*...*a_n.
        //std::cout << "factors: " << factors << std::endl;
        //return 0;

        /* check if the polynomial is 1-smooth. */
        bool is_one_smooth = true;
        for (int i=0; i<factors.length(); ++i) {
            auto factor = factors[i];
            NTL::ZZ_pX poly = factor.a;
            const long deg = NTL::deg(poly);
            if (deg > 1) {
                is_one_smooth = false;
                break;
            }
        }

        if (is_one_smooth) {
            std::cout << "A=" << A_offset << "; ONE-SMOOTH: " << factors << std::endl;
            //if (running_smallest_coeff
            for (int i=0; i<factors.length(); ++i) {
                auto factor = factors[i];
                NTL::ZZ_pX poly = factor.a;
                //std::cout << "\tConstant term: " << NTL::ConstTerm(poly) << std::endl;
                if (running_smallest_coeff == 0) {
                    NTL::conv(running_smallest_coeff, NTL::ConstTerm(poly));
                    std::cout << "\tset initial running_smallest_coeff to " << running_smallest_coeff << std::endl;
                } else {
                    NTL::conv(register_a, NTL::ConstTerm(poly));
                    if (register_a < running_smallest_coeff) {
                        running_smallest_coeff = register_a;
                        std::cout << "\tupdated running_smallest_coeff to " << running_smallest_coeff << std::endl;
                        NTL::conv(register_a, NTL::LeadCoeff(px_A));
                        if (register_a > running_smallest_coeff) {
                            std::cout << "\t**** the `a_i` (" << register_a << ") is greater than `running_smallest_coeff`!" << std::endl;
                            // then, we must factor a_i and check that each factor is lesser than `running_smallest_coeff`.
                        }
                    }
                }
            }
        }

        px_A = NTL::MulMod(px_A, a, F); // (px_A * a) `mod` F
        A_offset++;
    } while(A_offset <= block_sz);

    if (log >= 1) std::cout << (get_timestamp() - t0) / 1000.0L << "ms." << std::endl;
    if (log >= 1) std::cout << px_A << std::endl;
    return 0;
}

int main(int argc, char** argv) {
    ZZ_p::init(ZZ(NTL::INIT_VAL, g_p_str));    

ts_t t0 = get_timestamp();

    const int block_sz = 0xFFF;
    work_block(ZZ(0), ZZ(block_sz), 0);

ts_t t1 = get_timestamp();
const double ms = (t1 - t0) / 1000.0L;
std::cout << "searched " << block_sz << " in " << ms << "ms."
          << " (" << (65536/block_sz)*ms << "ms per standard 2^16 block size.)" << std::endl;

//ts_t t0 = get_timestamp();
//    ZZ_pX f(5, 1);
//    NTL::SetCoeff(f, 4, 12353);
//    NTL::SetCoeff(f, 3, 235723);
//    NTL::SetCoeff(f, 2, 38948);
//    NTL::SetCoeff(f, 1, 82391);
//    NTL::SetCoeff(f, 0, 289383);
//
//    std::cout << f << std::endl;
//
//    vec_pair_ZZ_pX_long factors;
//    for (int i=0; i<1000; ++i) {
//        NTL::CanZass(factors, f, 0);
//        //std::cout << "factored: " << factors << std::endl;
//    }
//ts_t t1 = get_timestamp();
//double ms = (t1 - t0) / 1000L;
//std::cout << ms << "ms." << std::endl;

    return 0;
}

/*NTL_CLIENT
long getSmoothness(const GF2X& f) {
	static vec_pair_GF2X_long r;
	DDF(r, f);
	return r[r.length()-1].b;
}

int main(int argc, char* argv[])
{
	ZZ a, inc;
	if (argc == 3) {
		// Use commandline specified a and inc
		a = conv<ZZ>(argv[1]);
		inc = conv<ZZ>(argv[2]);
	} else {
		// Use default a and inc
		a = ZZ(113851762);
		inc = ZZ(1000000000);
	}

	// Modulus Polynomial
	GF2X m = GF2X(607, 1);
	SetCoeff(m, 105);
	SetCoeff(m, 0);

	// Base Polynomial 
	GF2X b = GF2X(606, 1);
	SetCoeff(b, 165);
	SetCoeff(b, 0);

	GF2X r;
	long s;
	long smoothest = 606;
	int iterations = 0;

	for (;;) {
		PowerMod(r, b, a, m);
		s = getSmoothness(r);
		if (s < smoothest) {
			smoothest = s;
			std::cout << "HIT:" << s << ":" << a << std::endl;
		}
		a += inc;
		iterations++;
		if (iterations == 10000) {
			iterations = 0;
			std::cout << "ITR:" << a << std::endl;
		}
	}
}
*/

