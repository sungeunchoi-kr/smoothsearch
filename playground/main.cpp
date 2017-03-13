#include <iostream>
//#include <NTL/ZZ_p.h>
//#include <NTL/GF2X.h>
//#include <NTL/GF2XFactoring.h>
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

int main(int argc, char** argv) {
    std::string p_str("837583943092107483758343358937591");
    ZZ p(NTL::INIT_VAL, p_str.c_str());
    ZZ_p::init(ZZ(p));    

ts_t t0 = get_timestamp();
    ZZ_pX f(5, 1);
    NTL::SetCoeff(f, 4, 12353);
    NTL::SetCoeff(f, 3, 235723);
    NTL::SetCoeff(f, 2, 38948);
    NTL::SetCoeff(f, 1, 82391);
    NTL::SetCoeff(f, 0, 289383);

    std::cout << f << std::endl;

    vec_pair_ZZ_pX_long factors;
    for (int i=0; i<1; ++i) {
        NTL::CanZass(factors, f, 1);
        std::cout << "factored: " << factors << std::endl;
    }
ts_t t1 = get_timestamp();
double ms = (t1 - t0) / 1000L;
std::cout << ms << "ms." << std::endl;

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

