#ifndef TESTS_H
#define TESTS_H
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>

/* factor (x^15 + x^3 + 1) over GF2[x]. */
int changs_poly_factor_test() {
    std::string p_str("2");
    NTL::ZZ p(NTL::INIT_VAL, p_str.c_str());
    NTL::ZZ_p::init(NTL::ZZ(p));    

    // make x^15 + x^3 + 1
    NTL::ZZ_pX f(15, 1);
    NTL::SetCoeff(f, 3, 1);
    NTL::SetCoeff(f, 0, 1);

    std::cout << "factoring " << f << std::endl;

    NTL::vec_pair_ZZ_pX_long factors;
    NTL::CanZass(factors, f, 0);

    std::cout << "factored: " << factors << std::endl;
    return 0;
}

#endif

