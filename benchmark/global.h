#ifndef GLOBAL_H
#define GLOBAL_H

#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <pari/pari.h>

const NTL::ZZ g_A_seed(NTL::INIT_VAL, "11302398300000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
const NTL::ZZ g_p = NTL::ZZ(NTL::INIT_VAL, "837583943092107483758343358937591");
NTL::ZZ_pX g_p5;     /* x^5 + 2 */
NTL::ZZ_pX g_F;      /* x^6 + x - 44 */

/*
 initialize global resources and set NTL's modular context to `mod g_p`.
*/
void initialize_global_state(void) {
    /* set current modular context to `mod g_p`. */
    NTL::ZZ_p::init(g_p);    

    /* make a = x^5 + 2 */
    //NTL::ZZ_pX a(5, 1);
    NTL::SetCoeff(g_p5, 5, 1);
    NTL::SetCoeff(g_p5, 0, 2);

    /* make F = x^6 + x - 44; */
    //NTL::ZZ_pX F(6, 1);
    NTL::SetCoeff(g_F, 6, 1);
    NTL::SetCoeff(g_F, 1, 1);
    NTL::SetCoeff(g_F, 0, -44);

    /* initialize pari */
    pari_init(8000000, 500000);
    return;
}

#endif

