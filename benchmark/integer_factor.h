#ifndef INTEGER_FACTOR_H
#define INTEGER_FACTOR_H

#include <pari/pari.h>

char* largest_prime_factor(const char* n) {
    char* str;
    GEN z = gp_read_str(n);

    /* debug 
    str = GENtostr(nzz);
    std::cout << str << std::endl;
    free(str);
    */

    GEN factors = factor(z);

    /* debug
    str = GENtostr(factors);
    std::cout << str << std::endl;
    free(str);
    */

    const long rows_len = nbrows(factors);
    GEN entry = gmael(factors, 1, rows_len);
    str = GENtostr(entry);

    //std::cout << str << std::endl;

    //free(factors);
    //free(entry);
    //free(z);

    return str;
}

#endif

