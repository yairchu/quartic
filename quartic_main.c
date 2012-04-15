#include <stdio.h>
#include <stdlib.h>

#include "quartic.h"

int main(int argc, char** argv) {
    if (argc != 6) {
        printf(
            "To solve: a*x^4 + b*x^3 + c*x^2 + d*x + e\n"
            "Use: quartic <a> <b> <c> <d> <e>\n"
            "\n");
        return -1;
    }
    complex_t poly[5];
    poly[0].imag = poly[1].imag = poly[2].imag = poly[3].imag = poly[4].imag = 0;
    poly[4].real = atof(argv[1]);
    poly[3].real = atof(argv[2]);
    poly[2].real = atof(argv[3]);
    poly[1].real = atof(argv[4]);
    poly[0].real = atof(argv[5]);
    complex_t sols[4];
    const int num_sols = solve_poly(4, poly, sols);
    int i;
    for (i = 0; i < num_sols; ++i)
        printf("%f + %fi\n", sols[i].real, sols[i].imag);
    return 0;
}