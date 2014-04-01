#include <stdio.h>
#include <stdlib.h>

#include "quartic_real.h"

int main(int argc, char** argv) {
    if (argc != 6) {
        printf(
            "To solve: a*x^4 + b*x^3 + c*x^2 + d*x + e\n"
            "Use: quartic <a> <b> <c> <d> <e>\n"
            "\n");
        return -1;
    }
    double poly[5];
    poly[4] = atof(argv[1]);
    poly[3] = atof(argv[2]);
    poly[2] = atof(argv[3]);
    poly[1] = atof(argv[4]);
    poly[0] = atof(argv[5]);
    double sols[4];
    const int num_sols = solve_real_poly(4, poly, sols);
    int i;
    for (i = 0; i < num_sols; ++i)
        printf("%f\n", sols[i]);
    return 0;
}