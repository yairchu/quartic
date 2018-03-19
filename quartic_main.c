#include <stdio.h>
#include <stdlib.h>

#include "quartic.h"

complex_t parse_complex(char* str)
{
    complex_t result;
    char* comp_end = str;
    double comp = strtod(str, &comp_end);
    if (comp_end[0] == 'i')
    {
        result.real = 0.0;
        result.imag = comp;
        return result;
    }
    result.real = comp;
    result.imag = strtod(comp_end, NULL);
    return result;
}

int main(int argc, char** argv)
{
    if (argc < 2 || argc > 6)
    {
        printf(
            "To solve: a*x^4 + b*x^3 + c*x^2 + d*x + e\n"
            "Use: quartic <a> <b> <c> <d> <e>\n"
            "\n");
        return -1;
    }
    int degree = argc - 2;
    complex_t poly[5];
    for (int i = 0; i <= degree; ++i)
        poly[degree - i] = parse_complex(argv[1+i]);
    complex_t sols[4];
    const int num_sols = solve_poly(degree, poly, sols);
    int i;
    for (i = 0; i < num_sols; ++i)
        printf("%f + %fi\n", sols[i].real, sols[i].imag);
    return 0;
}