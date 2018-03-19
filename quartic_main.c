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
    if (argc != 6)
    {
        printf(
            "To solve: a*x^4 + b*x^3 + c*x^2 + d*x + e\n"
            "Use: quartic <a> <b> <c> <d> <e>\n"
            "\n");
        return -1;
    }
    complex_t poly[5];
    poly[4] = parse_complex(argv[1]);
    poly[3] = parse_complex(argv[2]);
    poly[2] = parse_complex(argv[3]);
    poly[1] = parse_complex(argv[4]);
    poly[0] = parse_complex(argv[5]);
    complex_t sols[4];
    const int num_sols = solve_poly(4, poly, sols);
    int i;
    for (i = 0; i < num_sols; ++i)
        printf("%f + %fi\n", sols[i].real, sols[i].imag);
    return 0;
}