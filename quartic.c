/* Solve Polynomials of up to the fourth degree. (over complex numbers)
 * Algorithms by Ferrari, Tartaglia, Cardano, et al. (16th century Italy)
 */

#include <math.h>

#include "quartic.h"

#define MAX_DEGREE 4
#define PI (3.141592653589793)

static double stableness_score(complex_t a, complex_t b);
static int solve_normalized_poly(int degree, const complex_t* poly, complex_t* results);
static void calc_shifted_coefs(complex_t shift, int degree, const complex_t* src, complex_t* dst);
static void calc_binomials(int num_binoms, int stride, double* dst);
static void calc_powers(complex_t x, int max_power, complex_t* dst);
static int solve_depressed_poly(int degree, const complex_t* poly, complex_t* results);
static int solve_depressed_quartic(const complex_t* poly, complex_t* results);
static int solve_depressed_cubic(const complex_t* poly, complex_t* results);
static int solve_depressed_quadratic(const complex_t* poly, complex_t* results);

static int complex_eq(complex_t a, complex_t b);
static complex_t complex_from_real(double x);
static complex_t complex_add(complex_t a, complex_t b);
static complex_t complex_mult(complex_t a, complex_t b);
static complex_t complex_mult_real(double c, complex_t x);
static complex_t complex_div(complex_t a, complex_t b);
static double complex_sqr_norm(complex_t x);
static complex_t complex_inverse(complex_t x);
static complex_t complex_negate(complex_t x);
static complex_t complex_sqrt(complex_t x);
static complex_t complex_pow_real(complex_t x, double power);

/* poly: pointer to coefficients array of size degree + 1.
 * results: pointer to results output array of size degree.
 */
int solve_poly(int degree, const complex_t* poly, complex_t* results) {
    complex_t normalized_poly[MAX_DEGREE + 1];
    int i;
    const complex_t a = poly[degree];
    if (complex_eq(a, complex_from_real(0.0)))
        return solve_poly(degree - 1, poly, results);
    if (degree > MAX_DEGREE)
        return -1;
    if (degree > 2 && stableness_score(poly[degree], poly[degree - 1]) > stableness_score(poly[0], poly[1])) {
        complex_t rev_poly[MAX_DEGREE + 1];
        int i, num_results;
        for (i = 0; i <= degree; ++i)
            rev_poly[i] = poly[degree - i];
        num_results = solve_poly(degree, rev_poly, results);
        for (i = 0; i < num_results; ++i)
            results[i] = complex_inverse(results[i]);
        return num_results;
    }
    for (i = 0; i < degree; ++i)
        normalized_poly[i] = complex_div(poly[i], a);
    normalized_poly[degree] = complex_from_real(1.0);
    return solve_normalized_poly(degree, normalized_poly, results);
}

static double stableness_score(complex_t a, complex_t b) {
    const double t = fabs(complex_sqr_norm(a) / complex_sqr_norm(b));
    return t + 1.0 / t;
}

/* Normalized polynomials have the form of
 *   x^n + a*x^(n-1) + ..
 * The coefficient for x^n is one.
 * solve_normalized_poly does expect to get this coefficient despite it being known.
 */
static int solve_normalized_poly(int degree, const complex_t* poly, complex_t* results) {
    const complex_t shift = complex_mult_real(-1.0 / (double) degree, poly[degree - 1]);
    complex_t shifted_coefs[MAX_DEGREE + 1];
	int i, num_results;
    calc_shifted_coefs(shift, degree, poly, shifted_coefs);
    num_results = solve_depressed_poly(degree, shifted_coefs, results);
    for (i = 0; i < num_results; ++i)
        results[i] = complex_add(results[i], shift);
    return num_results;
}

static void calc_shifted_coefs(complex_t shift, int degree, const complex_t* src, complex_t* dst) {
    double binomials[MAX_DEGREE + 1][MAX_DEGREE + 1];
    complex_t shift_powers[MAX_DEGREE + 1];
    int dst_i, src_i;
    for (dst_i = 0; dst_i <= degree; ++dst_i)
        dst[dst_i] = complex_from_real(0.0);
    calc_binomials(degree+1, sizeof(binomials[0]) / sizeof(binomials[0][0]), binomials[0]);
    calc_powers(shift, degree, shift_powers);
    for (src_i = 0; src_i <= degree; ++src_i)
        for (dst_i = 0; dst_i <= src_i; ++dst_i)
            dst[dst_i] = complex_add(dst[dst_i], complex_mult_real(binomials[src_i][dst_i], complex_mult(src[src_i], shift_powers[src_i - dst_i])));
}

static void calc_binomials(int num_binoms, int stride, double* dst) {
    int row;
    for (row = 0; row < num_binoms; ++row) {
        const int row_idx = row * stride;
        const int prev_row_idx = (row - 1) * stride;
        int col;
        dst[row_idx] = 1;
        for (col = 1; col < row; ++col) {
            dst[row_idx + col] = dst[prev_row_idx + col - 1] + dst[prev_row_idx + col];
        }
        dst[row_idx + row] = 1;
    }
}

static void calc_powers(complex_t x, int max_power, complex_t* dst) {
    int i;
    dst[0] = complex_from_real(1.0);
    if (max_power >= 1)
        dst[1] = x;
    for (i = 2; i <= max_power; ++i)
        dst[i] = complex_mult(x, dst[i - 1]);
}

/* Depressed polynomials have the form of:
 *   x^n + a*x^(n-2) + ..
 * The coefficient for x^n is 1 and for x^(n-1) is zero.
 * So it gets 3 coefficients for a depressed quartic polynom.
 */
static int solve_depressed_poly(int degree, const complex_t* poly, complex_t* results) {
    if (degree > 0 && complex_eq(poly[0], complex_from_real(0.0))) {
        results[0] = complex_from_real(0.0);
        return 1 + solve_depressed_poly(degree - 1, poly + 1, results + 1);
    }
    switch (degree) {
    case 4:
        return solve_depressed_quartic(poly, results);
    case 3:
        return solve_depressed_cubic(poly, results);
    case 2:
        return solve_depressed_quadratic(poly, results);
    case 1:
        results[0] = complex_from_real(0.0);
        return 1;
    case 0:
        return 0;
    default:
        return -1;
    }
}

/* Based on http://en.wikipedia.org/wiki/Quartic_function#Quick_and_memorable_solution_from_first_principles */
static int solve_depressed_quartic(const complex_t* poly, complex_t* results) {
    complex_t helper_cubic[4];
    complex_t helper_results[3];
    complex_t quadratic_factor[3];
	complex_t p, c_plus_p_sqr, d_div_p;
    const complex_t e = poly[0];
    const complex_t d = poly[1];
    const complex_t c = poly[2];
	int num_results;
    if (complex_eq(d, complex_from_real(0.0))) {
		int i, num_quad_results;
        complex_t quadratic[3];
        quadratic[0] = e;
        quadratic[1] = c;
        quadratic[2] = complex_from_real(1.0);
        complex_t quadratic_results[2];
        num_quad_results = solve_poly(2, quadratic, quadratic_results);
        for (i = 0; i < num_quad_results; ++i) {
            const complex_t s = complex_sqrt(quadratic_results[i]);
            results[2*i] = complex_negate(s);
            results[2*i + 1] = s;
        }
        return 2 * num_quad_results;
    }
    helper_cubic[0] = complex_negate(complex_mult(d, d));
    helper_cubic[1] = complex_add(complex_mult(c, c), complex_mult_real(-4.0, e));
    helper_cubic[2] = complex_mult_real(2.0, c);
    helper_cubic[3] = complex_from_real(1.0);
    if (solve_poly(3, helper_cubic, helper_results) < 1)
        return 0;
    p = complex_sqrt(helper_results[0]);
    c_plus_p_sqr = complex_add(c, complex_mult(p, p));
    d_div_p = complex_div(d, p);
    quadratic_factor[0] = complex_add(c_plus_p_sqr, complex_negate(d_div_p));
    quadratic_factor[1] = complex_mult_real(2.0, p);
    quadratic_factor[2] = complex_from_real(2.0);
    num_results = solve_poly(2, quadratic_factor, results);
    quadratic_factor[0] = complex_add(c_plus_p_sqr, d_div_p);
    quadratic_factor[1] = complex_negate(quadratic_factor[1]);
    return num_results + solve_poly(2, quadratic_factor, results + num_results);
}

/* Based on http://en.wikipedia.org/wiki/Cubic_equation#Cardano.27s_method */
static int solve_depressed_cubic(const complex_t* poly, complex_t* results) {
    const complex_t q = poly[0];
    const complex_t p = poly[1];
    complex_t t, u, cubic_root_of_unity;
    int i;
    if (complex_eq(p, complex_from_real(0.0))) {
        results[0] = complex_pow_real(complex_negate(q), 1.0/3.0);
        return 1;
    }
    t = complex_add(
        complex_mult_real(0.25, complex_mult(q, q)),
        complex_mult_real(1.0/27.0, complex_mult(p, complex_mult(p, p))));
    cubic_root_of_unity.real = -0.5;
    cubic_root_of_unity.imag = 0.5 * sqrt(3.0);
    for (i = 0; i < 3; ++i) {
        if (i == 0)
            u = complex_pow_real(complex_add(complex_mult_real(-0.5, q), complex_sqrt(t)), 1.0/3.0);
        else
            u = complex_mult(u, cubic_root_of_unity);
        results[i] = complex_add(u, complex_div(p, complex_mult_real(-3.0, u)));
    }
    return 3;
}

static int solve_depressed_quadratic(const complex_t* poly, complex_t* results) {
    const complex_t t = complex_sqrt(complex_negate(poly[0]));
    results[0] = complex_negate(t);
    results[1] = t;
    return 2;
}

static int complex_eq(complex_t a, complex_t b) {
    return a.real == b.real && a.imag == b.imag;
}

static complex_t complex_from_real(double x) {
    complex_t result;
    result.real = x;
    result.imag = 0;
    return result;
}

static complex_t complex_add(complex_t a, complex_t b) {
    complex_t result;
    result.real = a.real + b.real;
    result.imag = a.imag + b.imag;
    return result;
}

static complex_t complex_mult(complex_t a, complex_t b) {
    complex_t result;
    result.real = a.real * b.real - a.imag * b.imag;
    result.imag = a.real * b.imag + a.imag * b.real;
    return result;
}

static complex_t complex_mult_real(double c, complex_t x) {
    complex_t result;
    result.real = c * x.real;
    result.imag = c * x.imag;
    return result;
}

static complex_t complex_div(complex_t a, complex_t b) {
    return complex_mult(a, complex_inverse(b));
}

static double complex_sqr_norm(complex_t x) {
    return x.real * x.real + x.imag * x.imag;
}

static complex_t complex_inverse(complex_t x) {
    const double sqr_norm = complex_sqr_norm(x);
    complex_t result;
    result.real = x.real / sqr_norm;
    result.imag = -x.imag / sqr_norm;
    return result;
}

static complex_t complex_negate(complex_t x) {
    complex_t result;
    result.real = -x.real;
    result.imag = -x.imag;
    return result;
}

/* Based on http://en.wikipedia.org/wiki/Square_root#Algebraic_formula */
static complex_t complex_sqrt(complex_t x) {
    const double r = sqrt(complex_sqr_norm(x));
    complex_t result;
    result.real = sqrt(0.5 * (r + x.real));
    result.imag = sqrt(0.5 * (r - x.real));
    return result;
}

static complex_t complex_pow_real(complex_t x, double power) {
    const double norm = sqrt(pow(complex_sqr_norm(x), power));
    const double arg = atan2(x.imag, x.real) * power;
    complex_t result;
    result.real = norm * cos(arg);
    result.imag = norm * sin(arg);
    return result;
}
