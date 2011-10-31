/* Solve Polynomials of up to the fourth degree.
 * Algorithms by Ferrari, Tartaglia, Cardano, et al. (16th century Italy)
 */

#include <math.h>

#define MAX_DEGREE 4

int solve_poly(int degree, const double* poly, double* results);
static int solve_normalized_poly(int degree, const double* poly, double* results);
static void calc_shifted_coefs(double shift, int degree, const double* src, double* dst);
static void calc_binomials(int num_binoms, int stride, double* dst);
static void calc_powers(double x, int max_power, double* dst);
static int solve_depressed_poly(int degree, const double* poly, double* results);
static int solve_depressed_quartic(const double* poly, double* results);
static int solve_depressed_cubic(const double* poly, double* results);
static int solve_depressed_quadratic(const double* poly, double* results);

/* poly: pointer to coefficients array of size degree + 1.
 * results: pointer to results output array of size degree.
 */
int solve_poly(int degree, const double* poly, double* results) {
    const double a = poly[degree];
    if (a == 0)
        return solve_poly(degree - 1, poly, results);
    if (degree > MAX_DEGREE)
        return -1;
    double normalized_poly[MAX_DEGREE + 1];
    int i;
    for (i = 0; i < degree; ++i)
        normalized_poly[i] = poly[i] / a;
    normalized_poly[degree] = 1.0;
    return solve_normalized_poly(degree, normalized_poly, results);
}

/* Normalized polynomials have the form of
 *   x^n + a*x^(n-1) + ..
 * The coefficient for x^n is one.
 * solve_normalized_poly does expect to get this coefficient despite it being known.
 */
static int solve_normalized_poly(int degree, const double* poly, double* results) {
    const double shift = -poly[degree - 1] / degree;
    double shifted_coefs[MAX_DEGREE + 1];
    calc_shifted_coefs(shift, degree, poly, shifted_coefs);
    const int num_results = solve_depressed_poly(degree, shifted_coefs, results);
    int i;
    for (i = 0; i < num_results; ++i)
        results[i] += shift;
    return num_results;
}

static void calc_shifted_coefs(double shift, int degree, const double* src, double* dst) {
    int dst_i;
    for (dst_i = 0; dst_i <= degree; ++dst_i)
        dst[dst_i] = 0;
    double binomials[MAX_DEGREE + 1][MAX_DEGREE + 1];
    calc_binomials(degree, sizeof(binomials[0]) / sizeof(binomials[0][0]), binomials[0]);
    double shift_powers[MAX_DEGREE + 1];
    calc_powers(shift, degree, shift_powers);
    int src_i;
    for (src_i = 0; src_i <= degree; ++src_i)
        for (dst_i = 0; dst_i <= src_i; ++dst_i)
            dst[dst_i] += src[src_i] * shift_powers[src_i - dst_i] * binomials[src_i][dst_i];
}

static void calc_binomials(int num_binoms, int stride, double* dst) {
    int row;
    for (row = 0; row < num_binoms; ++row) {
        const int row_idx = row * stride;
        const int prev_row_idx = (row - 1) * stride;
        int col;
        for (col = 0; col <= row; ++col) {
            dst[row_idx + col] = dst[prev_row_idx + col] + dst[prev_row_idx + col - 1];
        }
        dst[row_idx + row] = 1;
    }
}

static void calc_powers(double x, int max_power, double* dst) {
    dst[0] = 1.0;
    if (max_power >= 1)
        dst[1] = x;
    int i;
    for (i = 2; i <= max_power; ++i)
        dst[i] = x * dst[i - 1];
}

/* Depressed polynomials have the form of:
 *   x^n + a*x^(n-2) + ..
 * The coefficient for x^n is 1 and for x^(n-1) is zero.
 * So it gets 3 coefficients for a depressed quartic polynom.
 */
static int solve_depressed_poly(int degree, const double* poly, double* results) {
    switch (degree) {
    case 4:
        return solve_depressed_quartic(poly, results);
    case 3:
        return solve_depressed_cubic(poly, results);
    case 2:
        return solve_depressed_quadratic(poly, results);
    case 1:
        results[0] = 0.0;
        return 1;
    default:
        return -1;
    }
}

/* Based on http://en.wikipedia.org/wiki/Quartic_function#Quick_and_memorable_solution_from_first_principles */
static int solve_depressed_quartic(const double* poly, double* results) {
    const double e = poly[0];
    const double d = poly[1];
    const double c = poly[2];
    double helper_cubic[4];
    helper_cubic[0] = -d*d;
    helper_cubic[1] = c*c - 4*e;
    helper_cubic[2] = 2*c;
    helper_cubic[3] = 1;
    double helper_results[3];
    if (solve_poly(3, helper_cubic, helper_results) < 1)
        return 0;
    const double p = sqrt(helper_results[0]);
    double quadratic_factor[3];
    quadratic_factor[0] = c + p*p - d/p;
    quadratic_factor[1] = 2*p;
    quadratic_factor[2] = 1;
    int num_results = solve_poly(2, quadratic_factor, results);
    quadratic_factor[0] = c + p*p + d/p;
    return num_results + solve_poly(2, quadratic_factor, results + num_results);
}

/* Based on http://en.wikipedia.org/wiki/Cubic_equation#Cardano.27s_method
 *
 * Only provides one real result out of 3 possibly complex results.
 */
static int solve_depressed_cubic(const double* poly, double* results) {
    const double q = poly[0];
    const double p = poly[1];
    const double t = q*q/4 + p*p*p/27;
    if (t >= 0) {
        const double u = pow(-q/2 - sqrt(t), 1.0 / 3.0);
        results[0] = u - p/3.0/u;
        return 1;
    }
    const double s_real = -q/2;
    const double s_abs = sqrt(s_real*s_real - t);
    const double s_phase = atan(sqrt(-t) / s_real);
    const double u_abs = pow(s_abs, 1.0 / 3.0);
    const double u_phase = s_phase / 3.0;
    const double u_real = u_abs * cos(u_phase);
    results[0] = 2 * u_real;
    return 1;
}

static int solve_depressed_quadratic(const double* poly, double* results) {
    const double t = sqrt(-poly[0]);
    results[0] = -t;
    results[1] = t;
    return 2;
}
