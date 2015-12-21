/* Solve Polynomials of up to the fourth degree. (over real numbers)
 * Algorithms by Ferrari, Tartaglia, Cardano, et al. (16th century Italy)
 */

#include <math.h>

#define MAX_DEGREE 4

static double stableness_score(double a, double b);
static int solve_normalized_poly(int degree, const double* poly, double* results);
static void calc_shifted_coefs(double shift, int degree, const double* src, double* dst);
static void calc_binomials(int num_binoms, int stride, double* dst);
static void calc_powers(double x, int max_power, double* dst);
static int solve_depressed_poly(int degree, const double* poly, double* results);
static int solve_depressed_quartic(const double* poly, double* results);
static int solve_depressed_cubic(const double* poly, double* results);
static double cubic_root(double x);
static int solve_depressed_quadratic(const double* poly, double* results);

/* poly: pointer to coefficients array of size degree + 1.
 * results: pointer to results output array of size degree.
 */
int solve_real_poly(int degree, const double* poly, double* results) {
    double normalized_poly[MAX_DEGREE + 1];
    int i;
    const double a = poly[degree];
    if (a == 0)
        return solve_real_poly(degree - 1, poly, results);
    if (degree > MAX_DEGREE)
        return -1;
    if (degree > 2 && stableness_score(poly[degree], poly[degree - 1]) > stableness_score(poly[0], poly[1])) {
        double rev_poly[MAX_DEGREE + 1];
        int i, num_results;
        for (i = 0; i <= degree; ++i)
            rev_poly[i] = poly[degree - i];
        num_results = solve_real_poly(degree, rev_poly, results);
        for (i = 0; i < num_results; ++i)
            results[i] = 1.0 / results[i];
        return num_results;
    }
    for (i = 0; i < degree; ++i)
        normalized_poly[i] = poly[i] / a;
    normalized_poly[degree] = 1.0;
    return solve_normalized_poly(degree, normalized_poly, results);
}

static double stableness_score(double a, double b) {
    const double t = fabs(a / b);
    return t + 1.0 / t;
}

/* Normalized polynomials have the form of
 *   x^n + a*x^(n-1) + ..
 * The coefficient for x^n is one.
 * solve_normalized_poly does expect to get this coefficient despite it being known.
 */
static int solve_normalized_poly(int degree, const double* poly, double* results) {
    const double shift = -poly[degree - 1] / (double) degree;
    double shifted_coefs[MAX_DEGREE + 1];
	int i, num_results;
    calc_shifted_coefs(shift, degree, poly, shifted_coefs);
    num_results = solve_depressed_poly(degree, shifted_coefs, results);
    for (i = 0; i < num_results; ++i)
        results[i] += shift;
    return num_results;
}

static void calc_shifted_coefs(double shift, int degree, const double* src, double* dst) {
    double binomials[MAX_DEGREE + 1][MAX_DEGREE + 1];
    double shift_powers[MAX_DEGREE + 1];
    int dst_i, src_i;
    for (dst_i = 0; dst_i <= degree; ++dst_i)
        dst[dst_i] = 0;
    calc_binomials(degree+1, sizeof(binomials[0]) / sizeof(binomials[0][0]), binomials[0]);
    calc_powers(shift, degree, shift_powers);
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
        dst[row_idx] = 1;
        for (col = 1; col < row; ++col) {
            dst[row_idx + col] = dst[prev_row_idx + col - 1] + dst[prev_row_idx + col];
        }
        dst[row_idx + row] = 1;
    }
}

static void calc_powers(double x, int max_power, double* dst) {
    int i;
    dst[0] = 1.0;
    if (max_power >= 1)
        dst[1] = x;
    for (i = 2; i <= max_power; ++i)
        dst[i] = x * dst[i - 1];
}

/* Depressed polynomials have the form of:
 *   x^n + a*x^(n-2) + ..
 * The coefficient for x^n is 1 and for x^(n-1) is zero.
 * So it gets 3 coefficients for a depressed quartic polynom.
 */
static int solve_depressed_poly(int degree, const double* poly, double* results) {
    if (degree > 0 && poly[0] == 0.0) {
        results[0] = 0;
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
        results[0] = 0.0;
        return 1;
    case 0:
        return 0;
    default:
        return -1;
    }
}

/* Based on http://en.wikipedia.org/wiki/Quartic_function#Quick_and_memorable_solution_from_first_principles */
static int solve_depressed_quartic(const double* poly, double* results) {
    double helper_cubic[4];
    double helper_results[3];
    double quadratic_factor[3];
	double p;
    const double e = poly[0];
    const double d = poly[1];
    const double c = poly[2];
	int num_results;
    if (d == 0) {
		int i, num_quad_results;
        double quadratic[3];
        double quadratic_results[2];
        quadratic[0] = e;
        quadratic[1] = c;
        quadratic[2] = 1;
        num_quad_results = solve_real_poly(2, quadratic, quadratic_results);
        for (i = 0; i < num_quad_results; ++i) {
            const double s = sqrt(quadratic_results[i]);
            results[2*i] = -s;
            results[2*i + 1] = s;
        }
        return 2 * num_quad_results;
    }
    helper_cubic[0] = -d*d;
    helper_cubic[1] = c*c - 4*e;
    helper_cubic[2] = 2*c;
    helper_cubic[3] = 1;
    if (solve_real_poly(3, helper_cubic, helper_results) < 1)
        return 0;
    p = sqrt(helper_results[0]);
    quadratic_factor[0] = c + p*p - d/p;
    quadratic_factor[1] = 2*p;
    quadratic_factor[2] = 2;
    num_results = solve_real_poly(2, quadratic_factor, results);
    quadratic_factor[0] = c + p*p + d/p;
    quadratic_factor[1] = -2*p;
    return num_results + solve_real_poly(2, quadratic_factor, results + num_results);
}

/* Based on http://en.wikipedia.org/wiki/Cubic_equation#Cardano.27s_method
 *
 * Only provides one real result out of 3 possibly complex results.
 */
static int solve_depressed_cubic(const double* poly, double* results) {
    const double q = poly[0];
    const double p = poly[1];
	double t, s_real, s_abs, s_phase, u_abs, u_phase, u_real;
    if (p == 0.0) {
        results[0] = cubic_root(-q);
        return 1;
    }
    t = q*q/4 + p*p*p/27;
    if (t >= 0.0) {
        const double u = cubic_root(-q/2 + sqrt(t));
        results[0] = u - p/3.0/u;
        return 1;
    }
    s_real = -q/2;
    s_abs = sqrt(s_real*s_real - t);
    s_phase = atan(sqrt(-t) / s_real) + (s_real >= 0 ? 0 : M_PI);
    u_abs = cubic_root(s_abs);
    u_phase = s_phase / 3.0;
    u_real = u_abs * cos(u_phase);
    results[0] = 2 * u_real;
    return 1;
}

static double cubic_root(double x) {
    const double t = pow(fabs(x), 1.0 / 3.0);
    return x >= 0.0 ? t : -t;
}

static int solve_depressed_quadratic(const double* poly, double* results) {
    const double t = sqrt(-poly[0]);
    results[0] = -t;
    results[1] = t;
    return 2;
}
