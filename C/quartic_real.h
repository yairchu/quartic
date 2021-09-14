/* Solve Polynomials of up to the fourth degree. (over real numbers)
 * Algorithms by Ferrari, Tartaglia, Cardano, et al. (16th century Italy)
 */

#ifndef QUARTIC_REAL_H
#define QUARTIC_REAL_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

int solve_real_poly(int degree, const double* poly, double* results);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* QUARTIC_REAL_H */
