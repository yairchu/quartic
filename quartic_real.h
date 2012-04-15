/* Solve Polynomials of up to the fourth degree. (over real numbers)
 * Algorithms by Ferrari, Tartaglia, Cardano, et al. (16th century Italy)
 */

#ifndef __QUARTIC_H__
#define __QUARTIC_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

int solve_real_poly(int degree, const double* poly, double* results);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __QUARTIC_H__ */