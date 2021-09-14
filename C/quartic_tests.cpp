#include <rapidcheck.h>

#include <complex>
#include <vector>
#include <algorithm>

#include "quartic.h"

using std::complex;
using std::vector;

// Missing support for std::complex in rapidcheck
// (https://github.com/emil-e/rapidcheck/issues/11)
template <typename T>
struct rc::Arbitrary<complex<T>> {
  static Gen<complex<T>> arbitrary() {
    return gen::map(
      gen::pair(
        Arbitrary<T>::arbitrary(),
        Arbitrary<T>::arbitrary()),
      [](std::pair<T, T> pair) { return complex<T> (pair.first, pair.second); }
      );
  }
};

template<typename T>
T evaluatePoly(const vector<T>& poly, T x)
{
    T y = 0, xE = 1;
    for (T c : poly)
    {
        y += c * xE;
        xE *= x;
    }
    return y;
}

void checkPoly(const vector<complex<double>>& poly)
{
    const int degree = poly.size()-1;
    complex<double> results[degree];
    const int numResults = solve_poly(degree, (const complex_t*) poly.data(), (complex_t*) results);
    if (numResults == 0)
    {
        for (int i = 1; i < poly.size(); ++i)
            RC_ASSERT(poly[i] == 0.0);
        return;
    }
    for (int i = 0; i < numResults; ++i)
        RC_ASSERT(std::abs(evaluatePoly(poly, results[i])) < 0.01);
}

template <typename T>
vector<T> convolve(vector<T> x, vector<T> y)
{
    vector<T> result(x.size()+y.size()-1, 0);
    for (int iX = 0; iX < x.size(); ++iX)
        for (int iY = 0; iY < y.size(); ++iY)
            result[iX+iY] += x[iX] * y[iY];
    return result;
}

template <typename T>
vector<T> polyFromRoots(vector<T> roots)
{
    vector<T> result = { 1 };
    for (T root : roots)
        result = convolve(result, {root, -1});
    return result;
}

void checks(const vector<complex<double>>& vals)
{
    // Avoid too large values due to accuracy issues
    for (const auto& x : vals)
        if (std::norm(x) > 1e20)
            return;

    checkPoly(vals);

    const vector<complex<double>> poly = polyFromRoots(vals);
    const int degree = poly.size()-1;
    complex<double> roots[degree];
    const int numRoots = solve_poly(degree, (const complex_t*) poly.data(), (complex_t*) roots);
    RC_ASSERT(numRoots >= 0);
    const double epsilon = 1e-5;
    for (int iR = 0; iR < numRoots; ++iR)
    {
        const complex<double> root = roots[iR];
        bool found = false;
        for (const auto& v : vals)
            if (std::norm(v - root) < epsilon)
            {
                found = true;
                break;
            }
        RC_ASSERT(found);
    }
    for (const auto& v : vals)
    {
        if (v == 0.0)
            continue;
        bool found = false;
        for (int iR = 0; iR < numRoots; ++iR)
            if (std::norm(v - roots[iR]) < epsilon)
            {
                found = true;
                break;
            }
        RC_ASSERT(found);
    }
}

void propsRealQuadratic(double a, double b, double c)
{
    const vector<complex<double>> poly = {a, b, c};
    complex<double> roots[2];
    const int numRoots = solve_poly(2, (const complex_t*) poly.data(), (complex_t*) roots);
    if (numRoots == 0)
    {
        RC_ASSERT(a == 0.0 && b == 0.0 && c == 0.0);
        return;
    }
    if (numRoots == 1)
    {
        RC_ASSERT(roots[0].imag() == 0.0);
        return;
    }
    RC_ASSERT(numRoots == 2);
    RC_ASSERT(
        (roots[0].imag() == 0.0 && roots[1].imag() == 0.0)
        ||
        (roots[0].real() == roots[1].real() && roots[0].imag() == -roots[1].imag())
        );
}

int main()
{
    typedef complex<double> C;
    rc::check("Real quadratic", propsRealQuadratic);
    rc::check("Degree=1", [](C a) { checks( {a} ); });
    rc::check("Degree=2", [](C a, C b) { checks( {a, b} ); });
    rc::check("Degree=3", [](C a, C b, C c) { checks( {a, b, c} ); });
    rc::check("Degree=4", [](C a, C b, C c, C d) { checks( {a, b, c, d} ); });
    return 0;
}
