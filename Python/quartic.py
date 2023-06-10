# Solve Polynomials of up to the fourth degree.
# Algorithms by Ferrari, Tartaglia, Cardano, et al. (16th century Italy)

import math

def stableSolvePoly(coefs):
    """Solve polynomial of floating-point coefficients.
    Will result in more accurate results than "solvePoly".
    """
    revCoefs = coefs[::-1]
    if stablenessScore(coefs) > stablenessScore(revCoefs):
        return [1 / x for x in solvePoly(revCoefs)]
    return solvePoly(coefs)

def stablenessScore(coefs):
    if len(coefs) < 2:
        return 1
    t = abs(coefs[-1] / coefs[-2])
    return t + 1 / t

def solvePoly(coefs):
    "Solve polynomial."
    if not coefs:
        return []
    if coefs[0] == 0:
        return [0]
    if coefs[-1] == 0:
        return solvePoly(coefs[:-1])
    return solveNormalizedPoly([x / coefs[-1] for x in coefs[:-1]])

def solveNormalizedPoly(coefs):
    """Normalized polynomials have the form of
        x^n + a*x^(n-1) + ..
    The coefficient for x^n is one.
    So it gets 4 coefficients for a normalized quartic polynomial.
    """
    degree = len(coefs)
    shift = -coefs[-1] / degree
    return [x + shift for x in solveDepressedPoly(shiftedCoefs(shift, coefs + [1])[:degree-1])]

def shiftedCoefs(shift, coefs):
    result = [0]*len(coefs)
    for coef, bin in zip(coefs, binomials()):
        x = 1
        for i, b in enumerate(bin):
            result[len(bin)-1-i] += coef * b * x
            x *= shift
    return result

def binomials():
    cur = [1]
    while True:
        yield cur
        cur = [x + y for x, y in zip([0]+cur, cur+[0])]

def solveDepressedPoly(coefs):
    """Depressed polynomials have the form of:
        x^n + a*x^(n-2) + ..
    The coefficient for x^n is 1 and for x^(n-1) is zero.
    So it gets 3 coefficients for a depressed quartic polynom.
    """
    if not coefs:
        # Poly is: x + 0 = 0
        return [0]
    if coefs[0] == 0:
        return solveDepressedPoly(coefs[1:])
    if len(coefs) == 1:
        # Quadratic
        return sqrts(-coefs[0])
    if len(coefs) == 2:
        return solveDepressedCubic(coefs[0], coefs[1])
    if len(coefs) == 3:
        return solveDepressedQuartic(coefs[0], coefs[1], coefs[2])
    raise ValueError("unsupported polynomial degree")

# Based on http://en.wikipedia.org/wiki/Quartic_function#Quick_and_memorable_solution_from_first_principles
def solveDepressedQuartic(e, d, c):
    if d == 0:
        return [s for x in solvePoly([e, c, 1]) for s in sqrts(x)]
    p = solvePoly([-d*d, c*c-4*e, 2*c, 1])[0]**0.5
    return solvePoly([c + p*p - d/p, 2*p, 2]) + solvePoly([c + p*p + d/p, -2*p, 2])

def sqrts(x):
    s = x**0.5
    return [-s, s]

# Based on http://en.wikipedia.org/wiki/Cubic_equation#Cardano.27s_method
def solveDepressedCubic(q, p):
    thirdRootUnity = math.e**(math.pi/3j)
    if p == 0:
        r = -q ** (1/3.)
    else:
        u = solvePoly([-p*p*p/27, q, 1])[0]**(1/3.)
        r = u - p/3/u
    return [r, r*thirdRootUnity, r*thirdRootUnity**2]
