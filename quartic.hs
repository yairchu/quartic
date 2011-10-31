-- Solve Polynomials of up to the fourth degree.
-- Algorithms by Ferrari, Tartaglia, Cardano, et al. (16th century Italy)

module Quartic (solvePoly) where

import Data.List (inits)

-- 'a' should be a type that supports 'sqrt (-1)'.
solvePoly :: Floating a => [a] -> [a]
solvePoly coefs
    | a == 0 = solvePoly $ init coefs
    | stablenessScore coefs > stablenessScore revCoefs = map (1 /) $ solvePoly revCoefs
    | otherwise = solveNormalizedPoly . map (/ a) $ init coefs
    where
        revCoefs = reverse coefs
        a = head revCoefs

stablenessScore poly :: Fractional a => [a] -> [a]
stablenessScore [] = 1
stablenessScore [_] = 1
stablenessScore xs =
    t + 1 / t
    where
        (a : b : _) = reverse xs
        t = abs $ a / b

-- Normalized polynomials have the form of
--   x^n + a*x^(n-1) + ..
-- The coefficient for x^n is one.
-- So it gets 4 coefficients for a normalized quartic polynom.
solveNormalizedPoly :: Floating a => [a] -> [a]
solveNormalizedPoly coefs =
    map (+ shift) . solveDepressedPoly . take (degree - 1) . shiftedCoefs shift $ coefs ++ [1]
    where
        shift = -last coefs / fromIntegral degree
        degree = length coefs

shiftedCoefs :: Num a => a -> [a] -> [a]
shiftedCoefs shift coefs =
    foldl1 (zipWithDefault 0 (+)) .
    zipWith (map . (*)) coefs .
    zipWith (zipWith (*)) binomials .
    map reverse . tail . inits $ iterate (* shift) 1

zipWithDefault :: a -> (a -> a -> b) -> [a] -> [a] -> [b]
zipWithDefault _ _ [] [] = []
zipWithDefault d f xs ys =
    f (mhead xs) (mhead ys) : zipWithDefault d f (mtail xs) (mtail ys)
    where
        mhead [] = d
        mhead (x:_) = x
        mtail [] = []
        mtail (_:rest) = rest

binomials :: Num a => [[a]]
binomials =
    iterate step [1]
    where
        step prev = zipWith (+) (0 : prev) (prev ++ [0])

-- Depressed polynomials have the form of:
--   x^n + a*x^(n-2) + ..
-- The coefficient for x^n is 1 and for x^(n-1) is zero.
-- So it gets 3 coefficients for a depressed quartic polynom.
solveDepressedPoly :: Floating a => [a] -> [a]
solveDepressedPoly coefs
    | null coefs = []
    | head coefs == 0 = 0 : solveDepressedPoly (tail coefs)
    | degree == 4 = solveDepressedQuartic coefs
    | degree == 3 = solveDepressedCubic coefs
    | degree == 2 = solveDepressedQuadratic coefs
    | degree == 1 = [0]
    | otherwise = error "unsupported polynomial degree"
    where
        degree = length coefs + 1

-- Based on http://en.wikipedia.org/wiki/Quartic_function#Quick_and_memorable_solution_from_first_principles
solveDepressedQuartic :: Floating a => [a] -> [a]
solveDepressedQuartic coefs
    | d == 0 = concatMap sqrts $ solvePoly [e, c, 1]
    | otherwise = solvePoly [c + p*p - d/p, 2*p, 2] ++ solvePoly [c + p*p + d/p, -2*p, 2]
    where
        p = sqrt . head $ solvePoly [-d*d, c*c-4*e, 2*c, 1]
        [e, d, c] = coefs

sqrts :: Floating a => a -> [a]
sqrts x =
    [-s, s]
    where
        s = sqrt x

-- Based on http://en.wikipedia.org/wiki/Cubic_equation#Cardano.27s_method
--
-- Currently only provides one solution out of three.
-- Providing the other two for Complex numbers would require using cis, or some typeclass providing 'primitive roots of unity'...
solveDepressedCubic :: Floating a => [a] -> [a]
solveDepressedCubic coefs
    | p == 0 = [cubicRoot (-q)]
    | otherwise = [u - p/3/u]
    where
        [q, p] = coefs
        u = cubicRoot . head $ solvePoly [-p*p*p/27, q, 1]

-- (** (1/3)) doesn't work for Doubles.
cubicRoot :: Floating a => a -> a
cubicRoot x
    | signum x == -1 = -(-x)**(1/3)
    | otherwise = x ** (1/3)

solveDepressedQuadratic :: Floating a => [a] -> [a]
solveDepressedQuadratic = sqrts . negate . head
