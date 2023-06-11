-- Solve Polynomials of up to the fourth degree.
-- Algorithms by Ferrari, Tartaglia, Cardano, et al. (16th century Italy)

module Quartic (solvePoly, stableSolvePoly) where

import Data.Complex (Complex, magnitude)
import Data.List (inits)

-- | Solve polynomial of floating-point coefficients.
-- Will result in more accurate results than @solvePoly@.
stableSolvePoly :: (RealFloat a, Ord a) => [Complex a] -> [Complex a]
stableSolvePoly coefs
    | stablenessScore coefs > stablenessScore revCoefs = map (1 /) $ solvePoly revCoefs
    | otherwise = solvePoly coefs
    where
        revCoefs = reverse coefs

stablenessScore :: RealFloat a => [Complex a] -> a
stablenessScore [] = 1
stablenessScore [_] = 1
stablenessScore xs =
    t + 1 / t
    where
        (a : b : _) = reverse xs
        t = magnitude $ a / b

-- | Solve polynomial.
-- @a@ should be a type that supports @sqrt (-1)@, i.e not @Double@ or @Float@.
solvePoly :: (Eq a, Floating a) => [a] -> [a]
solvePoly [] = []
solvePoly (0:_) = [0]
solvePoly coefs
    | a == 0 = solvePoly $ init coefs
    | otherwise = solveNormalizedPoly . map (/ a) $ init coefs
    where
        a = last coefs

-- Normalized polynomials have the form of
--   x^n + a*x^(n-1) + ..
-- The coefficient for x^n is one.
-- So it gets 4 coefficients for a normalized quartic polynom.
solveNormalizedPoly :: (Eq a, Floating a) => [a] -> [a]
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
    f (mhead xs) (mhead ys) : zipWithDefault d f (drop 1 xs) (drop 1 ys)
    where
        mhead [] = d
        mhead (x:_) = x

binomials :: Num a => [[a]]
binomials =
    iterate step [1]
    where
        step prev = zipWith (+) (0 : prev) (prev ++ [0])

-- Depressed polynomials have the form of:
--   x^n + a*x^(n-2) + ..
-- The coefficient for x^n is 1 and for x^(n-1) is zero.
-- So it gets 3 coefficients for a depressed quartic polynom.
solveDepressedPoly :: (Eq a, Floating a) => [a] -> [a]
solveDepressedPoly (0 : xs) = 0 : solveDepressedPoly xs
solveDepressedPoly [] = [0] -- Poly is: x + 0 = 0
solveDepressedPoly [c0] = sqrts (-c0) -- Quadratic
solveDepressedPoly [c0, c1] = solveDepressedCubic c0 c1
solveDepressedPoly [c0, c1, c2] = solveDepressedQuartic c0 c1 c2
solveDepressedPoly _ = error "unsupported polynomial degree"

-- Based on http://en.wikipedia.org/wiki/Quartic_function#Quick_and_memorable_solution_from_first_principles
solveDepressedQuartic :: (Eq a, Floating a) => a -> a -> a -> [a]
solveDepressedQuartic e 0 c = concatMap sqrts $ solvePoly [e, c, 1]
solveDepressedQuartic e d c =
    solvePoly [c + p*p - d/p, 2*p, 2] ++ solvePoly [c + p*p + d/p, -2*p, 2]
    where
        p = sqrt . head $ solvePoly [-d*d, c*c-4*e, 2*c, 1]

sqrts :: Floating a => a -> [a]
sqrts x =
    [-s, s]
    where
        s = sqrt x

-- Based on http://en.wikipedia.org/wiki/Cubic_equation#Cardano.27s_method
--
-- Currently only provides one solution out of three.
-- Providing the other two for Complex numbers would require using cis, or some typeclass providing 'primitive roots of unity'...
solveDepressedCubic :: (Eq a, Floating a) => a -> a -> [a]
solveDepressedCubic q 0 = [(-q)**(1/3)]
solveDepressedCubic q p =
    [u - p/3/u]
    where
        u = head (solvePoly [-p*p*p/27, q, 1]) ** (1/3)
