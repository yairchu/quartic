import Quartic

import Data.Complex
import Test.QuickCheck

toComplex :: RealFloat a => a -> Complex a
toComplex = (:+ 0)

checkPoly :: [Complex Double] -> Bool
checkPoly poly =
    checkRoots $ solvePoly poly
    where
        checkRoots [] = all (== 0) poly
        checkRoots roots = all ((< 0.01) . magnitude . evaluatePoly poly) roots

props_realQuadratic :: Double -> Double -> Double -> Bool
props_realQuadratic a b c =
    checkRoots $ solvePoly poly
    where
        poly = map toComplex [a, b, c]
        checkRoots [] =
            -- No roots only for a zero polynom. (or infinitly many roots..)
            poly == [0, 0, 0]
        checkRoots [root] = checkRoots [root, root]
        checkRoots [root0, root1] =
            -- Either both roots are real or they are conjugates.
            (isReal root0 && isReal root1) || root0 == conjugate root1
        checkRoots _ = error "bad resulting number of roots"
        isReal = (== 0) . imagPart

evaluatePoly :: Num a => [a] -> a -> a
evaluatePoly poly x = sum . zipWith (*) poly $ iterate (* x) 1

main :: IO ()
main = do
    quickCheck props_realQuadratic
    quickCheck $ \a b -> checkPoly [a, b]
    quickCheck $ \a b c -> checkPoly [a, b, c]
    quickCheck $ \a b c d -> checkPoly [a, b, c, d]
    quickCheck $ \a b c d e -> checkPoly [a, b, c, d, e]
