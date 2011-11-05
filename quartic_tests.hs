import Quartic

import Data.Complex
import Test.QuickCheck (quickCheck)

toComplex :: RealFloat a => a -> Complex a
toComplex = (:+ 0)

checks :: [Complex Double] -> Bool
checks vals =
    checkPoly vals && all (elemC vals) roots &&
    -- For Cubic equations we do not find all roots. See comment in Quartic.hs
    (length vals == 3 || all (elemC roots) vals)
    where
        poly = polyFromRoots vals
        roots = solvePoly poly
        elemC ys x = any ((< 0.01) . magnitude . (x -)) ys

checkPoly :: [Complex Double] -> Bool
checkPoly poly =
    checkRoots $ solvePoly poly
    where
        checkRoots [] = all (== 0) poly
        checkRoots roots = all ((< 0.01) . magnitude . evaluatePoly poly) roots

polyFromRoots :: Num a => [a] -> [a]
polyFromRoots = foldl convolve [1] . map (: [-1])

convolve :: Num a => [a] -> [a] -> [a]
convolve [] _ = []
convolve (x:xs) ys = zipWithDefault 0 (+) (map (* x) ys) $ 0 : convolve xs ys

zipWithDefault :: a -> (a -> a -> b) -> [a] -> [a] -> [b]
zipWithDefault _ _ [] [] = []
zipWithDefault d f xs ys =
    f (mhead xs) (mhead ys) : zipWithDefault d f (drop 1 xs) (drop 1 ys)
    where
        mhead [] = d
        mhead (x:_) = x

propsRealQuadratic :: Double -> Double -> Double -> Bool
propsRealQuadratic a b c =
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
    quickCheck propsRealQuadratic
    quickCheck $ \a b -> checks [a, b]
    quickCheck $ \a b c -> checks [a, b, c]
    quickCheck $ \a b c d -> checks [a, b, c, d]
    quickCheck $ \a b c d e -> checkPoly [a, b, c, d, e]
