{-# LANGUAGE TemplateHaskell #-}

module Tests where

import qualified Hedgehog.Gen   as Gen
import qualified Hedgehog.Range as Range
import Hedgehog
import Scintillans.BlockMatrix
import Scintillans.Solver
import Scintillans

-----------------
-- Hedgehog tests
-----------------
--
-- To run them, one can use `stack ghci` that exports all Scintillans modules into GHCi session,
-- and then run `tests` function.

fromListM21 [a00, a10] = M21 a00 a10

fromListM31 [a00, a10, a20] = M31 a00 a10 a20

fromListM22 [ a00, a01
            , a10, a11 ] = M22 a00 a01
                               a10 a11

fromListM33 [ a00, a01, a02
            , a10, a11, a12
            , a20, a21, a22 ] = M33 a00 a01 a02
                                    a10 a11 a12
                                    a20 a21 a22

-- random generator producing *n* ints in [-100, 100]
nRandomNumbers n = Gen.list (Range.singleton n) $ Gen.int (Range.constant (-100) 100)

-- random generator producing *n* doubles in [0, 1)
nRandomDoubles n = Gen.list (Range.singleton n) $ Gen.double (Range.constant 0 1)

-- trace (AB) = trace (BA)
prop_M22_trace :: Property
prop_M22_trace =
  property $ do
    as <- forAll $ nRandomNumbers 4
    bs <- forAll $ nRandomNumbers 4
    trace (mult' as bs) === trace (mult' bs as)
      where mult' :: [Int] -> [Int] -> Matrix22 Int
            mult' x y = mult (fromListM22 x) (fromListM22 y)
            trace ( M22 a _
                        _ b ) = a + b

-- trace (AB) = trace (BA)
prop_M33_trace :: Property
prop_M33_trace =
  property $ do
    as <- forAll $ nRandomNumbers 9
    bs <- forAll $ nRandomNumbers 9
    trace (mult' as bs) === trace (mult' bs as)
      where mult' :: [Int] -> [Int] -> Matrix33 Int
            mult' x y = mult (fromListM33 x) (fromListM33 y)
            trace ( M33 a _ _
                        _ b _
                        _ _ c ) = a + b + c

-- Norm || U x || = || x || for unitary matrix U and any vector x.
prop_M21_unitary :: Property
prop_M21_unitary =
  property $ do
    us <- forAll $ nRandomNumbers 2
    xs <- forAll $ nRandomNumbers 2
    -- `u / sqrt nu` is an unitary matrix
    norm (mult (u us) (x xs)) === (nu us) * norm (x xs)
      where u us = let [u0, u1] = us in M22 u0    u1
                                            (-u1) u0
            nu us = let [u0, u1] = us in u0 * u0 + u1 * u1
            x xs = fromListM21 xs
            norm (M21 x0 x1) = x0 * x0 + x1 * x1

-- Norm || U x || = || x || for unitary matrix U and any vector x.
-- We use Double matrices instead of Int matrices here for the sake of simplicity.
prop_M31_unitary :: Property
prop_M31_unitary =
  property $ do
    angles <- forAll $ nRandomDoubles 3
    xs <- forAll $ nRandomDoubles 3
    -- to get sophisticated unitary matrix, chained rotations are used
    let [phi, theta, psi] = angles

        a = M33 1 0         0
                0 (cos phi) (-sin phi)
                0 (sin phi) (cos phi)

        b = M33 (cos theta)  0 (sin theta)
                0            1 0
                (-sin theta) 0 (cos theta)

        c = M33 (cos psi) (-sin psi) 0
                (sin psi) (cos psi)  0
                0         0          1

        x = fromListM31 xs
        norm (M31 x y z) = x * x + y * y + z * z
    (abs ((norm $ a `mult` b `mult` c `mult` x) - (norm x)) < 1e-15) === True

-- Zipping and converting to a list, back and forth.
prop_zipM_unzipM :: Property
prop_zipM_unzipM =
  property $ do
    let n = 10
    arr <- forAll $ nRandomNumbers n
    brr <- forAll $ nRandomNumbers n
    crr <- forAll $ nRandomNumbers n
    let xs   = M11 arr
        xs'  = M21 arr brr
        xs'' = M31 arr brr crr
    (unzipM . toList . fromList . zipM) xs   === xs
    (unzipM . toList . fromList . zipM) xs'  === xs'
    (unzipM . toList . fromList . zipM) xs'' === xs''

tests :: IO Bool
tests = checkParallel $$(discover)
