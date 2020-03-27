{-# LANGUAGE FunctionalDependencies #-}

-- |This module provides helper functions to solve the following differential equation:
-- \[
-- \partial_t f = \hat A f,
-- \]
-- where \(\hat A\) is a constant matrix and \(f\) is a vector. The general
-- solution is \(f = \exp( \hat A t ) f_0\), with \(f_0\) the initial condition.

module Scintillans.Solver (
    Multable (mult)
  , mmultS
  , idMatrix
  , pow
  , exp) where

import Prelude hiding (exp)
import qualified Data.Array.Repa as R
import Data.Vector.Unboxed

-- |Class to express block matrices which can be multiplied. In general, \(\hat A\) and \(f\) above
-- are block matrices (a matrices of (sub-)matrices, see "Scintillans.BlockMatrix") and here we
-- prevent multiplication of matrices with inappropriate number of blocks. For instance,
-- multiplication of a matrix of 2x3 blocks with a matrix of 1x1 block wouldn't compile because
-- there is no instance of @Multable@ for them.
class Multable a b c | a b -> c where
  mult :: a -> b -> c

-- |Multiplication of Repa matrices (or block matrices represented as Repa matrices), as in the
-- linear algebra. The elements of the resulting matrix are computed sequentially.
mmultS :: (Unbox a, Unbox b, Unbox c, Num c, Multable a b c)
       => R.Array R.U R.DIM2 a -- ^\(\hat A\)
       -> R.Array R.U R.DIM2 b -- ^\(\hat B\)
       -> R.Array R.U R.DIM2 c -- ^\(\hat A \hat B\)
mmultS arr brr
  | wa /= hb  =
      error "Scintillans.Solver mmultS: width of the first matrix != height of the second one."
  | otherwise =
      arr `R.deepSeqArray`
      brr `R.deepSeqArray`
      ( R.computeS
      $ R.fromFunction (R.Z R.:. ha R.:. wb)
      $ \(R.Z R.:. i R.:. j) -> R.sumAllS
            ( R.zipWith mult
                        (R.slice arr (R.Any R.:. i R.:. R.All))
                        (R.slice brr (R.Any R.:. j))
            )
      )
    where (R.Z R.:. ha R.:. wa) = R.extent arr
          (R.Z R.:. hb R.:. wb) = R.extent brr

-- |@idMatrix n@ is the @nxn@ identity matrix.
idMatrix :: (Num a, Unbox a) => Int -> R.Array R.U R.DIM2 a
idMatrix n = R.computeS $ R.fromFunction (R.Z R.:. n R.:. n) $ \(R.Z R.:. i R.:. j) ->
  if i == j then fromInteger 1 else fromInteger 0

-- |\(m^n\), computed with exponentiation by squaring.
pow :: (Num a, Unbox a, Multable a a a) =>
       (R.Array R.U R.DIM2 a) -- ^m
    -> Int                    -- ^n \(\geqslant\) 0
    -> (R.Array R.U R.DIM2 a)
pow m n
    | n <  0         = error "Scintillans.Solver pow: call of @pow m n@ for n < 0."
    | n == 0         = idMatrix d
    | n == 1         = m
    | n `mod` 2 == 1 = mmultS m (pow m $ n - 1)
    | otherwise      = mmultS m' m'
      where m' = pow m (n `div` 2)
            (R.Z R.:. d R.:. _) = R.extent m

-- |\(\exp(\hat A t) \), computed approximately with finite-difference Euler's method, i.e. as
-- \[
-- \left( \hat 1 + \hat A \Delta t \right)^{[t / \Delta t]}.
-- \]
exp :: (Functor f, Num a, Num (f a), Unbox (f a), Multable (f a) (f a) (f a)) => 
       (R.Array R.U R.DIM2 (f a)) -- ^\(\hat A\)
    -> a                          -- ^dt
    -> Int                        -- ^\ [t / dt\ ]
    -> (R.Array R.U R.DIM2 (f a))
exp hatA dt n = pow hatA' n
  where hatA' = R.computeS $ (idMatrix h) R.+^ hatAdt
        R.Z R.:. h R.:. _ = R.extent hatA
        hatAdt = R.map (fmap ((*) dt)) hatA
