-- |This module provides helper functions to solve the following differential equation:
-- \[
-- \partial_t f = \hat A f,
-- \]
-- where \(\hat A\) is a constant matrix and \(f\) is a vector. The general
-- solution is \(f = \exp( \hat A t ) f_0\), with \(f_0\) the initial condition.

module Scintillans.Solver ( pow , exp) where

import Prelude hiding (exp)
import qualified Data.Array.Repa as R
import Data.Vector.Unboxed
import Scintillans.BlockMatrix

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
