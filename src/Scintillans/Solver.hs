module Scintillans.Solver where

import qualified Data.Array.Repa as R
import Data.Vector.Unboxed
import Scintillans.BlockMatrix

{--
This module provides helper functions to solve the following differential equation (LaTex notation
is used):

\partial_t f = \hat A f,

where $\hat A$ is a matrix anf $f$ is a vector. The general solution is f = exp( A t ) f_0, where
f_0 is the initial condition. Here, the equation is solved numerically with finite-difference
Euler's method, i.e. we compute the exponent of A t numerically.
--}

-- m^n, with exponentiation by squaring
pow :: (Num a, Unbox a, Multable a a a) =>
       (R.Array R.U R.DIM2 a) ->
       Int ->
       (R.Array R.U R.DIM2 a)
pow m n
    | n == 1         = m
    | n `mod` 2 == 1 = mmultS m (pow m $ n - 1)
    | otherwise      = let m' = pow m (n `div` 2) in mmultS m' m'

-- exp( hatA t ), computed approximately
exp :: (Functor f, Num a, Num (f a), Unbox (f a), Multable (f a) (f a) (f a)) => 
       (R.Array R.U R.DIM2 (f a)) -> -- hatA
       a ->                          -- dt
       Int ->                        -- n
       (R.Array R.U R.DIM2 (f a))
exp hatA dt n = pow hatA' n
  where hatA' = R.computeS $ (idMatrix h) R.+^ hatAdt
        R.Z R.:. h R.:. _ = R.extent hatA
        hatAdt = R.map (fmap ((*) dt)) hatA
