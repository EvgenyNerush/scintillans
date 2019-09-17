{-# LANGUAGE KindSignatures #-}

-- |This module provides a simple interface serving for simulation of synchrotron emission and pair
-- photoproduction in strong magnetic field. See "Scintillans.Synchrotron" for formulas and
-- normalization used here.
--
-- To find particle distribution \(f\) one should first prepare initial particle distribution,
-- \(f_0\).  Second, he/she should prepare matrix \(\hat S\) and apply it to \(f_0\),
-- \(f = \hat S f_0\). Then, one can extract electron, photon and positron spectra with pattern
-- matching. All this chain can be expressed as follows:
--
-- > M21 fe fph = unzipM . toList
-- >            $ mmultS (hatS b xa xb nx dt nt :: Block Matrix22)
-- >            $ fromList . zipM
-- >            $ f0
--
-- Here only photon emission is taken into account; @fe@ and @fph@ are the found distributions of
-- the electrons and the photons, respectively, @b@ is the magnetic field strength, @dt@ is the
-- time step and @nt@ is the number of time steps used to compute the solution. Note the type
-- specification @:: Block Matrix22@ that explicitly indicates two-component version of the matrix
-- \(\hat S\). For details on block matrix representation of \(f\) and \(\hat S\), see
-- "Scintillans.BlockMatrix".
--
-- The distribution functions are lists of values at @nx@ equidistant grid nodes placed from the
-- bottom energy boundary @xa@ to the top energy boundary @xb@ (both inclusive). @f0@ is the
-- initial particle distribution, e.g. a monoenergetic beam of electrons with energy of
-- @xb@ accompanied with no photons:
-- 
-- > f0 = M22 (diracDelta nx $ deltaX xa xb nx) zeros
--
-- See 'fromList'' for one more example on how to set the initial particle distribution.

module Scintillans
  ( deltaX
  , zeros
  , diracDelta
  , Zippable (zipM)
  , Unzippable (unzipM)
  , fromList'
  , fromList
  , toList
  , Block
  , S (hatS)
  , mmultS
  ) where

import qualified Data.Array.Repa   as R
import Scintillans.Solver          as S
import Data.Vector.Unboxed (Unbox)
import Scintillans.BlockMatrix
import Scintillans.Synchrotron

-- |Step of the energy grid, @deltaX xa xb nx = (xb - xa) / (nx - 1)@.
deltaX :: Double -> Double -> Int -> Double
deltaX xa xb nx = (xb - xa) / (fromIntegral $ nx - 1)

-- |@[0,0..]@
zeros :: Num a => [a]
zeros = 0:zeros

-- |Numerical representation of the Dirac delta function containing exactly one particle, at the
-- right energy boundary.
diracDelta :: Int    -- ^ number of the nodes
           -> Double -- ^ interval between the nodes
           -> [Double]
diracDelta nx dx = take (nx - 1) zeros ++ [1 / dx]

class Zippable (t :: * -> *) where
  -- |Zips a matrix of lists, e.g. @zipM (M21 [1, 2, 3] [4, 5, 6]) = [(M21 1 4), (M21 2 5), (M21 3
  -- 6)]@.
  zipM :: t [a] -> [t a]

instance Zippable Matrix11 where
  zipM (M11 xs) = map (\x -> M11 x) xs

instance Zippable Matrix21 where
  zipM (M21 [] _) = []
  zipM (M21 _ []) = []
  zipM (M21 (x:xs) (y:ys)) = (M21 x y):(zipM $ M21 xs ys)

instance Zippable Matrix31 where
  zipM (M31 [] _ _) = []
  zipM (M31 _ [] _) = []
  zipM (M31 _ _ []) = []
  zipM (M31 (x:xs) (y:ys) (z:zs)) = (M31 x y z):(zipM $ M31 xs ys zs)

class Unzippable (t :: * -> *) where
  -- |Unzips a list of matrices, e.g. @unzipM [(M21 1 2), (M21 3 4), (M21 5 6)] = M21 [1, 3, 5] [2,
  -- 4, 6]@.
  unzipM :: [t a] -> t [a]

instance Unzippable Matrix11 where
  unzipM = M11
         . map (\(M11 x) -> x)

instance Unzippable Matrix21 where
  unzipM [] = M21 [] []
  unzipM ((M21 x y):ms) = M21 (x:xs) (y:ys)
    where M21 xs ys = unzipM ms

instance Unzippable Matrix31 where
  unzipM [] = M31 [] [] []
  unzipM ((M31 x y z):ms) = M31 (x:xs) (y:ys) (z:zs)
    where M31 xs ys zs = unzipM ms

-- |Constructs \(f_0\) from an (infinite) list of blocks. E.g., for two-component distribution function (see
-- "Scintillans.Synchrotron#M2c") a rising electron distribution and an oscillating photon
-- distribution can be obtained with
--
-- > fromList' nx [M21 x (sin x) | x <- [xa, (xa + dx)..]]
fromList' :: (Unbox a) =>
            Int -- ^number of nodes in numerical approximation of \(f_e\)
         -> [a] -- ^values of \(f\) at nodes
         -> R.Array R.U R.DIM2 a
fromList' nx fvals = R.fromListUnboxed (R.Z R.:. nx R.:. (1 :: Int)) fvals

-- |Constructs \(f_0\) from a finite list of blocks.
fromList :: (Unbox a) =>
            [a] -- ^values of \(f\) at nodes
         -> R.Array R.U R.DIM2 a
fromList fvals = fromList' nx fvals
  where nx = length fvals

-- |Converts \(f\) from one representation to another.
toList :: (Unbox a) =>
          R.Array R.U R.DIM2 a -- ^\(f\) as a block matrix
       -> [a]                  -- ^\(f\) as a list
toList = R.toList

-- |Alias for block matrix type.
type Block (m :: * -> *) = R.Array R.U R.DIM2 (m Double)

class S m where
  -- |\(\hat S\), a matrix which propagates distribution functions in time from 0 to @t = nt * dt@.
  hatS :: Double -- ^magnetic field strength
       -> Double -- ^lower energy boundary
       -> Double -- ^upper energy boundary
       -> Int    -- ^dimension of the vector space which represents \(f_e\)
       -> Double -- ^time step, dt
       -> Int    -- ^number of timesteps, nt
       -> Block m

-- |Instance for single-component distribution function, see "Scintillans.Synchrotron#M1c".
instance S Matrix11 where
  hatS b xa xb nx dt nt = S.exp hatA dt nt
    where hatA = R.computeS $ R.map (\x -> M11 x) $ hatA00 b xa xb nx

-- |Instance for two-component distribution function, see "Scintillans.Synchrotron#M2c".
instance S Matrix22 where
  hatS b xa xb nx dt nt = S.exp hatA dt nt
    where hatA = R.computeS $ R.traverse2
            (hatA00 b xa xb nx)
            (hatA10 b xa xb nx)
            (\_ _ -> (R.Z R.:. nx R.:. nx))
            (\lf lf' (R.Z R.:. i R.:. j) -> M22 (lf  (R.Z R.:. i R.:. j))
                                                0
                                                (lf' (R.Z R.:. i R.:. j))
                                                0
            )
