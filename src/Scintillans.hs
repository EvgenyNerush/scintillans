{-# LANGUAGE KindSignatures #-}

-- |This module provides a simple interface serving for simulation of synchrotron emission and pair
-- photoproduction in strong magnetic field. See "Scintillans.Synchrotron" for formulas and
-- normalization used here.
--
-- To find particle distribution \(f\) one should first prepare initial particle distribution,
-- \(f_0\).  Second, he/she should prepare matrix \(\hat S\) with @hatS@ function and apply it to
-- \(f_0\), \(f = \hat S f_0\). Then, one can extract electron, photon and positron spectra with
-- pattern matching. All this chain can be expressed as follows:
--
-- > M21 fe fph = unzipM . toList
-- >            $ mmultS (hatS model b xa xb nx dt nt :: Block Matrix22)
-- >            $ fromList . zipM
-- >            $ f0
--
-- Here only photon emission is taken into account; @fe@ and @fph@ are the found distributions of
-- the electrons and the photons, respectively, @b@ is the magnetic field strength, @dt@ is the
-- time step and @nt@ is the number of time steps used to compute the solution. Note the type
-- specification @:: Block Matrix22@ that explicitly indicates two-component (for electrons and
-- photons) version of the matrix \(\hat S\). Functions 'zipM', 'fromList', 'toList' and 'unzipM'
-- are needed to convert distribution functions to (or from) internal /Scintillans/ representation,
-- which uses Repa arrays for the sake of performance (see "Scintillans.BlockMatrix" for details).
-- Before (after) the conversion to (from) the internal representations, block matrices are just
-- matrices that consist of blocks. For instance,
--
-- > M22 a00 a01 a10 a11 = | a00 a01 |
-- >                       | a10 a11 |,
--
-- whenewer @a00@ etc. are, or, for example @M12 [1, 2] [3, 4]@ is the matrix @| [1, 2] [3, 4] |@.
--
-- The distribution function @fe@ here is the list of values at @nx@ equidistant grid nodes placed
-- from the bottom energy boundary @xa@ to the top energy boundary @xb@ (both inclusive). For the
-- photons, @fph@ starts at @0@ and ends at @(xb - xa)@ (both inclusive) in the two-component
-- representation and boundaries of @fph@ are the same as for @fe@ in the three-component
-- representation. See "Scintillans.Synchrotron#M1c" for details.
--
-- @f0@ is the initial particle distribution, e.g. a monoenergetic beam of electrons with energy of
-- @xb@ accompanied with zero photon distribution function:
-- 
-- > f0 = M21 (diracDelta nx $ deltaX xa xb nx) zeros
--
-- See 'fromList'' for one more example on how to set the initial particle distribution.

module Scintillans
  ( deltaX
  , energyGrid
  , zeros
  , diracDelta
  , Zippable (zipM)
  , Unzippable (unzipM)
  , fromList'
  , fromList
  , toList
  , Model(..)
  , Block
  , S(..)
  , mmultS
  ) where

import qualified Data.Array.Repa   as R
import Scintillans.Solver          as Sol
import Data.Vector.Unboxed (Unbox)
import Scintillans.BlockMatrix
import Scintillans.Synchrotron
import Scintillans.SynchrotronClassical

-- |Step of the energy grid, @deltaX xa xb nx = (xb - xa) / (nx - 1)@.
deltaX :: Double -> Double -> Int -> Double
deltaX xa xb nx = (xb - xa) / (fromIntegral $ nx - 1)

-- |Electron energy grid, i.e. nodes where values of the spectrum \(f_e\) are given
energyGrid :: Double -> Double -> Int -> [Double]
energyGrid xa xb nx = [xa + dx * fromIntegral i | i <- [0..(nx -1)]]
  where dx = deltaX xa xb nx

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

-- |Type pointing which model to use, e.g. the Baier-Katkov-Strakhovenko quantum model (@BKS@), the
-- classical one where radiation recoil is neglected (@Classical@) or model based on
-- Landau-Lifshitz radiation force (@LL@) which neglects both radiation recoil and stochasticity.
-- Constructors @Model1@ - @Model4@ can be used for user-defined models.
data Model = BKS | Classical | LL | Model1 | Model2 | Model3 | Model4

-- |Alias for the block matrix type.
type Block (m :: * -> *) = R.Array R.U R.DIM2 (m Double)

class S m where
  -- |\(\hat S\), a matrix which propagates distribution functions in time from 0 to @t = nt * dt@.
  hatS :: Model  -- ^a model to use, e.g. @BKS@ or @Classical@
       -> Double -- ^magnetic field strength
       -> Double -- ^lower energy boundary
       -> Double -- ^upper energy boundary
       -> Int    -- ^dimension of the vector space which represents \(f_e\)
       -> Double -- ^time step, dt
       -> Int    -- ^number of timesteps, nt
       -> Block m

-- |Instance for single-component distribution function, see "Scintillans.Synchrotron#M1c".
instance S Matrix11 where
  hatS BKS       b xa xb nx dt nt = Sol.exp hatA dt nt
    where hatA = R.computeS $ R.map (\x -> M11 x) $ hatA00   b xa xb nx
  hatS Classical b xa xb nx dt nt = Sol.exp hatA dt nt
    where hatA = R.computeS $ R.map (\x -> M11 x) $ hatA00Cl b xa xb nx
  -- Here formulas from Chapter 76 of [L.D. Landau, E.M. Lifshitz, The Classical Theory of Fields]
  -- are used: 1 / epsilon(t) = 1 / epsilon(0) + (2/3) \alpha b t. Note that the method used here
  -- is accurate, i.e. it conserves the number of particles.
  hatS LL        b xa xb nx dt nt = R.computeS
    $ R.fromFunction (R.Z R.:. nx R.:. nx) f
    where f (R.Z R.:. i R.:. j) = if      k == i     then M11 rRight
                                  else if k == i - 1 then M11 rLeft
                                  else                    M11 0
            where x0  = xa + dx * fromIntegral j
                  x   = 1 / ( 1 / x0 + (2/3) * alpha * b * dt * fromIntegral nt )
                  dx = deltaX xa xb nx
                  xdx = (x - xa) / dx
                  k   = floor $ xdx
                  rLeft  = xdx - fromIntegral k
                  rRight = 1 - rLeft
     
-- |Instance for two-component distribution function based on the quantum formulas, see
-- "Scintillans.Synchrotron#M2c".
instance S Matrix22 where
  hatS BKS b xa xb nx dt nt = Sol.exp hatA dt nt
    where hatA = R.computeS $ R.traverse2
            (hatA00 b xa xb nx)
            (hatA10 b xa xb nx)
            (\_ _ -> (R.Z R.:. nx R.:. nx))
            (\lf lf' (R.Z R.:. i R.:. j) -> M22 (lf  (R.Z R.:. i R.:. j))
                                                0
                                                (lf' (R.Z R.:. i R.:. j))
                                                0
            )

