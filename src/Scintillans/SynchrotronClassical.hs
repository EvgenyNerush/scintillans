-- |This module provides alternative "hatS" realization which is based on the classical synchrotron
-- emission spectrum. See "Scintillans.Synchrotron" for details and for matrices which are based on
-- the quantum (Baier-Katkov-Strakhovenko) formula, and see "Scintillans" for BKS's @hatS@.
--
-- To use the classical @hatS@ instance one should use ...

module Scintillans.SynchrotronClassical
  ( wCl
  , hatWCl
  , hatUCl
  , hatA00Cl
  ) where

import qualified Data.Array.Repa as R
import Numeric.GSL.Special.Airy
import Scintillans.Synchrotron

-- |@wCl b@ is the probability rate of the photon emission in the magnetic field @b@, computed from
-- the classical formula for the synchrotron emission power. More precisely, @wCl b x y =@ \(
-- w_{Cl}(x \to y) \), where \( w_{cl}(x \to y) \Delta y \) is the probability per time unit for an
-- electron (or positron) with the energy \(x\) to emit a photon such that the resulting electron
-- energy is in the interval \( [y, y + \Delta y] \), where \(\Delta y\) is infinitesimal.  Thus,
-- the photon energy is \(x - y\). Note that this function is /unsafe/ and should be used only for
-- \(x > y \geqslant 0\).
wCl :: Double -- ^magnetic field strength
    -> Double -- ^initial electron energy
    -> Double -- ^final electron energy
    -> Double -- ^result
wCl b x y = if z > aiMaxArg then 0
            else -alpha / (b * x * x) * (intAi z + 2 / z * airy_Ai_deriv z PrecSingle)
  where z = ((x - y) / (b * x * x)) ** (2 / 3)

-- Matrices which are based on the classical synchrotron formula. For their quantum analog see
-- "Scintillans.Synchrotron".

hatWCl :: Double                    -- ^magnetic field strength
       -> Double                    -- ^lower bound of the considered energy interval
       -> Double                    -- ^upper bound of the considered energy interval
       -> Int                       -- ^dimension of the vector space which represents \(f_e\)
       -> R.Array R.D R.DIM2 Double -- ^result
hatWCl b xa xb n = R.fromFunction (R.Z R.:. n R.:. n)
  $ \(R.Z R.:. i R.:. j) ->
    if j > i then dx * wCl b (x j) (x i) else 0
  where x k = xa + dx * fromIntegral k
        dx = (xb - xa) / fromIntegral (n - 1)

hatUCl :: Double -> Double -> Double -> Int -> R.Array R.D R.DIM2 Double
hatUCl b xa xb n = R.fromFunction (R.Z R.:. n R.:. n)
  $ \(R.Z R.:. i R.:. j) -> if i /= j then 0 else
    R.sumAllS $ R.slice (hatWCl b xa xb n) (R.Any R.:. j)

hatA00Cl :: Double -> Double -> Double -> Int -> R.Array R.D R.DIM2 Double
hatA00Cl b xa xb n = hatWCl b xa xb n R.-^ hatUCl b xa xb n
