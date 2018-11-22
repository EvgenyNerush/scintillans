module Scintillans.Synchrotron where

import qualified Data.Array.Repa as R
import Numeric.GSL.Special.Airy

-- w x y means w(x -> y), where w(x -> y) dy is the probability per time unit from the particle
-- with energy x to produce particle with energy in the interval [y, y + dy], where dy is
-- infinitesimal.
-- u x means the full probability rate, i.e. u x = \int_0^x w(x -> y) dy.

---------------------------------------------------------------------------------------------------
-- Quantum synchrotron probabilities (Nikishov--Ritus or Baier--Katkov), see [Quantum
-- Electrodynamics, V. B. Berestetskii, E. M. Lifshitz and L. P. Pitaevskii, Pergamon, 1982].
---------------------------------------------------------------------------------------------------
--
-- Here we supose that an external constant magnetic field B is perpendicular to the electron
-- velocity. Time is normalized to the radiation formation time
-- t_{rf} = m c / (e B),
-- and energy is normalized to the rest-mass electron energy, mc^2.
-- The first parameter of the probabilities is the magnetic field strength normalized to
-- Sauter-Schwinger (critical) field B_S = m^2 c^3 / e hbar, where hbar = h / 2 pi, and h is the
-- Planck's constant. For these units chi = b x, where x is the normalized electron energy.
--
-- Note that the radiation time, i.e. time on that one photon is emitted in average, is about
-- t_rad ~ 1 / alpha
-- in the classical limit (chi << 1) and is
-- t_rad ~ chi^{1/3} / alpha
-- in the quantum limit (chi >> 1). In the Solver one should use dt << t_rad.

-- e^2 / hbar c
alpha = 1 / 137.04 :: Double

-- Ai(x) decay sharply with x and is about 1e-3 at x = 4. Ai(x) is already out of Double precision
-- approximately at x > 100.
aiMaxArg = 30 :: Double

-- Trapezoidal rule of numerical integration
trapRule :: (Double -> Double) -> Double -> Double -> Int -> Double
trapRule f a b n = 0.5 * dx * (f a + f b + 2 * (sum $ map f $ mid))
  where dx = (b - a) / fromIntegral n
        mid = [a + dx * (fromIntegral i) | i <- [1..(n - 1)]]

-- int_x^\infty Ai(y) dy, computed numerically with accuracy of about one percent.
intAi :: Double -> Double
intAi x = if x > aiMaxArg then 0 else trapRule (\y -> airy_Ai y PrecSingle) x x' 17
  where x' = if x < 25 / 16 then x + 4 else x + 5 / sqrt x

-- Probability rate for the photon emission... y>= x...
--w b x y = if y >= x || y <= 0 || z > aiMaxArg then 0 else w
w :: Double -> Double -> Double -> Double
--w b x y = if y >= x || y <= 0 || z > aiMaxArg then 0 else w
w b x y = if z > aiMaxArg then 0
          else -alpha / (b * x * x) * (intAi z + v * airy_Ai_deriv z PrecSingle)
  where chi = b * x
        z = ((x - y) / (y * chi)) ** (2 / 3)
        v = 2 / z + (x - y) / x * chi * sqrt z

-- Part of the full probability rate near the singularity..., int_(x-dx)^x...
-- here we assume dx << x
s :: Double -> Double -> Double -> Double
s b x dx = -alpha / (b * x * x) * (dx * intAi 0
  + 6 * airy_Ai_deriv 0 PrecSingle * ((b * x * x)**2 * dx)**(1 / 3))

-- Full probability of the photon emission... Computed numerically...
u :: Int -> Double -> Double -> Double
u n b x = s b x dx + trapRule (w b x) x' (x - dx) (n - 1)
  where x' = x / (1 + 8 * b * x) -- y < x' => arg of Ai >= 4
        dx = (x - x') / fromIntegral n

-- The emission power. There is no singularity in the power distribution, moreover,
-- w(x -> y) * (x - y) = 0 at y = x.
p n b x = trapRule dpdy x' x n
  where dpdy y = if y /= x then (x - y) * w b x y else 0
        x' = x / (1 + 8 * chi)
        chi = b * x

-----------------------------------------------------------------------
-- Right-hand-side of the Boltzmann's equation for synchrotron emission
-----------------------------------------------------------------------

-- E.g., Boltzmann's equation for electrons
-- \partial_t f_e(x) = - u(x) f_e(x) + \int_x^\infty f_e(y) w(y -> x) dy
-- in the matrix form can be written as
-- \partial_t f_e = -\hat U f_e + \hat W f_e

-- Note that the matrices here are of the R.D type and should be computed before use in the Solver

hatW :: Double -> Double -> Double -> Int -> R.Array R.D R.DIM2 Double
hatW b xa xb n = R.fromFunction (R.Z R.:. n R.:. n)
  $ \(R.Z R.:. i R.:. j) ->
    if j > i then dx * w b (x j) (x i) else 0
  where x k = xa + dx * fromIntegral k
        dx = (xb - xa) / fromIntegral (n - 1)

hatU :: Double -> Double -> Double -> Int -> R.Array R.D R.DIM2 Double
hatU b xa xb n = R.fromFunction (R.Z R.:. n R.:. n)
  $ \(R.Z R.:. i R.:. j) -> if i /= j then 0 else
    R.sumAllS $ R.slice (hatW b xa xb n) (R.Any R.:. j)

hatA :: Double -> Double -> Double -> Int -> R.Array R.D R.DIM2 Double
hatA b xa xb n = hatW b xa xb n R.-^ hatU b xa xb n
