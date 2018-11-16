module Scintillans.Probability where

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


-- Probability rate for the photon emission...
wqed :: Double -> Double -> Double -> Double
wqed b x y = if y >= x || y <= 0 || z > aiMaxArg then 0 else w
  where chi = b * x
        z = ((x - y) / (y * chi)) ** (2 / 3)
        v = 2 / z + (x - y) / x * chi * sqrt z
        w = - alpha / (b * x * x) * (intAi z + v * airy_Ai_deriv z PrecSingle)

-- Full probability of the photon emission... Computed numerically...
uqed :: Int -> Double -> Double -> Double
uqed n b x = if x == 0 then 0 else trapRule (wqed b x) x' (x - dx) (n - 1) + pec
  where x' = x / (1 + 8 * chi) -- y < x' => arg of Ai >= 4
        dx = (x - x') / fromIntegral n
        -- pec is the approximate value of int_(x - dx)^x (...)
        pec = -alpha / (b * x * x) * (
          dx * intAi 0 +
          airy_Ai_deriv 0 PrecSingle * (
            0.5 * ((dx / x)**5 * chi)**(1/3) +
            6 * ((x * chi)**2 * dx)**(1 / 3)))
        chi = b * x

-----------------------------------------------------------------------
-- Right-hand-side of the Boltzmann's equation for synchrotron emission
-----------------------------------------------------------------------

-- E.g., Boltzmann's equation for electrons
-- \partial_t f_e(x) = - u(x) f_e(x) + \int_x^\infty f_e(y) w(y -> x) dy
-- in the matrix form can be written as
-- \partial_t f_e = -\hat U f_e + \hat W f_e

-- Note that the matrices here are of the R.D type and should be computed before use in the Solver

hatUqed :: Int -> Double -> Double -> Double -> Int -> R.Array R.D R.DIM2 Double
hatUqed s b x1 xn n = R.fromFunction (R.Z R.:. n R.:. n) $ \(R.Z R.:. i R.:. j) -> if i == j then u i else 0
-- fromIntegral i else 0
-- u i else 0
  where u i = uqed s b $ x1 + fromIntegral i * dx
        dx = (xn - x1) / fromIntegral (n - 1)

hatWqed :: Double -> Double -> Double -> Int -> R.Array R.D R.DIM2 Double
hatWqed b x1 xn n = R.fromFunction (R.Z R.:. n R.:. n) $ \(R.Z R.:. i R.:. j) -> if j > i then dx * wqed b (x j) (x i) else 0
-- then 1 else 0
-- dx * wqed b (x j) (x i)
  where x k = x1 + fromIntegral k * dx
        dx = (xn - x1) / fromIntegral (n - 1)

hatAqed :: Int -> Double -> Double -> Double -> Int -> R.Array R.D R.DIM2 Double
hatAqed s b x1 xn n = hatWqed b x1 xn n R.-^ hatUqed s b x1 xn n
