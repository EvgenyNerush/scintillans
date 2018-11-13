module Scintillans.Probability where

import Numeric.GSL.Special.Airy

-- conventions:
-- w x y means w(x -> y), where w(x -> y) dy is the probability per time unit from the particle with energy x to produce particle
-- with energy in the interval [y, y + dy], where dy is infinitesimal.
-- u x means the full probability rate, i.e. u x = int_0^x w(x -> y) dy.

---------------------------------------------------------------------------------------------------
-- Quantum synchrotron probabilities (Nikishov--Ritus or Baier--Katkov), see [Quantum
-- Electrodynamics, V. B. Berestetskii, E. M. Lifshitz and L. P. Pitaevskii, Pergamon, 1982].
---------------------------------------------------------------------------------------------------
--
-- for emission in a constant mag field, perp
--
-- Here we supose that time is normalized to the radiation formation time
-- t_{rf} = m c / (e B),
-- and energy is normalized to the rest-mass electron energy, mc^2.
-- The first parameter of the probabilities is the magnetic field strength normalized to
-- Sauter-Schwinger (critical) field B_S = m^2 c^3 / e hbar, where hbar = h / 2 pi, and h is the
-- Planck's constant. For these units W is ... 137, chi is ...

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
uqed n b x = trapRule (wqed b x) x' (x - dx) (n - 1) + pec
  where x' = x / (1 + 8 * chi) -- y < x' => arg of Ai >= 4
        dx = (x - x') / fromIntegral n
        -- pec is the approximate value of int_(x - dx)^x (...)
        pec = -alpha / (b * x * x) * (
          dx * intAi 0 +
          airy_Ai_deriv 0 PrecSingle * (
            0.5 * ((dx / x)**5 * chi)**(1/3) +
            6 * ((x * chi)**2 * dx)**(1 / 3)))
            --}
        chi = b * x
