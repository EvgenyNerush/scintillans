-- |This module provides numerical approximations for probabilityes of the photon emission and pair
-- photoproduction in a constant magnetic field \(B\) which is perpendicular to the particle
-- velocity, in a synchrotron regime. See [Electromagnetic processes at high energies in oriented
-- single crystals, V. N. Baier, V.M. Katkov and V.M. Strakhovenko, World Scientific, Singapore,
-- 1998] and [Quantum Electrodynamics, V. B. Berestetskii, E. M.  Lifshitz and L. P. Pitaevskii,
-- Pergamon, 1982] for details.
--
-- Here time is normalized to the classical radiation formation time \[t_{rf} = m c / (e B),\] and
-- energy \(x\) is normalized to the rest-mass electron energy, \(mc^2\).  The first parameter of
-- the probabilities, \(b\), is the magnetic field strength \(B\) normalized to the
-- Sauter-Schwinger (critical) field \(B_S = m^2 c^3 / e \hbar\). For these units \(\chi = b x\),
-- where, again, \(x\) is the normalized electron energy.
--
-- Note that the radiation time, i.e. time on that one photon is emitted on average, is about
-- \[t_{rad} \sim 1 / \alpha\] in the classical limit ( \(\chi \ll 1\) ) and is \[t_{rad} \sim
-- \chi^{1/3} / \alpha\] in the quantum limit ( \(\chi \gg 1\) ). In the "Scintillans.Solver" one
-- should use the timestep \(\Delta t \ll t_{rad}\).

module Scintillans.Synchrotron
  ( -- *Probabilities
    alpha
  , trapRule
  , aiMaxArg
  , intAi
  , w
  , s
  , u
  , p
    -- *Matrices
    -- $Matrices
  , hatW
  , hatU
    -- **Matrices for electron distribution function
    -- $Matrices-1c
  , hatA00
    -- **Matrices for two-component distribution function
    -- $Matrices-2c
  , hatA10
  ) where

import qualified Data.Array.Repa as R
import Numeric.GSL.Special.Airy

-- |\(\alpha \equiv e^2 / \hbar c \approx 1 / 137\)
alpha = 1 / 137.04 :: Double

-- |Trapezoidal rule of numerical integration, \(\int_a^b f(x) \, dx \approx\) @trapRule f a b n@,
-- with @n@ the number of /intervals/ used for the integration.
trapRule :: (Double -> Double) -> Double -> Double -> Int -> Double
trapRule f a b n = 0.5 * dx * (f a + f b + 2 * (sum $ map f $ mid))
  where dx = (b - a) / fromIntegral n
        mid = [a + dx * (fromIntegral i) | i <- [1..(n - 1)]]

-- |We suppose \(\operatorname{Ai}(x) = 0\) and \(\operatorname{Ai}'(x) = 0\) for \(x >\)
-- @aiMaxArg@. Function \( \operatorname{Ai}(x) \) decays sharply with \(x\): it is is about 1e-3
-- at \(x = 4\) and is already out of Double precision approximately at \(x > 100\).
aiMaxArg = 30 :: Double

-- |@intAi x@ is the integral \(\int_x^\infty \operatorname{Ai}(y) \, dy \) computed numerically by
-- 'trapRule' with accuracy of about one percent.
intAi :: Double -> Double
intAi x = if x > aiMaxArg then 0 else trapRule (\y -> airy_Ai y PrecSingle) x x' 17
  where x' = if x < 25 / 16 then x + 4 else x + 5 / sqrt x

-- |@w b@ is the probability rate of the photon emission in the magnetic field @b@. Namely,
-- @w b x y =@ \( w(x \to y) \), where \( w(x \to y) \Delta y \) is the probability per time unit
-- for an electron (or positron) with energy \(x\) to emit a photon such that the resulting
-- electron energy is in the interval \( [y, y + \Delta y] \), where \(\Delta y\) is infinitesimal.
-- Thus, the photon energy is approximately \(x - y\). Note that this function is /unsafe/ and
-- should be used only for \(x > y \geqslant 0\).
w :: Double -- ^magnetic field strength
  -> Double -- ^initial electron energy
  -> Double -- ^final electron energy
  -> Double -- ^result
w b x y = if z > aiMaxArg then 0
          else -alpha / (b * x * x) * (intAi z + v * airy_Ai_deriv z PrecSingle)
  where chi = b * x
        z = ((x - y) / (y * chi)) ** (2 / 3)
        v = 2 / z + (x - y) / x * chi * sqrt z

-- |For magnetic field @b@ and small energy interval @dx@ \(\equiv \Delta x \ll x\), the integral
-- \(\int_{x - \Delta x}^x w(x \to y) \, dy \approx\) @s b x dx@. Therefore, @s b x dx@ is a
-- contribution of the vicinity of @x@ to the full probability rate.
s :: Double -> Double -> Double -> Double
s b x dx = -alpha / (b * x * x) * (dx * intAi 0
  + 6 * airy_Ai_deriv 0 PrecSingle * ((b * x * x)**2 * dx)**(1 / 3))

-- |@u n b x@ is the full probability rate of the photon emission \(\int_0^x w(x \to y) \, dy\) in
-- the field @b@ computed with trapezoidal rule on @n@ intervals, with @x@ the energy of the
-- electron emitting the photon.
u :: Int -> Double -> Double -> Double
u n b x = s b x dx + trapRule (w b x) x' (x - dx) (n - 1)
  where x' = x / (1 + 8 * b * x) -- y < x' => arg of Ai >= 4
        dx = (x - x') / fromIntegral n

-- |@p n b x@ is the power of the photon emission \(\int_0^x (x - y) w(x \to y) \, dy\)  in the
-- field @b@ computed with trapezoidal rule on @n@ intervals, with @x@ the initial energy of the
-- electron.  Note that opposite to \( w(x \to y) \), there is no singularity in the power
-- distribution, moreover, \( (x - y) w(x \to y) = 0 \) at \(y = x\).
p n b x = trapRule dpdy x' x n
  where dpdy y = if y /= x then (x - y) * w b x y else 0
        x' = x / (1 + 8 * chi)
        chi = b * x


-- $Matrices
-- /Scintillans/ uses matrix representation of the equations which describe electron, photon and
-- positron distribution functions \(f \equiv dN / d\epsilon\ \equiv dN / dx\) in a constant
-- magnetic field. Matrices 'hatW', 'hatU', 'hatA00' and 'hatA10' can be used to build the equation
-- for the two-component distribution function (distribution function of electrons and photons),
-- which is described in detail in [I I Artemenko et al 2019 Plasma Phys. Control. Fusion 61
-- 074003](https://doi.org/10.1088/1361-6587/ab1712) (or see freely available
-- [preprint](https://www.researchgate.net/publication/332283915_Global_constant_field_approximation_for_radiation_reaction_in_collision_of_high-intensity_laser_pulse_with_electron_beam)).
-- Note that the matrices below are /delayed/ (of
-- [Repa.D](http://hackage.haskell.org/package/repa-3.4.1.4/docs/Data-Array-Repa.html) type) and
-- should be computed before use with "Scintillans.Solver".  Also note that it is assumed that the
-- lower boundary of the considered energy interval \(x_a\) is such small that the photon emission
-- is negligible for electrons with \(x_a\), or there is almost no electrons with energy about
-- \(x_a\). Otherwise /Scintillans/ does not treat \(f_e\) correctly near \(x_a\).


-- |\(\hat W\), as described in /Artemenko 2019/, see links above for details.
hatW :: Double                    -- ^magnetic field strength
     -> Double                    -- ^lower bound of the considered energy interval
     -> Double                    -- ^upper bound of the considered energy interval
     -> Int                       -- ^dimension of the vector space which represents \(f_e\)
     -> R.Array R.D R.DIM2 Double -- ^result
hatW b xa xb n = R.fromFunction (R.Z R.:. n R.:. n)
  $ \(R.Z R.:. i R.:. j) ->
    if j > i then dx * w b (x j) (x i) else 0
  where x k = xa + dx * fromIntegral k
        dx = (xb - xa) / fromIntegral (n - 1)

-- |\(\hat U\), as described in /Artemenko 2019/. The arguments are the same as in 'hatW'.
hatU :: Double -> Double -> Double -> Int -> R.Array R.D R.DIM2 Double
hatU b xa xb n = R.fromFunction (R.Z R.:. n R.:. n)
  $ \(R.Z R.:. i R.:. j) -> if i /= j then 0 else
    R.sumAllS $ R.slice (hatW b xa xb n) (R.Any R.:. j)

-- $Matrices-1c
-- #M1c#
-- If one is interested only in the electron spectrum and pair photoproduction is not taken into
-- account, then Boltzmann\'s equation for single-component distribution function (i.e., electron
-- spectrum) can be used:
-- \[\partial_t f_e = \hat A_{00} f_e.\]

-- |\(\hat A_{00}\), as described in /Artemenko 2019/. The arguments are the same as in 'hatW'.
hatA00 :: Double -> Double -> Double -> Int -> R.Array R.D R.DIM2 Double
hatA00 b xa xb n = hatW b xa xb n R.-^ hatU b xa xb n

-- $Matrices-2c
-- #M2c#
-- In two-component representation of the distribution function 2x1 block matrix is used:
-- \[
-- f = \left( \begin{array}{c}
--   f_e \\
--   f_{ph}
-- \end{array} \right),
-- \]
-- with \(f_e\) represented by values at the nodes
-- \(x_a, x_a + \Delta x, x_a + 2 \Delta_x, ...  x_b\) and \(f_{ph}\) by values at the nodes
-- \(0, \Delta x, ...  x_b - x_a\). In two-component representation the right-hand-side of
-- Boltzmann\'s equation is a 2x2 block matrix:
-- \[
-- \partial_t f = \left( \begin{array}{cc} \hat A_{00} & \hat 0 \\
--                                         \hat A_{10} & \hat 0
--                       \end{array} \right) f.
-- \]

-- |\(\hat A_{10}\), as described in /Artemenko 2019/. The arguments are the same as in 'hatW'.
hatA10 :: Double -> Double -> Double -> Int -> R.Array R.D R.DIM2 Double
hatA10 b xa xb n = R.traverse
  (hatW b xa xb n)
  (\_ -> (R.Z R.:. n R.:. n))
  (\lf (R.Z R.:. i R.:. j) -> if j >= i
                              then lf (R.Z R.:. (j - i) R.:. j)
                              else 0)
