import qualified Data.Array.Repa         as R
import qualified Scintillans.Solver      as S
import Scintillans.Synchrotron
import Scintillans.BlockMatrix

-- Here we test the classical limit of QED synchrotron emission probability by comparison with the
-- results of the approach using the Landau--Lifshitz radiation reaction force. Namely, we compare
-- the mean particle energy with the theoretical value obtained with Landau--Lifshitz (LL) force. See
-- Eq. (A7) in Appendix A of Ref. [A. S. Samsonov et al., PRA, 2018,
-- https://arxiv.org/abs/1807.04071], that in the units used in Scintillans reads as
-- gamma = gamma0 / (1 + (2/3) * alpha * b * gamma0 * t)

-- Performance, use threadscope...

-- chi_0 = x0 * b. The classical formula for LL force overestimates the radiation losses more than
-- by 10% if chi > 0.02. Thus one should use really small chi if he uses LL force as a reference.
-- The difference between the mean electron energy and the theoretical value in this test is caused
-- mostly by this difference between quantum and classical formulas for the radiation losses.

-- The width of w(b, x0 -> y) is about x0 chi = x0^2 b = 0.015 x0. dx should be much smaller than
-- that.

-- The LL approach is applicable if an electron emits many photons on average, thus one should use
-- t >> 1 / alpha.

x0 = 100 -- the initial electron energy
b  = 1.5e-4 -- hence chi = 0.015
t  = 685 -- hence every electron emits about 5 photons during the interaction time
x1  = 0.95 * x0 -- the mean electron energy for these parameters from the LL approach
x2  = 0.91 * x0 -- the same for the interaction time of 2 * t

-- the width of the emission powed distribution dp(b, x0 -> y)/dy is about x0 chi = x0^2 b = 1.5,
-- thus dx should be << 0.015 x0. Note that due to the the singularity in w(b, x0 -> y) one
-- should use very small dx to reproduce w shape, however, it is not needed for f.
xl = 75 -- we solve BE on the interval [xl, x0]
n  = 100
dx = (x0 - xl) / fromIntegral (n - 1)
dt = 0.2 / alpha
nt = round $ t / dt

hatA :: R.Array R.U R.DIM2 (Matrix11 Double)
hatA = R.computeS $ R.map (\x -> M11 x) $ hatA00 b xl x0 n

sol = S.exp hatA dt nt

-- initially \int_0^x_0 f_0 (x) dx = 1
f0 = R.fromListUnboxed (R.Z R.:. n R.:. (1 :: Int))
  $ (replicate (n - 1) $ M11 0) ++ [M11 $ 1 / dx]

f1 = mmultS sol f0

f2 = mmultS sol f1

unwrap :: R.Array R.U R.DIM2 (Matrix11 Double) -> [Double]
unwrap h = R.toList $ R.map (\(M11 x) -> x) h

norm :: R.Array R.U R.DIM2 (Matrix11 Double) -> Double
norm h = dx * (sum . unwrap) h

meanX :: R.Array R.U R.DIM2 (Matrix11 Double) -> Double
meanX h = (*) (dx / norm h) $ sum $ zipWith (*) h' xs
  where h' = unwrap h
        xs = [xl + dx * fromIntegral i | i <- [0..(n - 1)]]

norm1  = norm  f1
norm2  = norm  f2
meanX1 = meanX f1
meanX2 = meanX f2

accuracy = 1e-2

v =  (abs $ norm1 - 1) < 1e-15
  && (abs $ norm2 - 1) < 1e-15
  && (abs $ meanX1 - x1) / x1 < accuracy
  && (abs $ meanX2 - x2) / x2 < accuracy

main =
  if v then
    putStrLn "SynchrotronTestClassical: \x1b[32mpassed\x1b[0m"
  else do
    putStrLn "SynchrotronTestClassical: \x1b[1;31mfailed\x1b[0m"
    print $ R.toList $ R.map (\(M11 x) -> x) $ f2
    print norm1
    print norm2
    print meanX1
    print meanX2
