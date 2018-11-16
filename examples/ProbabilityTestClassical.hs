import qualified Data.Array.Repa     as R
import qualified Scintillans.Solver  as S
import Scintillans.BlockMatrix
import Scintillans.Probability

-- Here we test the classical limit of QED synchrotron emission probability by comparison with the
-- results of the approach using the Landau--Lifshitz radiation reaction force. Namely, we compare
-- the mean particle energy with the theoretical value obtained with Landau--Lifshitz force. See
-- Appendix A in A. S. Samsonov et al., PRA, 2018, https://arxiv.org/abs/1807.04071

-- Performance, use threadscope...

-- ROUGH DRAFT !

epsilon0 = 100 -- 43.6
b        = 3e-2 -- 9.6e-2
t        = 137 * 1 -- 274
epsilon = epsilon0 / 5

n = 50
dx = epsilon0 / fromIntegral (n - 1)
dt = 0.1 --0.01 * alpha
nt = round $ t / dt

hatA :: R.Array R.U R.DIM2 (Matrix11 Double)
hatA = R.computeS $ R.map (\x -> M11 x) $ hatAqed 50 b 0 epsilon0 n

sol = S.exp hatA dt nt

-- initial condition is the beam wich energy per particle is epsilon_0, and
-- \int_0^\epsilon_0 f_0 (x) dx = 1
f0 = R.fromListUnboxed (R.Z R.:. n R.:. (1 :: Int))  $ (replicate (n - 1) $ M11 0) ++ [M11 $ 1 / dx]

transpose :: R.Array R.U R.DIM2 Double -> R.Array R.U R.DIM2 Double
transpose arr = R.computeS $ R.traverse arr (\_ -> R.extent arr) (\f (R.Z R.:. i R.:. j) -> f (R.Z R.:. j R.:. i))

f = mmultS sol f0

norm = (*) dx $ R.sumAllS $ R.map (\(M11 x) -> x) f
norm0 = (*) dx $ R.sumAllS $ R.map (\(M11 x) -> x) f0

main = do
  print $ R.toList $ R.map (\(M11 x) -> x) $ hatA
  --print $ R.toList $ R.map (\(M11 x) -> x) $ f0
  --print $ R.toList $ R.map (\(M11 x) -> x) $ f
  putStrLn " "
  print $ R.toList $ R.sumS $ transpose $ R.computeS $ R.map (\(M11 x) -> x) $ hatA
  putStrLn " "
  print norm0
  print norm
  putStrLn " "
  print $ uqed n b epsilon0
  print $ trapRule (wqed b epsilon0) 0 epsilon0 100
