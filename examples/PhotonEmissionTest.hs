import qualified Data.Array.Repa         as R
import qualified Scintillans.Solver      as S
import Scintillans.Synchrotron
import Scintillans.BlockMatrix

-- If the interaction time is small, every electron emits less than a single photon in the average,
-- and in the case of monoenergetic electron beam the emission spectrum f_{ph}(x) coincides with
-- the function w(x0 -> x0 - x). Here we compute the emission spectrum and then use it to compute
-- the emission power. The resulting value is compared with that found in the textbook (see
-- WTest.hs for more details). Also, the energy conservation is tested.

-- chi = x0 * b = 1 / 3
x0 = 1e3
b = 1 / (3 * x0)
t = 0.1 / alpha

nt = 10
dt = t / fromIntegral nt

n = 300
xa = 0.5 * x0
xb = x0
dx = (xb - xa) / (fromIntegral $ n - 1)

hatA :: R.Array R.U R.DIM2 (Matrix22 Double)
hatA = R.computeS $ R.traverse2
  (hatA00   b xa xb n)
  (hatA10 b xa xb n)
  (\_ _ -> (R.Z R.:. n R.:. n))
  (\lf lf' (R.Z R.:. i R.:. j) -> M22 (lf  (R.Z R.:. i R.:. j))
                                      0
                                      (lf' (R.Z R.:. i R.:. j))
                                      0
  )

sol = S.exp hatA dt nt

-- for electrons \int_0^\x_0 f_0 (x) dx = 1
f0 :: R.Array R.U R.DIM2 (Matrix21 Double)
f0 = R.fromListUnboxed
  (R.Z R.:. n R.:. (1 :: Int))
  $ (replicate (n - 1) $ M21 0 0) ++ [M21 (1 / dx) 0]

f = mmultS sol f0

-- electron distribution function, as list
fE = R.toList $ R.map (\(M21 x _) -> x) $ f

-- photon distribution function, as list
fPh = R.toList $ R.map (\(M21 _ x) -> x) $ f

xsE  = [xa + dx * fromIntegral i | i <- [0..(n - 1)]]
xsPh = [0  + dx * fromIntegral i | i <- [0..(n - 1)]]

electronEnergy = (*) dx $ sum $ zipWith (*) fE xsE
emittedEnergy  = (*) dx $ sum $ zipWith (*) fPh xsPh

power = emittedEnergy / t

powerCl = 2 / 3 * alpha * x0**2 * b

v = (abs $ power / powerCl - 0.375) / 0.375
v' = (abs $ x0 - electronEnergy - emittedEnergy) / emittedEnergy

main = 
  if v < 0.05 && v' < 1e-13 then
    putStrLn "PhotonEmissionTest: \x1b[32mpassed\x1b[0m"
  else do
    putStrLn "PhotonEmissionTest: \x1b[1;31mfailed\x1b[0m"
    print v
    print v'
