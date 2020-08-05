import qualified Data.Array.Repa         as R
import Scintillans
import Scintillans.BlockMatrix
import Scintillans.Synchrotron
import Scintillans.SynchrotronClassical

-- For parameters which lead to chi <<< 1, we compute the hatS matrix with classical and quantum
-- radiation reaction models, and compare the results. Then, for chi~0.01 we compare the classical
-- and quantum radiation power, and compare their difference with that is predicted by the theory.

x     = 2e3 -- the initial electron Lorentz factor
bLow  = 1e-7 -- hence chi = 2e-4
bHigh = 1e-5 -- hence chi = 0.02
xa bField = x * (1 - 10 * x * bField)
xLow  = xa bLow
xHigh = xa bHigh

------------------
-- Very low chi --
------------------

dt = 0.1 / alpha
t = 1 / alpha
nt = round $ t / dt
nx = 100 :: Int

hatSClLow   = hatS Classical bLow xLow x nx dt nt       :: Block Matrix11
hatSBKSLow  = hatS BKS       bLow xLow x nx dt nt       :: Block Matrix11
hatSBKSLow' = hatS BKS       bLow xLow x nx dt (2 * nt) :: Block Matrix11

unwrap = R.map (\(M11 x) -> x)

norm :: Double
norm = sqrt
     $ R.sumAllS
     $ R.map (\x -> x * x)
     $ R.zipWith (-) (unwrap hatSClLow) (unwrap hatSBKSLow')

diff :: Double
diff = sqrt
     $ R.sumAllS
     $ R.map (\x -> x * x)
     $ R.zipWith (-) (unwrap hatSClLow) (unwrap hatSBKSLow)

v =  diff < 2 * norm * bLow * x

-------------
-- Low chi --
-------------

-- If the magnetic field is weak (and quantum parameter $\chi$ is small), the emitted energy
-- computed with BKS formula is $ I \approx (1 - 6 \chi + 48 \chi^2) I_{cl} $, where $I_{cl}$ is
-- the emitted energy computed with the classical synchrotron formula. See [V. B. Berestetskii and
-- E.  M. Lifshitz and L. P.  Pitaevskii, Quantum Electrodynamics, Pergamon, New York, 1982]. Here
-- we check the relation between the classical and the quantum results.

radiationPower :: Model -> Double
radiationPower model
  = R.sumAllS
  $ R.traverse vec
               (\_ -> (R.Z R.:. nx))
               (\lf (R.Z R.:. i) -> fromIntegral (nx - 1 - i) * dx * lf (R.Z R.:. i))
  where vec = R.slice (w model bHigh xHigh x nx) (R.Any R.:. (nx - 1))
        w BKS       = hatW
        w Classical = hatWCl
        dx  = deltaX xHigh x nx

chi  = bHigh * x -- hence (1 - 6 chi + 48 chi^2) = (1 - 0.12 + 0.0192) = 0.8992
rpCl = radiationPower Classical
rpQ  = radiationPower BKS

v' = ( abs $ rpCl * (1 - 6 * chi + 48 * chi**2) - rpQ ) / rpQ < 0.01

main =
  if v && v' then
    putStrLn "SynchrotronClassicalTest: \x1b[32mpassed\x1b[0m"
  else do
    putStrLn "SynchrotronClassicalTest: \x1b[1;31mfailed\x1b[0m"
    print   diff
    print $ 2 * norm * bLow * x
    print $ rpCl * 0.8992
    print   rpQ
