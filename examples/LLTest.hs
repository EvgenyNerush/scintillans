-- Here we test what hatS for Landau-Lifshits (LL) radiation model (see src/Scintillans.hs) yield
-- right value of energy, conserves the number of particles and do not disperse the spectrum more
-- than it is expected

import Scintillans
import Scintillans.Synchrotron
import Scintillans.BlockMatrix

xa = 1e2
xb = 1e3
nx = 100

b = 2e-4

t  = 1e3
dt = 0.1 / alpha
nt = round $ t / dt

-- expected energy value, for these parameters about (very approximately) a half of xb
xc = 1 / ( 1 / xb + (2/3) * alpha * b * t )

dx = deltaX xa xb nx

-- delta-function at xb, containing one particle
f0 = M11 (diracDelta nx dx)

M11 fe = unzipM . toList
       $ mmultS (hatS LL b xa xb nx dt nt :: Block Matrix11)
       $ fromList . zipM
       $ f0

grid = energyGrid xa xb nx

numberOfParticles = dx * sum fe

meanEnergy = ( sum $ zipWith (*) grid fe ) * dx / numberOfParticles

sigma = sqrt $ ( sum $ zipWith g grid fe ) * dx / numberOfParticles
  where g x f = (x - meanEnergy)**2 * f

accuracy = 1e-3

v =  (abs $ numberOfParticles - 1) < accuracy
  && (abs $ xc - meanEnergy) / xc  < accuracy
  && sigma < dx / 2

main =
  if v then
    putStrLn "LLTest: \x1b[32mpassed\x1b[0m"
  else do
    putStrLn "LLTest: \x1b[1;31mfailed\x1b[0m"
    print numberOfParticles
    print xc
    print meanEnergy
    print $ sigma / dx
