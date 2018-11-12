import Scintillans.Probability

-- Here we test the dependence of the full probability rate *uqed* on *b* in classical (chi << 1)
-- and quantum (chi >> 1) limits. Also, we compute emission power for different chi values and
-- compare it with the values given in [V. B. Berestetskii and E. M. Lifshitz and L. P. Pitaevskii,
-- Quantum Electrodynamics, Pergamon, 1982].

n = 100 :: Int
x = 1e3 :: Double

u :: Double -> Double
u b = uqed n b x

-- In the classical limit the full probability rate is proportional to the magnetic field strength
-- *b*; however, being normalized to t_rf (that depend on time), *u* should not depend on b.
-- u(chi = 1e-2) / u(chi = 1e-3)
vcl = u 1e-5 / u 1e-6

-- In the quantum limit non-normalized *u* is proportional to b^{2/3}, and normalized *u* is
-- proportional to b^{-1/3}.
-- u(chi = 1e3) / u(chi = 1e2)
vq = u 1 / u 0.1

accuracy = 2e-2 -- for higher accuracy one should use higher/lower values of chi in the
                -- quantum/classical regimes

v =  (abs $ vcl - 1) < accuracy
  && (abs $ vq - 10**(-1/3)) / (10**(-1/3)) < accuracy

main =
  if v then
    putStrLn "WTest: \x1b[32mpassed\x1b[0m"
  else do
    putStrLn "WTest: \x1b[1;31mfailed\x1b[0m"
