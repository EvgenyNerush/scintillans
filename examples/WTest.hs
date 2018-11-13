import Scintillans.Probability

-- Here we test the dependence of the full probability rate *uqed* on *b* in classical (chi << 1)
-- and quantum (chi >> 1) limits. Also, we compute emission power for different chi values and
-- compare it with the values given in [V. B. Berestetskii and E. M. Lifshitz and L. P. Pitaevskii,
-- Quantum Electrodynamics, Pergamon, 1982].

n = 20 :: Int
x0 = 1e3 :: Double

u :: Double -> Double
u b = uqed n b x0

-------------
-- u scalings
-------------

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

--------
-- Power
--------

-- classical synchrotron emission power in units of mc^2 / t_{rf}, where t_{rf} = m c / e B
pcl b x = 2 / 3 * alpha * x**2 * b

-- the same method as in *uqed* from Scintillans.Probability is used
pqed b x = trapRule dpdy x' x n
  where dpdy y = (x - y) * wqed b x y
        x' = x / (1 + 8 * chi)
        chi = b * x

-- power for chi = 0.25 / 3 = 0.0833
p1cl  = pcl  (0.0833 / x0) x0
p1qed = pqed (0.0833 / x0) x0

-- power for chi = 1 / 3
p2cl  = pcl  (0.3333 / x0) x0
p2qed = pqed (0.3333 / x0) x0

-- power for chi = 0.5
p3cl  = pcl  (0.5 / x0) x0
p3qed = pqed (0.5 / x0) x0

accuracy' = 0.05 -- data to compare are not accurate because are taken from the plot (figure) in the
                 -- textbook

v' =  (abs $ p1qed / p1cl - 0.69 ) / 0.69  < accuracy'
   && (abs $ p2qed / p2cl - 0.375) / 0.375 < accuracy'
   && (abs $ p3qed / p3cl - 0.29 ) / 0.29  < accuracy'

main =
  if v && v' then
    putStrLn "WTest: \x1b[32mpassed\x1b[0m"
  else do
    putStrLn "WTest: \x1b[1;31mfailed\x1b[0m"
    print vcl
    print vq
    print (p1qed / p1cl)
    print (p2qed / p2cl)
    print (p3qed / p3cl)
