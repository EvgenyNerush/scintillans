import Scintillans.Probability

-- \int_x^\infty Ai(y) dy for different values of x
intAi_0   = 0.33333
intAi_0_5 = 0.18738
intAi_1   = 9.7016e-2
intAi_2   = 2.08006e-2
intAi_4   = 4.40688e-4
intAi_8   = 1.60908e-8
intAi_16  = 1.02749e-20

vals = [intAi_0, intAi_0_5, intAi_1, intAi_2, intAi_4, intAi_8, intAi_16]

xs :: [Double]
xs = 0:[0.5 * 2**(fromIntegral i) | i <- [0..5]]

ys = map intAi xs

relativeAccuracy = 0.01

f x y = abs(x - y) / x < relativeAccuracy

res = foldl (&&) True $ zipWith f vals ys

main =
  if res then do
    putStrLn "IntAiTest: \x1b[32mpassed\x1b[0m"
  else
    putStrLn "IntAiTest: \x1b[1;31mfailed\x1b[0m"
