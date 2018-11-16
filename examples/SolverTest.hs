import qualified Data.Array.Repa as R
import qualified Scintillans.Solver  as S
import Scintillans.BlockMatrix

{-- The exponent of the matrix

    | a   b |
A = |       |
    | 0   1 |

is the following:

        | x   y |
exp A = |       | ,
        | 0   e |

where x = exp a, y = b ((exp a) - e) / (a - 1), and e = exp 1 is the Euler constant.
Here we check these relations.

--}

a = 2 :: Double
b = 1 :: Double

-- the matrix A
hatA :: R.Array R.U R.DIM2 (Matrix11 Double)
hatA = R.fromListUnboxed (R.Z R.:. (2 :: Int) R.:. (2 :: Int)) [M11 a, M11 b, M11 0, M11 1]

-- obviously, n * dt should be equal to 1
n = 1000 :: Int
dt = (1 / fromIntegral n) :: Double

[x, y, u, v] = R.toList $ R.map (\(M11 x) -> x) $ S.exp hatA dt n

accuracy = 15 / fromIntegral n -- Euler's absolute accuracy is proportional to 1 / n
res = (abs $ x - exp a) < accuracy
   && (abs $ y - b * (exp a - exp 1) / (a - 1)) < accuracy
   && u == 0
   && (abs $ v - exp 1) < accuracy

main =
  if res then
    putStrLn "SolverTest: \x1b[32mpassed\x1b[0m"
  else
    putStrLn "SolverTest: \x1b[1;31mfailed\x1b[0m"
