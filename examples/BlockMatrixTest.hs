import qualified Data.Array.Repa as R
import Scintillans.BlockMatrix

-- This test program multiplies two long 1D vectors (1, 0, 1, 0, ...) implemented as block
-- matrices, and also can multiply the same vectors implemented as Repa arrays. Use *threadscope*
-- to see that block matrices bring little overhead.

-- This test takes about 200 MB of memory and 5 seconds of CPU time.

-----------------
-- Block matrices
-----------------

-- number of elements in a block
n = 3000000 :: Int

-- helper  function
f i = sin $ 0.5 * pi * fromIntegral i

-- matrices to multiply

m1 :: R.Array R.U R.DIM2 (Matrix13 Double)
m1 = R.computeS $ R.fromFunction
  (R.Z R.:. (1 :: Int) R.:. n)
  $ \(R.Z R.:. _ R.:. j) ->
    (M13
      (f $ 3 * j)
      (f $ 3 * j + 1)
      (f $ 3 * j + 2)
    )

m2 :: R.Array R.U R.DIM2 (Matrix31 Double)
m2 = R.computeS $ R.fromFunction
  (R.Z R.:. n R.:. (1 :: Int))
  $ \(R.Z R.:. i R.:. _) ->
    (M31
      (f $ 3 * i)
      (f $ 3 * i + 1)
      (f $ 3 * i + 2)
    )

-- resulting value
v :: Double
v = R.sumAllS $ R.map (\(M11 x) -> x) $ mmultS m1 m2

-----------------------------------
-- The same things with Repa arrays
-----------------------------------

a1 :: R.Array R.U R.DIM2 Double
a1 = R.computeS $ R.fromFunction
  (R.Z R.:. (1 :: Int) R.:. (3 * n))
  $ \(R.Z R.:. _ R.:. j) -> (f $ fromIntegral j)

a2 :: R.Array R.U R.DIM2 Double
a2 = R.computeS $ R.fromFunction
  (R.Z R.:. (3 * n) R.:. (1 :: Int))
  $ \(R.Z R.:. i R.:. _) -> (f $ fromIntegral i)

a2' :: R.Array R.U R.DIM2 Double
a2' = R.computeS $ R.reshape (R.extent a1) a2

-- resulting value
w :: Double
w = R.sumAllS $ a1 R.*^ a2'

main =
  -- use *w* instead of *v* below to compare Block matrices with Repa arrays
  if v == (1.5 * fromIntegral n) then
    putStrLn "BlockMatrixTest: \x1b[32mpassed\x1b[0m"
  else
    putStrLn "BlockMatrixTest: \x1b[1;31mfailed\x1b[0m"
