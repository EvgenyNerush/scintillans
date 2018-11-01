{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies, TypeFamilies, BangPatterns #-}

module Scintillans.BlockMatrix where

import qualified Data.Array.Repa             as R
import qualified Data.Vector.Generic         as G
import qualified Data.Vector.Generic.Mutable as M
import Data.Vector.Unboxed
import Control.Monad ( liftM )

-- 1x1 matrix
data Matrix11 a = M11 !a

-- 1x2 matrix, i.e. M12 x y= | x y |
data Matrix12 a = M12 !a !a

-- 2x1 matrix, i.e. M21 x y = | x |
--                            | y |
data Matrix21 a = M21 !a !a

-- M22 a11 a12 a21 a22 = | a11 a12 |
--                       | a21 a22 |
data Matrix22 a = M22 !a !a !a !a

data Matrix13 a = M13 !a !a !a

data Matrix31 a = M31 !a !a !a

data Matrix33 a = M33 !a !a !a !a !a !a !a !a !a

-- class to express matrices which can be multiplied
class Multable a b c | a b -> c where
  mult :: a -> b -> c

-- E.g., multiplication of 1x3 matrix to 3x1 matrix results 1x1 matrix
instance Num a => Multable (Matrix13 a) (Matrix31 a) (Matrix11 a) where
  mult (M13 u v w) (M31 x y z)  = M11 $ u * x + v * y + w * z

instance Num a => Multable (Matrix11 a) (Matrix11 a) (Matrix11 a) where
  mult (M11 x) (M11 y) = M11 (x * y)

instance Functor Matrix11 where
    fmap f (M11 x) = M11 $ f x

instance Num a => Num (Matrix11 a) where
    (+) (M11 x) (M11 y) = M11 (x + y)
    fromInteger x = M11 $ fromInteger x

-- identity matrix
idMatrix :: (Num a, Unbox a) => Int -> R.Array R.U R.DIM2 a
idMatrix h = R.computeS $ R.fromFunction (R.Z R.:. h R.:. h) $ \(R.Z R.:. i R.:. j) ->
  if i == j then fromInteger 1 else fromInteger 0

-- Matrix multiplication, as in linear algebra.
mmultS :: (Unbox a, Unbox b, Unbox c, Num c, Multable a b c) => R.Array R.U R.DIM2 a -> R.Array R.U R.DIM2 b -> R.Array R.U R.DIM2 c
mmultS arr brr
  | wa /= hb  = error "Firefly.BlockMatrix mmultS: width of the first matrix != height of the second one."
  | otherwise = arr `R.deepSeqArray` brr `R.deepSeqArray` (R.computeS $ R.fromFunction (R.Z R.:. ha R.:. wb) $ \(R.Z R.:. i R.:. j) -> R.sumAllS $ R.zipWith mult
    (R.slice arr (R.Any R.:. i R.:. R.All)) (R.slice brr (R.Any R.:. j)))
    where (R.Z R.:. ha R.:. wa) = R.extent arr
          (R.Z R.:. hb R.:. wb) = R.extent brr

-----------
-- Matrix11
-----------

newtype instance MVector s (Matrix11 a) = MV_Matrix11 (MVector s a)
newtype instance Vector    (Matrix11 a) = V_Matrix11  (Vector    a)

instance (Unbox a) => Unbox (Matrix11 a)

instance (Unbox a) => M.MVector MVector (Matrix11 a) where
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicOverlaps #-}
  {-# INLINE basicUnsafeNew #-}
  {-# INLINE basicInitialize #-}
  {-# INLINE basicUnsafeReplicate #-}
  {-# INLINE basicUnsafeRead #-}
  {-# INLINE basicUnsafeWrite #-}
  {-# INLINE basicClear #-}
  {-# INLINE basicSet #-}
  {-# INLINE basicUnsafeCopy #-}
  {-# INLINE basicUnsafeGrow #-}
  basicLength (MV_Matrix11 v) = M.basicLength v
  basicUnsafeSlice i n (MV_Matrix11 v) = MV_Matrix11 $ M.basicUnsafeSlice i n v
  basicOverlaps (MV_Matrix11 v1) (MV_Matrix11 v2) = M.basicOverlaps v1 v2
  basicUnsafeNew n = MV_Matrix11 `liftM` M.basicUnsafeNew n
  basicInitialize (MV_Matrix11 v) = M.basicInitialize v
  basicUnsafeReplicate n (M11 x) = MV_Matrix11 `liftM` M.basicUnsafeReplicate n x
  basicUnsafeRead (MV_Matrix11 v) i = M11 `liftM` M.basicUnsafeRead v i
  basicUnsafeWrite (MV_Matrix11 v) i (M11 x) = M.basicUnsafeWrite v i x
  basicClear (MV_Matrix11 v) = M.basicClear v
  basicSet (MV_Matrix11 v) (M11 x) = M.basicSet v x
  basicUnsafeCopy (MV_Matrix11 v1) (MV_Matrix11 v2) = M.basicUnsafeCopy v1 v2
  basicUnsafeMove (MV_Matrix11 v1) (MV_Matrix11 v2) = M.basicUnsafeMove v1 v2
  basicUnsafeGrow (MV_Matrix11 v) n = MV_Matrix11 `liftM` M.basicUnsafeGrow v n

instance (Unbox a) => G.Vector Vector (Matrix11 a) where
  {-# INLINE basicUnsafeFreeze #-}
  {-# INLINE basicUnsafeThaw #-}
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicUnsafeIndexM #-}
  {-# INLINE elemseq #-}
  basicUnsafeFreeze (MV_Matrix11 v) = V_Matrix11 `liftM` G.basicUnsafeFreeze v
  basicUnsafeThaw (V_Matrix11 v) = MV_Matrix11 `liftM` G.basicUnsafeThaw v
  basicLength (V_Matrix11 v) = G.basicLength v
  basicUnsafeSlice i n (V_Matrix11 v) = V_Matrix11 $ G.basicUnsafeSlice i n v
  basicUnsafeIndexM (V_Matrix11 v) i
                = M11 `liftM` G.basicUnsafeIndexM v i
  basicUnsafeCopy (MV_Matrix11 mv) (V_Matrix11 v)
                = G.basicUnsafeCopy mv v
  elemseq _ (M11 x) z = seq x z

-----------
-- Matrix13
-----------

newtype instance MVector s (Matrix13 a) = MV_Matrix13 (MVector s (a, a, a))
newtype instance Vector    (Matrix13 a) = V_Matrix13  (Vector    (a, a, a))

instance (Unbox a) => Unbox (Matrix13 a)

instance (Unbox a) => M.MVector MVector (Matrix13 a) where
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicOverlaps #-}
  {-# INLINE basicUnsafeNew #-}
  {-# INLINE basicInitialize #-}
  {-# INLINE basicUnsafeReplicate #-}
  {-# INLINE basicUnsafeRead #-}
  {-# INLINE basicUnsafeWrite #-}
  {-# INLINE basicClear #-}
  {-# INLINE basicSet #-}
  {-# INLINE basicUnsafeCopy #-}
  {-# INLINE basicUnsafeGrow #-}
  basicLength (MV_Matrix13 v) = M.basicLength v
  basicUnsafeSlice i n (MV_Matrix13 v) = MV_Matrix13 $ M.basicUnsafeSlice i n v
  basicOverlaps (MV_Matrix13 v1) (MV_Matrix13 v2) = M.basicOverlaps v1 v2
  basicUnsafeNew n = MV_Matrix13 `liftM` M.basicUnsafeNew n
  basicInitialize (MV_Matrix13 v) = M.basicInitialize v
  basicUnsafeReplicate n (M13 x y z) = MV_Matrix13 `liftM` M.basicUnsafeReplicate n (x, y, z)
  basicUnsafeRead (MV_Matrix13 v) i = (\(x, y, z) -> M13 x y z) `liftM` M.basicUnsafeRead v i
  basicUnsafeWrite (MV_Matrix13 v) i (M13 x y z) = M.basicUnsafeWrite v i (x, y, z)
  basicClear (MV_Matrix13 v) = M.basicClear v
  basicSet (MV_Matrix13 v) (M13 x y z) = M.basicSet v (x, y, z)
  basicUnsafeCopy (MV_Matrix13 v1) (MV_Matrix13 v2) = M.basicUnsafeCopy v1 v2
  basicUnsafeMove (MV_Matrix13 v1) (MV_Matrix13 v2) = M.basicUnsafeMove v1 v2
  basicUnsafeGrow (MV_Matrix13 v) n = MV_Matrix13 `liftM` M.basicUnsafeGrow v n

instance (Unbox a) => G.Vector Vector (Matrix13 a) where
  {-# INLINE basicUnsafeFreeze #-}
  {-# INLINE basicUnsafeThaw #-}
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicUnsafeIndexM #-}
  {-# INLINE elemseq #-}
  basicUnsafeFreeze (MV_Matrix13 v) = V_Matrix13 `liftM` G.basicUnsafeFreeze v
  basicUnsafeThaw (V_Matrix13 v) = MV_Matrix13 `liftM` G.basicUnsafeThaw v
  basicLength (V_Matrix13 v) = G.basicLength v
  basicUnsafeSlice i n (V_Matrix13 v) = V_Matrix13 $ G.basicUnsafeSlice i n v
  basicUnsafeIndexM (V_Matrix13 v) i
                = (\(x, y, z) -> M13 x y z) `liftM` G.basicUnsafeIndexM v i
  basicUnsafeCopy (MV_Matrix13 mv) (V_Matrix13 v)
                = G.basicUnsafeCopy mv v
  elemseq _ (M13 x y z) t = x `seq` y `seq` z `seq` t

-----------
-- Matrix31
-----------

newtype instance MVector s (Matrix31 a) = MV_Matrix31 (MVector s (a, a, a))
newtype instance Vector    (Matrix31 a) = V_Matrix31  (Vector    (a, a, a))

instance (Unbox a) => Unbox (Matrix31 a)

instance (Unbox a) => M.MVector MVector (Matrix31 a) where
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicOverlaps #-}
  {-# INLINE basicUnsafeNew #-}
  {-# INLINE basicInitialize #-}
  {-# INLINE basicUnsafeReplicate #-}
  {-# INLINE basicUnsafeRead #-}
  {-# INLINE basicUnsafeWrite #-}
  {-# INLINE basicClear #-}
  {-# INLINE basicSet #-}
  {-# INLINE basicUnsafeCopy #-}
  {-# INLINE basicUnsafeGrow #-}
  basicLength (MV_Matrix31 v) = M.basicLength v
  basicUnsafeSlice i n (MV_Matrix31 v) = MV_Matrix31 $ M.basicUnsafeSlice i n v
  basicOverlaps (MV_Matrix31 v1) (MV_Matrix31 v2) = M.basicOverlaps v1 v2
  basicUnsafeNew n = MV_Matrix31 `liftM` M.basicUnsafeNew n
  basicInitialize (MV_Matrix31 v) = M.basicInitialize v
  basicUnsafeReplicate n (M31 x y z) = MV_Matrix31 `liftM` M.basicUnsafeReplicate n (x, y, z)
  basicUnsafeRead (MV_Matrix31 v) i = (\(x, y, z) -> M31 x y z) `liftM` M.basicUnsafeRead v i
  basicUnsafeWrite (MV_Matrix31 v) i (M31 x y z) = M.basicUnsafeWrite v i (x, y, z)
  basicClear (MV_Matrix31 v) = M.basicClear v
  basicSet (MV_Matrix31 v) (M31 x y z) = M.basicSet v (x, y, z)
  basicUnsafeCopy (MV_Matrix31 v1) (MV_Matrix31 v2) = M.basicUnsafeCopy v1 v2
  basicUnsafeMove (MV_Matrix31 v1) (MV_Matrix31 v2) = M.basicUnsafeMove v1 v2
  basicUnsafeGrow (MV_Matrix31 v) n = MV_Matrix31 `liftM` M.basicUnsafeGrow v n

instance (Unbox a) => G.Vector Vector (Matrix31 a) where
  {-# INLINE basicUnsafeFreeze #-}
  {-# INLINE basicUnsafeThaw #-}
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicUnsafeIndexM #-}
  {-# INLINE elemseq #-}
  basicUnsafeFreeze (MV_Matrix31 v) = V_Matrix31 `liftM` G.basicUnsafeFreeze v
  basicUnsafeThaw (V_Matrix31 v) = MV_Matrix31 `liftM` G.basicUnsafeThaw v
  basicLength (V_Matrix31 v) = G.basicLength v
  basicUnsafeSlice i n (V_Matrix31 v) = V_Matrix31 $ G.basicUnsafeSlice i n v
  basicUnsafeIndexM (V_Matrix31 v) i
                = (\(x, y, z) -> M31 x y z) `liftM` G.basicUnsafeIndexM v i
  basicUnsafeCopy (MV_Matrix31 mv) (V_Matrix31 v)
                = G.basicUnsafeCopy mv v
  elemseq _ (M31 x y z) t = x `seq` y `seq` z `seq` t
