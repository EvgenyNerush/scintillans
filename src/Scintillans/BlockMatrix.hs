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

-- M22 a00 a01 a10 a11 = | a00 a01 |
--                       | a10 a11 |
data Matrix22 a = M22 !a !a !a !a

data Matrix13 a = M13 !a !a !a

data Matrix31 a = M31 !a !a !a

data Matrix33 a = M33 !a !a !a !a !a !a !a !a !a

------------
-- Multables
------------

-- class to express block matrices which can be multiplied
class Multable a b c | a b -> c where
  mult :: a -> b -> c

-- E.g., multiplication of 1x3 matrix to 3x1 matrix results 1x1 matrix
instance Num a => Multable (Matrix13 a) (Matrix31 a) (Matrix11 a) where
  mult (M13 u v w) (M31 x y z)  = M11 $ u * x + v * y + w * z

instance Num a => Multable (Matrix11 a) (Matrix11 a) (Matrix11 a) where
  mult (M11 x) (M11 y) = M11 (x * y)

instance Num a => Multable (Matrix22 a) (Matrix21 a) (Matrix21 a) where
  mult (M22 x y z t) (M21 u v) = M21 (x * u + y * v) (z * u + t * v)

instance Num a => Multable (Matrix22 a) (Matrix22 a) (Matrix22 a) where
  mult (M22 x y z t) (M22 u v w s) =
    M22 (x * u + y * w) (x * v + y * s) (z * u + t * w) (z * v + t * s)

instance Num a => Multable (Matrix33 a) (Matrix31 a) (Matrix31 a) where
  mult (M33 a00 a01 a02 a10 a11 a12 a20 a21 a22) (M31 b00 b10 b20) =
    M31 (a00 * b00 + a01 * b10 + a02 * b20)
        (a10 * b00 + a11 * b10 + a12 * b20)
        (a20 * b00 + a21 * b10 + a22 * b20) 

instance Num a => Multable (Matrix33 a) (Matrix33 a) (Matrix33 a) where
  mult (M33 a00 a01 a02 a10 a11 a12 a20 a21 a22)
       (M33 b00 b01 b02 b10 b11 b12 b20 b21 b22) =
    M33 (a00*b00 + a01*b10 + a02*b20) (a00*b01 + a01*b11 + a02*b21) (a00*b02 + a01*b12 + a02*b22) 
        (a10*b00 + a11*b10 + a12*b20) (a10*b01 + a11*b11 + a12*b21) (a10*b02 + a11*b12 + a12*b22) 
        (a20*b00 + a21*b10 + a22*b20) (a20*b01 + a21*b11 + a22*b21) (a20*b02 + a21*b12 + a22*b22) 

-- identity matrix, smth should be added to its type because fromInteger should be called for
-- square matrices only
idMatrix :: (Num a, Unbox a) => Int -> R.Array R.U R.DIM2 a
idMatrix n = R.computeS $ R.fromFunction (R.Z R.:. n R.:. n) $ \(R.Z R.:. i R.:. j) ->
  if i == j then fromInteger 1 else fromInteger 0

-- Multiplication of Repa matrices, as in linear algebra. The type contained in the resulting
-- matrix should be an instance of Num.
mmultS :: (Unbox a, Unbox b, Unbox c, Num c, Multable a b c) => R.Array R.U R.DIM2 a -> R.Array R.U R.DIM2 b -> R.Array R.U R.DIM2 c
mmultS arr brr
  | wa /= hb  = error "Scintillans.BlockMatrix mmultS: width of the first matrix != height of the second one."
  | otherwise = arr `R.deepSeqArray` brr `R.deepSeqArray` (R.computeS $ R.fromFunction (R.Z R.:. ha R.:. wb) $ \(R.Z R.:. i R.:. j) -> R.sumAllS $ R.zipWith mult
    (R.slice arr (R.Any R.:. i R.:. R.All)) (R.slice brr (R.Any R.:. j)))
    where (R.Z R.:. ha R.:. wa) = R.extent arr
          (R.Z R.:. hb R.:. wb) = R.extent brr

-----------
-- Matrix11
-----------

instance Functor Matrix11 where
    fmap f (M11 x) = M11 $ f x

instance Num a => Num (Matrix11 a) where
    (+) (M11 x) (M11 y) = M11 (x + y)
    fromInteger x = M11 $ fromInteger x

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

-----------
-- Matrix12
-----------

newtype instance MVector s (Matrix12 a) = MV_Matrix12 (MVector s (a, a))
newtype instance Vector    (Matrix12 a) = V_Matrix12  (Vector    (a, a))

instance (Unbox a) => Unbox (Matrix12 a)

instance (Unbox a) => M.MVector MVector (Matrix12 a) where
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
  basicLength (MV_Matrix12 v) = M.basicLength v
  basicUnsafeSlice i n (MV_Matrix12 v) = MV_Matrix12 $ M.basicUnsafeSlice i n v
  basicOverlaps (MV_Matrix12 v1) (MV_Matrix12 v2) = M.basicOverlaps v1 v2
  basicUnsafeNew n = MV_Matrix12 `liftM` M.basicUnsafeNew n
  basicInitialize (MV_Matrix12 v) = M.basicInitialize v
  basicUnsafeReplicate n (M12 x y) = MV_Matrix12 `liftM` M.basicUnsafeReplicate n (x, y)
  basicUnsafeRead (MV_Matrix12 v) i = (\(x, y) -> M12 x y) `liftM` M.basicUnsafeRead v i
  basicUnsafeWrite (MV_Matrix12 v) i (M12 x y) = M.basicUnsafeWrite v i (x, y)
  basicClear (MV_Matrix12 v) = M.basicClear v
  basicSet (MV_Matrix12 v) (M12 x y) = M.basicSet v (x, y)
  basicUnsafeCopy (MV_Matrix12 v1) (MV_Matrix12 v2) = M.basicUnsafeCopy v1 v2
  basicUnsafeMove (MV_Matrix12 v1) (MV_Matrix12 v2) = M.basicUnsafeMove v1 v2
  basicUnsafeGrow (MV_Matrix12 v) n = MV_Matrix12 `liftM` M.basicUnsafeGrow v n

instance (Unbox a) => G.Vector Vector (Matrix12 a) where
  {-# INLINE basicUnsafeFreeze #-}
  {-# INLINE basicUnsafeThaw #-}
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicUnsafeIndexM #-}
  {-# INLINE elemseq #-}
  basicUnsafeFreeze (MV_Matrix12 v) = V_Matrix12 `liftM` G.basicUnsafeFreeze v
  basicUnsafeThaw (V_Matrix12 v) = MV_Matrix12 `liftM` G.basicUnsafeThaw v
  basicLength (V_Matrix12 v) = G.basicLength v
  basicUnsafeSlice i n (V_Matrix12 v) = V_Matrix12 $ G.basicUnsafeSlice i n v
  basicUnsafeIndexM (V_Matrix12 v) i
                = (\(x, y) -> M12 x y) `liftM` G.basicUnsafeIndexM v i
  basicUnsafeCopy (MV_Matrix12 mv) (V_Matrix12 v)
                = G.basicUnsafeCopy mv v
  elemseq _ (M12 x y) z = x `seq` y `seq` z

-----------
-- Matrix21
-----------

instance Num a => Num (Matrix21 a) where
    (+) (M21 x y) (M21 z t) = M21 (x + z) (y + t)
    fromInteger x = M21 (fromInteger x) (fromInteger x)

newtype instance MVector s (Matrix21 a) = MV_Matrix21 (MVector s (a, a))
newtype instance Vector    (Matrix21 a) = V_Matrix21  (Vector    (a, a))

instance (Unbox a) => Unbox (Matrix21 a)

instance (Unbox a) => M.MVector MVector (Matrix21 a) where
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
  basicLength (MV_Matrix21 v) = M.basicLength v
  basicUnsafeSlice i n (MV_Matrix21 v) = MV_Matrix21 $ M.basicUnsafeSlice i n v
  basicOverlaps (MV_Matrix21 v1) (MV_Matrix21 v2) = M.basicOverlaps v1 v2
  basicUnsafeNew n = MV_Matrix21 `liftM` M.basicUnsafeNew n
  basicInitialize (MV_Matrix21 v) = M.basicInitialize v
  basicUnsafeReplicate n (M21 x y) = MV_Matrix21 `liftM` M.basicUnsafeReplicate n (x, y)
  basicUnsafeRead (MV_Matrix21 v) i = (\(x, y) -> M21 x y) `liftM` M.basicUnsafeRead v i
  basicUnsafeWrite (MV_Matrix21 v) i (M21 x y) = M.basicUnsafeWrite v i (x, y)
  basicClear (MV_Matrix21 v) = M.basicClear v
  basicSet (MV_Matrix21 v) (M21 x y) = M.basicSet v (x, y)
  basicUnsafeCopy (MV_Matrix21 v1) (MV_Matrix21 v2) = M.basicUnsafeCopy v1 v2
  basicUnsafeMove (MV_Matrix21 v1) (MV_Matrix21 v2) = M.basicUnsafeMove v1 v2
  basicUnsafeGrow (MV_Matrix21 v) n = MV_Matrix21 `liftM` M.basicUnsafeGrow v n

instance (Unbox a) => G.Vector Vector (Matrix21 a) where
  {-# INLINE basicUnsafeFreeze #-}
  {-# INLINE basicUnsafeThaw #-}
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicUnsafeIndexM #-}
  {-# INLINE elemseq #-}
  basicUnsafeFreeze (MV_Matrix21 v) = V_Matrix21 `liftM` G.basicUnsafeFreeze v
  basicUnsafeThaw (V_Matrix21 v) = MV_Matrix21 `liftM` G.basicUnsafeThaw v
  basicLength (V_Matrix21 v) = G.basicLength v
  basicUnsafeSlice i n (V_Matrix21 v) = V_Matrix21 $ G.basicUnsafeSlice i n v
  basicUnsafeIndexM (V_Matrix21 v) i
                = (\(x, y) -> M21 x y) `liftM` G.basicUnsafeIndexM v i
  basicUnsafeCopy (MV_Matrix21 mv) (V_Matrix21 v)
                = G.basicUnsafeCopy mv v
  elemseq _ (M21 x y) z = x `seq` y `seq` z

-----------
-- Matrix22
-----------

instance Functor Matrix22 where
    fmap f (M22 x y z t) = M22 (f x) (f y) (f z) (f t)

instance Num a => Num (Matrix22 a) where
    (+) (M22 x y z t) (M22 u v w s) = M22 (x + u) (y + v) (z + w) (t + s)
    fromInteger x = M22 (fromInteger x) (fromInteger 0) (fromInteger 0) (fromInteger x)

newtype instance MVector s (Matrix22 a) = MV_Matrix22 (MVector s (a, a, a, a))
newtype instance Vector    (Matrix22 a) = V_Matrix22  (Vector    (a, a, a, a))

instance (Unbox a) => Unbox (Matrix22 a)

instance (Unbox a) => M.MVector MVector (Matrix22 a) where
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
  basicLength (MV_Matrix22 v) = M.basicLength v
  basicUnsafeSlice i n (MV_Matrix22 v) = MV_Matrix22 $ M.basicUnsafeSlice i n v
  basicOverlaps (MV_Matrix22 v1) (MV_Matrix22 v2) = M.basicOverlaps v1 v2
  basicUnsafeNew n = MV_Matrix22 `liftM` M.basicUnsafeNew n
  basicInitialize (MV_Matrix22 v) = M.basicInitialize v
  basicUnsafeReplicate n (M22 x y z t) = MV_Matrix22 `liftM` M.basicUnsafeReplicate n (x, y, z, t)
  basicUnsafeRead (MV_Matrix22 v) i = (\(x, y, z, t) -> M22 x y z t) `liftM` M.basicUnsafeRead v i
  basicUnsafeWrite (MV_Matrix22 v) i (M22 x y z t) = M.basicUnsafeWrite v i (x, y, z, t)
  basicClear (MV_Matrix22 v) = M.basicClear v
  basicSet (MV_Matrix22 v) (M22 x y z t) = M.basicSet v (x, y, z, t)
  basicUnsafeCopy (MV_Matrix22 v1) (MV_Matrix22 v2) = M.basicUnsafeCopy v1 v2
  basicUnsafeMove (MV_Matrix22 v1) (MV_Matrix22 v2) = M.basicUnsafeMove v1 v2
  basicUnsafeGrow (MV_Matrix22 v) n = MV_Matrix22 `liftM` M.basicUnsafeGrow v n

instance (Unbox a) => G.Vector Vector (Matrix22 a) where
  {-# INLINE basicUnsafeFreeze #-}
  {-# INLINE basicUnsafeThaw #-}
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicUnsafeIndexM #-}
  {-# INLINE elemseq #-}
  basicUnsafeFreeze (MV_Matrix22 v) = V_Matrix22 `liftM` G.basicUnsafeFreeze v
  basicUnsafeThaw (V_Matrix22 v) = MV_Matrix22 `liftM` G.basicUnsafeThaw v
  basicLength (V_Matrix22 v) = G.basicLength v
  basicUnsafeSlice i n (V_Matrix22 v) = V_Matrix22 $ G.basicUnsafeSlice i n v
  basicUnsafeIndexM (V_Matrix22 v) i
                = (\(x, y, z, t) -> M22 x y z t) `liftM` G.basicUnsafeIndexM v i
  basicUnsafeCopy (MV_Matrix22 mv) (V_Matrix22 v)
                = G.basicUnsafeCopy mv v
  elemseq _ (M22 x y z t) u = x `seq` y `seq` z `seq` t `seq` u

-----------
-- Matrix33
-----------

instance Functor Matrix33 where
    fmap f (M33 a00 a01 a02 a10 a11 a12 a20 a21 a22) =
      M33 (f a00) (f a01) (f a02)
          (f a10) (f a11) (f a12)
          (f a20) (f a21) (f a22)

newtype instance MVector s (Matrix33 a) = MV_Matrix33 (MVector s ( (a, a, a)
                                                                 , (a, a, a)
                                                                 , (a, a, a) ))
newtype instance Vector    (Matrix33 a) = V_Matrix33  (Vector    ( (a, a, a)
                                                                 , (a, a, a)
                                                                 , (a, a, a) ))

instance (Unbox a) => Unbox (Matrix33 a)

instance (Unbox a) => M.MVector MVector (Matrix33 a) where
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
  basicLength (MV_Matrix33 v) = M.basicLength v
  basicUnsafeSlice i n (MV_Matrix33 v) = MV_Matrix33 $ M.basicUnsafeSlice i n v
  basicOverlaps (MV_Matrix33 v1) (MV_Matrix33 v2) = M.basicOverlaps v1 v2
  basicUnsafeNew n = MV_Matrix33 `liftM` M.basicUnsafeNew n
  basicInitialize (MV_Matrix33 v) = M.basicInitialize v
  basicUnsafeReplicate n (M33 a00 a01 a02
                              a10 a11 a12
                              a20 a21 a22)
    = MV_Matrix33 `liftM` M.basicUnsafeReplicate n ( (a00, a01, a02)
                                                   , (a10, a11, a12)
                                                   , (a20, a21, a22) )
  basicUnsafeRead (MV_Matrix33 v) i
    = (\( (a00, a01, a02)
        , (a10, a11, a12)
        , (a20, a21, a22) ) -> (M33 a00 a01 a02
                                    a10 a11 a12
                                    a20 a21 a22)) `liftM` M.basicUnsafeRead v i
  basicUnsafeWrite (MV_Matrix33 v) i (M33 a00 a01 a02
                                          a10 a11 a12
                                          a20 a21 a22)
    = M.basicUnsafeWrite v i ( (a00, a01, a02)
                             , (a10, a11, a12)
                             , (a20, a21, a22) )
  basicClear (MV_Matrix33 v) = M.basicClear v
  basicSet (MV_Matrix33 v) (M33 a00 a01 a02
                                a10 a11 a12
                                a20 a21 a22)
    = M.basicSet v ( (a00, a01, a02)
                   , (a10, a11, a12)
                   , (a20, a21, a22) )

  basicUnsafeCopy (MV_Matrix33 v1) (MV_Matrix33 v2) = M.basicUnsafeCopy v1 v2
  basicUnsafeMove (MV_Matrix33 v1) (MV_Matrix33 v2) = M.basicUnsafeMove v1 v2
  basicUnsafeGrow (MV_Matrix33 v) n = MV_Matrix33 `liftM` M.basicUnsafeGrow v n

instance (Unbox a) => G.Vector Vector (Matrix33 a) where
  {-# INLINE basicUnsafeFreeze #-}
  {-# INLINE basicUnsafeThaw #-}
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicUnsafeIndexM #-}
  {-# INLINE elemseq #-}
  basicUnsafeFreeze (MV_Matrix33 v) = V_Matrix33 `liftM` G.basicUnsafeFreeze v
  basicUnsafeThaw (V_Matrix33 v) = MV_Matrix33 `liftM` G.basicUnsafeThaw v
  basicLength (V_Matrix33 v) = G.basicLength v
  basicUnsafeSlice i n (V_Matrix33 v) = V_Matrix33 $ G.basicUnsafeSlice i n v
  basicUnsafeIndexM (V_Matrix33 v) i
    = (\( (a00, a01, a02)
        , (a10, a11, a12)
        , (a20, a21, a22) ) -> M33 a00 a01 a02
                                   a10 a11 a12
                                   a20 a21 a22) `liftM` G.basicUnsafeIndexM v i
  basicUnsafeCopy (MV_Matrix33 mv) (V_Matrix33 v)
                = G.basicUnsafeCopy mv v
  elemseq _ (M33 a00 a01 a02
                 a10 a11 a12
                 a20 a21 a22) t = a00 `seq` a01 `seq` a02 `seq` 
                                  a10 `seq` a11 `seq` a12 `seq` 
                                  a20 `seq` a21 `seq` a22 `seq` t

