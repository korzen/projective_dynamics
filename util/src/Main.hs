--TODO: this is overengineered and complex for such simple thing

{-# LANGUAGE FlexibleInstances, OverloadedStrings #-}
module Main where


import Control.Applicative (many)
import qualified Data.Aeson as A
import qualified Data.ByteString.Lazy as BL
import Data.List (nub, sort, subsequences)
import qualified Data.Vector as V
import Linear.Metric (norm)
import Linear.V3
import Linear.V4
import Text.Megaparsec (parse, try)
import Text.Megaparsec.Char (anyChar, char, newline, space, string)
import Text.Megaparsec.Combinator (endBy, manyTill, sepBy, sepEndBy)
import Text.Megaparsec.String (Parser)
import Text.Megaparsec.Lexer (decimal, number, signed)
import System.Environment (getArgs)


type Vertex      = V3 Double
type Triangle    = V3 Integer
type Tetrahedron = V4 Integer


data Mesh = TetrahedralMesh { vertices   :: [Vertex]
                            , indices    :: [Triangle]
                            , tetrahedra :: [Tetrahedron]
                            }
        deriving Show


data Constraint = Attachment { i :: Integer, position :: Vertex }
                | Spring     { is :: (Integer, Integer), restLength :: Double }
        deriving Show


instance A.ToJSON Mesh where
        toJSON (TetrahedralMesh vs is ts) =
                A.object ["objects" A..= A.object ["data" A..= d]]
                where d = A.object ["attachments" A..= as
                                   ,"indices"     A..= is
                                   ,"springs"     A..= ss
                                   ,"vertices"    A..= vs
                                   ]
                      ss = mkSprings vs ts
                      as = [] :: [Int]


instance A.ToJSON Vertex where
        toJSON (V3 x y z) = A.toJSON [x, y, z]


instance A.ToJSON Triangle where
        toJSON (V3 a b c) = A.toJSON [a, b, c]


instance A.ToJSON Tetrahedron where
        toJSON (V4 a b c d) = A.toJSON [a, b, c, d]


instance A.ToJSON Constraint where
        toJSON x@(Attachment {}) =
                A.object ["i" A..= i x, "position" A..= position x]
        toJSON x@(Spring {}) =
                A.object ["i" A..= is x, "rest_length" A..= restLength x]


mkSprings :: [Vertex] -> [Tetrahedron] -> [Constraint]
mkSprings vs ts = mkSpring <$> (nub . sort . concatMap (mkEdges . f) $ ts)
        where f (V4 a b c d) = [a, b, c, d]
              mkEdges = filter (\x -> length x == 2) . subsequences
              mkSpring [a,b] = let s = if a < b then (a, b) else (b, a)
                               in Spring s (norm (vs' V.! (fromIntegral a) - vs' V.! (fromIntegral b)))
              vs' = V.fromList vs



list :: String -> (Parser a) -> Parser [a]
list name f = do
        string name
        newline
        decimal
        newline
        xs <- many (f >>= (\x -> newline >> return x))
        return xs


vertex :: Parser Vertex
vertex = do
        (x:y:z:_) <- signedDouble `sepBy` char ' '
        return (V3 x y z)


triangle :: Parser Triangle
triangle = do
        (i0:i1:i2:_) <- decimal `sepBy` char ' '
        return (V3 i0 i1 i2 - pure 1)


tetrahedron :: Parser Tetrahedron
tetrahedron = do
        (i0:i1:i2:i3:_) <- decimal `sepBy` char ' '
        return (V4 i0 i1 i2 i3 - pure 1)


signedDouble :: Parser Double
signedDouble = either fromIntegral id <$> signed space number


mesh :: Parser Mesh
mesh = do
        -- TODO
        manyTill anyChar newline
        manyTill anyChar newline
        vs <- list "Vertices" vertex
        is <- list "Triangles" triangle
        ts <- list "Tetrahedra" tetrahedron
        return (TetrahedralMesh vs is ts)


main :: IO ()
main = do
        [input, output] <- getArgs
        m <- fmap (parse mesh "") . readFile $ input
        case m of
                Left x  -> print x
                Right x -> BL.writeFile output . A.encode $ x
