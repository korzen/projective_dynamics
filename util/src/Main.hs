module Main where

import Control.Applicative (many)
import Text.Megaparsec (parse, try)
import Text.Megaparsec.Char (anyChar, char, newline, space, string)
import Text.Megaparsec.Combinator (endBy, manyTill, sepBy, sepEndBy)
import Text.Megaparsec.String (Parser)
import Text.Megaparsec.Lexer (decimal, number, signed)
import System.Environment (getArgs)


data V3 = V3 Double Double Double
        deriving Show
data I3 = I3 Integer Integer Integer
        deriving Show
data I4 = I4 Integer Integer Integer Integer
        deriving Show


data Mesh = TetrahedralMesh { vertices   :: [V3]
                            , indices    :: [I3]
                            , tetrahedra :: [I4]
                            }
        deriving Show


list :: String -> (Parser a) -> Parser [a]
list name f = do
        string name
        newline
        decimal
        newline
        xs <- many (f >>= (\x -> newline >> return x))
        return xs


vertex :: Parser V3
vertex = do
        (x:y:z:_) <- signedDouble `sepBy` char ' '
        return (V3 x y z)


triangle :: Parser I3
triangle = do
        (i0:i1:i2:_) <- decimal `sepBy` char ' '
        return (I3 i0 i1 i2)


tetrahedron :: Parser I4
tetrahedron = do
        (i0:i1:i2:i3:_) <- decimal `sepBy` char ' '
        return (I4 i0 i1 i2 i3)


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
        print =<< (fmap (parse mesh "") . readFile $ input)
        return ()
