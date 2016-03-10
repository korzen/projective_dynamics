{-# LANGUAGE RecordWildCards, OverloadedStrings #-}

module Main where

import qualified Data.Aeson as A
import qualified Data.ByteString.Lazy as BL
import GHC.Generics
import Graphics.Rendering.Chart.Easy
import Graphics.Rendering.Chart.Backend.Cairo
import System.Environment (getArgs)


data Unit = Ms
        deriving Show

instance A.FromJSON Unit where
        parseJSON (A.String "ms") = pure Ms


data Run = Run { nX, nY    :: Int
               , localCMA  :: Double
               , globalCMA :: Double
               , unit      :: Unit
               }
        deriving Show

instance A.FromJSON Run where
        parseJSON (A.Object v) = do
                nX <- v A..: "n_x"
                nY <- v A..: "n_y"
                localCMA <- v A..: "local_cma"
                globalCMA <- v A..: "global_cma"
                unit <- v A..: "unit"
                return (Run {..})


parseRuns :: BL.ByteString -> Maybe [Run]
parseRuns = A.decode


plotRuns :: FilePath -> [Run] -> IO ()
plotRuns output xs = toFile def output $ do
        layout_title .= "Mesh Grid Simulation"
        layout_x_axis . laxis_title .= "#vertices"
        layout_y_axis . laxis_title .= "time (ms)"
--        plot (line "global/local" [fmap (f (\x -> globalCMA x/localCMA x)) xs])
        plot (line "global" [fmap (f globalCMA) xs])
        plot (line "local" [fmap (f localCMA) xs])
        where f :: (Run -> Double) -> Run -> (Double, Double)
              f g x = (fromIntegral (nX x*nY x), g x)


main :: IO ()
main = do
        [input, output] <- getArgs
        xs <- BL.readFile input
        case parseRuns xs of
                Just runs -> plotRuns output runs
                Nothing   -> putStrLn "Failed to parse file"
