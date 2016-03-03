{-# LANGUAGE RecordWildCards, OverloadedStrings #-}

module Main where

import qualified Data.Aeson as A
import qualified Data.ByteString.Lazy as BL
import GHC.Generics
import Graphics.Rendering.Chart.Easy
import Graphics.Rendering.Chart.Backend.Cairo


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


plotRuns :: [Run] -> IO ()
plotRuns xs = toFile def "grid.png" $ do
        layout_title .= "Mesh Grid Simulation Time Global/Local Ratio"
        layout_x_axis . laxis_title .= "#vertices"
        layout_y_axis . laxis_title .= "ratio"
        plot (line "global/local" [run xs])

run :: [Run] -> [(Double, Double)]
run xs = [(fromIntegral (nX x*nY x), globalCMA x/localCMA x) | x <- xs]


main :: IO ()
main = do
        xs <- BL.readFile "../benchmark/grid.json"
        case parseRuns xs of
                Just runs -> plotRuns runs
                Nothing   -> putStrLn "Failed to parse file"
