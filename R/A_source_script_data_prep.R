library(here)
library(dplyr)
library(sf)
library(purrr)
library(future)
library(furrr)

## Configuration file ✅
source(here::here("R/00_Configuration.R"))

## Get data from sql database ✅
# source(here::here("R/01_Get_data.R"))

## Preprocess data ✅
## (handle cells not sampled, introduced sp., etc)
# source(here::here("R/02_Preprocess_data.R"))

## Calculate logRatio & Jaccard ✅
# source(here::here("R/03_Calculate_atlas_variables.R"))

## Predictors H1: Mean probability to co occur ✅
# source(here::here("R/04_Predictors_mean_prob_cooccurr.R"))

## Predictors H3: Diversity Metrics ✅
# source(here::here("R/05_Predictors_diversity_metrics.R"))

## Predictors H1: Avonet traits ✅
# source(here::here("R/06_Predictors_avonet.R"))

## Predictors H1: Phylogenetic distinctiveness ✅
# source(here::here("R/07_Predictors_phylogenetic_metrics.R))

## Predictors H1: IUCN Threat status ✅
# source(here::here("R/08_Predictors_IUCN.R))

## Predictors H1: Global climate niche breadth ✅
# source(here::here("R/09_Predictors_climate_niches.R"))

## Predictors H2: Geometry features ✅
# source(here::here("R/10_Predictors_geometry.R"))

## Predictors H2: Spatial autocorrelation (Morans I, Join count) ✅
source(here::here("R/11_Predictors_spatial_autocorrelation.R"))

## Predictors H2: Lacunarity
# source(here::here("R/12a_Predictors_Lacunarity_rasterize_ranges.R")) ✅
# source(here::here("R/12b_Predictors_Lacunarity_calculation.R"))


## Merge predictors:
source(here::here("R/13_Merge_predictors.R"))

