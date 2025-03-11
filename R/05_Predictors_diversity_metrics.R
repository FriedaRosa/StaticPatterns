#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#               05_Predictors_diversity_metrics.R
#                
#
#                    Friederike Wölke 
#                        2025
#
#----------------------------------------------------------#

source(here::here("R/00_Configuration.R"))
lapply(package_list, require, character = TRUE)


#----------------------------------------------------------#
# Load data  -----
#----------------------------------------------------------#

occ_data_final <- readRDS(here::here("Data/output/1_data/1_occ_data_final.rds"))

#----------------------------------------------------------#
# Calculate diversity metrics -----
#----------------------------------------------------------#

div_metrics <- occ_data_final %>%
  ungroup() %>%
  #filter(scalingID == 1) %>%
  dplyr::select(datasetID, samplingPeriodID, siteID, verbatimIdentification) %>%
  distinct() %>%

  group_by(datasetID, samplingPeriodID) %>%
  dplyr::mutate(GammaSR = sum(n_distinct(verbatimIdentification))) %>%
  ungroup() %>%

  group_by(datasetID, samplingPeriodID, siteID) %>%
  dplyr::mutate(AlphaSR = sum(n_distinct(verbatimIdentification))) %>%
  dplyr::mutate(BetaSR = GammaSR / AlphaSR) %>%
  ungroup() %>%

  group_by(datasetID, samplingPeriodID, verbatimIdentification) %>%
  dplyr::mutate(AlphaSR_sp = mean(AlphaSR)) %>%
  dplyr::mutate(BetaSR_sp = GammaSR / AlphaSR_sp) %>%
  dplyr::select(datasetID, samplingPeriodID, verbatimIdentification, AlphaSR_sp, BetaSR_sp, GammaSR) %>%
  distinct() %>%
  ungroup()

#----------------------------------------------------------#
# Save results to .rds  -----
#----------------------------------------------------------#

saveRDS(div_metrics, here::here("Data/output/1_data/single_predictors/2_div_metrics.rds"))
