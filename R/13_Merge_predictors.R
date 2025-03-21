#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#                      13_Merge_predictors.R
#                
#
#                    Friederike Wölke 
#                        2025
#
#----------------------------------------------------------#
rm(list=ls())
gc()
#----------------------------------------------------------#
# Install and load libraries -----
#----------------------------------------------------------#
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(skimr)
library(sf)
#fix_windows_histograms()

source(here::here("R/00_Configuration.R"))


summary_fun <- function(my_data){

  my_tab <- my_data %>%
    group_by(datasetID, samplingPeriodID) %>%
    rstatix::get_summary_stats(type = "mean_sd")
  return(my_tab)
}

# Raw species data (including columns to filter just to check if other species were lost)

dta0 <- readRDS(here("Data/output/1_data/1_data_sf.rds")) %>%
  st_drop_geometry() %>%
  filter(scalingID == 1 & cells_keep == 1) %>%
  distinct(
    datasetID, samplingPeriodID, verbatimIdentification, scientificName,
    introduced, sp_remove_expert, sp_sampling_repeats, species_keep
  ) %>%
  filter(!is.na(verbatimIdentification)) %>%
  distinct(datasetID, samplingPeriodID, verbatimIdentification, species_keep, .keep_all = TRUE)


#----------------------------------------------------------#
# Load predictors -----
#----------------------------------------------------------#

# birdlife <-
#   readRDS(here("Data/output/1_data/single_predictors/2_all_preds_BirdLife.rds"))

# avonet <-
#   readRDS(here("Data/output/1_data/single_predictors/2_Avonet.rds")) # Avonet data

# iucn <-
#   readRDS(here("Data/output/1_data/single_predictors/2_IUCN_20250225.rds")) # IUCN data

# phylo_dist <-
#   readRDS(here("Data/output/1_data/single_predictors/2_phylo_distinct.rds")) # Phylogenetic distinctiveness

# species_predictors <-
#   full_join(birdlife, avonet) %>%
#   full_join(iucn) %>%
#   full_join(phylo_dist) %>%
#   distinct()

# saveRDS(species_predictors, here("Data/output/1_data/single_predictors/2_predictorsSpecies.rds"))

species_predictors <-
  readRDS(here("Data/output/1_data/single_predictors/2_predictorsSpecies.rds")) %>% # climate niche, range size, pd, avonet
  distinct()

setdiff(dta0 %>% filter(species_keep == 1) %>% pull(verbatimIdentification), species_predictors$verbatimIdentification)

species_predictors2 <- dta0 %>%
  filter() %>%
  left_join(species_predictors) %>%
  rename("IUCN" = "code")

#-----------------------------------------------------#

big_table <-
  readRDS(here("Data/output/1_data/single_predictors/2_big_table.rds"))

summary_fun(big_table)


co_occurrence <-
  readRDS(here("Data/output/1_data/single_predictors/2_cooccurrence.rds")) %>%
  distinct(samplingPeriodID, datasetID, verbatimIdentification, mean_prob_cooccur)

geometry <-
  readRDS(here("Data/output/1_data/single_predictors/2_range_geometries.rds")) # Species ranges, Atlas geometry

sac_metrics <-
  readRDS(here("Data/output/1_data/single_predictors/2_spatial_auto.rds"))

diversity <-
  readRDS(here("Data/output/1_data/single_predictors/2_div_metrics.rds"))

lacunarity <-
  readRDS(here("Data/output/1_data/single_predictors/2_lacunarity_terra_focal.rds")) %>%
  ungroup() %>%
  select(-name) %>%
  mutate(samplingPeriodID = as.numeric(as.character(samplingPeriodID)),
         datasetID = as.numeric(as.character(datasetID))) %>%
  mutate(verbatimIdentification = gsub("_", " ", verbatimIdentification)) %>%
  as.data.frame()
names(lacunarity) <- c("r", "ln(r)", "lac", "ln(lac)", "datasetID", "samplingPeriodID", "verbatimIdentification")


# Calculate mean across increasing window sizes
mean_lac <- lacunarity %>%
  group_by(datasetID, samplingPeriodID, verbatimIdentification) %>%
  summarize(mean_lnLac = mean(`ln(lac)`, na.rm = TRUE))

# quick check on lacunarity data
mean_lac %>%
  group_by(datasetID, samplingPeriodID) %>%
  summarize(n_sp = n_distinct(verbatimIdentification))


#----------------------------------------------------------#
# Merge predictors -----
#----------------------------------------------------------#

predictors <-
  species_predictors2 %>%
  full_join(big_table) %>%
  full_join(co_occurrence) %>%
  full_join(sac_metrics) %>%
  full_join(diversity) %>%
  full_join(geometry) %>%
  full_join(mean_lac) %>%
  distinct(datasetID, verbatimIdentification, samplingPeriodID,
           .keep_all = TRUE) %>%
  mutate(
    across(
      where(is.character) & !matches("verbatimIdentification") & !matches("scientificName"),
      as.factor)) %>%
  mutate(
    across(c("datasetID","samplingPeriodID",
             "Habitat", "IUCN", "Habitat.Density",
             "Migration", "Primary.Lifestyle",
             "Trophic.Niche", "Trophic.Level",
             "Family1", "Order1",
             "introduced", "sp_remove_expert", "sp_sampling_repeats", "species_keep"),
           as.factor)) %>%
  distinct(datasetID, samplingPeriodID, verbatimIdentification, species_keep,
           .keep_all = TRUE)

#----------------------------------------------------------#
# Check predictors -----
#----------------------------------------------------------#

predictors %>%
  glimpse()

str(predictors)
predictors %>%
  is.na() %>%
  colSums()

predictors %>%
  filter(samplingPeriodID == 1) %>%
  skim() %>%
  to_long() %>%
  setNames(c("variable_class", "variable", "metric", "value"))

predictors %>% filter(is.na(species_keep))

predictors %>%
  filter(samplingPeriodID == 1) %>%
  group_by(datasetID) %>%
  skim() %>%
  as_tibble() %>%
  write.csv(here("Documentation/META_predictors_skim_summary.csv"))

names(predictors$D_AOO_a) <- NULL
names(predictors$morans_I) <- NULL
names(predictors$morans_I_p) <- NULL
names(predictors$Lac) <- NULL
str(predictors)


final_predictors <- predictors %>%
  filter(species_keep == 1) %>%
  select(-sp_remove_expert, -sp_sampling_repeats, -species_keep, -introduced) %>%
  distinct()

saveRDS(final_predictors, here("Data/output/1_data/2_predictors.rds"))


final_predictors %>% is.na() %>% colSums()

summary_fun(final_predictors)

skimr::skim(final_predictors)

