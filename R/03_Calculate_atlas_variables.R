#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#                 03_Calculate_atlas_variables.R
#                
#
#                    Friederike Wölke 
#                        2025
#
#----------------------------------------------------------#


# Start with clean environment
rm(list = ls())
gc()

# Source 00_Configuration.R
source(here::here("R/00_Configuration.R"))
lapply(package_list, require, character = TRUE)

# Enable s2 processing for sf package with projected data
sf_use_s2(TRUE)


#----------------------------------------------------------#
# Load data -----
#----------------------------------------------------------#

# Read spatial dataset
data_sf <- readRDS(here("Data/output/1_data/1_data_sf.rds"))

#----------------------------------------------------------#
# Calculate grid information  -----
#----------------------------------------------------------#

# Custom function to get total and sampled area/cells from datasets
calculate_grid_info <- function(data_sf) {
  grid_info <-
    full_join(
      data_sf %>%
      st_drop_geometry() %>%
      distinct(datasetID, scalingID, siteID, cell_sampling_repeats, croppedArea) %>%
      group_by(datasetID, scalingID) %>%
      summarise(
        Total_Ncells = n_distinct(siteID)),
    data_sf %>%
      st_drop_geometry() %>%
      distinct(datasetID, scalingID, siteID, cell_sampling_repeats, croppedArea) %>%
      filter(cell_sampling_repeats == 2) %>%
      group_by(datasetID, scalingID) %>%
      summarise(
        Total_Ncells_samp = n_distinct(siteID),
        Total_area_samp = sum(croppedArea, na.rm = TRUE))
  )
  return(grid_info)
}

#--------------------------------------------------#

# Apply function to compute grid information
atlas_areas <-
  calculate_grid_info(data_sf)
#--------------------------------------------------#


# Read processed species data & merge with atlas_areas
pres_dat_final <- readRDS(here("Data/output/1_data/1_data_filtered.rds")) %>%
  left_join(atlas_areas) %>%
  mutate(scalingID = factor(scalingID, levels = vars$desired_levels)) %>%
  filter(!is.na(verbatimIdentification)) # Remove unsampled cells



#----------------------------------------------------------#

# Overview
atlas_areas %>%
  kable()

# Summary statistics
atlas_areas %>%
  skim()

# Save results
write.csv(atlas_areas,
          paste0(vars$Documentation, "atlas_areas_METADATA.csv"))


#----------------------------------------------------------#
# Calculate jaccard index   -----
#----------------------------------------------------------#

# Custom function to calculate jaccard across sites
jaccard <- function(set1, set2) {
  a <- length(intersect(set1, set2))
  b <- length(setdiff(set2, set1))
  c <- length(setdiff(set1, set2))
  n_union  <- length(union(set1,set2))
  return(list(jaccard_index = a / (a + b + c), a = a, b = b, c = c, n_union = n_union))
}


#----------------------------------------------------------#


# Calculate jaccard for all datasets
Jaccard_df <- pres_dat_final %>%
  filter(scalingID == 1) %>%
  group_by(datasetID, verbatimIdentification) %>%
  mutate(
    jaccard_results = list(jaccard(
      filter(pick(everything()), samplingPeriodID == 1)$siteID,
      filter(pick(everything()), samplingPeriodID == 2)$siteID
    ))
  ) %>%
  unnest_wider(jaccard_results) %>%
  group_by(verbatimIdentification, datasetID) %>%
  mutate(d = Total_Ncells_samp - n_union) %>%
  rename(Jaccard_sim = jaccard_index) %>%
  mutate(Jaccard_dissim = 1 - Jaccard_sim) %>%
  ungroup() %>%
  distinct(datasetID, verbatimIdentification, Jaccard_dissim, Jaccard_sim, a, b, c, d, Total_Ncells_samp)

# Save jaccard results
write.csv(Jaccard_df, paste0(vars$Documentation, "Jaccard_df.csv"))


#----------------------------------------------------------#

# Merge:
pres_dat_final_v2 <- pres_dat_final %>%
  left_join(Jaccard_df)

head(pres_dat_final_v2)

#----------------------------------------------------------#

rm(pres_dat_final, Jaccard_df, jaccard, atlas_areas, calculate_grid_info)


#----------------------------------------------------------#
# Calculate area of occupancy   -----
#----------------------------------------------------------#

occ_data_final <-
  pres_dat_final_v2 %>%
  ungroup() %>%
  distinct(datasetID, samplingPeriodID, scalingID,
           verbatimIdentification, siteID,
           .keep_all = TRUE) %>%
  group_by(datasetID, samplingPeriodID, scalingID, verbatimIdentification) %>%
  # Calculate AOO:
  dplyr::summarise(
    mean_area = mean(croppedArea, na.rm = TRUE),
    AOO = sum(croppedArea, na.rm = TRUE),
    occ_Ncells = n_distinct(siteID)) %>%
  left_join(pres_dat_final_v2 %>%
              mutate(scalingID = as.factor(as.character(scalingID)))) %>%
  # Calculate relative Occupancy:
  dplyr::mutate(
    rel_AOO = AOO / Total_area_samp,
    rel_occ_Ncells = occ_Ncells / Total_Ncells_samp) %>% # Prevalence
  # Remove duplicated rows:
  distinct()

#----------------------------------------------------------#

# Save AOO data to .rds
saveRDS(occ_data_final, here::here("Data", "output", "1_data", "1_occ_data_final.rds"))
rm(pres_dat_final_v2)

#----------------------------------------------------------#

# Create species-level data
species_data <-
  occ_data_final %>%
  dplyr::select(
    -siteID, -cell_sampling_repeats, -croppedArea, -sp_sampling_repeats
  ) %>%
  distinct(datasetID, samplingPeriodID, scalingID, verbatimIdentification,
           .keep_all = TRUE)

#----------------------------------------------------------#
# Calculate occupancy-area-relationship   -----
#----------------------------------------------------------#

# Custom function to calculate OAR and Fractal dimension
calculate_OAR <- function(species_data) {
  datasets <-
    unique(species_data$datasetID)
  sampling_periods <-
    unique(species_data$samplingPeriodID)

  sp_dta <-
    species_data %>%
    select(datasetID, samplingPeriodID, scalingID,
         verbatimIdentification,
         rel_AOO, rel_occ_Ncells, mean_area, AOO)

# Expand the grid for all combinations of datasetID and samplingPeriodID
species_data_new <-
  expand_grid(
  datasetID = datasets,
  samplingPeriodID = sampling_periods) %>%
  inner_join(sp_dta, by = c("datasetID", "samplingPeriodID")) %>%
  # Group by identifiers
  group_by(datasetID, samplingPeriodID, verbatimIdentification) %>%
  # Get available scales where relative occupancy is not saturated (< 1)
  summarise(
    available_scales = n(),
    mean_relAOO = mean(rel_AOO, na.rm = TRUE),
    exclude_sp_OAR = if_else(available_scales < 2, 1, 0),
    .groups = "drop"
  ) %>%
  # remove those where the range is saturated at the smallest resolution
  mutate(
    exclude_sp_OAR = if_else(available_scales == 0, 1, exclude_sp_OAR),
    mean_relAOO = if_else(available_scales == 0, 1, mean_relAOO)
  ) %>%
  full_join(sp_dta, by = c("datasetID", "samplingPeriodID", "verbatimIdentification")) %>%
  filter(
    exclude_sp_OAR == 0,
    rel_occ_Ncells < 1
  ) %>%
  distinct() %>%
  filter_at(vars(scalingID, AOO, mean_area), any_vars(!is.na(.))) %>%
  ungroup()

# Fit models using purrr::map
species_data_new_v2 <-
  species_data_new %>%
  group_by(datasetID, samplingPeriodID, verbatimIdentification) %>%
  nest(data = c(scalingID, AOO, mean_area, rel_AOO, rel_occ_Ncells)) %>%
  mutate(
    coefficients = map(
      data,
      ~ .x %>%
        filter(
          !is.na(AOO) & AOO > 0,
          !is.na(mean_area) & mean_area > 0
        ) %>%
        lm(log(AOO) ~ log(mean_area), data = .) %>%
        coef() %>%
        {
          tibble(
            m_AOO_a = .[2],
            b_AOO_a = .[1],
            D_AOO_a = -2 * .[2] + 2
          )
        }
    )
  ) %>%
  unnest(coefficients) %>%
  ungroup() %>%
  # Final cleanup of results
  select(datasetID, samplingPeriodID, verbatimIdentification, m_AOO_a, b_AOO_a, D_AOO_a) %>%
  distinct() %>%
  # Merge with the original dataset
  full_join(species_data) %>%
  filter(scalingID == 1) %>%
  select(
    datasetID, verbatimIdentification, samplingPeriodID, Total_area_samp,
    Total_Ncells, Total_Ncells_samp, AOO, occ_Ncells, rel_occ_Ncells,
    rel_AOO, Jaccard_dissim, m_AOO_a, b_AOO_a, D_AOO_a
  )

return(species_data_new_v2)
}


#----------------------------------------------------------#

# Apply OAR-function:
species_data_new <- calculate_OAR(species_data) %>%
  distinct(datasetID, samplingPeriodID, verbatimIdentification, .keep_all = TRUE)   %>%
  mutate(
    D_AOO_a = case_when(is.na(D_AOO_a) ~ 2,
                        .default = D_AOO_a)) %>%
      left_join(species_data %>% filter(scalingID == 1)) %>%
      # reduce columns
      select(
        datasetID, verbatimIdentification, samplingPeriodID,
        Total_area_samp, Total_Ncells, Total_Ncells_samp,
        AOO, occ_Ncells, rel_occ_Ncells, rel_AOO,
        Jaccard_dissim, a,b,c,d, D_AOO_a
      ) %>%
      distinct(datasetID, samplingPeriodID, verbatimIdentification,
               .keep_all = TRUE)


#----------------------------------------------------------#
# Calculate log ratio AOO   -----
#----------------------------------------------------------#

time_periods <- vars$time_periods

# Custom function to transform data to wide-format
transform_to_wide <- function(species_data_new, time_periods = c(1, 2)) {
  # Create a list to store wide data for each time period
  wide_dfs <- list()

  for (i in seq_along(time_periods)) {
    wide_dfs[[i]] <- species_data_new %>%
      distinct(datasetID, samplingPeriodID, verbatimIdentification, AOO) %>%
      group_by(datasetID, samplingPeriodID, verbatimIdentification) %>%
      filter(samplingPeriodID == time_periods[i]) %>%
      setNames(paste0("samplingPeriodID", i, "_", names(.))) %>%
      ungroup() %>%
      select(-c(paste0("samplingPeriodID", i, "_samplingPeriodID"))) %>%
      dplyr::rename(
        verbatimIdentification = paste0("samplingPeriodID", i, "_verbatimIdentification"),
        datasetID = paste0("samplingPeriodID", i, "_datasetID")
      )
  }

  # Merge the wide data frames sequentially
  sp_dat_wide <- reduce(wide_dfs, full_join, by = c("verbatimIdentification", "datasetID"))

  cat("NA counts in wide data after processing:\n")
  print(colSums(is.na(sp_dat_wide)))

  cat("Preview of wide data:\n")
  print(head(sp_dat_wide))

  return(sp_dat_wide)
}


#----------------------------------------------------------#


# Apply function:
sp_dat_wide <-
  transform_to_wide(species_data_new, time_periods) %>%
  na.omit() # drop species lost or gained completely

#----------------------------------------------------------#

# Calculate log ratio of AOO
logRatio <-
  sp_dat_wide %>%
  mutate(
    log_R2_1 = log(samplingPeriodID2_AOO / samplingPeriodID1_AOO),
    ratio_R2_1 = samplingPeriodID2_AOO / samplingPeriodID1_AOO) %>%
  select(-samplingPeriodID1_AOO, -samplingPeriodID2_AOO)

# save final predictor table from grids/atlases
big_table <-
  full_join(species_data_new, logRatio) %>%
  distinct(datasetID, samplingPeriodID, verbatimIdentification,
           .keep_all = T) %>%
  mutate_if(is.numeric, round, 3)


#----------------------------------------------------------#
# Save results to .rds  -----
#----------------------------------------------------------#

saveRDS(species_data_new, here("Data/output/1_data/1_species_data.rds"))
saveRDS(big_table, here("Data/output/1_data/2_big_table.rds"))
