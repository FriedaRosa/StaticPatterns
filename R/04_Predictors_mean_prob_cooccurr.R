#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#                      04_Predictors_mean_prob_cooccurr.R
#                
#
#                    Friederike Wöle 
#                        2025
#
#----------------------------------------------------------#

source(here::here("R/00_Configuration.R"))
lapply(package_list, require, character = TRUE)

library(fossil)
library(cooccur)

#----------------------------------------------------------#
# Function -----
#----------------------------------------------------------#

process_co_occurrence <- function(species_data, atlas_names) {

  #browser()

  # Helper function to create species-by-site matrices
  create_comm_matrices <- function(data, time_values) {
    lapply(time_values, function(time_val) {
      fossil::create.matrix(
        data,
        tax.name = "verbatimIdentification",
        locality = "siteID",
        time.col = "samplingPeriodID",
        time = time_val,
        abund = FALSE
      )
    })
  }

  # Helper function to calculate co-occurrence metrics
  calculate_cooccurrence <- function(comm_matrix) {
    cooccur::cooccur(comm_matrix, spp_names = TRUE)$results %>%
      group_by(sp1_name) %>%
      dplyr::summarise(
        mean_prob_cooccur = mean(prob_cooccur, na.rm = TRUE) # mean per species
      )
  }

  # Process each atlas and compute co-occurrence results
  co_occ_list <- lapply(atlas_names, function(atlas_name) {

    #browser()

    # Filter data for the current atlas
    comm_dat <- species_data %>%
      filter(datasetID == atlas_name) %>%
      distinct(samplingPeriodID, verbatimIdentification, siteID) %>%
      as.data.frame()

    # Create community matrices for time periods
    comm_matrices <- create_comm_matrices(comm_dat, time_values = c("1", "2"))

    # Calculate co-occurrence results for both time periods
    co_occ_results_l <- lapply(
      seq_along(comm_matrices), function(tp_idx) {

        calculate_cooccurrence(comm_matrices[[tp_idx]]) %>%
          mutate(samplingPeriodID = as.integer(tp_idx),
                 datasetID = atlas_name)

      })

    saveRDS(co_occ_results_l, here::here(paste0("Data/output/1_data/single_predictors/co_occ_res/co_occ_res_list_", atlas_name, ".rds" )))

    co_occ_results <- co_occ_results_l %>%
      bind_rows() %>%
      dplyr::rename(verbatimIdentification = sp1_name)

    co_occ_results$datasetID <- as.integer(as.character(co_occ_results$datasetID))

    # Combine with original species data
    species_data %>%
      filter(datasetID == atlas_name) %>%
      distinct() %>%
      ungroup() %>%
      left_join(co_occ_results)

  })

  # Combine results for all atlases into a single dataframe
  final_df <- bind_rows(co_occ_list)
  return(final_df)

}

#----------------------------------------------------------#
# Load data  -----
#----------------------------------------------------------#

# Load data
species_data <- readRDS(here::here("Data/output/1_data/1_occ_data_final.rds")) %>%
  filter(scalingID == 1) %>%
  distinct(datasetID, samplingPeriodID, verbatimIdentification, siteID)

#----------------------------------------------------------#
# Calculate co-occurrence metrics -----
#----------------------------------------------------------#

# Calculate co-occurrence metrics
mean_cooccurrence_df <- process_co_occurrence(
  species_data = species_data,
  atlas_names = vars$atlas_names)

#----------------------------------------------------------#
# Save results to .rds  -----
#----------------------------------------------------------#

saveRDS(mean_cooccurrence_df,
        here::here("Data/output/1_data/single_predictors/2_cooccurrence.rds"))
