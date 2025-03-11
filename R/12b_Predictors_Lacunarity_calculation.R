#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#                      12b_Predictors_Lacunarity_calculation.R
#                
#
#                    Friederike Wölke 
#                        2025
#
#----------------------------------------------------------#



#----------------------------------------------------------#
# Install and load libraries -----
#----------------------------------------------------------#
rm(list=ls())


# Load necessary packages
library(here)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(furrr)  # Parallel processing for rasterization
library(terra)


# Source 00_Configuration.R
source(here::here("R/00_Configuration.R"))

# Source lacunarity function for matrices:
source(here("anxiliary_code/fun_lacunarity_matrix.R"))

# Load required packages
list_p <- c("here", "sf", "tidyverse", "purrr", "furrr", "tictoc", "terra")
lapply(list_p, require, character.only = TRUE)
sf_use_s2(FALSE)


# Read spatial data
all_sf <- readRDS("Data/output/1_data/1_data_sf.rds") %>%
  filter(scalingID == 1)


#----------------------------------------------------------#
# Load rasters -----
#----------------------------------------------------------#

pattern <- "^5_"
# Function to Load and Process Rasters
load_rasters <- function(pattern) {
  files <- list.files(pattern = pattern, here("Data/input/species_ranges_tiff/"), full.names = TRUE)
  rast_obj <- rast(files)
  names(rast_obj) <- gsub(pattern = ".tif", "", basename(files))
  return(rast_obj)
}

# Load datasets
ranges_cz <- load_rasters("^5_")
ranges_ny <- load_rasters("^6_")
ranges_jp <- load_rasters("^13_")
ranges_eu <- load_rasters("^26_")

# Store all datasets in a list
ranges_list <- list(ranges_cz, ranges_ny, ranges_jp, ranges_eu)


#----------------------------------------------------------#
# Example use: lacunarty_R() -------
#----------------------------------------------------------#

r_vec <- c(3, 5, 9, 17, 33, 65)

## parallel:
library(terra)
library(furrr)
library(tictoc)

plan(sequential)
# Set up parallel processing with increased memory limit
options(future.globals.maxSize = 3 * 1024^3)  # 3 GB limit
plan(multisession, workers = 4)

library(furrr)
library(tictoc)
library(purrr)
library(dplyr)

# **Preallocate results list with correct length**
list_lac_res <- vector("list", length(ranges_list))

# Define r_vec outside the loop (avoids repeated allocation)
r_vec <- c(3, 5, 9, 17, 33, 65)


# **Loop through datasets for automation**
for (data_id in seq_along(ranges_list)) {
  print(paste("Processing dataset", data_id))

  # **Extract current dataset**
  current_data <- ranges_list[[data_id]]  # Extract single dataset
  current_ranges_list <- as.list(current_data)  # Convert to list of single-layer rasters
  ranges_matrices <- map(current_ranges_list, ~ as.matrix(.x, wide = TRUE))  # Convert each to matrix
  species_names <- names(current_data)  # Extract species names

  # **Run Parallel Lacunarity Computation**
  tic("Parallel Lacunarity Computation")
  res_list <- future_map2(
    ranges_matrices, species_names,  # Pass one matrix & name at a time
    ~ {
      species_label <- .y  # Explicitly extract the correct species name

      # Run lacunarity function
      res <- lacunarity_R(.x, r_vec, progress = FALSE, ncores = 1L, save_plot = FALSE, plot = FALSE)

      # Add species name as a column
      res <- res %>%
        mutate(species_name = species_label)  # Assign the correct species name


      return(res)
    },
    .progress = FALSE  # No progress bar for cleaner output
  )
  toc()

  # **Bind results into a single dataframe**
  data_lacunarity <- bind_rows(res_list, .id = "species_name")

  saveRDS(res_list,
          paste0(here("Data/output/1_data/single_predictors/2_lacunarity_"),
                 paste0(data_id, ".rds")))

  # **Store results in preallocated list**
  list_lac_res[[data_id]] <- data_lacunarity
}

# Print summary
print("All datasets processed successfully!")


saveRDS(results_list, "Data/output/1_data/single_predictors/3_lacunarity_CZ.rds")

















### Without loop: example for CZ

# **Convert Each Layer to a Matrix Before Parallel Processing**
ranges_cz_list <- as.list(ranges_cz)  # Converts multi-layer SpatRaster to list of single-layer SpatRasters
ranges_cz_matrices <- map(ranges_cz_list, ~ as.matrix(.x, wide = TRUE))  # Convert each layer to matrix

#plan(sequential)

# **Prepare Data for Parallel Execution**
species_names <- names(ranges_cz)  # Extract species names

# **Run Parallel Lacunarity Computation**
tic("Parallel Lacunarity Computation")
res_list <- future_map2(
  ranges_cz_matrices, species_names,  # Pass one matrix & name at a time
  ~ {
    # Run lacunarity function on the matrix
    res <- lacunarity_R(.x, r_vec, r_max, progress = FALSE, ncores = 1L, save_plot = FALSE, plot = FALSE)

    # Add species name as a column
    res <- res %>%
      mutate(name = .y)

    return(res)
  },
  .progress = TRUE  # Show progress bar
)
toc()

# close parallel session
plan(sequential)

# Print summary
print(res_list)


# Assign names properly
names(res_list) <- names(ranges_cz)  # Preserve original names of raster layers

data_lacunarity <- res_list %>% bind_rows()
# Print summary
print(data_lacunarity)

saveRDS(data_lacunarity, here("Data/output/1_data/3_lacunarity_CZ.rds"))

#----------------------------------------------------------#
