#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#                      12a_Predictors_Lacunarity_rasterize_ranges.R
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

# Load required packages
list_p <- c("here", "sf", "tidyverse", "purrr", "furrr", "tictoc", "terra")
lapply(list_p, require, character.only = TRUE)
sf_use_s2(FALSE)


# Read spatial data
all_sf <- readRDS("Data/output/1_data/1_data_sf.rds") %>%
  filter(scalingID == 1)



# Start Loop for rasterizing sf objects :

for (data_id in unique(all_sf$datasetID)){
print(data_id)
data_id <- 26
#----------------------------------------------------------#
# Load and Filter Data -----
#----------------------------------------------------------#

# Read spatial data
data_sf <- all_sf %>%
  filter(datasetID == data_id) %>%
  st_as_sf()

# Read grid
sf_grid <- readRDS(here("Data/input/grid.rds")) %>%
  filter(scalingID == 1, datasetID == data_id)

#----------------------------------------------------------#
# Prepare Species Data -----
#----------------------------------------------------------#

# Handle unsampled sites (expand in a memory-efficient way)
unsampled_sites_expanded <- data_sf %>%
  filter(cell_sampling_repeats == 0) %>%
  distinct(datasetID, siteID, geometry) %>%
  tidyr::expand_grid(samplingPeriodID = c(1, 2)) %>%
  mutate(verbatimIdentification = NA)

# Combine with sampled data
data_sf_with_unsampled <- bind_rows(
  data_sf %>% filter(cell_sampling_repeats == 2 | !is.na(verbatimIdentification)),
  unsampled_sites_expanded
)

# Get unique species list (excluding NA)
species_list <- unique(data_sf_with_unsampled$verbatimIdentification) %>%
  na.omit()

#----------------------------------------------------------#
# Create Raster Templates -----
#----------------------------------------------------------#

# Compute resolutions efficiently
calculate_resolution <- function(sf_obj) {
  sf_obj %>%
    st_make_valid() %>%
    group_split(datasetID) %>%
    map(~ {
      bbox <- st_bbox(.x)
      list(
        datasetID = unique(.x$datasetID),
        resolution = c(
          abs(bbox$xmax - bbox$xmin) / (length(unique(st_coordinates(.x)[, 1])) - 1),
          abs(bbox$ymax - bbox$ymin) / (length(unique(st_coordinates(.x)[, 2])) - 1)
        ),
        bbox = bbox
      )
    })
}

#----------------------------------------------------------#

resolutions <- calculate_resolution(sf_grid)

# Create masked raster templates
template_list <- lapply(resolutions, function(res) {
  bbox <- res$bbox
  rast(ext(bbox$xmin, bbox$xmax, bbox$ymin, bbox$ymax), resolution = res$resolution, crs = st_crs(sf_grid)$wkt)
})

masked_templates <- lapply(seq_along(template_list), function(i) {
  template <- template_list[[i]]
  values(template) <- 0
  country_boundary <- sf_grid %>%
    filter(datasetID == unique(sf_grid$datasetID)[i]) %>%
    st_union() %>%
    vect()
  mask(template, country_boundary)
})

#----------------------------------------------------------#
# Parallelized Rasterization & Saving -----
#----------------------------------------------------------#

# Directory for output rasters
output_dir <- "Data/input/species_ranges_tiff/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Select raster template (Assuming one dataset)
template_i <- masked_templates[[1]]

# **Parallel Rasterization Function**
save_raster <- function(sampling_period, species_name) {
  # Extract current species-period data
  current_sf <- data_sf_with_unsampled %>%
    filter(samplingPeriodID == sampling_period & verbatimIdentification == species_name) %>%
    st_as_sf() %>%
    mutate(presence = 1)

  # Skip if there's no data for this combination
  if (nrow(current_sf) == 0) return(NULL)

  # **Generate valid filename**
  safe_species_name <- gsub("[^A-Za-z0-9_-]", "_", species_name)
  filename <- paste0(output_dir, data_id, "_", sampling_period, "_", safe_species_name, ".tif")

  # **Rasterize**
  r <- rasterize(current_sf, template_i, field = "presence", update = TRUE)

  # **Save raster**
  writeRaster(r, filename = filename, overwrite = TRUE)

  print(paste("Saved:", filename))
}

map2(
  rep(c(1, 2), each = length(species_list)),  # Expand time periods
  rep(species_list, times = 2),  # Expand species list
  save_raster,
  .progress = TRUE
)
}


#----------------------------------------------------------#
# Load rasters -----
#----------------------------------------------------------#


#----------------------------------------------------------#
files <- list.files(pattern = "5_",
                    here("Data/input/species_ranges_tiff/"),
                    full.names = T)
ranges_cz <- rast(files)
names(ranges_cz) <- gsub(pattern = ".tif", "", basename(files))
# plot(ranges_cz[[1]])
# ranges_cz[[1]] %>% as.matrix(wide =T)
#----------------------------------------------------------#
files <- list.files(pattern = "^6_",
                    path = here("Data/input/species_ranges_tiff/"),
                    full.names = TRUE)
ranges_ny <- rast(files)
names(ranges_ny) <- gsub(pattern = ".tif", "", basename(files))
plot(ranges_ny[[1]])

#----------------------------------------------------------#
files <- list.files(pattern = "13_",
                    here("Data/input/species_ranges_tiff/"),
                    full.names = T)
ranges_jp <- rast(files)
names(ranges_jp) <- gsub(pattern = ".tif", "", basename(files))
plot(ranges_jp[[1]])
#----------------------------------------------------------#
files <- list.files(pattern = "26_",
                    here("Data/input/species_ranges_tiff/"),
                    full.names = T)
ranges_eu <- rast(files)
names(ranges_eu) <- gsub(pattern = ".tif", "", basename(files))
plot(ranges_eu[[1]])

#----------------------------------------------------------#
ranges_list <- list(ranges_cz, ranges_ny, ranges_jp, ranges_eu)
#----------------------------------------------------------#

source(here("anxiliary_code/fun_lacunarity_matrix.R"))

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



# Try with loop through datasets for automatization

list_lac_res <- replicate(4, list())

for (data_id in seq_along(ranges_list)){
  print(data_id)

  current_data <- ranges_list[[data_id]]
  ranges_list <- as.list(current_data)
  ranges_matrices <- map(ranges_list, ~ as.matrix(.x, wide = TRUE))

  species_names <- names(current_data)

  # **Run Parallel Lacunarity Computation**
  tic("Parallel Lacunarity Computation")
  res_list <- future_map2(
    ranges_matrices, species_names,  # Pass one matrix & name at a time
    ~ {
      # Run lacunarity function on the matrix
      res <- lacunarity_R(.x, r_vec, r_max, progress = FALSE, ncores = 1L, save_plot = FALSE, plot = FALSE)

      # Add species name as a column
      res <- res %>%
        mutate(name = .y)

      return(res)
    },
    .progress = FALSE  # Show progress bar
  )
  toc()

  # Assign names properly
  names(res_list) <- species_names  # Preserve original names of raster layers

  # bind to dataframe
  data_lacunarity <- res_list %>% bind_rows()

  list_lac_res[[data_id]] <- data_lacunarity

}


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
