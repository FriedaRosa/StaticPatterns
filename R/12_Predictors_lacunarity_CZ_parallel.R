#----------------------------------------------------------#
#
#  12_Predictors_lacunarity_CZ_parallel.R
#
#----------------------------------------------------------#

# Start with clean environment
rm(list = ls())
gc()

#----------------------------------------------------------#
# Load required libraries -----
#----------------------------------------------------------#

# Load necessary packages
library(here)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(furrr)  # Parallel processing for rasterization
library(terra)

sf_use_s2(FALSE)  # Avoid unnecessary computations

# Set parallelization plan (adjust worker count based on system)
plan(multisession, workers = 4)

data_id <- 5

#----------------------------------------------------------#
# Load data -----
#----------------------------------------------------------#

# Read spatial data
data_sf <- readRDS("Data/output/1_data/1_data_sf.rds") %>%
  filter(scalingID == 1, datasetID == data_id)

# Read grid
sf_grid <- readRDS(here("Data/input/grid.rds")) %>%
  filter(scalingID == 1, datasetID == data_id)

#----------------------------------------------------------#
# Prepare species data -----
#----------------------------------------------------------#

# Handle unsampled sites
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

# Extract unique sites
sites <- data_sf_with_unsampled %>%
  st_drop_geometry() %>%
  distinct(datasetID, samplingPeriodID, siteID, cell_sampling_repeats)

# Extract all species per dataset and site
all_species <- data_sf_with_unsampled %>%
  st_drop_geometry() %>%
  distinct(datasetID, siteID, samplingPeriodID, verbatimIdentification) %>%
  mutate(
    presence = if_else(is.na(verbatimIdentification), NA_real_, 1)
  )

#----------------------------------------------------------#
# Expand Combinations (Sequential for Stability) -----
#----------------------------------------------------------#

expand_combinations_chunked <- function(data_chunk, species_data) {
  data_sites <- data_chunk %>%
    st_drop_geometry() %>%
    distinct(datasetID, samplingPeriodID, siteID, cell_sampling_repeats)

  data_species <- species_data %>%
    filter(
      datasetID == unique(data_chunk$datasetID),
      samplingPeriodID == unique(data_chunk$samplingPeriodID)
    ) %>%
    mutate(presence = if_else(is.na(verbatimIdentification), "no_data", "1"))

  expanded_chunk <- data_sites %>%
    mutate(cell_sampling_repeats = replace_na(cell_sampling_repeats, 0)) %>%
    group_split(siteID) %>%
    map_dfr(~ tidyr::expand_grid(.x, verbatimIdentification = data_species$verbatimIdentification)) %>%
    left_join(data_species) %>%
    mutate(
      presence = case_when(
        cell_sampling_repeats == 0 ~ "no_data",
        is.na(presence) ~ "0",
        TRUE ~ as.character(presence)
      ),
      presence = as.numeric(if_else(presence == "no_data", NA_character_, presence))
    )

  return(expanded_chunk)
}

# Run function sequentially for stability
expanded <- data_sf_with_unsampled %>%
  group_by(datasetID, samplingPeriodID) %>%
  group_split() %>%
  map(~ expand_combinations_chunked(.x, all_species))

expanded_sf <- bind_rows(expanded)

# Re-attach geometry
expanded_sf <- expanded_sf %>%
  right_join(
    data_sf_with_unsampled %>%
      distinct(datasetID, siteID, samplingPeriodID, geometry)
  )

# Split into list by dataset, species, and sampling period
species_list <- expanded_sf %>%
  group_split(datasetID, samplingPeriodID, verbatimIdentification)

names(species_list) <- expanded_sf %>%
  distinct(datasetID, samplingPeriodID, verbatimIdentification) %>%
  transmute(name = paste(datasetID, samplingPeriodID, verbatimIdentification, sep = "_")) %>%
  pull(name)

#----------------------------------------------------------#
# Create raster templates -----
#----------------------------------------------------------#

calculate_resolution <- function(sf_obj) {
  sf_obj %>%
    group_split(datasetID) %>%
    map(function(dataset) {
      bbox <- st_bbox(dataset)
      x_coords <- st_coordinates(dataset)[, 1]
      y_coords <- st_coordinates(dataset)[, 2]
      x_res <- abs(bbox$xmax - bbox$xmin) / (length(unique(x_coords)) - 1)
      y_res <- abs(bbox$ymax - bbox$ymin) / (length(unique(y_coords)) - 1)

      list(
        datasetID = unique(dataset$datasetID),
        resolution = c(x_res, y_res),
        bbox = bbox
      )
    })
}

resolutions <- sf_grid %>% calculate_resolution()

# Create templates
template_list <- lapply(resolutions, function(res) {
  bbox <- res$bbox
  rast(ext(bbox$xmin, bbox$xmax, bbox$ymin, bbox$ymax), resolution = res$resolution, crs = crs(sf_grid))
})

# Create masked templates
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
# Parallelized Raster Processing (Only This Step) -----
#----------------------------------------------------------#

rasterize_species <- function(species_sf, template) {
  sf_data <- st_as_sf(species_sf)
  rasterize(sf_data, template, field = "presence")
}

# Run rasterization in parallel (ONLY HERE)
masked_rasters <- future_map2(
  species_list, list(masked_templates[[1]]),
  rasterize_species,
  .progress = TRUE
)

# Convert to raster stack
masked_raster_stack <- rast(masked_rasters)

#----------------------------------------------------------#
# Save the raster stack -----
#----------------------------------------------------------#

writeRaster(masked_raster_stack,
            filename = "Data/output/CZ_species_raster_stack.tif",
            overwrite = TRUE)

#----------------------------------------------------------#
