#----------------------------------------------------------#
#
#
#                     Static Patterns
#
#                   12_Predictors_lacunarity.R
#                
#
#                     Friederike Wölke
#                          2025
#
#----------------------------------------------------------#


# Start with clean environment
rm(list = ls())
gc()


#----------------------------------------------------------#
# Install and load libraries -----
#----------------------------------------------------------#

# Source 00_Configuration.R
source(here::here("R/00_Configuration.R"))

# subset of packages to load to avoid conflicts:
list_p <- c("here", "sf" ,"tidyverse", "purrr", "tictoc", "terra")
lapply(list_p, require, character = TRUE)
sf_use_s2(FALSE)

#----------------------------------------------------------#
# Custom functions -----
#----------------------------------------------------------#

# Function to expand combinations for a subset of data
expand_combinations_chunked <- function(data_chunk, species_data) {
  # Extract unique site and species combinations for the current chunk
  data_sites <- data_chunk %>%
    st_drop_geometry() %>%
    dplyr::select(datasetID, samplingPeriodID, siteID, cell_sampling_repeats) %>%
    distinct()

  print(paste0("sites_data, datasetID = ", unique(data_sites$datasetID)))
  print(paste0("sites_data, samplingPeriodID = ", unique(data_sites$samplingPeriodID)))

  data_species <- species_data %>%
    filter(datasetID == unique(data_chunk$datasetID)) %>%
    filter(samplingPeriodID == unique(data_chunk$samplingPeriodID)) %>%
    mutate(presence = as.character(presence)) %>%
    mutate(presence = if_else(is.na(verbatimIdentification), "no_data", "1"))

  print(paste0("species_data, datasetID = ", unique(data_species$datasetID)))
  print(paste0("species_data, samplingPeriodID = ", unique(data_species$samplingPeriodID)))

  # Expand combinations within this chunk
  expanded_chunk <- data_sites %>%
    mutate(
      cell_sampling_repeats = case_when(is.na(cell_sampling_repeats) ~ 0,
                                        .default = cell_sampling_repeats)) %>%
    group_split(siteID) %>% # Process each site separately
    purrr::map_dfr(~ tidyr::expand_grid(
      .x,
      verbatimIdentification = data_species$verbatimIdentification
    )) %>%
    left_join(data_species) %>%
    mutate(presence = case_when(cell_sampling_repeats == 0 ~ "no_data",
                                is.na(presence) ~ "0",
                                .default = as.character(presence)
    )) %>%
    mutate(
      presence = as.numeric(
        if_else(presence == "no_data", NA, presence)))

  print(paste0("expanded presences = ", unique(expanded_chunk$presence)))



  return(expanded_chunk)
}

#----------------------------------------------------------#
# Load data -----
#----------------------------------------------------------#

# Read data and filter
data_sf <-
  readRDS("Data/output/1_data/1_data_sf.rds") %>%
  ungroup() %>%
  filter(scalingID == 1)

#----------------------------------------------------------#
# 1. Prepare species data for rasterizing -----
#----------------------------------------------------------#

# Handle unsampled sites
unsampled_sites_expanded <- data_sf %>%
  filter(cell_sampling_repeats == 0) %>%
  distinct(datasetID, siteID, geometry) %>%
  tidyr::expand_grid(samplingPeriodID = c(1, 2)) %>%
  mutate(verbatimIdentification = NA)

# Combine with sampled data
data_sf_with_unsampled <- data_sf %>%
  filter(cell_sampling_repeats == 2 | !is.na(verbatimIdentification)) %>%
  bind_rows(unsampled_sites_expanded)

#----------------------------------------------------------#

# Extract unique site and sampling period combinations
sites <- data_sf_with_unsampled %>%
  st_drop_geometry() %>%
  distinct(datasetID, samplingPeriodID, siteID, cell_sampling_repeats)

# Extract all species per dataset and site
all_species <- data_sf_with_unsampled %>%
  st_drop_geometry() %>%
  distinct(datasetID, siteID, samplingPeriodID, verbatimIdentification) %>%
  mutate(
    presence = if_else(is.na(verbatimIdentification), NA_real_, 1)
    # Assign presence = NA for unsampled sites
  )

#----------------------------------------------------------#

# run expand-function for dataset~samplingperiodID (i.e., split data)
expanded <- data_sf_with_unsampled %>%
  group_by(datasetID, samplingPeriodID) %>%
  group_split() %>%
  purrr::map(~ expand_combinations_chunked(.x, all_species))

# Re-attach geometry from the original data
expanded_sf <- expanded %>%
  right_join(
    data_sf_with_unsampled %>%
      distinct(datasetID, siteID, geometry))

# Split into list by dataset, species, and sampling period
species_list <- expanded_sf %>%
  group_split(datasetID, samplingPeriodID, verbatimIdentification)

names(species_list) <- expanded_sf %>%
  distinct(datasetID, samplingPeriodID, verbatimIdentification) %>%
  transmute(name = paste(datasetID, samplingPeriodID, verbatimIdentification, sep = "_")) %>%
  pull(name)


#----------------------------------------------------------#
# 2. Convert range maps from sf to raster -----
#----------------------------------------------------------#

## step 1: convert species distribution sf objects to raster using terra

templates_list <- list()
sf_grid <- readRDS(here("Data/input/grid.rds")) %>%
  filter(scalingID == 1)

list_ext <- sf_grid %>%
  group_split(datasetID) %>%
  map(~ ext(.))

my_crs <- crs(sf_grid)

calculate_resolution <- function(sf_obj) {
  sf_obj %>%
    group_split(datasetID) %>%
    map(function(dataset) {
      bbox <- st_bbox(dataset)

      # Infer resolution based on number of unique x and y coordinates
      x_coords <- dataset %>% st_coordinates() %>% .[, 1]
      y_coords <- dataset %>% st_coordinates() %>% .[, 2]
      x_res <- abs(bbox$xmax - bbox$xmin) / (length(unique(x_coords)) - 1)
      y_res <- abs(bbox$ymax - bbox$ymin) / (length(unique(y_coords)) - 1)

      list(
        datasetID = unique(dataset$datasetID),
        resolution = c(x_res, y_res),
        bbox = bbox
      )
    })
}


# Calculate resolution from the sf object
resolutions <- sf_grid %>% calculate_resolution()

# make templates
template_list <- lapply(resolutions, function(res) {
  bbox <- res$bbox
  ext <- ext(bbox$xmin, bbox$xmax, bbox$ymin, bbox$ymax)
  # Create a SpatRaster with the calculated resolution
  rast(ext, resolution = res$resolution, crs = st_crs(sf_grid)$wkt)
})

# Inspect the SpatRaster templates
template_list

# Create layers with values = 0 within the mask
masked_templates <- lapply(seq_along(template_list), function(i) {
  # Start with the template
  template <- template_list[[i]]

  # Populate all cells with 0
  values(template) <- 0

  # Create a mask using the dataset's boundaries
  country_boundary <- sf_grid %>%
    filter(datasetID == unique(sf_grid$datasetID)[i]) %>%
    st_union() %>%
    vect() # Convert sf object to terra-compatible SpatVector

  # Apply the mask (keep NA outside the boundary)
  masked_layer <- mask(template, country_boundary)

  return(masked_layer)
})

# Inspect the masked layers
masked_templates
crs(masked_templates[[1]]) == my_crs





# Convert each group to raster and calculate lacunarity
masked_rasters <- list()
res_list <- list()

for (i in seq_along(species_list)) {
  current_sf <- species_list[[i]]
  atlas_i <- current_sf$datasetID[1]

  # Convert sf to raster
  r <- rasterize(current_sf, template_CZ, field = "presence")
  r_masked <- mask(r, as_Spatial(country_boundary))

  # Calculate lacunarity
  lacunarity_results <- calculate_lacunarity(
    as.matrix(r_masked),
    box_sizes
  )

  # Store results
  res <- data.frame(
    samplingPeriodID = current_sf$samplingPeriodID[1],
    verbatimIdentification = current_sf$verbatimIdentification[1],
    lacunarity_results
  )
  res_list[[i]] <- res

  # Store masked raster
  masked_rasters[[paste(current_sf$samplingPeriodID[1], current_sf$verbatimIdentification[1], sep = "_")]] <- r_masked
}

# Combine results
res_df <- bind_rows(res_list)

# Combine masked rasters into a RasterStack
masked_raster_stack <- stack(masked_rasters)
