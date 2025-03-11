#----------------------------------------------------------#
#
# Optimized 12_Predictors_lacunarity_NY.R without Grid Expansion
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

# Load required packages
list_p <- c("here", "sf", "tidyverse", "purrr", "tictoc", "terra")
lapply(list_p, require, character.only = TRUE)
sf_use_s2(FALSE)

data_id <- 6

#----------------------------------------------------------#
# Load and Filter Data -----
#----------------------------------------------------------#

# Read spatial data
data_sf <- readRDS("Data/output/1_data/1_data_sf.rds") %>%
  filter(scalingID == 1, datasetID == data_id)

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

# # Extract unique site-sampling period combinations
# sites <- data_sf_with_unsampled %>%
#   st_drop_geometry() %>%
#   distinct(datasetID, samplingPeriodID, siteID, cell_sampling_repeats)

# Extract all species per dataset and site
all_species <- data_sf_with_unsampled %>%
  st_drop_geometry() %>%
  distinct(datasetID, siteID, samplingPeriodID, verbatimIdentification) %>%
  mutate(presence = if_else(is.na(verbatimIdentification), NA_real_, 1))

#----------------------------------------------------------#
# Create Raster Templates -----
#----------------------------------------------------------#

# Compute resolutions efficiently
calculate_resolution <- function(sf_obj) {
  sf_obj %>%
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

plot(masked_templates[[1]]) # there is only a single element since we are doing datasets separately

#----------------------------------------------------------#
# Memory-Efficient & Organized Raster Saving -----
#----------------------------------------------------------#

# Directory for output rasters
output_dir <- "Data/input/species_ranges_tiff/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Select raster template (Assuming one dataset)
template_i <- masked_templates[[1]]
species_list <- unique(data_sf_with_unsampled$verbatimIdentification) %>%
  na.omit()

# **Loop through species and periods, saving each as a separate raster**
for (tp_i in seq_along(c(1, 2))) {
  for (sp_i in species_list) {
    # Extract current data
    current_sf <- current_sf <- data_sf_with_unsampled %>%
      filter(samplingPeriodID == tp_i & verbatimIdentification == sp_i) %>%
      st_as_sf() %>%
      mutate(presence = 1)

    # Extract names for the file
    dataset_id <- unique(current_sf$datasetID)
    sampling_period <- unique(current_sf$samplingPeriodID)
    species_name <- unique(current_sf$verbatimIdentification)

    # **Ensure valid filenames (remove special characters)**
    safe_species_name <- gsub("[^A-Za-z0-9_-]", "_", species_name)
    filename <- paste0(output_dir, dataset_id, "_", sampling_period, "_", safe_species_name, ".tif")
    print(filename)
    # **Rasterize and save**
    r <- rasterize(current_sf, template_i, field = "presence", update = TRUE)

    # Save the raster
    writeRaster(r, filename = filename, overwrite = TRUE)

    print(paste("Saved:", filename)) # Log progress
  }
}

#----------------------------------------------------------#




