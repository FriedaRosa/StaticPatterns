#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#             12a_Predictors_Lacunarity_rasterize_ranges.R
#         Script to save species ranges (sf) as individual .tiff rasters
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

# Source 00_Configuration.R
source(here::here("R/00_Configuration.R"))

# Load required packages
list_p <-
  c("here", "sf", "tidyverse", "purrr", "furrr", "tictoc", "terra")
lapply(list_p, require, character.only = TRUE)
sf_use_s2(FALSE) # not lat long data

#----------------------------------------------------------#
# Function to calculate resolution -----
#----------------------------------------------------------#

calculate_resolution <- function(sf_obj) {
  sf_obj %>%
    st_make_valid() %>%
    st_cast("MULTIPOLYGON") %>%
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
# Function to rasterize ranges -----
#----------------------------------------------------------#

save_raster <- function(sampling_period, species_name) {
  # Extract current species-period data
  current_sf <- data_sf_with_unsampled %>%
    filter(samplingPeriodID == sampling_period & verbatimIdentification == species_name) %>%
    st_as_sf() %>%
    st_transform(crs = vars$crs[data_nr]) %>%
    mutate(presence = 1)

  # Skip if there's no data for this combination
  if (nrow(current_sf) == 0) return(NULL)

  # **Generate valid filename**
  safe_species_name <- gsub("[^A-Za-z0-9_-]", "_", species_name)
  filename <- paste0(output_dir, data_id, "_", sampling_period, "_", safe_species_name, ".tif")

  # **Rasterize**
  # #, datatype = "INT1U" --> save as integer to save memory
  r <- rasterize(current_sf, template_i, field = "presence", update = TRUE, datatype = "INT1U")

  # **Save raster**
  writeRaster(r, filename = filename, overwrite = TRUE)
  gc()
  rm(list = c("current_sf", "r", "safe_species_name"))
  print(paste("Saved:", filename))
}


#----------------------------------------------------------#
# Read spatial data and set variables for single atlas -----
#----------------------------------------------------------#

all_sf <-
  readRDS("Data/output/1_data/1_data_sf.rds") %>%
  filter(scalingID == 1)


# Read grid
grid_sf_all <-
  readRDS(here("Data/input/grid.rds")) %>%
  filter(scalingID == 1)

## Note: to not overflow the memory, set data_id manually and skip the outer loop
# data_id <- 6; data_nr <- 2
# data_id <- 5; data_nr <- 1
# data_id <- 13; data_nr <- 3
data_id <- 26; data_nr <- 4

#----------------------------------------------------------#
# Filter Data for one dataset -----
#----------------------------------------------------------#
for (data_id in unique(all_sf$datasetID)){
  if (data_id == 5) data_nr <- 1 else
    if (data_id == 6) data_nr <- 2 else
      if (data_id == 13) data_nr <- 3 else
        if (data_id == 26) data_nr <- 4 else
          stop("data_id not found")
  if(data_id == 26) sf_use_s2(FALSE)

  # Filter spatial data to dataset
data_sf <-
  all_sf %>%
  filter(datasetID == data_id) %>%
  st_as_sf() %>%
  st_transform(crs = vars$crs[data_nr])

# Filter grid to dataset
grid_sf <-
  grid_sf_all %>%
  filter(datasetID == data_id) %>%
  st_transform(crs = vars$crs[data_nr])

#----------------------------------------------------------#
# Create Raster Templates -----
#----------------------------------------------------------#

resolutions <-
  calculate_resolution(grid_sf)

# Create masked raster templates
template_list <- lapply(resolutions, function(res) {
  bbox <- res$bbox
  rast(ext(bbox$xmin, bbox$xmax, bbox$ymin, bbox$ymax),
       resolution = res$resolution,
       crs = st_crs(grid_sf)$wkt)
})

masked_templates <- lapply(seq_along(template_list), function(i) {
  template <- template_list[[i]]
  values(template) <- 0
  country_boundary <- grid_sf %>%
    filter(datasetID == unique(grid_sf$datasetID)[i]) %>%
    st_union() %>%
    vect()
  mask(template, country_boundary)
})

#----------------------------------------------------------#
# Prepare Species Data -----
#----------------------------------------------------------#

# Handle unsampled sites (expand in a memory-efficient way)
unsampled_sites_expanded <-
  data_sf %>%
  filter(cell_sampling_repeats == 0) %>%
  distinct(datasetID, siteID, geometry) %>%
  tidyr::expand_grid(samplingPeriodID = c(1, 2)) %>%
  mutate(verbatimIdentification = NA)

# Combine with sampled data
data_sf_with_unsampled <-
  bind_rows(
  data_sf %>% filter(cell_sampling_repeats == 2 | !is.na(verbatimIdentification)),
  unsampled_sites_expanded
)

# Get unique species list (excluding NA)
species_list <-
  unique(data_sf_with_unsampled$verbatimIdentification) %>%
  na.omit()


#----------------------------------------------------------#
# Rasterization & Saving -----
#----------------------------------------------------------#

# Directory for output rasters
output_dir <- "Data/input/species_ranges_tiff/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Select raster template (Assuming one dataset)
template_i <- masked_templates[[1]]


map2(
  rep(c(1, 2), each = length(species_list)),  # Expand time periods
  rep(species_list, times = 2),  # Expand species list
  save_raster,
  .progress = TRUE
)

rm(list = c("grid_sf", "data_sf", "resolutions", "template_list",
            "masked_templates", "unsampled_sites_expanded",
            "data_sf_with_unsampled", "species_list", "template_i", "output_dir"))
gc()


#----------------------------------------------------------#
}
