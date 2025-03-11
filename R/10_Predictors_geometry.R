#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#                   10_Predictors_geometry.R
#                
#
#              Friederike Wölke, Gabriel Ortega-Solís  
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
list_p <- c("here", "sf" ,"tidyverse", "purrr", "tictoc", "terra", "broom", "geosphere")
lapply(list_p, require, character = TRUE)


# Source poly-attribute function
source(here("R/src/poly_attr.R"))

poly_attr <- possibly(poly_attr, otherwise = NA_real_, quiet = FALSE)
#----------------------------------------------------------#

#----------------------------------------------------------#

sf_use_s2(FALSE) # we use planar projections so this can be turned off

#----------------------------------------------------------#
# Load data -----
#----------------------------------------------------------#

atlas_sf <-
  readRDS("Data/output/1_data/1_data_sf.rds") %>%
  filter(scalingID == 1 & cell_sampling_repeats == 2) %>%
  group_split(datasetID)

#----------------------------------------------------------#
# Atlas geometries  -----
#----------------------------------------------------------#

# Define dataset-specific CRS mapping
dataset_crs_map <-
  tibble(
    datasetID = c(5, 6, 13, 26),
    crs = vars$crs[1:4]  # Extract corresponding CRS values
  )

#----------------------------------------------------------#

# Define function for processing a single atlas
process_atlas <- function(atlas, crs_value) {

  # Filter and transform dataset
  sf_current_atlas <- atlas %>%
    group_by(datasetID) %>%
    select(datasetID, geometry) %>%
    st_as_sf() %>%
    summarise(geometry = st_union(geometry), .groups = "drop") %>%
    st_transform(crs = crs_value) %>%
    st_make_valid()

  atlasID <- unique(sf_current_atlas$datasetID)

  message("Processing atlas ID = ", atlasID)

  # Convert to terra SpatVector
  terra_current_atlas <- vect(sf_current_atlas)

  # Store border lines
  atlas_border_lines <- terra::as.lines(terra_current_atlas)

  # Bounding box metrics
  bbox <- terra::ext(terra_current_atlas)

  # Calculate predictors
  data_atlas_res <- data.frame(
    datasetID = atlasID,
    atlas_xmin = bbox[1],
    atlas_xmax = bbox[2],
    atlas_xhalf = bbox[1] + (bbox[2] - bbox[1]) / 2,
    atlas_ymin = bbox[3],
    atlas_ymax = bbox[4],
    atlas_yhalf = bbox[3] + (bbox[4] - bbox[3]) / 2,
    atlas_nsDist = poly_attr(terra_current_atlas, "nsDist"),
    atlas_ewDist = poly_attr(terra_current_atlas, "ewDist"),
    atlas_maxDist = poly_attr(terra_current_atlas, "maxDist"),
    atlas_lengthMinRect = poly_attr(terra_current_atlas, "lengthMinRect"),
    atlas_widthMinRect = poly_attr(terra_current_atlas, "widthMinRect"),
    atlas_elonMinRect = poly_attr(terra_current_atlas, "elonMinRect"),
    atlas_elonRatio = poly_attr(terra_current_atlas, "elonRatio"),
    atlas_circ = poly_attr(terra_current_atlas, "circ"),
    atlas_circNorm = poly_attr(terra_current_atlas, "circNorm"),
    atlas_relCirc = poly_attr(terra_current_atlas, "relCirc"),
    atlas_lin = poly_attr(terra_current_atlas, "lin"),
    atlas_bearingMinRect = poly_attr(terra_current_atlas, "bearingMinRect"),
    atlas_bearing = poly_attr(terra_current_atlas, "bearing"),
    atlas_perimeter = terra_current_atlas %>% project("epsg:4326") %>% perim()
  )

  return(
    list(
      geometry_metrics = data_atlas_res,
      borders = atlas_border_lines))
}

#----------------------------------------------------------#

# Apply function to each atlas using map2()
results <- map2(atlas_sf, vars$crs, process_atlas)

#----------------------------------------------------------#

# Extract results into separate lists
data_atlas_geometry <- map_dfr(results, "geometry_metrics")

# Combine results into a single data frame
row.names(data_atlas_geometry) <- NULL

# Extract borders
list_borders <- map(results, "borders")


# works until here.
# Next step is to calculate range geometries for each species


#----------------------------------------------------------#
# Single species range geometries  -----
#----------------------------------------------------------#

# Not parallel #


list_range_geometries <-
  replicate(4, list())

#----------------------------------------------------------#

tictoc::tic()
for (atlas_i in seq_along(data_atlas_geometry$datasetID)){

  atlas <- data_atlas_geometry[[atlas_i]]
  border_lines <- list_borders[[atlas_i]]
  atlas_centroid <- centroids(border_lines)

  grouped_list <- atlas_sf[[atlas_i]] %>%
    st_as_sf() %>%
    group_by(samplingPeriodID, verbatimIdentification) %>%
    group_split()

  #----------------------------------------#
    range_geometries <- map_dfr(grouped_list, function(dta) {

    #message("Processing atlas nr = ", atlas_i)
    print(paste0("processing = ", unique(dta$verbatimIdentification),
                 ", datasetID = ", unique(dta$datasetID),
                 ", tp = ", unique(dta$samplingPeriodID) ))


    sp_range <- dta %>%
      dplyr::select(geometry, siteID, datasetID) %>%
      summarise() %>%
      terra::vect()

    sp_range_centroid <- centroids(sp_range)

    #----------------------------------------#

    # Compute range attributes
    tibble(
      datasetID = unique(dta$datasetID),
      samplingPeriodID = unique(dta$samplingPeriodID),
      verbatimIdentification = unique(dta$verbatimIdentification),
      nsDist = poly_attr(sp_range, "nsDist"),
      ewDist = poly_attr(sp_range, "ewDist"),
      maxDist = poly_attr(sp_range, "maxDist"),
      lengthMinRect = poly_attr(sp_range, "lengthMinRect"),
      widthMinRect = poly_attr(sp_range, "widthMinRect"),
      elonMinRect = poly_attr(sp_range, "elonMinRect"),
      elonRatio = poly_attr(sp_range, "elonRatio"),
      circ = poly_attr(sp_range, "circ"), # based on area-perimeter ratio (C=P^2/A) --> scale-dependent ()
      circNorm = poly_attr(sp_range, "circNorm"), # based on area-perimeter ratio (Cn=P^2/4πA) --> scale independent ( = 1 if it is a circle,  > 1 if it deviates from a circle)
      relCirc = poly_attr(sp_range, "relCirc"), # like minRect but with a circle (difference between circle and range; )
      lin = poly_attr(sp_range, "lin"),
      bearingMinRect = poly_attr(sp_range, "bearingMinRect"),
      bearing = poly_attr(sp_range, "bearing"),
      Dist_centroid_to_COG = distance(crds(sp_range_centroid), crds(atlas_centroid), lonlat = F)[,1], # in meters
      minDist_toBorder_centr = min(distance(crds(sp_range_centroid), crds(border_lines), lonlat = F)),
      maxDist_toBorder_centr = max(distance(crds(sp_range_centroid), crds(border_lines), lonlat = F)),
      minDist_toBorder_border = min(distance(crds(sp_range), crds(border_lines), lonlat = F)),
      maxDist_toBorder_border = max(distance(crds(sp_range), crds(border_lines), lonlat = F)),
      range_centr_long = crds(sp_range_centroid)[1],
      range_centr_lat = crds(sp_range_centroid)[2]
    )
  })


  list_range_geometries[[atlas_i]] <- range_geometries

}
tictoc::toc() # 297.64 sec elapsed

#----------------------------------------------------------#

# Bind back together:
data_range_geometries <- list_range_geometries %>%
  bind_rows()

#----------------------------------------------------------#

data_final <- data_range_geometries %>%
  left_join(data_atlas_geometry, by = "datasetID") %>%
  group_by(datasetID, samplingPeriodID, verbatimIdentification) %>%
  mutate(
    southerness = 1 - ((range_centr_lat - atlas_ymin) / ((atlas_ymax - atlas_ymin))),
    westernnes = 1 - ((range_centr_long - atlas_xmin) / ((atlas_xmax - atlas_xmin))),
    rel_maxDist = maxDist / atlas_maxDist,
    rel_ewDist = ewDist / atlas_ewDist,
    rel_nsDist = nsDist / atlas_nsDist,
    rel_elonRatio = elonRatio / atlas_elonRatio,
    rel_lin = lin / atlas_lin,
    rel_circNorm = circNorm / atlas_circNorm
  ) %>%
  ungroup() %>%
  distinct(datasetID, samplingPeriodID, verbatimIdentification, .keep_all = TRUE)

#----------------------------------------------------------#

colSums(is.na(data_final)) #relCirc is NA.

#----------------------------------------------------------#
# Save data to .rds  -----
#----------------------------------------------------------#


saveRDS(data_final, here("Data/output/1_data/2_range_geometries.rds"))
