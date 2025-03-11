#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#                      02_Preprocess_data.R
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

#----------------------------------------------------------#
# Load data -----
#----------------------------------------------------------#

grids <-
  readRDS(vars$grid)

data <-
  readRDS(vars$data)

data_sf <-
  readRDS(vars$data_sf)

meta <-
  read.csv(paste0(vars$Documentation, "METADATA_datasets.csv"))

#----------------------------------------------------------#
# 1. Data wrangling -----
#----------------------------------------------------------#


#--------------------------------------------------#
# 1.1. Cell-indicator for repeated sampling -----
#--------------------------------------------------#

data_sf2 <-
  data_sf %>%
  group_by(datasetID, scalingID, siteID, samplingPeriodID) %>%
  mutate(cell_sampled = if_else(is.na(verbatimIdentification), 0, 1)) %>%
  ungroup()

#--------------------------------------------------#

cells_rep <-
  data_sf2 %>%
  st_drop_geometry() %>%
  distinct(datasetID, scalingID, siteID, samplingPeriodID, cell_sampled) %>%
  group_by(datasetID, scalingID, siteID) %>%
  dplyr::summarise(
    cell_sampling_repeats = sum(cell_sampled), .groups = "keep"
  )

#--------------------------------------------------#

# Merge indicators and clean data
data_sf3 <-
  data_sf2 %>%
  left_join(cells_rep) %>%
  select(
    datasetID, scalingID, siteID,
    cell_sampling_repeats, samplingPeriodID,
    verbatimIdentification, scientificName,
    centroidDecimalLongitude, centroidDecimalLatitude,
    croppedArea) %>%
  unique()

# Check numbers
glimpse(data_sf3)

# Drop geometry to get data:
data <-
  data_sf3 %>%
  st_drop_geometry()


#--------------------------------------------------#
# 1.2. Native vs. Introduced -----
#--------------------------------------------------#

sf_use_s2(FALSE)

# Load range maps for introduced species
BirdLife_introduced <-
  st_read(here("Data/input/shp_introduced/"))

# Read and preprocess the data
countries <-
  readRDS(here("Data/input/grid.rds")) %>%
  filter((datasetID != 5 | scalingID == 64) & (!(datasetID %in% c(6, 13, 26)) | scalingID == 128)) %>%
  select(datasetID, geometry) %>%
  unique() %>%
  st_make_valid()


#--------------------------------------------------#

# Handle datasetID == 6 separately
countries_6 <-
  countries %>%
  filter(datasetID == 6) %>%
  summarize(
    geometry = st_union(geometry),
    .groups = "keep"
  ) # Combine all geometries for datasetID == 6

countries_6$datasetID <- 6

# Handle all other datasetIDs
countries_others <-
  countries %>%
  filter(datasetID != 6)

# Combine back into a single object
countries_final <-
  bind_rows(countries_6, countries_others) %>%
  st_make_valid()

#--------------------------------------------------#

# Spatial join: Countries and introduced ranges
introduced_sp <-
  st_join(countries_final, BirdLife_introduced,
          join = st_intersects)

# Check numbers
introduced_sp %>%
  st_drop_geometry() %>%
  group_by(datasetID) %>%
  summarize(
    n_sp = n_distinct(sci_name),
    .groups = "keep")

# Save list of introduced species to .csv for documentation
introduced_sp %>%
  st_drop_geometry() %>%
  group_by(datasetID) %>%
  write.csv(here("Documentation/introduced_species.csv"))

#--------------------------------------------------#

# Remove introduced species spatially:

data_sf4 <- introduced_sp %>%
  st_drop_geometry() %>%
  select(datasetID, sci_name) %>%
  rename("scientificName" = "sci_name") %>%
  unique() %>%
  mutate(introduced = 1) %>%
  right_join(data_sf3) %>%
  mutate(
    introduced = case_when(is.na(introduced) ~ 0, .default = introduced)) %>%
  filter(introduced == 0) # remove them

#--------------------------------------------------#

# Save list of removed species:

introduced_sp %>%
  st_drop_geometry() %>%
  select(datasetID, sci_name) %>%
  rename("scientificName" = "sci_name") %>%
  unique() %>%
  mutate(introduced = 1) %>%
  right_join(data_sf3) %>%
  mutate(
    introduced = case_when(is.na(introduced) ~ 0, .default = introduced)) %>%
  filter(introduced == 1) %>%
  st_drop_geometry() %>%
  distinct(verbatimIdentification, datasetID) %>%
  write.csv(paste0(vars$Documentation, "META_removed_introduced_sp.csv"))


#--------------------------------------------------#
# 1.3. Species lost or gained  -----
#--------------------------------------------------#

# Filter species data:

all_sp_raw <-
  data %>%
  filter(scalingID == 1) %>%
  distinct(datasetID, samplingPeriodID, verbatimIdentification) %>%
  na.omit()


# Get species lost or gained completely for each dataset:

list_lost_gained_species <- lapply(vars$atlas_names, function(atlas) {

  df <-
    all_sp_raw %>%
    filter(datasetID == atlas)

  sp1 <-
    df %>%
    filter(samplingPeriodID == 1) %>%
    pull(verbatimIdentification)

  sp2 <-
    df %>%
    filter(samplingPeriodID == 2) %>%
    pull(verbatimIdentification)

  list(
    lost = setdiff(sp1, sp2),
    gained = setdiff(sp2, sp1))

})

names(list_lost_gained_species) <-
  names(vars$atlas_names)

#--------------------------------------------------#

# Filter species sampled twice (in cells sampled twice)
common_sp <-
  data %>%
  filter(cell_sampling_repeats == 2, scalingID == 1) %>%
  group_by(datasetID, verbatimIdentification) %>%
  dplyr::summarise(sp_sampling_repeats = n_distinct(samplingPeriodID), .groups = "drop")

# those sampled only once
excluded_sp <-
  data %>%
  filter(cell_sampling_repeats == 2, scalingID == 1) %>%
  group_by(datasetID, verbatimIdentification) %>%
  mutate(sp_sampling_repeats = n_distinct(samplingPeriodID)) %>%
  filter(sp_sampling_repeats == 1) %>%
  ungroup() %>%
  select(datasetID, verbatimIdentification, samplingPeriodID, sp_sampling_repeats) %>%
  unique()

excluded_sp

#--------------------------------------------------#


# Summarize dropped species:
excluded_sp %>%
  group_by(datasetID, samplingPeriodID) %>%
  summarise(n_sp = n_distinct(verbatimIdentification)) %>%
  mutate(
    lost = case_when(samplingPeriodID == 1 ~ n_sp),
    gained = case_when(samplingPeriodID == 2 ~ n_sp)
  ) %>%
  group_by(datasetID) %>%
  summarise(
    lost = sum(lost, na.rm = TRUE),
    gained = sum(gained, na.rm = TRUE),
    .groups = "drop"
  )


#----------------------------------------------------------#
# Apply data filters  -----
#----------------------------------------------------------#

presence_data_filt <-
  data %>%
  full_join(common_sp) %>%
  filter(cell_sampling_repeats == 2, sp_sampling_repeats == 2) %>%
  mutate(cell_sampling_repeats = as.factor(cell_sampling_repeats)) %>%
  unique()

# final sf
data_sf5 <-
  data_sf4 %>%
  left_join(common_sp)


#----------------------------------------------------------#
# Remove underrepresented species from Japan  -----
#----------------------------------------------------------#
jp_sp_remove <-
  read.csv(here("Documentation/META_removed_sp_Japan_expert_knowledge.csv"),
           header = FALSE,
           strip.white = TRUE) %>%
  pull(V1)

jp_sp_remove <-
  gsub("[\u00A0\\s]+$",
       "",
       jp_sp_remove,
       perl = TRUE)

presence_data_filt2 <-
  presence_data_filt %>%
  filter(
    !c(verbatimIdentification %in% jp_sp_remove &
         datasetID == 13))

data_sf6 <-
  data_sf5 %>%
  filter(
  !c(verbatimIdentification %in% jp_sp_remove &
       datasetID == 13))

#----------------------------------------------------------#
# Save data to .rds  -----
#----------------------------------------------------------#

saveRDS(data_sf6, here::here("Data", "output", "1_data", "1_data_sf.rds"))
saveRDS(presence_data_filt2, here::here("Data", "output",  "1_data", "1_data_filtered.rds"))
