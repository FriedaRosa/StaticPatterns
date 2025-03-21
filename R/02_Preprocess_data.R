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
sf_use_s2(TRUE)

#----------------------------------------------------------#
# Load data -----
#----------------------------------------------------------#

grids <-
  readRDS(vars$grid)

data_sf <-
  readRDS(vars$data_sf)

#----------------------------------------------------------#
# 1. Spatial filter -----
#----------------------------------------------------------#

# get list of cells and sampling repeats
# get list of species and filters (or whether to include them)

#--------------------------------------------------#
# 1.1. Cell-indicator for repeated sampling -----
#--------------------------------------------------#

cells <-
  data_sf %>%
  st_drop_geometry() %>%
  group_by(datasetID, scalingID, siteID, samplingPeriodID) %>%
  mutate(cell_sampled = if_else(is.na(verbatimIdentification), 0, 1)) %>%
  ungroup() %>%
  distinct(datasetID, scalingID, siteID, samplingPeriodID, cell_sampled, croppedArea) %>%
  group_by(datasetID, scalingID, siteID, croppedArea) %>%
  dplyr::summarise(
    cell_sampling_repeats = n_distinct(samplingPeriodID, na.rm = TRUE),
    .groups = "keep") %>%
  mutate(cells_keep =
           case_when(
             cell_sampling_repeats == 2 & !is.na(croppedArea) ~ 1,
             cell_sampling_repeats %in% c(0,1) | is.na(croppedArea) ~ 0,
            .default = NA)) %>%
  unique()

# Checks:
colSums(is.na(cells))

#--------------------------------------------------#
# 2. Species filter:
#--------------------------------------------------#

species_df <-
  data_sf %>%
  st_drop_geometry() %>%
  filter(scalingID == 1) %>%
  distinct(datasetID, verbatimIdentification, scientificName, samplingPeriodID) %>%
  filter(!is.na(verbatimIdentification))

# Checks:
colSums(is.na(species_df))

#--------------------------------------------------#
# 2.1 Under-sampled species from Japan  -----
#--------------------------------------------------#

jp_sp_remove <-
  read.csv(here("Documentation/META_removed_sp_Japan_expert_knowledge.csv"),
           header = FALSE,
           strip.white = TRUE) %>%
  pull(V1)
#--------------------------------------------------#

species_df2 <- species_df %>%
  mutate(
    sp_remove_expert =
      case_when(
        verbatimIdentification %in% jp_sp_remove & datasetID == 13 ~ 1,
        TRUE ~ 0))

#--------------------------------------------------#
# 2.2. Native vs. Introduced -----
#--------------------------------------------------#
sf_use_s2(FALSE)

# Read grid shapefiles at largest resolution(1-cell-aggregations)
countries <-
  grids %>%
  group_by(datasetID) %>%
  filter(scalingID == max(scalingID)) %>%
  select(datasetID, geometry) %>%
  summarize(
    geometry = st_union(geometry),
    .groups = "keep"
  ) %>%
  st_make_valid()

# Load range maps for introduced species
BirdLife_introduced <-
  st_read(here("Data/input/shp_introduced/")) %>%
  st_transform(crs = st_crs(countries))

# check if they are identical now:
st_crs(countries) == st_crs(BirdLife_introduced)

# Spatial join: Countries and introduced ranges
introduced_sp <-
  st_join(countries, BirdLife_introduced,
          join = st_intersects) %>%
  rename("scientificName" = "sci_name") %>%
  mutate(introduced = 1) %>%
  st_drop_geometry() %>%
  select(datasetID, scientificName, introduced) %>%
  arrange(datasetID)

#--------------------------------------------------#

species_df3 <-
  species_df2 %>%
  left_join(introduced_sp) %>%
  mutate(introduced = case_when(is.na(introduced) ~ 0,
                                .default = introduced))


#--------------------------------------------------#
# 2.3. Species lost or gained  -----
#--------------------------------------------------#

common_sp <-
  data_sf %>%
  st_drop_geometry() %>%
  filter(scalingID == 1) %>%
  left_join(species_df3) %>%
  left_join(cells %>% filter(scalingID == 1)) %>%
  na.omit() %>%
  filter(cells_keep == 1 & introduced == 0 & sp_remove_expert == 0) %>%
  group_by(datasetID, verbatimIdentification) %>%
  dplyr::summarise(sp_sampling_repeats = n_distinct(samplingPeriodID),
                   .groups = "drop")


# checks
common_sp %>%
  group_by(datasetID, sp_sampling_repeats) %>%
  dplyr::summarise(n_sp = n_distinct(verbatimIdentification),
                   .groups = "keep")

#--------------------------------------------------#

species_df4 <-
  species_df3 %>%
  left_join(common_sp) %>%
  mutate(
    species_keep =
      case_when(
        sp_remove_expert == 0 & sp_sampling_repeats == 2 & introduced == 0 ~ 1,
        TRUE ~ 0))



#----------------------------------------------------------#
# Apply data filters  -----
#----------------------------------------------------------#

data_sf2 <-
  data_sf %>%
  left_join(cells) %>%
  left_join(species_df4)

data_filt <-
  data_sf2 %>%
  st_drop_geometry() %>%
  filter(!is.na(verbatimIdentification) & species_keep == 1 & cells_keep == 1)

# Checks:
data_filt %>%
  group_by(datasetID, samplingPeriodID) %>%
  skimr::skim()

#----------------------------------------------------------#
# Save data to .rds  -----
#----------------------------------------------------------#

saveRDS(data_sf2, here::here("Data", "output", "1_data", "1_data_sf.rds"))
saveRDS(data_filt, here::here("Data", "output",  "1_data", "1_data_filtered.rds"))

