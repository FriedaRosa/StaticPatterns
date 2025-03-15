#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#                 Sensitivity_samplingPeriod_3.R
#                
#
#                    Friederike Wölke 
#                        2025
#
#----------------------------------------------------------#




# Sensitivity : Check predictions with third atlas replications.
# only for CZ & Japan available.



# Start with clean environment
rm(list = ls())
gc()

# Source 00_Configuration.R
source(here::here("R/00_Configuration.R"))
lapply(package_list, require, character = TRUE)


#----------------------------------------------------------#
# Get atlas data -----
#----------------------------------------------------------#

# Connect to the database
con <- dbConnect(Postgres(),
                 dbname = "MOBI_atlases_v1",
                 host = "localhost",
                 port = 5432,
                 user = "frieda",
                 password = Sys.getenv('PASSWORD_SERVER')
)
dbListTables(con)
#--------------------------------------------------#

# Get filtered data:
sql_query_3 <-
  "SELECT
  \"datasetID\", \"scalingID\", \"siteID\", \"startYear\",
  \"croppedArea\", \"verbatimIdentification\", \"scientificName\",
  \"centroidDecimalLongitude\", \"centroidDecimalLatitude\"
  FROM
  \"MOBI_vw_FINAL_presence_records\"
  WHERE
  \"croppedArea\" IS NOT NULL
  AND
  \"datasetID\" IN (5,13)"


data <-
  tbl(con, sql(sql_query_3)) %>%
  collect() %>%
  mutate(samplingPeriodID = case_when(
    datasetID == 5 & startYear == 1985 ~ 1,
    datasetID == 5 & startYear == 2001 ~ 2,
    datasetID == 13 & startYear == 1974 ~ 1,
    datasetID == 13 & startYear == 1997 ~ 2,
    TRUE ~ 3)) %>%
  distinct(datasetID, scalingID, siteID, samplingPeriodID, verbatimIdentification,
           .keep_all = TRUE) %>%
  filter(samplingPeriodID %in% c(2,3))

#--------------------------------------------------#

# Get original grids:
grids <-
  st_read(con, query = vars$sql_query_grid) %>%
  filter(datasetID %in% c(5,13))


#--------------------------------------------------#

# Create sf that includes cells not sampled at all
data_sf <- grids %>%
  left_join(data)

library(ggplot2)
ggplot(data_sf %>%
         filter(samplingPeriodID == 3 & scalingID == 1 & verbatimIdentification == "Accipiter gentilis" | is.na(verbatimIdentification)))+
  geom_sf(aes(fill = verbatimIdentification), show.legend = FALSE)+
  theme_minimal()

# disconnect from MOBI db
dbDisconnect(con)

# ---------------------------------------------#


# filtered by sampling period & croppedArea & recordFilter

saveRDS(data, here("Data/input/data_3.rds"))

# filtered by sampling period & croppedArea & recordFilter
#   but with unsampled cells as NA in verbatimIdentification
saveRDS(data_sf,  here("Data/input/data_sf_3.rds"))

#----------------------------------------------------------#
# preprocess data -----
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
countries_final <-
  readRDS(here("Data/input/grid.rds")) %>%
  filter(datasetID %in% c(5,13)) %>%
  filter((datasetID != 5 | scalingID == 64) & (!(datasetID %in% c(13)) | scalingID == 128)) %>%
  select(datasetID, geometry) %>%
  unique() %>%
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
  distinct(verbatimIdentification, samplingPeriodID, datasetID) %>%
  filter(samplingPeriodID == 3)



# Filter species data:

all_sp_raw <-
  data %>%
  filter(scalingID == 1) %>%
  distinct(datasetID, samplingPeriodID, verbatimIdentification) %>%
  na.omit()


# Get species lost or gained completely for each dataset:

list_lost_gained_species <- lapply(vars$atlas_names[c(1,3)], function(atlas) {

  df <-
    all_sp_raw %>%
    filter(datasetID == atlas)

  sp1 <-
    df %>%
    filter(samplingPeriodID == 2) %>%
    pull(verbatimIdentification)

  sp2 <-
    df %>%
    filter(samplingPeriodID == 3) %>%
    pull(verbatimIdentification)

  list(
    lost = setdiff(sp1, sp2),
    gained = setdiff(sp2, sp1))

})

names(list_lost_gained_species) <-
  names(vars$atlas_names[c(1,3)])

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
  filter(sp_sampling_repeats %in% c(1)) %>%
  ungroup() %>%
  select(datasetID, verbatimIdentification, samplingPeriodID, sp_sampling_repeats) %>%
  unique()

excluded_sp %>%
  group_by(datasetID, samplingPeriodID) %>%
    summarize(n_sp = n_distinct(verbatimIdentification))

#--------------------------------------------------#


# Summarize dropped species:
excluded_sp %>%
  group_by(datasetID, samplingPeriodID) %>%
  summarise(n_sp = n_distinct(verbatimIdentification)) %>%
  mutate(
    lost = case_when(samplingPeriodID == 2 ~ n_sp),
    gained = case_when(samplingPeriodID == 3 ~ n_sp)
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

saveRDS(data_sf6, here::here("Data", "output", "1_data", "1_data_sf_3.rds"))
saveRDS(presence_data_filt2, here::here("Data", "output",  "1_data", "1_data_filtered_3.rds"))


#### Calculate atlas variables:



# Start with clean environment
rm(list = ls())
gc()

# Source 00_Configuration.R
source(here::here("R/00_Configuration.R"))
lapply(package_list, require, character = TRUE)

# Enable s2 processing for sf package with projected data
sf_use_s2(TRUE)


#----------------------------------------------------------#
# Load data -----
#----------------------------------------------------------#

# Read spatial dataset
data_sf <- readRDS(here("Data/output/1_data/1_data_sf_3.rds"))

#----------------------------------------------------------#
# Calculate grid information  -----
#----------------------------------------------------------#

# Custom function to get total and sampled area/cells from datasets
calculate_grid_info <- function(data_sf) {
  grid_info <-
    full_join(
      data_sf %>%
        st_drop_geometry() %>%
        distinct(datasetID, scalingID, siteID, cell_sampling_repeats, croppedArea) %>%
        group_by(datasetID, scalingID) %>%
        summarise(
          Total_Ncells = n_distinct(siteID)),
      data_sf %>%
        st_drop_geometry() %>%
        distinct(datasetID, scalingID, siteID, cell_sampling_repeats, croppedArea) %>%
        filter(cell_sampling_repeats == 2) %>%
        group_by(datasetID, scalingID) %>%
        summarise(
          Total_Ncells_samp = n_distinct(siteID),
          Total_area_samp = sum(croppedArea, na.rm = TRUE))
    )
  return(grid_info)
}

#--------------------------------------------------#

# Apply function to compute grid information
atlas_areas <-
  calculate_grid_info(data_sf)
#--------------------------------------------------#


# Read processed species data & merge with atlas_areas
pres_dat_final <- readRDS(here("Data/output/1_data/1_data_filtered_3.rds")) %>%
  left_join(atlas_areas) %>%
  mutate(scalingID = factor(scalingID, levels = vars$desired_levels)) %>%
  filter(!is.na(verbatimIdentification)) # Remove unsampled cells



#----------------------------------------------------------#

# Overview
atlas_areas %>%
  kable()

# Summary statistics
atlas_areas %>%
  skim()

# Save results
write.csv(atlas_areas,
          paste0(vars$Documentation, "atlas_areas_METADATA.csv"))


#----------------------------------------------------------#
# Calculate jaccard index   -----
#----------------------------------------------------------#

# Custom function to calculate jaccard across sites
jaccard <- function(set1, set2) {
  a <- length(intersect(set1, set2))
  b <- length(setdiff(set2, set1))
  c <- length(setdiff(set1, set2))
  n_union  <- length(union(set1,set2))
  return(list(jaccard_index = a / (a + b + c), a = a, b = b, c = c, n_union = n_union))
}

#----------------------------------------------------------#


# Calculate jaccard for all datasets
Jaccard_df <- pres_dat_final %>%
  filter(scalingID == 1) %>%
  group_by(datasetID, verbatimIdentification) %>%
  mutate(
    jaccard_results = list(jaccard(
      filter(pick(everything()), samplingPeriodID == 2)$siteID,
      filter(pick(everything()), samplingPeriodID == 3)$siteID
    ))
  ) %>%
  unnest_wider(jaccard_results) %>%
  group_by(verbatimIdentification, datasetID) %>%
  mutate(d = Total_Ncells_samp - n_union) %>%
  rename(Jaccard_sim = jaccard_index) %>%
  mutate(Jaccard_dissim = 1 - Jaccard_sim) %>%
  distinct(datasetID, verbatimIdentification, Jaccard_dissim, Jaccard_sim, a, b, c, d, Total_Ncells_samp)

# Save jaccard results
write.csv(Jaccard_df, paste0(vars$Documentation, "Jaccard_df_3.csv"))


#----------------------------------------------------------#

# Merge:
pres_dat_final_v2 <- pres_dat_final %>%
  left_join(Jaccard_df)

head(pres_dat_final_v2)

#----------------------------------------------------------#

rm(pres_dat_final, Jaccard_df, jaccard, atlas_areas, calculate_grid_info)


#----------------------------------------------------------#
# Calculate area of occupancy   -----
#----------------------------------------------------------#

occ_data_final <-
  pres_dat_final_v2 %>%
  ungroup() %>%
  distinct(datasetID, samplingPeriodID, scalingID,
           verbatimIdentification, siteID,
           .keep_all = TRUE) %>%
  group_by(datasetID, samplingPeriodID, scalingID, verbatimIdentification) %>%
  # Calculate AOO:
  dplyr::summarise(
    mean_area = mean(croppedArea, na.rm = TRUE),
    AOO = sum(croppedArea, na.rm = TRUE),
    occ_Ncells = n_distinct(siteID)) %>%
  left_join(pres_dat_final_v2 %>%
              mutate(scalingID = as.factor(as.character(scalingID)))) %>%
  # Calculate relative Occupancy:
  dplyr::mutate(
    rel_AOO = AOO / Total_area_samp,
    rel_occ_Ncells = occ_Ncells / Total_Ncells_samp) %>% # Prevalence
  # Remove duplicated rows:
  distinct()

#----------------------------------------------------------#

# Save AOO data to .rds
saveRDS(occ_data_final, here::here("Data", "output", "1_data", "1_occ_data_final_3.rds"))
rm(pres_dat_final_v2)

#----------------------------------------------------------#

# Create species-level data
species_data <-
  occ_data_final %>%
  dplyr::select(
    -siteID, -cell_sampling_repeats, -croppedArea, -sp_sampling_repeats
  ) %>%
  distinct(datasetID, samplingPeriodID, scalingID, verbatimIdentification,
           .keep_all = TRUE)

#----------------------------------------------------------#
# Calculate occupancy-area-relationship   -----
#----------------------------------------------------------#

# Custom function to calculate OAR and Fractal dimension
calculate_OAR <- function(species_data) {
  datasets <-
    unique(species_data$datasetID)
  sampling_periods <-
    unique(species_data$samplingPeriodID)

  sp_dta <-
    species_data %>%
    select(datasetID, samplingPeriodID, scalingID,
           verbatimIdentification,
           rel_AOO, rel_occ_Ncells, mean_area, AOO)

  # Expand the grid for all combinations of datasetID and samplingPeriodID
  species_data_new <-
    expand_grid(
      datasetID = datasets,
      samplingPeriodID = sampling_periods) %>%
    inner_join(sp_dta, by = c("datasetID", "samplingPeriodID")) %>%
    # Group by identifiers
    group_by(datasetID, samplingPeriodID, verbatimIdentification) %>%
    # Get available scales where relative occupancy is not saturated (< 1)
    summarise(
      available_scales = n(),
      mean_relAOO = mean(rel_AOO, na.rm = TRUE),
      exclude_sp_OAR = if_else(available_scales < 2, 1, 0),
      .groups = "drop"
    ) %>%
    # remove those where the range is saturated at the smallest resolution
    mutate(
      exclude_sp_OAR = if_else(available_scales == 0, 1, exclude_sp_OAR),
      mean_relAOO = if_else(available_scales == 0, 1, mean_relAOO)
    ) %>%
    full_join(sp_dta, by = c("datasetID", "samplingPeriodID", "verbatimIdentification")) %>%
    filter(
      exclude_sp_OAR == 0,
      rel_occ_Ncells < 1
    ) %>%
    distinct() %>%
    filter_at(vars(scalingID, AOO, mean_area), any_vars(!is.na(.))) %>%
    ungroup()

  # Fit models using purrr::map
  species_data_new_v2 <-
    species_data_new %>%
    group_by(datasetID, samplingPeriodID, verbatimIdentification) %>%
    nest(data = c(scalingID, AOO, mean_area, rel_AOO, rel_occ_Ncells)) %>%
    mutate(
      coefficients = map(
        data,
        ~ .x %>%
          filter(
            !is.na(AOO) & AOO > 0,
            !is.na(mean_area) & mean_area > 0
          ) %>%
          lm(log(AOO) ~ log(mean_area), data = .) %>%
          coef() %>%
          {
            tibble(
              m_AOO_a = .[2],
              b_AOO_a = .[1],
              D_AOO_a = -2 * .[2] + 2
            )
          }
      )
    ) %>%
    unnest(coefficients) %>%
    ungroup() %>%
    # Final cleanup of results
    select(datasetID, samplingPeriodID, verbatimIdentification, m_AOO_a, b_AOO_a, D_AOO_a) %>%
    distinct() %>%
    # Merge with the original dataset
    full_join(species_data) %>%
    filter(scalingID == 1) %>%
    select(
      datasetID, verbatimIdentification, samplingPeriodID, Total_area_samp,
      Total_Ncells, Total_Ncells_samp, AOO, occ_Ncells, rel_occ_Ncells,
      rel_AOO, Jaccard_dissim, m_AOO_a, b_AOO_a, D_AOO_a
    )

  return(species_data_new_v2)
}


#----------------------------------------------------------#

# Apply OAR-function:
species_data_new <- calculate_OAR(species_data) %>%
  distinct(datasetID, samplingPeriodID, verbatimIdentification, .keep_all = TRUE)   %>%
  mutate(
    D_AOO_a = case_when(is.na(D_AOO_a) ~ 2,
                        .default = D_AOO_a)) %>%
  left_join(species_data %>% filter(scalingID == 1)) %>%
  # reduce columns
  select(
    datasetID, verbatimIdentification, samplingPeriodID,
    Total_area_samp, Total_Ncells, Total_Ncells_samp,
    AOO, occ_Ncells, rel_occ_Ncells, rel_AOO,
    Jaccard_dissim, a,b,c,d, D_AOO_a
  ) %>%
  distinct(datasetID, samplingPeriodID, verbatimIdentification,
           .keep_all = TRUE)


#----------------------------------------------------------#
# Calculate log ratio AOO   -----
#----------------------------------------------------------#

time_periods <- c(2,3)

# Custom function to transform data to wide-format
transform_to_wide <- function(species_data_new, time_periods = c(2,3)) {
  # Create a list to store wide data for each time period
  wide_dfs <- list()

  for (i in seq_along(time_periods)) {
    wide_dfs[[i]] <- species_data_new %>%
      distinct(datasetID, samplingPeriodID, verbatimIdentification, AOO) %>%
      group_by(datasetID, samplingPeriodID, verbatimIdentification) %>%
      filter(samplingPeriodID == time_periods[i]) %>%
      setNames(paste0("samplingPeriodID", time_periods[i], "_", names(.))) %>%
      ungroup() %>%
      select(-c(paste0("samplingPeriodID", time_periods[i], "_samplingPeriodID"))) %>%
      dplyr::rename(
        verbatimIdentification = paste0("samplingPeriodID", time_periods[i], "_verbatimIdentification"),
        datasetID = paste0("samplingPeriodID", time_periods[i], "_datasetID")
      )
  }

  # Merge the wide data frames sequentially
  sp_dat_wide <- reduce(wide_dfs, full_join, by = c("verbatimIdentification", "datasetID"))

  cat("NA counts in wide data after processing:\n")
  print(colSums(is.na(sp_dat_wide)))

  cat("Preview of wide data:\n")
  print(head(sp_dat_wide))

  return(sp_dat_wide)
}


#----------------------------------------------------------#


# Apply function:
sp_dat_wide <-
  transform_to_wide(species_data_new, time_periods) %>%
  na.omit() # drop species lost or gained completely

#----------------------------------------------------------#

# Calculate log ratio of AOO
logRatio <-
  sp_dat_wide %>%
  mutate(
    log_R3_2 = log(samplingPeriodID3_AOO / samplingPeriodID2_AOO),
    ratio_R3_2 = samplingPeriodID3_AOO / samplingPeriodID2_AOO) %>%
  select(-samplingPeriodID2_AOO, -samplingPeriodID3_AOO)

# save final predictor table from grids/atlases
big_table_3 <-
  full_join(species_data_new, logRatio) %>%
  distinct(datasetID, samplingPeriodID, verbatimIdentification,
           .keep_all = T) %>%
  mutate_if(is.numeric, round, 3)


#----------------------------------------------------------#
# Save results to .rds  -----
#----------------------------------------------------------#

saveRDS(species_data_new, here("Data/output/1_data/1_species_data_3.rds"))
saveRDS(big_table_3, here("Data/output/1_data/2_big_table_3.rds"))
