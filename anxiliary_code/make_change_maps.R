#----------------------------------------------------#
# Start clean
#----------------------------------------------------#

rm(list = ls())

#----------------------------------------------------#
# Load libraries
#----------------------------------------------------#

library(sf)
library(dplyr)
library(here)
library(skimr)
library(ggplot2)
library(purrr)


#----------------------------------------------------#
# Part A: Create spatial change data
#----------------------------------------------------#

#----------------------------------------------------#
# Load data
#----------------------------------------------------#

grids <-
  readRDS(here("Data/input/grid.rds")) %>%
  filter(scalingID == 1)

dat <-
  readRDS(here("Data/output/1_data/1_data_filtered.rds")) %>%
  filter(scalingID == 1) %>%
  filter(cell_sampling_repeats == 2 & sp_sampling_repeats == 2)

#----------------------------------------------------#
# Merge grids & data
#----------------------------------------------------#

data_sf <-
  full_join(grids, dat)

#----------------------------------------------------#
# Extract all sites
#----------------------------------------------------#

all_cells <-
  data_sf %>%
  st_drop_geometry() %>%
  distinct(datasetID,
           siteID,
           cell_sampling_repeats) %>%
  mutate(cell_sampling_repeats = case_when(
    is.na(cell_sampling_repeats) ~ factor(0),
    .default = cell_sampling_repeats))

#----------------------------------------------------#
# Split data by groups into a list
#----------------------------------------------------#

list_dat <-
  dat %>%
  group_by(datasetID) %>%
  group_split()

#----------------------------------------------------#
# Function to get spatial change in occupied sites for 2 atlas periods
#----------------------------------------------------#

get_occupancy_change <- function(data_i, all_cells) {

  # dataset identificator
  atlas_id <- unique(data_i$datasetID)

  # get species list for this atlas
  sp_list <- na.omit(unique(data_i$verbatimIdentification))

  # filter all cells to atlas
  atlas_cells <- all_cells %>% filter(datasetID == atlas_id)

  # Process each species in parallel using map
  map_dfr(sp_list, function(this_species) {
    # subset data to species
    data_species <- filter(data_i, verbatimIdentification == this_species)

    # get occpied sites from period 1 and 2
    sites1 <- unique(data_species %>% filter(samplingPeriodID == 1) %>% pull(siteID))
    sites2 <- unique(data_species %>% filter(samplingPeriodID == 2) %>% pull(siteID))

    # get occupancy status for all cells
    stable <- intersect(sites1, sites2)
    colonized <- setdiff(sites2, sites1)
    extirpated <- setdiff(sites1, sites2)

    # get not sampled cells
    not_sampled <- atlas_cells %>%
      filter(cell_sampling_repeats == 0) %>%
      pull(siteID) %>%
      unique()

    # get cells not occupied
    not_occupied <- setdiff(atlas_cells$siteID, c(stable, colonized, extirpated, not_sampled))

    # bind to dataframe
    atlas_cells %>%
      mutate(
        verbatimIdentification = this_species,
        change = factor(case_when(
          siteID %in% stable ~ "stable occupancy",
          siteID %in% colonized ~ "newly colonized",
          siteID %in% extirpated ~ "newly extirpated",
          siteID %in% not_sampled ~ "not sampled",
          siteID %in% not_occupied ~ "not occupied",
          TRUE ~ NA_character_
        ))
      )
  })
}

#----------------------------------------------------#
# Apply function to each dataset
#----------------------------------------------------#
list_change_species <-
  map(list_dat, ~get_occupancy_change(.x, all_cells))

#----------------------------------------------------#
# make sf data with change info
#----------------------------------------------------#
change_sf <-
  list_change_species %>%
  bind_rows() %>%
  full_join(grids) %>%
  st_as_sf()
#----------------------------------------------------#
# check NAs
#----------------------------------------------------#
change_sf %>%
  is.na() %>%
  colSums()
#----------------------------------------------------#
# split change data by dataset into lists
#----------------------------------------------------#
list_change_sf <-
  change_sf %>%
  group_by(datasetID) %>%
  group_split()

#----------------------------------------------------#
# Get country boundaries for plotting:
#----------------------------------------------------#


library(rnaturalearth)
sf_use_s2(TRUE)

#----------------------------------------------------#
# Europe:
#----------------------------------------------------#
my_crs_eu <-
  grids %>%
  st_drop_geometry() %>%
  filter(datasetID == 26) %>%
  distinct(verbatimFootprintSRS) %>%
  pull()

# crop europe
buffered_bbox <-
  st_buffer(
    grids %>%
      filter(datasetID == 26) %>%
      st_transform(crs = my_crs_eu) %>%
      select(siteID) %>%
      unique(),
    dist = 15)

europe <-
  ne_countries(continent = "europe",
               scale = 10) %>%
  st_transform(crs = my_crs_eu) %>%
  st_crop(buffered_bbox) %>%
  select(name)

#----------------------------------------------------#
# Czechia:
#----------------------------------------------------#
my_crs_cz <-
  grids %>%
  st_drop_geometry() %>%
  filter(datasetID == 5) %>%
  distinct(verbatimFootprintSRS) %>%
  pull()

czechia <-
  ne_countries(country = "czechia", scale = 10) %>%
  st_transform(crs = my_crs_cz) %>%
  select(name)

#----------------------------------------------------#
# Japan:
#----------------------------------------------------#
my_crs_jp <-
  grids %>%
  st_drop_geometry() %>%
  filter(datasetID == 13) %>%
  distinct(verbatimFootprintSRS) %>%
  pull()

japan <-
  ne_countries(country = "japan", scale = 10) %>%
  st_transform(crs = my_crs_jp ) %>%
  select(name)

#----------------------------------------------------#
# New York state:
#----------------------------------------------------#
my_crs_ny <-
  grids %>%
  st_drop_geometry() %>%
  filter(datasetID == 6) %>%
  distinct(verbatimFootprintSRS) %>%
  pull()

new_york <-
  ne_states(country = "United States of America") %>%
  filter(name == "New York") %>%
  st_transform(crs = my_crs_ny) %>%
  select(name)



#----------------------------------------------------#
# Part B: Plotting maps
#----------------------------------------------------#

#----------------------------------------------------#
# Transform change data to local CRS
#----------------------------------------------------#

list_country_borders <-
  list(czechia, new_york, japan, europe)

list_crs <-
  list(my_crs_cz,
       my_crs_ny,
       my_crs_jp,
       my_crs_eu)

list_change_sf_transformed <-
  map2(list_change_sf, list_crs, ~st_transform(.x, .y))

#----------------------------------------------------#
# Colors:
#----------------------------------------------------#

colors <-
  setNames(
    c("#574f7d", "#998ec3", "#f1a340",
      "#F7F7F7","lightgrey"),
    c("stable occupancy", "newly colonized", "newly extirpated",
      "not occupied", "not sampled"))

#----------------------------------------------------#
# Theme:
#----------------------------------------------------#

my_theme <-
  ggthemes::theme_map() +  # Use theme_map() as a base
  theme(
    legend.position = "right",   # Move legend outside
    legend.box = "vertical",     # Arrange legend items vertically
    legend.margin = margin(10, 10, 10, 10),  # Add space around legend
    plot.margin = margin(10, 50, 10, 10)  # Expand right margin to fit legend
  )

#----------------------------------------------------#
# Map:
#----------------------------------------------------#

ggplot()+
  geom_sf(
    data =list_change_sf[[1]] %>%
      filter(verbatimIdentification == "Cinclus cinclus"),
          aes(fill = change),
          color = "lightgrey")+
  geom_sf(
    data = list_country_borders[[1]],
          fill = NA,
          col = "black")+
  scale_fill_manual(values = colors)+
  my_theme




#----------------------------------------------------#
# Map all species via loop
#----------------------------------------------------#
out_paths <- list(here("figures/3_maps/czechia//"),
                  here("figures/3_maps/new_york//"),
                  here("figures/3_maps/japan//"),
                  here("figures/3_maps/europe//"))

for(atlas_i in seq_along(list_change_sf)){
  this_atlas <- list_change_sf[[atlas_i]]
  atlas_id <- unique(this_atlas$datasetID)
  sp_list <- unique(this_atlas$verbatimIdentification)

  for(sp_i in seq_along(sp_list)){

    this_species <- sp_list[sp_i]
    data_species <- this_atlas %>%
      filter(verbatimIdentification == this_species)

    species_expr <- bquote(italic(.(this_species)))
    atlas_expr <- paste0("dataset = ", atlas_id)
    title_expr <- bquote(.(species_expr) * " " * .(atlas_expr))


    map <-
      ggplot()+
      my_theme+
      geom_sf(data = data_species,
              aes(fill = change),
              col = "lightgrey")+
      geom_sf(data = list_country_borders[[atlas_i]],
              fill = NA,
              col = "black")+
      scale_fill_manual(values = colors)


    # Export to powerpoint vector file
    export::graph2ppt(map,
                      width=9,
                      height=9,
                      file=paste0(out_paths[[atlas_i]],
                                 gsub(" ", "_", this_species),
                                  "_change_map.pptx"))


    }

  }
