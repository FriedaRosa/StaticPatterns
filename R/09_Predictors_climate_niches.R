#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#               09_Predictors_climate_niches.R
#                
#
#                    Friederike Wölke 
#                        2025
#
#----------------------------------------------------------#


source(here::here("R/00_Configuration.R"))
lapply(package_list, require, character = TRUE)

#if(!"rasterSp" %in% installed.packages()[,"Package"]) remotes::install_github("RS-eco/rasterSp", build_vignettes = T)
#if(!"climateNiche" %in% installed.packages()[,"Package"]) remotes::install_github("RS-eco/climateNiche", build_vignettes = T)

library(rasterSp)
library(climateNiche)
library(factoextra)

tax_path <-
  here("Documentation/Tax_lookup.csv")

# planar CRS
Mollweide_CRS <- 'PROJCS["ProjWiz_Custom_Mollweide",
 GEOGCS["GCS_WGS_1984",
  DATUM["D_WGS_1984",
   SPHEROID["WGS_1984",6378137.0,298.257223563]],
  PRIMEM["Greenwich",0.0],
  UNIT["Degree",0.0174532925199433]],
 PROJECTION["Mollweide"],
 PARAMETER["False_Easting",0.0],
 PARAMETER["False_Northing",0.0],
 PARAMETER["Central_Meridian",0],
 UNIT["Meter",1.0]]'

sf_use_s2(FALSE)

#----------------------------------------------------------#
# Load data  -----
#----------------------------------------------------------#

# Global ranges
BirdLife <- st_read("Data/input/shp_global/")

# Read Climate data
climate_stack_agg <- readRDS(here::here("Data/input/climate_stack.rds"))

# Taxonomic match table
Tax <-
  read.csv(tax_path) %>%
  select(verbatimIdentification, scientificName)


#----------------------------------------------------------#
# Rasterize ranges individually   -----
#----------------------------------------------------------#

rasterizeRange(dsn = here::here("Data/input/shp_global"),
               id = "sci_name", touches=TRUE, save = TRUE, resolution = 1,
               path = here::here("Data/input/species_rasters//"))


#----------------------------------------------------------#
# Match ranges and climate    -----
#----------------------------------------------------------#


# Read individual species .tifs back in:
r_birds <- rast(
  list.files(here::here("Data/input/species_rasters/"),
             pattern =".tif",
             full.names = TRUE)
)


species_names_BL <- list.files(
  here::here("Data/input/species_rasters/"),
  full.names = FALSE
) %>%
  stringr::str_remove_all("\\.tif$") %>%  # Remove ".tif" extension
  stringr::str_replace_all("_1$", "") %>% # Remove "_1" at the end
  stringr::str_replace_all("_", " ")      # Replace underscores with spaces

names(r_birds) <- species_names_BL

#--------------------------------------------------#


# Extract species coordinates
try({
  coords_sp <- lapply(seq_along(names(r_birds)), function(i) {
    r_b <- r_birds[[i]]

    # Match resolution, extent, and projection of climate raster to species raster
    r_b_matched <- project(r_b, climate_stack_agg, method = "average")

    # Convert to data frame and add species column
    r_b2 <- as.data.frame(r_b_matched, xy = TRUE) %>%
      dplyr::rename(Presence = names(r_b)) %>%
      dplyr::mutate(species = names(r_birds)[i])

    return(r_b2)
  })

  coords_df <- do.call(rbind, coords_sp)
})


#----------------------------------------------------------#
# Construct climate space    -----
#----------------------------------------------------------#

# Climate PCA (construct global climate space)
pca_res <- prcomp(climate_stack_agg, scale. = T)


# extract loadings from PCA for documentation
pca_res$rotation %>%
  write.csv(here::here("Documentation/META_PCA_loadings_climNiche.csv"))


#--------------------------------------------------#


# Visualize contributions to PC1 and PC2
p_pc1 <- fviz_contrib(pca_res, choice = "var", axes = 1, top = 11) +
  ggthemes::theme_base() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = ))

ggsave("figures/1_data/PC1_var_contributions.pdf", p_pc1)

p_pc2 <- fviz_contrib(pca_res, choice = "var", axes = 2, top = 11) +
  ggthemes::theme_base() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/1_data/PC2_var_contributions.pdf", p_pc2)


#----------------------------------------------------------#
# Project bird ranges into climate space    -----
#----------------------------------------------------------#

# Unique species
sp_v <- unique(coords_df$species)

# Process climate values and PCA projections
results <- lapply(sp_v, function(species_name) {
  # Filter coordinates for the current species
  sp <- coords_df %>%
    filter(species == species_name) %>%
    dplyr::select(-species, -Presence)

  # Extract climate values for the species
  clim_vals <- extract(climate_stack_agg, sp, xy = TRUE, df = TRUE) %>%
    as.data.frame() %>%
    mutate(species = species_name)

  # Perform PCA projections
  sp_pca <- predict(pca_res, clim_vals[, 2:14])[, 1:2] %>%
    as.data.frame() %>%
    mutate(species = species_name)

  # Return both climate values and PCA results
  list(clim_vals = clim_vals, sp_pca = sp_pca)
})

#--------------------------------------------------#


# Combine results
clim_vals_all <- do.call(rbind, lapply(results, `[[`, "clim_vals"))
pca_all <- do.call(rbind, lapply(results, `[[`, "sp_pca")) %>%
  data.frame()


#----------------------------------------------------------#
# Calculate climate niche breadth    -----
#----------------------------------------------------------#

Niches_df <-
  pca_all %>%
  group_by(species) %>%
  mutate(PC1_sd = sd(PC1),
         PC2_sd = sd(PC2)) %>%
  distinct(species, PC1_sd, PC2_sd) %>%
  dplyr::rename("scientificName" = "species") %>%
  right_join(Tax) %>%
  distinct(verbatimIdentification, scientificName, .keep_all = TRUE) %>%
  ungroup()

skim(Niches_df)

#----------------------------------------------------------#
# Calculate global range size    -----
#----------------------------------------------------------#

BirdLife_list <- BirdLife %>%
  group_by(sci_name) %>%   # Group by species
  group_split()
sf_use_s2(FALSE)

#--------------------------------------------------#

# Set up parallel session:
library(future.apply)

# Set up parallel processing
plan(multisession, workers = parallel::detectCores() - 4)  # Leave 4 cores free

#--------------------------------------------------#

RangeSizes_l <- future_lapply(BirdLife_list, function(sp_data) {
  sp_data %>%
    group_by(sci_name) %>%
    st_transform(crs = Mollweide_CRS) %>%
    st_make_valid() %>%
    summarise() %>%
    mutate(GlobRangeSize = as.numeric(st_area(.))) %>%  # Convert to km²
    st_drop_geometry() %>%
    group_by(sci_name) %>%
    summarise(
      GlobRangeSize_km2 = sum(GlobRangeSize, na.rm = TRUE) / 1e6,  # Sum range sizes
      .groups = "keep"
    )
}, future.seed = TRUE)  # Ensure reproducibility


# Combine lists
RangeSizes <- bind_rows(RangeSizes_l) %>%
  rename(scientificName = sci_name) %>%
  right_join(Tax, by = "scientificName") %>%
  distinct(scientificName, verbatimIdentification, .keep_all = TRUE) %>%
  mutate(GlobRangeSize_km2 = ifelse(GlobRangeSize_km2==0, NA, GlobRangeSize_km2)) %>%
  ungroup()


#----------------------------------------------------------#
# Save results to .rds  -----
#----------------------------------------------------------#


BirdLife_Predictors <- Niches_df %>%
  full_join(RangeSizes) %>%
  distinct(scientificName,verbatimIdentification, .keep_all=TRUE) %>%
  ungroup()

skim(BirdLife_Predictors)
saveRDS(BirdLife_Predictors, here::here("Data/output/1_data/2_all_preds_BirdLife.rds"))
