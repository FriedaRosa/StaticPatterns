#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#               06_Predictors_avonet.R
#                
#
#                    Friederike Wölke 
#                        2025
#
#----------------------------------------------------------#

source(here::here("R/00_Configuration.R"))
lapply(package_list, require, character = TRUE)

tax_path <-
  here("Documentation/Tax_lookup.csv")

#----------------------------------------------------------#
# Load data  -----
#----------------------------------------------------------#

sci_name <-
  readRDS(here::here("Data/output/1_data/1_data_filtered.rds")) %>%
  distinct(verbatimIdentification, scientificName) # BirdLife 2024 taxonomy

Tax <- read.csv(tax_path)

traits <-
  read_excel("C:/Users/wolke/OneDrive - CZU v Praze/Dokumenty/PhD_Projects/Project-folders/StaticPredictors/Data/AVONET/AVONET Supplementary dataset 1.xlsx",
              sheet = "AVONET1_BirdLife") %>%
  dplyr::select(
    Species1, Family1, Order1, # tax info
    Habitat, Migration,Trophic.Level, Trophic.Niche,Primary.Lifestyle, # char
    Beak.Width,'Hand-Wing.Index', Mass, Habitat.Density # dbl
  ) %>%
  mutate(
    across(c(Habitat, Migration, Trophic.Level,
             Trophic.Niche, Primary.Lifestyle,
             Habitat.Density),
           as.factor)) %>%
  rename("ScientificName2018" = "Species1") %>%
  right_join(Tax %>%
               select(verbatimIdentification, scientificName, ScientificName2018)) %>%
  select(-ScientificName2018)

#----------------------------------------------------------#
# Save results to .rds  -----
#----------------------------------------------------------#

saveRDS(traits, here::here("Data/output/1_data/single_predictors/2_Avonet.rds"))
