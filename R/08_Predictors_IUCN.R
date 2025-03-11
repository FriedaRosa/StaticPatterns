#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#               08_Predictors_IUCN.R
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

Tax <-
  read_csv(tax_path) %>%
  distinct(verbatimIdentification, scientificName)

#----------------------------------------------------------#
# Get IUCN status  -----
#----------------------------------------------------------#


IUCN.list <-
  iucn_summary(unique(Tax$scientificName),
               distr.detail = F,
               key = Sys.getenv('IUCN_REDLIST_KEY'))

# Extract codes
IUCN <- lapply(names(IUCN.list), function(name) {
  item <- IUCN.list[[name]]

  # Check if the element is a list and contains `red_list_category`
  if (is.list(item) && !is.null(item$red_list_category) && !is.null(item$red_list_category$code)) {
    return(data.frame(name = name, code = item$red_list_category$code, stringsAsFactors = FALSE))
  } else {
    return(data.frame(name = name, code = NA, stringsAsFactors = FALSE)) # Fill missing with NA
  }
}) %>% do.call(rbind, .)

# Check NA species
species_with_NA <-
  IUCN %>% filter(is.na(code))


IUCN.list2 <-
  iucn_summary(unique(species_with_NA$name),
               distr.detail = F,
               key = Sys.getenv('IUCN_REDLIST_KEY'))

IUCN2 <- lapply(names(IUCN.list2), function(name) {
  item <- IUCN.list2[[name]]

  # Check if the element is a list and contains `red_list_category`
  if (is.list(item) && !is.null(item$red_list_category) && !is.null(item$red_list_category$code)) {
    return(data.frame(name = name, code = item$red_list_category$code, stringsAsFactors = FALSE))
  } else {
    return(data.frame(name = name, code = NA, stringsAsFactors = FALSE)) # Fill missing with NA
  }
}) %>% do.call(rbind, .)

IUCN_merged <- rbind(IUCN, IUCN2) %>% unique() %>% na.omit()
# Reshape
IUCN_df <-
  as.data.frame(IUCN_merged,
                row.names = c(IUCN_merged$name)) %>%
  rownames_to_column(var = "scientificName") %>%
  right_join(Tax) %>%
  distinct(verbatimIdentification, scientificName, code) %>%
  mutate(
    code = as.factor(case_when(scientificName == "Nannopterum auritus" ~ "LC",
                     .default = code))
  )

#----------------------------------------------------------#
# Save results to .rds  -----
#----------------------------------------------------------#

saveRDS(IUCN_df,
        here::here("Data/output/1_data/single_predictors/2_IUCN_20250225.rds"))
