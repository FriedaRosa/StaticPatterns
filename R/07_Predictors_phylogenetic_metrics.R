#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#               07_Predictors_phylogenetic_metrics.R
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
  read_csv(tax_path)

tree <-
  ape::read.tree(
  here::here("Data/input/MRC_consensus_BirdTree.tre"))

#----------------------------------------------------------#
# Calculate phylogenetic distinctiveness -----
#----------------------------------------------------------#

pd <-
  phyloregion::evol_distinct(
  tree,
  type = "fair.proportion",
  scale = FALSE,
  use.branch.lengths = TRUE) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "tip.label") %>%
  rename(pd = ".") %>%
  right_join(Tax) %>%
  select(verbatimIdentification, scientificName, pd)

#----------------------------------------------------------#
# Save results to .rds  -----
#----------------------------------------------------------#

saveRDS(pd,
        here::here("Data/output/1_data/single_predictors/2_phylo_distinct.rds"))
