## Taxonomic Stuff here ##

#################################################################
#################################################################
#
# Every dataset uses it's own taxonomy.
# In order to assess as many species based on field data as possibe,
# we decided to keep the original taxonomies for all analyses.
# That means that same species that are now
# recognized to belong to different clades,
# will be assessed based on their old taxonomy.
# Subspecies and hybrids will also be based on the old taxonomy
# although more updated nomenclature is available.
#
#################################################################
#
# There are 3 data sources that may lead to information loss if
# species are not standardized to match the taxonomy in those datastets.
# Those are:
#
# 1. BirdLife global range maps for the calculation of global climate niche breadth
# (version: 2024)
#
# 2. BirdTree phylogeny for the calculation of phylogenetic distinctiveness and
# assessment of phylogenetic autocorrelation in the model residuals.
# (version: 2012)
#
# 3. AVONET trait database for the extraction of species traits
# (version 2022)
#
# (4. IUCN but this is part of BirdLife2024 and has the same taxonomy)
#
#################################################################
#################################################################


# 1. BirdLife 2024 + dataset verbatim names
# The data was matched to BirdLife 2024 taxonomy (there was only one change in scientific name)
tax_lookup <- readRDS(here::here("Data/output/1_data_filtered.rds")) %>%
  distinct(verbatimIdentification, scientificName) %>% na.omit() # BirdLife 2024 taxonomy
# should match BirdLife range maps (2024) and IUCN data

length(unique(tax_lookup$scientificName)) # 727 standardized species
length(unique(tax_lookup$verbatimIdentification)) # 765 from raw data



# Crosswalk for BL2010-BL2018
crosswalk <- readxl::read_excel("C:/Users/wolke/OneDrive - CZU v Praze/Dokumenty/PhD_Projects/Project-folders/StaticPredictors/Carmen/BL_taxonomy_crosswalk.xlsx",
                                sheet = "Crosswalk") %>%
  dplyr::select(1,41) %>% # 2018 & 2010 taxonomy
  setNames(c("ScientificName2018", "ScientificName2010"))

skim(crosswalk)
crosswalk %>% filter(ScientificName2018 == "Not Recognised") %>% kable()
crosswalk %>% filter(ScientificName2010 == "Not Recognised") %>% kable()


###############################################
## Differences between BL2024 & BL2018
nm <- setdiff(tax_lookup$scientificName, crosswalk$ScientificName2018) # 21 not matched
setdiff(tax_lookup$verbatimIdentification, crosswalk$ScientificName2018) #100 not matched
setdiff(tax_lookup$verbatimIdentification, crosswalk$ScientificName2010) #95 not matched

setdiff(tax_lookup %>% filter(scientificName %in% nm) %>% pull(verbatimIdentification),
        crosswalk$ScientificName2010) # 8 out of 21 not matched (but 13 more matched)
###############################################

## step 1: merge scientificName & ScientificName2018
## step 2: merge verbatimIdentification & ScientificName2010


tax_lookup2 <- tax_lookup %>%
  left_join(
    crosswalk,
    by = join_by(scientificName == ScientificName2018), keep = TRUE) %>%
  left_join(
    crosswalk,
    by = join_by(verbatimIdentification == ScientificName2010),
    keep = TRUE) %>%
  mutate(ScientificName2010 = coalesce(ScientificName2010.x, ScientificName2010.y)) %>%
  select(-ScientificName2010.x, -ScientificName2010.y) %>%
  mutate(ScientificName2018 = coalesce(ScientificName2018.x, ScientificName2018.y)) %>%
  select(-ScientificName2018.x, -ScientificName2018.y) %>%
  mutate(ScientificName2010 = case_when(ScientificName2010 == "Not Recognised" ~ NA,.default = ScientificName2010),
         ScientificName2018 = case_when(ScientificName2018 == "Not Recognised" ~ NA,.default = ScientificName2018))

skim(tax_lookup2)

tax_lookup2 %>% filter(is.na(ScientificName2010)) %>% kable()
tax_lookup2 %>% filter(is.na(ScientificName2018)) %>% kable()

# Notes: there are NAs in the ScientificName2018, 2010 columns
# We have to be careful not to match NAs to NAs which would add wrong data.
#################
#################

# Avonet crosswalk: Birdlife / BirdTree (Species 1 = ScientificName2018)
BL_BT_crosswalk <- read_excel("C:/Users/wolke/OneDrive - CZU v Praze/Dokumenty/PhD_Projects/Project-folders/StaticPredictors/Data/AVONET/AVONET Supplementary dataset 1.xlsx",
                              sheet = "BirdLife–BirdTree crosswalk")

dd_BL_notMatched <- setdiff(tax_lookup2$ScientificName2018 %>% na.omit(), BL_BT_crosswalk$Species1 %>% na.omit()) # 1 sp not matched

intersect(dd_BL_notMatched, BL_BT_crosswalk$Species3) # missing species found in Species3 ("Psittacula krameri")

tax_lookup3 <- tax_lookup2 %>%
  select(ScientificName2018) %>% unique() %>%
  na.omit() %>%
  left_join(BL_BT_crosswalk %>% filter(!is.na(Species1)) %>%
              select(-Match.type, -Match.notes),
            by = join_by("ScientificName2018" == "Species1"),
                         keep = FALSE) %>%
  setNames(c("ScientificName2018", "BirdTree")) %>%
  # left_join(BL_BT_crosswalk %>%
  #             filter(Species3 == "Psittacula krameri") %>% select(Species1, Species3),
  #           by = join_by("ScientificName2018" == "Species3"), keep = TRUE) %>%
  right_join(tax_lookup2) %>%
  select(verbatimIdentification, scientificName, ScientificName2010, ScientificName2018, BirdTree)

skim(tax_lookup3)

# check NAs:
tax_lookup3 %>% filter(if_any(everything(), is.na))

# 10 species missing from phylo
# 8 species missing from AVONET (9 in ScientificName2018)

# BirdTree consensus (based on recommendations from: https://doi.org/10.1093/czoolo/61.6.959)
## created in python
library(ape)
tips <- ape::read.tree(
  here::here("Data/input/MRC_consensus_BirdTree.tre"))$tip.label
tree_sp <- data.frame(tip.label = tips)
tree_sp$tip_label <- gsub("_", " ", tree_sp$tip.label)

setdiff(unique(tax_lookup3$BirdTree)[-520], tree_sp$tip_label)
# all matched (we dropped the NA from BirdTree column)


# Merge with actual bird tree tip labels (not just the ones from Avonet)
tax_lookup4 <- tax_lookup3 %>%
  left_join(tree_sp, by = join_by("BirdTree" == "tip_label"), keep = FALSE)



write.csv(tax_lookup4, here::here("Documentation/Tax_lookup.csv"))

tax_lookup4 %>% filter(is.na(BirdTree)) %>% kable() #11
