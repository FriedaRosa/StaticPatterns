library(here)
library(dplyr)
library(skimr)


introduced <- read.csv(here::here("Documentation/introduced_species.csv")) %>%
  rename("scientificName" = "sci_name") %>%
  .[2:3] %>%
  mutate(introduced = 1)
Jacc <- read.csv(here::here("Documentation/Jaccard_df.csv"))
Tax <- read.csv(here::here("Documentation/Tax_lookup.csv")) %>% .[2:3] %>% unique()

introduced

Jacc %>%
  .[2:8] %>%
  left_join(Tax) %>%
  left_join(introduced, by = join_by(datasetID == datasetID, scientificName == scientificName)) %>%
  filter(is.na(scientificName))



tax_sp <- readRDS(here("Data/input/data.rds")) %>%
  select(datasetID, verbatimIdentification, scientificName) %>% distinct()

tax_sp %>%
  filter(startsWith(scientificName, "Myiopsitta"))

tax_sp %>%
  filter(startsWith(scientificName, "Euplectes"))


Jacc_introduced <- Jacc %>%
  left_join(tax_sp) %>%
  left_join(introduced) %>%
  mutate(introduced = case_when(is.na(introduced) ~ 0,
                                !is.na(introduced) ~ introduced))

library(ggplot2)

Jacc_introduced %>%
  ggplot(aes(x = introduced, y = Jaccard_dissim, fill = factor(introduced) ))+
  geom_boxplot()


Jacc_introduced %>%
  group_by(introduced) %>%
  rstatix::get_summary_stats(Jaccard_dissim, type = "full")

t.test(Jaccard_dissim ~ introduced, data = Jacc_introduced)


# temporal autocorrelation in AOO:
predictors <- readRDS(here("Data/output/1_data/2_predictors.rds"))
library(dplyr)
library(tidyr)

df_wide <- predictors %>%
  ungroup() %>%
  select(datasetID, samplingPeriodID, AOO) %>%
  group_by(datasetID, samplingPeriodID) %>%
  mutate(row_id = row_number()) %>%  # Create an index to align values across periods
  pivot_wider(names_from = samplingPeriodID, values_from = AOO) %>%
  select(-row_id)  # Remove temporary row ID

# Compute correlation for each datasetID
df_correlation <- df_wide %>%
  group_by(datasetID) %>%
  summarize(cor = cor(`1`, `2`, method = "spearman", use = "pairwise.complete.obs"), .groups = "drop")

df_correlation



# All three vs. Jaccard
GGally::ggpairs(Jacc_introduced, columns = c("a", "b", "c", "Jaccard_dissim"),
                ggplot2::aes(col = factor(datasetID)))



# A vs Jaccard
ggplot(Jacc_introduced, aes(x = a, y = Jaccard_dissim, col = datasetID)) +
  geom_point(alpha = 0.5) +
  labs(title = "Jaccard Index vs. a", x = "a (intersection size)", y = "Jaccard Index")

ggplot(Jacc_introduced, aes(x = b, y = Jaccard_dissim, col = datasetID)) +
  geom_point(alpha = 0.5) +
  labs(title = "Jaccard Index vs. b", x = "b (set2 unique elements)", y = "Jaccard Index")


ggplot(Jacc_introduced, aes(x = c, y = Jaccard_dissim, col = datasetID)) +
  geom_point(alpha = 0.5) +
  labs(title = "Jaccard Index vs. c", x = "c (set1 unique elements)", y = "Jaccard Index")



ggplot(Jacc_introduced, aes(x = a, y = b + c, fill = Jaccard_dissim)) +
  geom_bin2d(bins = 30) +
  scale_fill_viridis_c() +
  labs(title = "Heatmap of Jaccard Index by a and (b+c)", x = "a (intersection size)", y = "b + c (set differences)")



results %>%
  group_by(datasetID) %>%
  count(a, b, c, Jaccard_dissim, sort = TRUE) %>%
  arrange(desc(a))


library(GGally)

results %>%
  select(a, b, c, Jaccard_dissim, datasetID) %>%
  ggpairs(columns = 1:4, ggplot2::aes(col = factor(datasetID), alpha = 0.05))

results %>%
  tidyr::pivot_longer(cols = c(a,b,c)) %>%
  group_by(datasetID) %>%
  ggplot(aes(y = Jaccard_dissim, x = value, col = name))+
  geom_smooth()



predictors %>%
  right_join(Jacc_introduced %>% mutate(datasetID = as.factor(as.character(datasetID)))) %>%
  select(Jaccard_dissim, datasetID, a, b, c, Total_Ncells_samp) %>%
  tidyr::pivot_longer(cols = c(a,b,c)) %>%
  mutate(prop_value = case_when(datasetID == 5 ~ value / 628,
                                datasetID == 6 ~ value / 5319,
                                datasetID == 13 ~ value / 1184,
                                datasetID == 26 ~ value / 2821)) %>%
  ggplot(aes(y = Jaccard_dissim, x = prop_value, col = datasetID))+
  geom_smooth()+
  #geom_point()+
  facet_wrap(~ name)

