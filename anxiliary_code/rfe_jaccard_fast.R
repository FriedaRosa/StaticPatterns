rm(list = ls())
library(caret)
library(dplyr)
library(tidymodels)

## get data
model_dd <- readRDS(here::here("Data/2_predictors.rds")) %>%
  filter(samplingPeriodID == 1) %>%
  select(
    #log_R2_1,
    Jaccard_dissim, # y
    verbatimIdentification, # speciesID
    datasetID,
    Total_area_samp,
    AOO,
    D_AOO_a,
    mean_lnLac,
    rel_occ_Ncells,
    Hand.Wing.Index,
    Mass,
    Beak.Width,
    Migration,
    Range.Size,
    Habitat,
    Habitat.Density,
    Trophic.Level,
    Trophic.Niche,
    Primary.Lifestyle,
    Range.Size,
    pd,
    PC1_sd,
    PC2_sd,
    IUCN,
    morans_I,
    joincount_delta,
    AlphaSR_sp,
    BetaSR_sp,
    GammaSR,
    circ,
    circNorm,
    Dist_centroid_to_COG,
    minDist_toBorder_border,
    southerness,
    westernnes,
    rel_circNorm,
    rel_maxDist,
    rel_elonRatio,
    rel_lin,
    atlas_elonMinRect,
    atlas_circNorm,
    atlas_bearing
  )

n_pred <- ncol(model_dd) - 2



ranger_funcs <- list(
  summary = defaultSummary,
  fit = function(x,y, first, last, ...){
    set.seed(123)
    loadNamespace("ranger")
    data = data.frame(cbind(x,y))
    ranger::ranger(y ~ .,
                   data = data,
                   num.trees = 5000, #new
                   oob.error = T,
                   importance="permutation",
                   respect.unordered.factors = TRUE)
  },
  pred = function(object, x)  {
    predict(object, data = x)$predictions
  },
  rank =
    function(object, x, y) {
      vimp <- data.frame(Overall = object$variable.importance)#ranger
      # vimp <- caret::varImp(object, type = 1, scale = TRUE) #randomForest
      # vimp <- varImp(object) # default setting
      if (is.factor(y)) {
        if (all(levels(y) %in% colnames(vimp))) {
          avImp <- apply(vimp[, levels(y), drop = TRUE], 1, mean)
          vimp$Overall <- avImp
        }
      }
      vimp <- vimp[order(vimp$Overall, decreasing = TRUE), , drop = FALSE]
      if (ncol(x) == 1) {
        vimp$var <- colnames(x)
      } else {
        vimp$var <- rownames(vimp)
      }
      vimp
    },

  selectSize =
    function(x, metric, maximize) {
      best <- if (maximize)
        which.max(x[, metric])
      else which.min(x[, metric])
      min(x[best, "Variables"])
    },

  selectVar = function(y, size) {
    library(dplyr)
    finalImp <- y %>%
      # keep only the needed columns
      select(Overall, var) %>%
      # group by var
      group_by(var) %>%
      # summarize
      summarise(Overall = mean(Overall, na.rm = TRUE)) %>%
      # arrange descending
      arrange(desc(Overall))

    # Now pick out the top variables
    as.character(finalImp$var[1:size])
  }
)

## Splitting
set.seed(123)
split_info <- initial_split(model_dd, strata = datasetID, prop = 0.8)
train <- training(split_info)
test <- testing(split_info)


## Recipe
set.seed(123)
mod_rec <- recipe(Jaccard_dissim ~ . , data = train) %>% # Jaccard
  update_role(verbatimIdentification,
              new_role = "speciesID",
              old_role = "predictor") %>%
  step_impute_knn(all_predictors(), neighbors = 10) %>%
  step_corr(all_numeric_predictors(), threshold = 0.99) # remove those that are super highly correlated (new step)

## rfe
set.seed(123)
rfe_ctrl <- rfeControl(
  functions = ranger_funcs,
  method = "repeatedcv",
  number = 5,
  repeats = 3,
  returnResamp = "all", # we need all resamples
  verbose = TRUE,
  saveDetails = TRUE)

sizes <- seq(from = n_pred, to = 1, by = -1)

tictoc::tic()
set.seed(123)
rf_filter <- rfe(
  form = Jaccard_dissim ~ .,
  na.action = "na.omit",
  mod_rec,
  data = train,
  sizes = sizes,
  rfeControl = rfe_ctrl,
  metric = "RMSE",
  maximze = FALSE
)
tictoc::toc() #3609.2 sec elapsed

rf_filter$optVariables

saveRDS(rf_filter, here::here("Data/3_ml/rfe_jaccard.rds")) # 17 vars // 31 vars
rf_filter <- readRDS(here::here("Data/3_ml/rfe_jaccard.rds")) # 37 vars

