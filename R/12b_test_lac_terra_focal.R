#----------------------------------------------------------#
#
#            Static Patterns - Lacunarity Calculation
#
#                12b_Predictors_Lacunarity_calculation.R
#
#               Friederike Wölke, 2025
#
#----------------------------------------------------------#


#----------------------------------------------------------#
# Install and Load Libraries -----
#----------------------------------------------------------#
rm(list=ls())  # Clear workspace
library(terra)
library(here)
library(ggplot2)
library(dplyr)

#----------------------------------------------------------#
# Function to run terra::focal() for multiple window-sizes -----
#----------------------------------------------------------#
apply_focal_multiple_w <- function(r, w_vec) {

  results <- lapply(w_vec, function(w) {
    f <- focal(r, w = w, fun = sum, na.policy="omit", na.rm=TRUE)

    # n_S_r = number of filled cells in the box (of size r)
    # Q_S_r = normalized frequency distribution of S for box size r

    box_masses_int <- as.integer(values(f, na.rm = T))
    max_value <- max(box_masses_int, na.rm = TRUE)

    # initialize empty vector (add +1 to enable indexing with 0 values)
    n_S_r <- integer(max_value + 1L)
    for (val in box_masses_int) {
      if (val >= 0 && val <= max_value) {
        n_S_r[val + 1L] <- n_S_r[val + 1L] + 1
      }
    }
    Q_S_r <- n_S_r / length(box_masses_int[!is.na(box_masses_int)])
    S_vals <- seq_along(Q_S_r) - 1
    Z_1 <- sum(S_vals * Q_S_r) #first moment = mean
    Z_2 <- sum((S_vals^2) * Q_S_r) # second moment = variance (mean of squared values)

    if (Z_1 == 0) {
      return(NA_real_)
    }
    lac <- Z_2 / (Z_1 * Z_1)

    return(lac)
  })

  names(results) <- paste0("w=", w_vec)
  r_name <- gsub(".tif", "", basename(sources(r)))
  message(r_name, ": done")
  # lac_res <- data.frame(lac = do.call(rbind, results), w = w_vec, species = r_name)
  # return(lac_res)


  this_out <- dplyr::tibble(name = rep(r_name, length(results)),
                            r = w_vec, "ln(r)" = log(w_vec),
                            Lac = results %>% do.call(rbind,.) %>% .[,1], "ln(Lac)" = log(results %>% do.call(rbind,.) %>% .[,1]))
  return(this_out)
}

#----------------------------------------------------------#
# Run focal function across multiple window-sizes -----
#----------------------------------------------------------#

w_vec <- c(3L, 5L, 9L, 17L, 33L)
all_pattern <- ".tif"
files <- list.files(pattern = all_pattern,
                    here("Data/input/species_ranges_tiff"),
                    full.names = TRUE)


tictoc::tic()
res_list <- list()
for (file_i in seq_along(files)) {
  r <- rast(files[file_i])
  test_res <- apply_focal_multiple_w(r, w_vec)
# add atlas and species info
  name_string <- strsplit(test_res$name, "_")[[1]]
  test_res$datasetID <- name_string[1]
  test_res$samplingPeriodID <- name_string[2]
  sp_name <- paste0(name_string[-c(1,2)], collapse = "_")
  test_res$verbatimIdentification <- sp_name

  res_list[[file_i]] <- test_res
}
tictoc::toc()

df_res <- do.call(rbind, res_list) %>%
  mutate(verbatimIdentification = gsub("_", " ", verbatimIdentification))


df_res <- df_res %>%
  mutate(verbatimIdentification =
           case_when(
             verbatimIdentification == "Anas platyrhynchos x A  rubripes" ~ "Anas platyrhynchos x A. rubripes",
             verbatimIdentification == "Vermivora chrysoptera x V  pinus" ~ "Vermivora chrysoptera x V. pinus",
             verbatimIdentification == "Columba livia f  domestica" ~ "Columba livia f. domestica",
             verbatimIdentification == "Vermivora pinus x V  chrysoptera" ~ "Vermivora pinus x V. chrysoptera",
             .default = verbatimIdentification
           ))

saveRDS(df_res, here("Data/output/1_data/single_predictors/2_lacunarity_terra_focal.rds"))





