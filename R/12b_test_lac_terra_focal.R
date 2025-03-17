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
  res_list[[file_i]] <- test_res
}
tictoc::toc()

df_res <- do.call(rbind, res_list)

saveRDS(df_res, here("Data/output/1_data/3_lacunarity_terra_focal.rds"))
df_res <- readRDS(here("Data/output/1_data/3_lacunarity_terra_focal.rds"))



df_res2 <- df_res %>%
  group_by(name) %>%
  mutate(datasetID = as.factor(strsplit(name, "_")[[1]][1]),
         samplingPeriodID = as.factor(strsplit(name, "_")[[1]][2]),
         verbatimIdentification = paste0(strsplit(name, "_")[[1]][3], "_", strsplit(name, "_")[[1]][4]))


ggplot(df_res2, aes(y = Lac, x = r, col = verbatimIdentification))+
  geom_line(aes(group = name), alpha = 0.4, show.legend = FALSE)+
  scale_y_log10()+
  scale_x_log10()+
  facet_wrap(~datasetID)





library(ggplot2)
library(dplyr)
library(terra)





# Help: https://r.geocompx.org/spatial-operations

w <- matrix(1, nrow = 3, ncol = 3)
r <- rast(files[1])
f <- focal(r, w = w, fun = "sum", na.rm = T, na.policy = "omit", expand = TRUE, fillvalue = NaN)

plot(r)
plot(f)

r_df <- as.data.frame(r, xy = TRUE)
f_df <- as.data.frame(f, xy = TRUE)

ggplot(r_df, aes(x=x, y=y, fill = as.integer(presence)))+
  geom_tile(aes(fill = as.integer(presence)))+
  scale_fill_viridis_c() +
  theme_minimal()

ggplot(f_df, aes(x=x, y=y, fill = as.integer(focal_sum)))+
  geom_tile(aes(fill = (focal_sum)))+
  scale_fill_viridis_c() +
  theme_minimal()
