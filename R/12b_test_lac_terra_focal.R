library(terra)
library(here)


w_vec <- c(3L, 5L, 9L, 17L, 33L)
apply_focal_multiple_w <- function(r, w_vec) {

  results <- lapply(w_vec, function(w) {
    f <- focal(r, w = w, fun = sum, na.rm = TRUE)
    # n_S_r = number of filled cells in the box (of size r)
    # Q_S_r = normalized frequency distribution of S for box size r
    box_masses_int <- as.integer(values(f, na.rm = T))
    max_value <- max(box_masses_int, na.rm = TRUE)

    n_S_r <- integer(max_value + 1L)
    for (val in box_masses_int) {
      if (val >= 0 && val <= max_value) {
        n_S_r[val + 1L] <- n_S_r[val + 1L] + 1
      }
    }
    Q_S_r <- n_S_r / length(box_masses_int[!is.na(box_masses_int)])
    S_vals <- seq_along(Q_S_r) - 1
    Z_1 <- sum(S_vals * Q_S_r)
    Z_2 <- sum((S_vals^2) * Q_S_r)

    if (Z_1 == 0) {
      return(NA_real_)
    }
    lac <- Z_2 / (Z_1 * Z_1)

    return(lac)
  })

  names(results) <- paste0("w=", w_vec)
  r_name <- gsub(".tif", "", basename(sources(r)))
  message(r_name)
  lac_res <- data.frame(lac = do.call(rbind, results), w = w_vec, species = r_name)
  return(lac_res)
}


# Test function :
w_vec <- c(3L, 5L, 9L, 17L, 33L)
files <- list.files(pattern = ".tif", here("Data/input/species_ranges_tiff"), full.names = TRUE)


tictoc::tic()

res_list <- list()

for (file_i in seq_along(files)) {
  r <- rast(files[file_i])
  test_res <- apply_focal_multiple_w(r, w_vec)
  res_list[[file_i]] <- test_res
}
tictoc::toc()

df_res <- do.call(rbind(res_list))

saveRDS(df_res, here("Data/output/1_data/3_lacunarity_terra_focal.rds"))
