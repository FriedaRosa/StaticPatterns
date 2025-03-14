rm(list=ls())

# Source 00_Configuration.R
source(here::here("R/00_Configuration.R"))


# Source lacunarity function for matrices:



source(here("R/src/fun_lacunarity_matrix.R"))

# Load required packages
list_p <- c("here", "sf", "tidyverse", "purrr", "furrr", "tictoc", "terra")
lapply(list_p, require, character.only = TRUE)
sf_use_s2(FALSE)

r_vec = c(3, 5, 9, 17, 33)

pattern_CZ <- "^5_"
pattern_NY <- "^6_"
pattern_JP <- "^13_"
pattern_EU <- "^26_"

 files <- list.files(pattern = pattern_NY, here("Data/input/species_ranges_tiff"), full.names = TRUE)
 res_list <- vector("list", length(files))
 for(file_i in seq_along(files)){
    file_i <- 1
    this_file <- files[file_i]
    sp_mat <- as.matrix(rast(this_file), wide = T)
    sp_name <- gsub(pattern = ".tif", "", basename(this_file))

    print(sp_name)
    names(sp_mat) <- sp_name

    # **Run lacunarity function**
    lac_result <- lacunarity_R(sp_mat,
        r_vec = r_vec,
        progress = FALSE,
        ncores = 1L,
        save_plot = FALSE,
        plot = FALSE
    )
    res_list[[file_i]] <- lac_result


 }

