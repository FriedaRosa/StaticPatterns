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

 files <- list.files(pattern = pattern_CZ, here("Data/input/species_ranges_tiff"), full.names = TRUE)
 res_list <- vector("list", length(files))
 for(file_i in seq_along(files)){
    #file_i <- 1
    this_file <- files[file_i]
    sp_mat <- as.matrix(rast(this_file), wide = T)
    sp_name <- gsub(pattern = ".tif", "", basename(this_file))

    print(sp_name)
    names(sp_mat) <- sp_name

    # **Run lacunarity function**
    lac_result <- lacunarity_R(sp_mat,
        r_vec = c(3, 5, 9, 17, 33),
        ncores = 1L,
        save_plot = FALSE,
        plot = FALSE
    )
    res_list[[file_i]] <- lac_result

    gc()
    rm(list = c("this_file", "sp_mat", "sp_name", "lac_result"))


 }








 ###########
 ###########

 # Clean Environment
 rm(list = ls())

 # Load Required Libraries
 library(here)
 library(sf)
 library(tidyverse)
 library(purrr)
 #library(furrr)
 library(tictoc)
 library(terra)

 # Disable spherical processing for `sf`
 sf_use_s2(FALSE)

 # Load Configurations & Functions
 source(here::here("R/00_Configuration.R"))
 source(here("R/src/fun_lacunarity_matrix.R"))

 # Define Analysis Parameters
 r_vec <- c(3, 5, 9, 17, 33)
 #pattern_CZ <- "^5_"

 # Get List of Relevant Raster Files
 files <- list.files(path = here("Data/input/species_ranges_tiff"),
                     pattern = ".tiff", full.names = TRUE)

 # Define Function to Process Each File
 process_file <- function(file_path) {
     sp_mat <- as.matrix(rast(file_path), wide = TRUE)  # Convert raster to matrix
     sp_name <- gsub(pattern = ".tif", "", basename(file_path))

     message("Processing: ", sp_name)

     # Run lacunarity function
     # lac_result <- lacunarity_R(sp_mat, r_vec = r_vec, ncores = parallel::detectCores() - 1L)
     lac_result <- lacunarity_R(sp_mat, r_vec = r_vec, ncores = 1L)


     gc()  # Clean memory
     return(lac_result)
 }

 # Enable Parallel Processing
 # plan(multisession, workers = parallel::detectCores() - 1L)

 # Process Files in Parallel
 tic("Processing All Files")
 # res_list <- future_map(files, process_file)
 res_list <- map(files, process_file)
 toc()

 # Disable Parallel Backend
 # plan(sequential)

 # Save Results
 saveRDS(res_list, here("Data/output/lacunarity_results_16_03_25.rds"))


 ###########
