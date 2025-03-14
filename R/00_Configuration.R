#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#                      00_Configuration.R
#                
#
#                    Friederike Wölke 
#                        2025
#
#----------------------------------------------------------#


#----------------------------------------------------------#
# Install and load libraries -----
#----------------------------------------------------------#
# make list of packages to install
package_list <- c("RPostgres", "here",

                  "tidyverse", "tidyr", "dplyr",

                  "ggplot2", "ggthemes", "ggfortify", "plotly",

                  "spdep", "geosphere", "geodata",

                  "purrr", "broom",

                  "sf", "tictoc", "skimr", "kableExtra",

                  "readxl", "terra",  "taxize",

                  "ape", "phyloregion",

                  "doParallel")




# install packages
# lapply(
#   package_list, utils::install.packages
# )

# lapply(
#   package_list, require, character = TRUE
# )

# Check if everything worked:
if (
  isTRUE(
    all(
      package_list %in%
        as.data.frame(
          utils::installed.packages()
        )[, 1]
    )
  )
) {
  cat("Everything is good to go")
} else {
  warning("All required packages are not installed")
  warning(
    paste(
      "Please install the following packages:",
      paste(
        setdiff(
          package_list,
          as.data.frame(
            utils::installed.packages()
          )[, 1]
        ),
        collapse = ", "
      )
    )
  )
}


#----------------------------------------------------------#
# Define variables and paths  -----
#----------------------------------------------------------#

vars <-
  list(
    predictors =
      "c:/Users/wolke/OneDrive - CZU v Praze/Dokumenty/PhD_Projects/Project-folders/StaticPredictors/Data/",
    out =
      "c:/Users/wolke/OneDrive - CZU v Praze/Dokumenty/Documents/GitHub/Project-folders/StaticPatterns/Data/output/1_data/",
    data_sf =
      "c:/Users/wolke/OneDrive - CZU v Praze/Dokumenty/Documents/GitHub/Project-folders/StaticPatterns/Data/input/data_sf.rds",
    data =
      "c:/Users/wolke/OneDrive - CZU v Praze/Dokumenty/Documents/GitHub/Project-folders/StaticPatterns/Data/input/data.rds",
    grid =
      "c:/Users/wolke/OneDrive - CZU v Praze/Dokumenty/Documents/GitHub/Project-folders/StaticPatterns/Data/input/grid.rds",
    atlas_names =
      setNames(
        factor(c(5, 6, 13, 26)),
        c("Czechia", "NewYork", "Japan", "Europe")
      ),
    time_periods =
      c(1, 2),
    desired_levels =
      factor(c("1", "2", "4", "8", "16", "32", "64", "128"),
        ordered = T,
        levels = c("1", "2", "4", "8", "16", "32", "64", "128")
      ),
    crs =
      c(
        "Czechia" = "epsg:5514",
        "NewYork" = "epsg:32118",
        "Japan" = "epsg:6684",
        "Europe" = "epsg:3035"
      ),
    # SQL queries
    sql_query_grid =
      "SELECT * FROM \"MOBI_vw_FINAL_site_metrics\" WHERE \"datasetID\" IN (5,6,13,26)",
    sql_query_data =
      "SELECT
  \"datasetID\", \"scalingID\", \"siteID\", \"startYear\",
  \"croppedArea\", \"verbatimIdentification\", \"scientificName\",
  \"centroidDecimalLongitude\", \"centroidDecimalLatitude\"
  FROM
  \"MOBI_vw_FINAL_presence_records\"
  WHERE
  \"croppedArea\" IS NOT NULL
  AND
  (\"datasetID\" <> 26 OR \"recordFilter\" IN (1, 2))
  AND
  (
    (\"datasetID\" = 5 AND \"startYear\" IN (1985, 2001))
    OR (\"datasetID\" = 6 AND \"startYear\" IN (1980, 2000))
    OR (\"datasetID\" = 13 AND \"startYear\" IN (1974, 1997))
    OR (\"datasetID\" = 26 AND \"startYear\" IN (1972, 2013))
  )",
    Documentation =
      "c:/Users/wolke/OneDrive - CZU v Praze/Dokumenty/Documents/GitHub/Project-folders/StaticPatterns/Documentation/"
  )
