#----------------------------------------------------------#
#
#
#                     Static Patterns 
#
#                      01_Get_data.R
#                
#
#                    Friederike Wölke 
#                        2025
#
#----------------------------------------------------------#


# Start with clean environment
rm(list = ls())
gc()

# Source 00_Configuration.R
source(here::here("R/00_Configuration.R"))
lapply(package_list, require, character = TRUE)

#----------------------------------------------------------#
# Get atlas data -----
#----------------------------------------------------------#

# Connect to the database
con <- dbConnect(Postgres(),
                 dbname = "MOBI_atlases_v1",
                 host = "localhost",
                 port = 5432,
                 user = "frieda",
                 password = Sys.getenv('PASSWORD_SERVER')
)
dbListTables(con)
#--------------------------------------------------#

# Get filtered data:
tic()
data <- tbl(con, sql(vars$sql_query_data)) %>%
  collect() %>%
  mutate(samplingPeriodID = case_when(
    datasetID == 5 & startYear == 1985 ~ 1,
    datasetID == 5 & startYear == 2001 ~ 2,
    datasetID == 6 & startYear == 1980 ~ 1,
    datasetID == 6 & startYear == 2000 ~ 2,
    datasetID == 13 & startYear == 1974 ~ 1,
    datasetID == 13 & startYear == 1997 ~ 2,
    datasetID == 26 & startYear == 1972 ~ 1,
    datasetID == 26 & startYear == 2013 ~ 2,
    TRUE ~ NA_integer_ # Default case: NA if no match
  )) %>%
  distinct(datasetID, scalingID, siteID, samplingPeriodID, verbatimIdentification, .keep_all = TRUE)
toc() # 6.78 seconds

#--------------------------------------------------#

# Get original grids:
grids <- st_read(con, query = vars$sql_query_grid)


#--------------------------------------------------#

# Create sf that includes cells not sampled at all
data_sf <- grids %>%
  left_join(data)

#--------------------------------------------------#


# Get metadata for datasets
meta <- tbl(con, "MOBI_dataset") %>%
  filter(datasetID %in% c(5, 6, 13, 26)) %>%
  collect()

# disconnect from MOBI db
dbDisconnect(con)

#----------------------------------------------------------#
# Save data -----
#----------------------------------------------------------#

# save METADATA for atlases to documentation file:
write.csv(meta, here("Documentation","METADATA_datasets.csv"))

# write to input data folder
saveRDS(grids, vars$grid)

# filtered by sampling period & croppedArea & recordFilter
saveRDS(data, vars$data)

# filtered by sampling period & croppedArea & recordFilter
#   but with unsampled cells as NA in verbatimIdentification
saveRDS(data_sf, vars$data_sf)



#----------------------------------------------------------#
# Get 'BirdLife 2024' data -----
#----------------------------------------------------------#

# Connect to the database
con <-
  dbConnect(Postgres(),
            dbname = "Birds_of_the_World",
            host = "localhost",
            port = 5432,
            user = "frieda",
            password = Sys.getenv('PASSWORD_SERVER')
)

#--------------------------------------------------#


dbListTables(con)

# Get species names (BirdLife 2024 taxonomy - some not matched)
name_vector_sql <- unique(data$scientificName)

# Ensure name_vector_sql is correctly formatted as an SQL-compatible string
name_vector_sql <- paste0("'", paste(name_vector_sql, collapse = "', '"), "'")

# Construct the query
bl_query <-
  paste0(
    "SELECT \"sci_name\", \"presence\", \"origin\", \"seasonal\",\"geometry\" ",
    "FROM \"MOBI_botw_multipolygon_2024\" ",
    "WHERE \"sci_name\" IN (", name_vector_sql, ") ",
    # Extant, probably extant, possibly extant (not possibly extinct or extinct)
    "AND \"presence\" IN (1, 2, 3)"
)

# For global range maps we need:
# origin = c(1) (native)
# seasonal = c(1,2) (breeding and resident)

# For invasive species we need:
# origin = 3 (introduced)

#--------------------------------------------------#

tic()
BirdLife <- st_read(con, query = bl_query)
toc() # 106.88 sec

# disconnect from BirdLife db
dbDisconnect(con)

#----------------------------------------------------------#
# Save BirdLife24 to shp  -----
#----------------------------------------------------------#
BirdLife_global <- BirdLife %>% filter(origin == 1 & seasonal %in% c(1,2))
BirdLife_introduced <- BirdLife %>% filter(origin == 3)

#--------------------------------------------------#

# save to shp
st_write(BirdLife_global, here::here("Data/input/shp_global/BirdLife_global.shp"), append = F)
st_write(BirdLife_introduced, here::here("Data/input/shp_introduced/BirdLife_introduced.shp"), append = F)
