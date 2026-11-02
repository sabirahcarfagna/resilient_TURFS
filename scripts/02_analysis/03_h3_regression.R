library(terra)
library(sf)
library(dplyr)
library(readr)
library(fixest)

# H3: EQUATORWARD RANGE POSITION AND HABITAT SUITABILITY CHANGE ---------------


# LOAD DATA -------------------------------------------------------------------

#SDM rasters needed to know wherfe each spp currently has SH
raster_dir <- "data/processed/sdm_rasters"

#Species specific cutoff to know which cells count as suitable
cutoffs <- readr::read_csv(
  "data/processed/sdm_cutoffs.csv",
  show_col_types = FALSE
)

turf_metrics <- readr::read_csv(
  "data/processed/turf_metrics/turf_species_metrics_exact.csv",
  show_col_types = FALSE
)


# EQUATORWARD REFERENCE LATITUDE FOR EACH SPECIES ------------------------------

#function to identify latitude of equatorward end of each species' suitable habitat

# go through every AphiaID in the cutoff CSV
# for each species, use its AphiaID to find its current SDM and its corresponding cutoff.
species_edges <- lapply(cutoffs$aphia_id, function(aphia) {
  
  #find the cutoff associated with this ID in the cutoff CSV
  threshold <- cutoffs$cutoff[cutoffs$aphia_id == aphia]
  
  #find and open the current SDM file for this ID ("aphiaID_current.tif")
  r <- rast(file.path(raster_dir, paste0(aphia, "_current.tif")))
  
  #compare each cell's HSI to the species-specific cutoff
  #create a raster where T/F is suitable/unsuitable.
  suitable <- r > threshold
  
  #convert suitability raster into table containing each cell's 
  #long, lat, and suitability status
  #remove cells with NA values
  #keep only TRUE cells (suitable)
  cells <- as.data.frame(suitable, xy = TRUE, na.rm = TRUE) %>%
    filter(current == TRUE)
  
  #if table has no suitable cells (like purple sea urchin)
  #return a row with ID and NA 
  if (nrow(cells) == 0) {
    return(data.frame(
      aphia_id = aphia,
      edge_latitude = NA_real_,
      crosses_equator = NA
    ))
  }
  
  #find southernmost and northernmost latitudes where specie has suitable habitat
  min_lat <- min(cells$y)
  max_lat <- max(cells$y)
  
  #if min latitude is negative and max is positive
  #the specie crosses  the equator, so "crosses" is TRUE
  crosses <- min_lat < 0 && max_lat > 0
  
  #detemrine which lat to use as the species equatorward edge 
  
  #if specie crosses the equator, use 0 deg (equator)
  edge <- if (crosses) {
    0
    
  #otherwise use 5th percentile of suitable latitudes
  } else if (min_lat >= 0) {
    quantile(cells$y, 0.05)
  }
  
  #create table with one row for current specie
  #save aphiaID, reference latitude, and whether suitable habitat crosses equator
  data.frame(
    aphia_id = aphia,
    edge_latitude = as.numeric(edge),
    crosses_equator = crosses
  )
  
  #finish loop through all AphiaIDs and combine in one table
  #called species_edges
}) %>%
  bind_rows()


# TURF COORDINATES ------------------------------------------------------------

#open geopkg with mex turf boundaries
turfs <- st_read("data/turfs/mex_turfs_combined.gpkg")

#create table containing one location for each turf 
turf_coords <- turfs %>%
  #keep only subid and geom
  select(sub_id) %>%
  #project into ESPG6933
  st_transform(6933) %>%
  #group all polygons belonging to same TURF (some turfs have multiple subid)
  group_by(sub_id) %>%
  summarise() %>%
  #find center of turf 
  st_centroid() %>%
  #convert the coords back to long and lat
  st_transform(4326) %>%
  mutate(
    longitude = st_coordinates(.)[, 1],
    latitude = st_coordinates(.)[, 2]
  ) %>%
  #remove geometry, leaving just TurfID and coordinates
  st_drop_geometry()


# DISTANCE FROM EQUATORWARD REFERENCE -----------------------------------------


# keep the TURF-species pairs under current conditions
h3_distances <- turf_metrics %>%
  filter(scenario == "current") %>%
  select(sub_id, aphia_id) %>%
  #add each turfs centroid latitude
  left_join(
    turf_coords %>% select(sub_id, latitude),
    by = "sub_id"
  ) %>%
  #add each species equatorward reference latitude
  left_join(species_edges, by = "aphia_id") %>%
  # subtract species' reference latitude from the TURF's latitude.
  # the difference is in degrees of latitude, so multiply by 111.2
  # to convert it into approximate kilometers
  # +ive values mean TURF is farther north than the reference
  # -ive values mean its farther south.
  mutate(
    distance_km = (latitude - edge_latitude) * 111.2
  )


# CHANGE IN MEAN HSI ----------------------------------------------------------

# create temporary table containing current mean HSI for each
# TURF-species pair. Rename mean_hsi to current_hsi so we can distinguish
# the starting HSI from the future HSI when calculating changes.
h3_baseline <- turf_metrics %>%
  filter(scenario == "current") %>%
  select(sub_id, aphia_id, current_hsi = mean_hsi)

#create dataset containing future observations and their HSI changes
#start by keeping only future scenarios from turf_metrics.
h3_data <- turf_metrics %>%
  filter(scenario != "current") %>%
  # add current HSI from h3_baseline to each future observation,
  # matching observations belonging to the same TURF and species.
  # so I can compare future HSI with the starting HSI.
  left_join(h3_baseline, by = c("sub_id", "aphia_id")) %>%
  left_join(
    h3_distances %>% select(sub_id, aphia_id, distance_km),
    by = c("sub_id", "aphia_id")
  ) %>%
  #current HSI - Future HSI 
  #divide distance by 100 so regression measure effect 
  #of being 100km farther from references 
  mutate(
    change_hsi = mean_hsi - current_hsi,
    distance_100km = distance_km / 100
  )


# PREPARE REGRESSION DATA ------------------------------------------------------

#create final dataset 
#start with h3_data (contains future HSI changes and distance from reference latitude)
h3_regression <- h3_data %>%
  #add long and lat of each TURF centroid 
  #towe need these for Conley SEs 
  left_join(
    turf_coords %>% select(sub_id, longitude, latitude),
    by = "sub_id"
  ) %>%
  #remove missing values 
  filter(
    !is.na(change_hsi),
    !is.na(distance_100km),
    !is.na(latitude),
    !is.na(longitude)
  ) %>%
  #extract everything before the underscore in the scenario name
  #rcp85_2100 becomes "2100"
  #set 2050 as reference period for regression 
  mutate(
    period = factor(
      sub(".*_", "", scenario),
      levels = c("2050", "2100")
    ),
    #rcp85_2100 becomes "rcp85"
    #set rcp26 as reference emission scenario
    rcp = factor(
      sub("_.*", "", scenario),
      levels = c("rcp26", "rcp45", "rcp85")
    )
  )

# H3 REGRESSION ---------------------------------------------------------------

# Estimate one common relationship between distance from the equatorward
# reference and change in mean HSI, accounting for differences between
# species, future periods, and emissions scenarios.
h3_model <- feols(
  change_hsi ~ distance_100km | aphia_id^period + rcp,
  data = h3_regression,
  vcov = conley(100) ~ latitude + longitude
)

summary(h3_model)

# H3 REGRESSION RESULTS TABLE -------------------------------------------------

h3_results <- broom::tidy(h3_model, conf.int = TRUE) %>%
  mutate(
    Metric = "Change in mean HSI",
    N = nobs(h3_model)
  )

View(h3_results)

# Compare the northernmost and southernmost Panulirus interruptus TURFs.

h3_regression %>%
  filter(aphia_id == 382895) %>%
  distinct(sub_id, latitude) %>%
  arrange(latitude) %>%
  summarise(
    southern_turf = first(sub_id),
    northern_turf = last(sub_id),
    southern_latitude = first(latitude),
    northern_latitude = last(latitude),
    distance_km = (last(latitude) - first(latitude)) * 111.2,
    predicted_hsi_difference = distance_km / 100 * 5.853506
  )