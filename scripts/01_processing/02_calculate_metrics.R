# WORKFLOW OVERVIEW ------------------------------------------------------------
#
# This script calculates habitat-suitability metrics for target species
# within Mexican TURFs.
#
# 1. Helper functions identify the species and scenario represented by each
#    raster and find the TURFs that target that species.
#
# 2. Metric functions calculate, for each TURF x species x scenario:
#    - mean HSI of suitable cells
#    - percent of the TURF that is suitable
#    - species presence based on the persistence threshold
#
# 3. process_raster_metrics() brings these functions together and calculates
#    the metrics for one species x scenario raster.
#
# 4. All raster files are processed and combined into one species-level
#    metrics table.
#
# 5. Species richness is then calculated from the species-level results,
#    producing one richness value for each TURF x scenario.
#
# Final outputs:
# - all_metrics_combined = TURF x species x scenario metrics
# - species_richness     = TURF x scenario species richness
# -----------------------------------------------------------------------------

# PACKAGES ---------------------------------------------------------------------

library(terra)
library(sf)
library(dplyr)

# FILE PATHS -------------------------------------------------------------------

# folder containing processed SDM rasters
raster_folder <- "data/processed/sdm_rasters"

cutoff_file <- "data/processed/sdm_cutoffs.csv"

# file of cleaned ere TURFs and species 
turf_file <- "data/turfs/mex_turfs_combined.gpkg"

# file where the calculated TURF metrics will be saved.
metrics_output_file <- "data/processed/turf_metrics.csv"

# LOAD TURF DATA ---------------------------------------------------------------

# Load  cleaned TURF polygons and species info
turfs <- st_read(turf_file)

# Load species-specific AquaX suitability thresholds

species_cutoffs <- readr::read_csv(
  cutoff_file,
  show_col_types = FALSE
)
# LIST RASTERS -----------------------------------------------------------------

raster_files <- list.files(
  path = raster_folder,
  pattern = "\\.tif$",
  full.names = TRUE
)

# FUNCTION: get_raster_aphia_id -----------------------------------------------
# gets the aphia id 

get_raster_aphia_id <- function(file) {
  
  sub(
    "_.*",
    "",
    basename(file)
  )
  
}

# FUNCTION: get_raster_scenario -----------------------------------------------
#gets the scenario 

get_raster_scenario <- function(file) {
  
  sub(
    "^[0-9]+_",
    "",
    tools::file_path_sans_ext(
      basename(file)
    )
  )
  
}

# FUNCTION: get_species_turfs -----------------------------------------------
#given an aphiaID, which turfs target that specie

get_species_turfs <- function(target_aphia_id) {
  
  turfs |>
    dplyr::filter(
      aphia_id == target_aphia_id
    )
  
}

# FUNCTION: get_species_cutoff -----------------------------------------------

get_species_cutoff <- function(target_aphia_id) {
  
  species_cutoffs |>
    dplyr::filter(
      as.character(aphia_id) == target_aphia_id
    ) |>
    dplyr::pull(cutoff)
  
}

# FUNCTION: load_matching_raster_data -----------------------------------------
#prepares everything needed for the calculations
#converts the coordinates 

#one function that takes a single raster file
#(one specie, one scenario)
load_matching_raster_data <- function(file) {
  
  aphia_id <- get_raster_aphia_id(file)
  
  cutoff <- get_species_cutoff(aphia_id)
  
  scenario <- get_raster_scenario(file)
  
  raster <- rast(file)
  
  species_turfs <- get_species_turfs(aphia_id)
  if (nrow(species_turfs) == 0) {
    return(NULL)
  }
  
  #transform CRS 
  species_turfs <- st_transform(
    species_turfs,
    crs = crs(raster)
  )
  
  #for each TURF polygon what are the HSI values
  # of every raster cell that falls inside it
  extracted_values <- terra::extract(
    raster,
    terra::vect(species_turfs)
  )
  
  #give the raster-value column the same name
  #for every scenario
  names(extracted_values)[2] <- "hsi"
  
  #starting to create a lookup table
  turf_lookup <- species_turfs |>
    st_drop_geometry() |>
    # Create temporary ID that matches the numbering
    # assigned by terra::extract() so extracted raster
    # values can be linked back to the correct TURFs.
    dplyr::mutate(
      ID = dplyr::row_number()
    ) %>% 
    # Keep only the identifiers needed
    # for the metric calculations
    dplyr::select(
      ID,
      sub_id,
      turf_id,
      aphia_id
    )
  # Add the corresponding TURF identifiers
  # to each extracted raster value
  extracted_values <- extracted_values |>
    dplyr::left_join(
      turf_lookup,
      by = "ID"
    )
  # Bundle all prepared data into a single object
  return(
    list(
      aphia_id = aphia_id,
      scenario = scenario,
      raster = raster,
      species_turfs = species_turfs,
      cutoff = cutoff,
      extracted_values = extracted_values
    )
  )
  
}

# FUNCTION: calculate_mean_hsi -----------------------------------------------

calculate_mean_hsi <- function(extracted_values, cutoff) {
  
  # Create one row per TURF so that TURFs with
  # no suitable cells are still kept in the final output
  all_turfs <- extracted_values |>
    dplyr::distinct(
      sub_id,
      turf_id,
      aphia_id
    )
  
  # Keep only cells above the suitability threshold
  mean_hsi <- extracted_values |>
    dplyr::filter(
      hsi > cutoff 
    ) |>
    # Group the cells by TURF and species
    dplyr::group_by(
      sub_id,
      turf_id,
      aphia_id
    ) |>
    # Calculate average HSI of suitable cells
    # within each TURF, ignoring missing values
    dplyr::summarise(
      mean_hsi = mean(hsi, na.rm = TRUE),
      .groups = "drop"
    )
  # Join calculated means back to the full TURF list
  # so TURFs with no suitable habitat remain as NA 
  # instead of disappearing
  all_turfs |>
    dplyr::left_join(
      mean_hsi,
      by = c("sub_id", "turf_id", "aphia_id")
    )
  
}

# FUNCTION: calculate_percent_suitable ---------------------------------------

calculate_percent_suitable <- function(extracted_values, cutoff) {
  
  extracted_values |>
    dplyr::group_by(
      sub_id,
      turf_id,
      aphia_id
    ) |>
    dplyr::summarise(
      # Count total number of cells
      # with valid habitat suitability values
      percent_suitable =
        100 *
        sum(hsi > cutoff, na.rm = TRUE) /
        sum(!is.na(hsi)),
      
      .groups = "drop"
    )
  
}

# FUNCTION: calculate_species_presence ----------------------------------------

calculate_species_presence <- function(percent_suitable) {
  
  percent_suitable |>
    dplyr::mutate(
      present = percent_suitable > 10
    )
  
}

# FUNCTION: calculate_species_richness ----------------------------------------

calculate_species_richness <- function(species_presence) {
  
  species_presence |>
    dplyr::group_by(
      sub_id,
      turf_id,
      scenario
    ) |>
    dplyr::summarise(
      species_richness = sum(present, na.rm = TRUE),
      .groups = "drop"
    )
  
}

# FUNCTION: process_raster_metrics --------------------------------------------

# Process one species-scenario raster from start to finish 
# and calculate habitat suitability metrics 
# for all TURFs that target that species
process_raster_metrics <- function(file) {
  
  # Prepare the raster, matching TURFs, scenario,
  # and extracted HSI values for this file
  data <- load_matching_raster_data(file)
  
  if (is.null(data)) {
    return(NULL)
  }
  
  #calculate mean HSI
  mean_hsi <- calculate_mean_hsi(
    data$extracted_values,
    data$cutoff
  )
  
  #caluclate percent of suitable habitat
  percent_suitable <- calculate_percent_suitable(
    data$extracted_values,
    data$cutoff
  )
  
  #combine the two metrics tables
  #using the TURF and species identifiers
  metrics <- mean_hsi |>
    dplyr::left_join(
      percent_suitable,
      by = c("sub_id", "turf_id", "aphia_id")
    )
  
  # use percent suitable habitat to classify whether
  # the species is present in each TURF.
  metrics <- calculate_species_presence(
    metrics
  )
  
  #add scenario associated to this raster
  metrics <- metrics |>
    dplyr::mutate(
      scenario = data$scenario
    )
  
  #return completed metric table 
  return(metrics)
  
}

# PROCESS ALL RASTERS ----------------------------------------------------------

all_metrics <- lapply(
  raster_files,
  process_raster_metrics
)
# COMBINE ALL RASTER METRICS --------------------------------------------------

all_metrics_combined <- dplyr::bind_rows(
  all_metrics
)

# CALCULATE SPECIES RICHNESS --------------------------------------------------

species_richness <- calculate_species_richness(
  all_metrics_combined
)

# SAVE OUTPUT TABLES ----------------------------------------------------------

readr::write_csv(
  all_metrics_combined,
  "data/processed/turf_metrics/turf_species_metrics.csv"
)

readr::write_csv(
  species_richness,
  "data/processed/turf_metrics/turf_species_richness.csv"
)

