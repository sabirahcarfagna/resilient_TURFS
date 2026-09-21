# WORKFLOW OVERVIEW ------------------------------------------------------------
#
# This script calculates habitat-suitability metrics for target species
# within Mexican TURFs.
#
# 1. "Helper" functions identify the species and scenario represented by each
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

# Load species-specific AquaX HSI cutoffs
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

# FUNCTION: get_modeled_turf_areas --------------------------------------------
# Creates one unique polygon for each spatial sub-ID included
# in the existing AquaX analysis.
#
# Unlike get_species_turfs(), these polygons are NOT filtered
# according to which species the TURF currently targets.
#
# These areas will later be used to evaluate every modeled
# AquaX species inside every modeled sub-ID.

get_modeled_turf_areas <- function() {
  
  # Aphia IDs for which we actually have AquaX rasters
  modeled_aphia_ids <- unique(
    vapply(
      raster_files,
      get_raster_aphia_id,
      character(1)
    )
  )
  
  # Identify the sub-IDs already represented in the
  # target-species AquaX analysis
  modeled_sub_ids <- turfs |>
    st_drop_geometry() |>
    dplyr::filter(
      as.character(aphia_id) %in% modeled_aphia_ids
    ) |>
    dplyr::distinct(sub_id) |>
    dplyr::pull(sub_id)
  
  # Keep one spatial geometry per modeled sub-ID.
  #
  # A sub-ID can occur on multiple rows because it can target
  # multiple species, so we group and union its geometry.
  modeled_turf_areas <- turfs |>
    dplyr::filter(
      sub_id %in% modeled_sub_ids
    ) |>
    dplyr::group_by(
      sub_id
    ) |>
    dplyr::summarise(
      turf_id = dplyr::first(turf_id),
      .groups = "drop"
    )
  
  return(modeled_turf_areas)
  
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

# FUNCTION: load_all_turf_raster_data -----------------------------------------
# Prepares one species x scenario raster for ALL modeled TURF areas,
# regardless of whether each TURF currently targets that species.
#
# Allows to detect potential suitable habitat for species
# that are not currently targeted by a TURF.

load_all_turf_raster_data <- function(file) {
  
  aphia_id <- get_raster_aphia_id(file)
  
  cutoff <- get_species_cutoff(aphia_id)
  
  scenario <- get_raster_scenario(file)
  
  raster <- rast(file)
  
  # Use all 36 modeled spatial sub-IDs rather than only
  # the TURFs that currently target this species
  all_turfs <- get_modeled_turf_areas()
  
  # Transform TURFs to match raster CRS
  all_turfs <- st_transform(
    all_turfs,
    crs = crs(raster)
  )
  
  # Extract HSI values from this species raster
  # within every modeled TURF area
  extracted_values <- terra::extract(
    raster,
    terra::vect(all_turfs)
  )
  
  # Give raster-value column a consistent name
  names(extracted_values)[2] <- "hsi"
  
  # Create lookup table connecting terra extraction IDs
  # back to the correct spatial sub-ID
  turf_lookup <- all_turfs |>
    st_drop_geometry() |>
    dplyr::mutate(
      ID = dplyr::row_number()
    ) |>
    dplyr::select(
      ID,
      sub_id,
      turf_id
    )
  
  extracted_values <- extracted_values |>
    dplyr::left_join(
      turf_lookup,
      by = "ID"
    ) |>
    dplyr::mutate(
      aphia_id = aphia_id
    )
  
  return(
    list(
      aphia_id = aphia_id,
      scenario = scenario,
      raster = raster,
      all_turfs = all_turfs,
      cutoff = cutoff,
      extracted_values = extracted_values
    )
  )
  
}

# FUNCTION: calculate_mean_hsi -----------------------------------------------

calculate_mean_hsi <- function(extracted_values, cutoff) {
  
  extracted_values |>
    dplyr::group_by(
      sub_id,
      turf_id,
      aphia_id
    ) |>
    dplyr::summarise(
      
      # Number of cells with an AquaX HSI prediction
      n_valid_cells = sum(!is.na(hsi)),
      
      # Mean HSI of suitable cells
      mean_hsi = {
        
        if (n_valid_cells == 0) {
          
          # No numeric AquaX prediction is available
          NA_real_
          
        } else if (sum(hsi > cutoff, na.rm = TRUE) == 0) {
          
          # AquaX predictions exist, but no cells are suitable
          0
          
        } else {
          
          # At least one suitable cell exists
          mean(
            hsi[hsi > cutoff],
            na.rm = TRUE
          )
          
        }
        
      },
      
      .groups = "drop"
    ) |>
    dplyr::select(
      -n_valid_cells
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
      
      n_valid_cells = sum(!is.na(hsi)),
      
      percent_suitable = {
        
        if (n_valid_cells == 0) {
          
          NA_real_
          
        } else {
          
          100 *
            sum(hsi > cutoff, na.rm = TRUE) /
            n_valid_cells
          
        }
        
      },
      
      .groups = "drop"
    ) |>
    dplyr::select(
      -n_valid_cells
    )
}

# FUNCTION: calculate_species_presence ----------------------------------------

calculate_species_presence <- function(percent_suitable) {
  
  percent_suitable |>
    dplyr::mutate(
      present = percent_suitable > 0
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

# FUNCTION: process_all_turf_metrics ------------------------------------------
# Process one species-scenario raster across ALL modeled TURF areas,
# regardless of whether each TURF currently targets that species.
#
# This produces percent suitable habitat for every
# modeled TURF x modeled species x scenario combination.

process_all_turf_metrics <- function(file) {
  
  # Prepare this species raster for all 36 modeled TURFs
  data <- load_all_turf_raster_data(file)
  
  # Calculate percent suitable habitat
  percent_suitable <- calculate_percent_suitable(
    data$extracted_values,
    data$cutoff
  )
  
  # Classify species as having suitable habitat
  # whenever percent suitable > 0
  metrics <- calculate_species_presence(
    percent_suitable
  )
  
  # Add scenario
  metrics <- metrics |>
    dplyr::mutate(
      scenario = data$scenario
    )
  
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

# CALCULATE POTENTIAL SPECIES REDISTRIBUTION ---------------------------------
#
# Evaluate every modeled AquaX species within every modeled TURF sub-ID,
# regardless of whether that species is currently targeted there.
#
# 36 sub-IDs x 20 modeled species x 7 scenarios = 5,040 observations.

potential_metrics <- lapply(
  raster_files,
  process_all_turf_metrics
)

potential_metrics_combined <- dplyr::bind_rows(
  potential_metrics
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

