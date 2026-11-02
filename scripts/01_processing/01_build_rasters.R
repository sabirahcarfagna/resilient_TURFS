# BUILD AQUAX RASTERS -----------------------------------------------------------

# Goal:
# Convert AquaX SDM tables into GeoTIFF rasters covering the selected
# Latin American EEZs, including Mexico, Central America, and South America.

# PACKAGES ---------------------------------------------------------------------

library(terra) #handles rasters
library(sf) #handles vector data

# FILE PATHS -------------------------------------------------------------------

# GOAL: 
# defining locations of input and output files 

sdm_folder <- "data/raw/2026_Sabira_JC_VILLASENOR_2/SDM"

mexico_eez_file <- "data/eez_v12.gpkg"

#where rasters will be saved
raster_output_folder <- "data/processed/sdm_rasters"

cutoff_output_file <- "data/processed/sdm_cutoffs.csv"

# LOAD STUDY AREA EEZS ---------------------------------------------------------

all_eez <- st_read(
  mexico_eez_file,
  layer = "eez_v12"
)

countries <- c(
  "Mexico",
  "Belize",
  "Guatemala",
  "El Salvador",
  "Honduras",
  "Nicaragua",
  "Costa Rica",
  "Panama",
  "Colombia",
  "Venezuela",
  "Guyana",
  "Suriname",
  "Brazil",
  "Ecuador",
  "Peru",
  "Chile",
  "Argentina",
  "Uruguay"
)

study_eez <- all_eez |>
  dplyr::filter(
    SOVEREIGN1 %in% countries |
      SOVEREIGN2 %in% countries |
      SOVEREIGN3 %in% countries
  )

# LIST SDM FILES ---------------------------------------------------------------

#GOAL: 
# Create a list of all SDM files to be processed later 
# keep files ending in .Rdata
sdm_files <- list.files(
  path = sdm_folder,
  pattern = "\\.Rdata$",
  full.names = TRUE #store complete file path 
)

# FUNCTION: get_aphia_id -------------------------------------------------------

#GOAL:
# extract the aphiaID from the filename to use it later when saving
get_aphia_id <- function(file) {
  
  sub(
    "FINAL_EMSDM_EMMEAN_SP_([0-9]+)\\.Rdata",
    "\\1",
    basename(file)
  )
  
}

# FUNCTION: load_species_sdm ---------------------------------------------------

#GOAL: 
# load a specie SDM, standardize variable names 

load_species_sdm <- function(file) { #function(input)
  
  load(file) 
  
  sdm <- FINALEMMEAN #now finalemmean is saved as sdm 
  
  rm(FINALEMMEAN) #remove finalemmean from env. 
  
  #renaming all variables to lowercase
  # standardize variable names
  sdm <- janitor::clean_names(sdm)
  
  return(sdm)
}

# FUNCTION: build_scenario_raster ----------------------------------------------

build_scenario_raster <- function(sdm, scenario) {
  
  #keep only x, y and HS values for selected scenario
  scenario_data <- sdm[, c("x", "y", scenario)]
  
  #convert x,y and HS values into SpatRaster
  scenario_raster <- rast(
    scenario_data,#defined above
    type = "xyz",
    crs = "EPSG:4326"
  )
  
  return(scenario_raster)
}

# FUNCTION: crop_to_eez --------------------------------------------------------
crop_to_eez <- function(raster, eez) {
  
  eez_vector <- vect(eez)
  
  cropped_raster <- crop(
    raster,
    eez_vector
  )
  
  return(cropped_raster)
}

# FUNCTION: save_scenario_raster -----------------------------------------------

save_scenario_raster <- function(raster, filename) {
  
  output_file <- file.path(
    raster_output_folder,
    filename
  )
  
  writeRaster(
    raster,
    output_file,
    overwrite = TRUE
  )
  
}

# FUNCTION: process_species ----------------------------------------------------

process_species <- function(file) {
  
  aphia_id <- get_aphia_id(file) #takes full file name and extracts aphiaID
  
  sdm <- load_species_sdm(file) #loads file, renames cols and returns it as sdm
  
  cutoff <- unique(sdm$cutoff)
  
  scenarios <- c( #create list of scenarios 
    "current",
    "rcp26_2050",
    "rcp26_2100",
    "rcp45_2050",
    "rcp45_2100",
    "rcp85_2050",
    "rcp85_2100"
  )
  
  for (scenario in scenarios) { #loop over scenarios
    
    scenario_raster <- build_scenario_raster( #build the raster
      sdm,
      scenario
    )
    
    study_raster <- crop_to_eez( #crop to eez's
      scenario_raster,
      study_eez
    )
    
    filename <- paste0( #build file name with aphia and scenario 
      aphia_id,
      "_",
      scenario,
      ".tif"
    )
    
    save_scenario_raster(
      study_raster,
      filename
    )
    
  }
  
  return(
    data.frame(
      aphia_id = aphia_id,
      cutoff = cutoff
    )
  )
  
}

  # PROCESS ALL SPECIES ----------------------------------------------------------
  
  species_cutoffs <- lapply(
    sdm_files,
    process_species
  )
  
  species_cutoffs <- dplyr::bind_rows(
    species_cutoffs
  )
  
  readr::write_csv(
    species_cutoffs,
    cutoff_output_file
  )
  
  