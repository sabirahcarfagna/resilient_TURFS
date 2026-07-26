# BUILD AQUAX RASTERS -----------------------------------------------------------

# Goal:
# Convert AquaX SDM tables into GeoTIFF rasters cropped to Mexico's EEZ.

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

# LOAD MEXICO EEZ --------------------------------------------------------------

#GOAL:
#load global EEZ polygon then only keep mexico

#read EEZ layer from geopackage file 
all_eez <- st_read(
  mexico_eez_file,
  layer = "eez_v12"
)

#keep only mexico layers
mexico_eez <- all_eez |>
  dplyr::filter(
    SOVEREIGN1 == "Mexico" |
      SOVEREIGN2 == "Mexico" |
      SOVEREIGN3 == "Mexico"
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
  names(sdm)[names(sdm) == "Current"] <- "current"
  names(sdm)[names(sdm) == "RCP26_2050"] <- "rcp26_2050"
  names(sdm)[names(sdm) == "RCP26_2100"] <- "rcp26_2100"
  names(sdm)[names(sdm) == "RCP45_2050"] <- "rcp45_2050"
  names(sdm)[names(sdm) == "RCP45_2100"] <- "rcp45_2100"
  names(sdm)[names(sdm) == "RCP85_2050"] <- "rcp85_2050"
  names(sdm)[names(sdm) == "RCP85_2100"] <- "rcp85_2100"
  
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
  
  cropped_raster <- crop(raster, eez_vector)
  
  masked_raster <- mask(cropped_raster, eez_vector)
  
  return(masked_raster)
  
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
    
    mexico_raster <- crop_to_eez( #crop to mexico eez
      scenario_raster,
      mexico_eez
    )
    
    filename <- paste0( #build file name with aphia and scenario 
      aphia_id,
      "_",
      scenario,
      ".tif"
    )
    
    save_scenario_raster( #save it to outout folder
      mexico_raster,
      filename
    )
    
  } 
  
}

# PROCESS ALL SPECIES ----------------------------------------------------------

for (file in sdm_files) { #now loop over species 
  
  process_species(file)
  
}