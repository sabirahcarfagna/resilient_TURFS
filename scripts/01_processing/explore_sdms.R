################################################################################
# Explore Gabriel SDMs
################################################################################

# Goal:
# One by one, load each species SDM file, build a current suitability raster,
# crop it to Mexico's EEZ, and plot it with all TURFs.
# This is only to check whether Gabriel's SDM data make spatial sense.

# SET UP -----------------------------------------------------------------------

library(here)
library(tidyverse)
library(terra)
library(sf)

# Load all TURFs
all_turfs <- read_sf("https://github.com/mex-fisheries/mex-TURFs/raw/refs/heads/master/data/output/mex_turfs_combined.gpkg") %>%
  st_transform(crs = "EPSG:4326")

#List all 20 SDM files
#list.files to look inside folder and list the files found there
sdm_files <- list.files(
  path = here("data/raw/2026_Sabira_JC_VILLASENOR_2/SDM"),
  pattern = "\\.Rdata$",
  full.names = TRUE
)

sdm_files #each element is a path to one file 
load(sdm_files[1]) #now create object from what was saved in file path 1
names(FINALEMMEAN) #see all columns in FINALEMMEAN
sdm <- FINALEMMEAN #working copy of FINALEMMEAN
#every time I run a new sdm_files[], it overwrites the previous
#I analyze it
#save the results 
#clear the workbench
#put the next species on 
