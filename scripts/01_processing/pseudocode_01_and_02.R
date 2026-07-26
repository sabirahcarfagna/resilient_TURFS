# PSEUDOCODE FOR MAIN SCRIPT ------------------------------------------

#maybe start with a script for raster makiing and one for metrics 


# GOAL -------------------------------------------------------------------------
# 
# Calculate climate-driven changes in three TURF-level metrics:
# - Mean HSI
# - Percent suitable habitat area
# - Species richness / persistence
# For each:
# - species
# - TURF
# - climate scenario
# - time period
#
# Compare:
# - future conditions (mid- and end-century) to current conditions
# - climate scenarios within each time period

# WORKFLOW ---------------------------------------------------------------------

# 1. Load final cleaned TURF dataset:
#    - data/output/mex_turfs_combined.gpkg
#    - one row per TURF polygon × species combination
#    - includes sub_id, species names, Aphia IDs, and polygon geometry

#create a list of species aphia ids so we do for (aphia_id in species_list)
# sdm <- load_species_sdm(aphia_id)
# species_turfs <- get_species_turfs(turfs, aphia_id)
# ...

# 2. Load AquaX SDM files:
#    - one SDM file per species
#    - includes HSI values for all climate scenarios and time periods
#    - includes the species-specific HSI cutoff

# 3. For each species with an available AquaX SDM:
#    - load the corresponding SDM
#    - extract the corresponding TURF polygons
#    - extract the species-specific HSI cutoff
#    - build rasters for all scenarios and time periods
#    - calculate mean HSI
#    - calculate suitable habitat area and percent suitable habitat
#    - classify presence/absence using τ = 10%

# 4. Combine outputs across all species into a single results table.

# 5. Compare:
#    - future values to current values
#    - climate scenarios within each time period

# FUNCTION 1: load_species_sdm() -----------------------------------------------

# Purpose:
# Load the AquaX SDM for the one species being processed 

# Input:
# - aphia_id 

# Steps:
# 1. Find the SDM file corresponding to the Aphia ID.
#    Eg:
#    - aphia_id = 382891
#    - FINAL_EMSDM_EMMEAN_SP_382891.Rdata

# 2. Check whether the SDM file exists.
#    - if not, stop the function
#    - return message saying no SDM file was found for that aphiaID

# 3. Load the SDM file.
#    - loading the file creates the object FINALEMMEAN

# 4. Verify that FINALEMMEAN contains the variables needed for later analyses
#    - x
#    - y
#    - cutoff
#    - Current
#    - RCP26_2050
#    - RCP45_2050
#    - RCP85_2050
#    - RCP26_2100
#    - RCP45_2100
#    - RCP85_2100
#rename things Current= current, no need for NR and other stuff

# 5. Return the loaded SDM table.

# FUNCTION 2: get_species_turfs() -----------------------------------------------

# Purpose:
# Extract the filtered TURF dataset for the species currently being processed

# Input:
# - turfs:
#   - data/output/mex_turfs_combined.gpkg
# - aphia_id:
#   - Aphia ID for the species currently being processed

# Steps:
# 1. Filter the cleaned TURF dataset to rows where aphia_id matches the
#    species currently being processed

# 2. Transform selected TURFs to the same CRS as the SDM data.
#    - EPSG:4326
#fix my dataset maybe

# 3. Return the filtered TURF dataset

# FUNCTION 3: build_scenario_raster() -------------------------------------------
# export as TIFs, crop them to extent of mex eez 
#create a single thing in the nevironemnt, a single raster 
# Purpose:
# Build one HSI raster for the species, scenario, and time period currently
# being processed

# Inputs:
# - sdm:
#   - SDM table for the species currently being processed
#   - returned by load_species_sdm()
#
# - species_turfs:
#   - filtered TURF dataset for the species currently being processed
#   - returned by get_species_turfs()
#
# - hsi_column: OR 1. period and scenario
#   - HSI column currently being processed (being converted into raster)
#   - eg: Current, RCP45_2050, RCP85_2100 etc

# Steps:
# 1. Start with the loaded SDM table.

# 2. Keep only cells where the selected HSI column is greater than cutoff.

# 3. Keep only:
#    - x
#    - y
#    - selected HSI column

# 4. Rename the selected HSI column to hsi.
#    - to standardize the column name so the same code can be used
#      for every scenario and time period 

# 5. Convert the x-y-hsi table into a raster.

# 6. Crop the raster to the filtered TURF dataset.

# 7. Return the cropped HSI raster.

# FUNCTION 4: calculate_mean_hsi() ----------------------------------------------
#make extraction to polygon level 

# Purpose:
# Calculate the mean HSI within each TURF for one scenario and time period.

# Inputs:
# - scenario_raster:
#   - HSI raster 
#   - returned by build_scenario_raster()

# - species_turfs:
#   - filtered TURF dataset for the species currently being processed
#   - returned by get_species_turfs()

# Steps:
# 1. Extract HSI values from the scenario raster within each 
#    TURF polygon in species_turfs

# 2. Calculate the mean HSI for each TURF.
#    - ignore NA values

# 3. Return a table containing the mean HSI for each TURF.
#could create one function (get metrics) instead of 3 functions
# FUNCTION 5: calculate_suitable_area() -----------------------------------------

# Purpose:
# Calculate the area of suitable habitat within each TURF for one scenario
# and time period.

# Inputs:
# - scenario_raster:
#   - HSI raster 
#   - returned by build_scenario_raster()

# - species_turfs:
#   - filtered TURF dataset
#   - returned by get_species_turfs()

# Steps:
# 1. Calculate area of each raster cell.

# 2. Keep cell-area values only where scenario_raster has suitable habitat.

# 3. Sum suitable cell areas within each TURF in species_turfs.
#    - ignore NA values

# 4. Return a table containing suitable habitat area for each TURF.

#same here, 5 and 6 can be one function 
# FUNCTION 6: calculate_percent_suitable() --------------------------------------

# Purpose:
# Calculate the percentage of each TURF that is suitable habitat for one
# scenario and time period.

# Inputs:
# - suitable_area:
#   - table returned by calculate_suitable_area()

# - species_turfs:
#   - filtered TURF dataset returned by get_species_turfs()

# Steps:
# 1. Calculate the total area of each TURF in species_turfs

# 2. Divide suitable_area by total TURF area.

# 3. Multiply by 100 to get percentage of suitable habitat.

# 4. Return a table containing the percentage of suitable habitat
#   for each TURF.

# FUNCTION 7: classify_presence() ---------------------------------------------

# Purpose:
# Classify each species as present or absent within each TURF for one
# scenario and time period.

# Inputs:
# - percent_suitable:
#   - table returned by calculate_percent_suitable()

# - threshold:
#   - minimum percentage of suitable habitat required for presence
#   - tau = 10% ? 

# Steps:
# 1. Compare percent_suitable to threshold.

# 2. Classify each TURF as:
#    - Present if percent_suitable ≥ threshold
#    - Absent if percent_suitable < threshold

# 3. Return a table containing the presence/absence classification for each TURF.