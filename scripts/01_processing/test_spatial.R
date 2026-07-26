################################################################################
# title
################################################################################
#
# You Name Here
# Your Email Here
# date
#
# Description of what this script will do
#
################################################################################
# Goal:
# For each target species, use SDMs to calculate:
# 1. Mean habitat suitability inside each TURF
# 2. Area of suitable habitat inside each TURF
# 3. Change from current to future scenarios
#
# General workflow:
# 1. Load packages needed for file paths, data cleaning, rasters, and spatial polygons.
# 2. Load Gabriel's SDM data for one species.
# 3. Rename the SDM object from FINALEMMEAN -> species specific object
# 4. Load the combined Mexican TURF polygons.
# 5. Filter TURFs to the species being analyzed.
# 6. Transform TURFs to the same CRS as SDM data.
# 7. Build the current HSI raster:
#    a. Start from the species SDM table.
#    b. Keep only cells where Current HSI is greater than the cutoff.
#    c. Keep x, y, and Current HSI.
#    d. Rename Current to current.
#    e. Convert the x-y-current table into a raster.
#    f. Crop the raster to the species TURF area.
# 8. Build the mid-century HSI raster:
#    a. Start from the same species SDM table.
#    b. Keep only cells where RCP45_2050 HSI is greater than the cutoff.
#    c. Keep x, y, and RCP45_2050 HSI.
#    d. Rename RCP45_2050 to mid.
#    e. Convert the x-y-mid table into a raster.
#    f. Crop the raster to the species TURF area.
# 9. Plot current HSI raster with TURF outlines.
# 10. Plot mid-century HSI raster with TURF outlines.
# 11. Plot the change raster, mid minus current, with TURF outlines.
# 12. Extract mean current HSI inside each TURF.
# 13. Extract mean mid-century HSI inside each TURF.
# 14. Join current and mid-century mean HSI values into one table.
# 15. Convert the current raster into a cell-area raster.
# 16. Remove area values where current habitat is unsuitable or missing.
# 17. Sum suitable cell areas inside each TURF to calculate current suitable habitat area.
# 18. Repeat the same area workflow for mid-century.
# 19. Later, repeat the whole workflow for all species using a function.

# SET UP #######################################################################

## Load packages ---------------------------------------------------------------
library(here)
library(tidyverse) #filter,select,left_join, %>% 
library(janitor) 
library(terra) # for raster spatial data: crop,extract,cellsize
library(sf) # for vector spatial read_sf, st_transform

## Load data -------------------------------------------------------------------

#load Gabriel's lobster data
load(here("data/raw/2026_Sabira_JC_VILLASENOR_2/SDM/FINAL_EMSDM_EMMEAN_SP_382891.Rdata"))
#FINALEMMEAN is the name for SDM table for the species
#the file is FINAL_EMDSM_etc and the object inside is FINALEMMEAN
#rename FINALEMMEAN to lobster
lobster<-FINALEMMEAN
#now the next specie I load will replace FINALEMMEAN

#load Ere's turf data
#downloading directly from github using URL
#create object "turfs" which contains all TURFs in Mex
turfs<-read_sf("https://github.com/mex-fisheries/mex-TURFs/raw/refs/heads/master/data/output/mex_turfs_combined.gpkg")%>%
  #take full TURF dataset and keep only rows where sci_name is panulirus
  #now "turfs" contains only lobster 
  filter(scientific_name=="Panulirus argus") %>% 
  #same CRS as Gabriel's data
  st_transform(crs="EPSG:4326")

# PROCESSING ###################################################################

## building raster  -------------------------------------------------------------------

#BUILD CURRENT HSI RASTER FOR LOBSTER

#Take lobster SDM table, keep only currently suitable cells
#keep coords and current HSI, convert into raster map 
#and crop that map to area around the lobster turfs 

#take object called lobster, send it thru next steps
#and save final results as current_rast
#at the start, lobster is still a table, not a raster yet
current_rast<-lobster %>% #take stuff on the left and pass it into next function
  #keep rows where current HSI is > cutoff (both are Gabriel data from "lobster" table)
  filter(Current>cutoff) %>% 
  #keep only columns with x, y and current HSI 
  #rename Current as current top standardize name since we will also use mid, future. 
  select(x,y,current=Current) %>% 
  #tell rast() that the data is in x-y-value (x-y-current) format
  #and that x and y coords are long and lat
  rast(type="xyz",crs="EPSG:4326") %>% 
  #keep only the part of the raster that overlaps with the general bounding area of TURF polygons
  #extend allows the crop area to be slightly expanded if needed
  crop(y=turfs, extend=TRUE)
#now current_rast is a current suitale habitat raster for 
#lobster, cropped to relevant TURF area. 

#BUILD MID CENTURY(2050) HSI RASTER FOR LOBSTER UNDER RCP45 SCENARIO

mid_rast<-lobster%>%
  filter(RCP45_2050>cutoff) %>% 
  #rename RCP45_2050 as mid (although it should be mid45)
  select(x,y,mid=RCP45_2050) %>% 
  rast(type="xyz",crs="EPSG:4326") %>% 
  crop(y=turfs, extend=TRUE)

# VISUALIZE ####################################################################

## Another step ----------------------------------------------------------------

#plot current HSI raster
#showing where is lobster suitable right now after applying the cutoff
plot(current_rast)
#overlay the TURF polygons on existing raster plot 
plot(turfs[,1],add=TRUE, col="transparent")

#plot mid HSI raster
plot(mid_rast)
#overlay the TURF polygons
plot(turfs[,1],add=TRUE, col="transparent")

#plot difference/change from present to mid 
plot(mid_rast - current_rast)
plot(turfs[,1],add=TRUE, col="transparent")

# ANALYSIS #####################################################################

## Almost last step ------------------------------------------------------------

#for each TURF polygon, look at raster cells inside it (which have HSI values)
#and calculate mean current HSI

#x and y and arguments of the extract function 
#x = raster I'm extracting values from (current rast)
#y = spatial object used to decide where to extract (TURFs)
#fun = function, which is mean now
#na.rm = remove NAs
hsi_i_current<-extract(x=current_rast, y=turfs, fun=mean, na.rm=TRUE)

#same but for mid century raster now 
hsi_i_mid<-extract(x=mid_rast, y=turfs, fun=mean, na.rm=TRUE)

#combine mid and current mean HSI values in a single table
#left_join= join two tables together, keeping all rows from left table
#x = left table, which has mean current HSI per TURF
#y = right table, which has mean mid HSI per TURF
#match by using ID column, which extract() automatically adds 
#when loopng thru the polygons
#ID is the row number of each polygon in the "turfs" object
hsi_i<-left_join(x=hsi_i_current,y=hsi_i_mid, by="ID")

## Calculate area --------------------------------------------------------------

#cellSize() = calculate physical area of each raster area
#not every cell is the same area

#current_rast= each cell value = current HSI
#area_current = each cell value = area of that cell 
area_current<-cellSize(current_rast)
#only count area where current_rast is suitable habitat
#wherever current_rast has no suitable habitat value,
#erase the area value too.
area_current[is.na(current_rast)]<-NA

#for each turf polygon(y), find all suitable area cells inside(x)
#and sum them (fn=sum)
area_i_current<- extract(x=area_current, y=turfs, fun=sum, na.rm=TRUE)

area_mid<- cellSize(mid_rast)
area_mid[is.na(mid_rast)]<-NA
area_i_mid<- extract(x=area_mid, y=turfs, fun=sum, na.rm=TRUE)

# EXPORT #######################################################################

## The final step --------------------------------------------------------------

#NOTES
#think on making this more efficient 
#for all spp.
#make function, standardize 
#write pseudo code first 
