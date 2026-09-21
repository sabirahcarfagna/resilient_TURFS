# H1: HABITAT SUITABILITY THROUGH TIME ----------------------------------------

# PACKAGES --------------------------------------------------------------------

library(dplyr)
library(tidyr)

# LOAD DATA -------------------------------------------------------------------

#turf_metrics will read CSV file into R
turf_metrics <- readr::read_csv(
  "data/processed/turf_metrics/turf_species_metrics.csv",
  show_col_types = FALSE
)

# EXTRACT TIME PERIOD FROM SCENARIO ----------------------------------------------------

#H1 cares about projection period not emissions severity 
#adds new column called period 

# turf_metrics <- turf_metrics %>% 
#   mutate( #create "period" column
#     period = case_when( 
#       #if scenario equals "current" then put "current" in the new period column 
#       scenario == "current" ~ "current",
#       #for each row look in the scenario coln, if it contains 2050 then assign to period
#       grepl("2050", scenario) ~ "2050",
#       grepl("2100", scenario) ~ "2100"
#     )
#   )
# 
# # SET CURRENT AS REFERENCE PERIOD ---------------------------------------------
# 
# turf_metrics$period <- relevel(
#   factor(turf_metrics$period), #treat period as categorical variable
#   ref = "current" #set current as the reference category
#   
# )

# RESHAPE MEAN HSI FOR CHANGE CALCULATION -------------------------------------

# Current structure:
# sub_id   scenario       mean_hsi
# A        current        0.80
# A        rcp26_2050     0.75
# A        rcp26_2100     0.70
#
# Reshape it so the scenarios become separate cols:
# sub_id   current   rcp26_2050   rcp26_2100
# A        0.80      0.75          0.70

# Reshape mean HSI so each scenario has its own column
hsi_change <- turf_metrics |>
  select(sub_id, turf_id, aphia_id, scenario, mean_hsi) |>
  pivot_wider(
    names_from = scenario, #take values inside scenario col and turn into new cols names 
    values_from = mean_hsi #and fill them with HSI values
  )

# CALCULATE CHANGE IN MEAN HSI ------------------------------------------------

hsi_change <- hsi_change %>%
  mutate(
    delta_rcp26_2050 = rcp26_2050 - current,
    delta_rcp26_2100 = rcp26_2100 - current,
    
    delta_rcp45_2050 = rcp45_2050 - current,
    delta_rcp45_2100 = rcp45_2100 - current,
    
    delta_rcp85_2050 = rcp85_2050 - current,
    delta_rcp85_2100 = rcp85_2100 - current
  )

# RESHAPE HSI CHANGE FOR REGRESSION -------------------------------------------

# Current structure:
# sub_id   delta_rcp26_2050   delta_rcp26_2100
# A             -8.76               -6.40
#
# Reshape it so each delta becomes its own row:
# sub_id   projection          delta_hsi
# A        delta_rcp26_2050      -8.76
# A        delta_rcp26_2100      -6.40

hsi_regression <- hsi_change %>%
  # turn the six delta columns into rows
  pivot_longer(
    cols = starts_with("delta_"),#put old delta column names into a new column called "projection"
    names_to = "projection", #put delta HSI into a new column called "delta_hsi"
    values_to = "delta_hsi"
  )

# EXTRACT RCP AND PROJECTION PERIOD -------------------------------------------

#separate at each _ and put those pieces into three coloumns 

hsi_regression <- hsi_regression %>%
  separate(
    projection,
    into = c("delta", "rcp", "period"),
    sep = "_"
  )

# REMOVE UNNECESSARY DELTA LABEL COLUMN ---------------------------------------

hsi_regression <- hsi_regression %>%
  select(-delta)