# H1: HABITAT SUITABILITY THROUGH TIME ----------------------------------------

# PACKAGES --------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(fixest)
library(sf)

# LOAD DATA -------------------------------------------------------------------

turf_metrics <- readr::read_csv(
  "data/processed/turf_metrics/turf_species_metrics_exact.csv",
  show_col_types = FALSE
)

richness <- readr::read_csv(
  "data/processed/turf_metrics/turf_species_richness.csv",
  show_col_types = FALSE
)

# TURF COORDS FOR CONLEY STANDARD ERRORS --------------------------------

turfs <- st_read("data/turfs/mex_turfs_combined.gpkg")

turf_coords <- turfs %>%
  select(sub_id) %>%
  st_transform(6933) %>%
  group_by(sub_id) %>%
  summarise() %>%
  st_point_on_surface() %>%
  st_transform(4326) %>%
  mutate(
    longitude = st_coordinates(.)[, 1],
    latitude  = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry()

# PREPARE H1 REGRESSION DATA ------------------------------------------------

# Calculate changes relative to current conditions and reshape for regression

h1_regression <- turf_metrics %>%
  group_by(sub_id, aphia_id) %>%
  mutate(
    delta_hsi = mean_hsi - mean_hsi[scenario == "current"],
    delta_percent = percent_suitable - percent_suitable[scenario == "current"]
  ) %>%
  filter(scenario != "current") %>%
  pivot_longer(
    cols = c(delta_hsi, delta_percent),
    names_to = "metric",
    values_to = "change"
  ) %>%
  mutate(
    metric = sub("delta_", "", metric),
    period = factor(sub(".*_", "", scenario)),
    rcp = sub("_.*", "", scenario)
  ) %>%
  left_join(turf_coords, by = "sub_id") %>%
  ungroup()

# MEAN HSI REGRESSION -------------------------------------------------------

hsi_regression <- h1_regression %>%
  filter(metric == "hsi")

hsi_model <- feols(
  change ~ 0 + period,
  data = hsi_regression,
  vcov = conley(100) ~ latitude + longitude
)

summary(hsi_model)
confint(hsi_model)

# PERCENT SUITABLE HABITAT REGRESSION ---------------------------------------

percent_regression <- h1_regression %>%
  filter(metric == "percent")

percent_model <- feols(
  change ~ 0 + period,
  data = percent_regression,
  vcov = conley(100) ~ latitude + longitude
)

summary(percent_model)
confint(percent_model)

# SPECIES RICHNESS REGRESSION -----------------------------------------------

# Calculate changes relative to current richness and prepare regression data

richness_regression <- richness %>%
  group_by(sub_id) %>%
  mutate(
    change = species_richness - species_richness[scenario == "current"]
  ) %>%
  filter(scenario != "current") %>%
  mutate(
    period = factor(sub(".*_", "", scenario)),
    rcp = sub("_.*", "", scenario)
  ) %>%
  left_join(turf_coords, by = "sub_id") %>%
  ungroup()

# Regression with Conley spatial standard errors

richness_model <- feols(
  change ~ 0 + period,
  data = richness_regression,
  vcov = conley(100) ~ latitude + longitude
)

summary(richness_model)
confint(richness_model)

# H1 REGRESSION RESULTS TABLE -----------------------------------------------
# (here I promted AI "write a function that combines the results of all
# 3 regression models into 1 summary table, including estimates,
# SEs, CIs, p vals, and sample sizes)

tidy_h1 <- function(model, metric) {
  ci <- confint(model)
  ct <- summary(model)$coeftable
  
  data.frame(
    Metric = metric,
    Period = names(coef(model)),
    Estimate = unname(coef(model)),
    SE = ct[, "Std. Error"],
    CI_lower = ci[, 1],
    CI_upper = ci[, 2],
    p_value = ct[, "Pr(>|t|)"],
    N = nobs(model),
    row.names = NULL
  )
}

h1_results <- bind_rows(
  tidy_h1(hsi_model, "Mean HSI"),
  tidy_h1(percent_model, "Percent suitable"),
  tidy_h1(richness_model, "Species richness")
)

View(h1_results)