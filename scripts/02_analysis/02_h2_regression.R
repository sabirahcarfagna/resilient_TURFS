# H2: EFFECT OF CLIMATE SCENARIO ON HABITAT SUITABILITY ----------------------

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

# TURF COORDINATES FOR CONLEY STANDARD ERRORS ---------------------------------

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
    latitude = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry()

# PREPARE H2 REGRESSION DATA --------------------------------------------------

h2_regression <- turf_metrics %>%
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
    rcp = factor(
      sub("_.*", "", scenario),
      levels = c("rcp26", "rcp45", "rcp85")
    )
  ) %>%
  left_join(turf_coords, by = "sub_id") %>%
  ungroup()

# MEAN HSI REGRESSION ---------------------------------------------------------

hsi_regression <- h2_regression %>%
  filter(metric == "hsi")

hsi_h2_model <- feols(
  change ~ 0 + period + period:rcp,
  data = hsi_regression,
  vcov = conley(100) ~ latitude + longitude
)

summary(hsi_h2_model)
confint(hsi_h2_model)
# 
# PERCENT SUITABLE HABITAT REGRESSION ----------------------------------------

percent_regression <- h2_regression %>%
  filter(metric == "percent")

percent_h2_model <- feols(
  change ~ 0 + period + period:rcp,
  data = percent_regression,
  vcov = conley(100) ~ latitude + longitude
)

summary(percent_h2_model)
confint(percent_h2_model)

# SPECIES RICHNESS REGRESSION -------------------------------------------------

richness_regression <- richness %>%
  group_by(sub_id) %>%
  mutate(
    change = species_richness - species_richness[scenario == "current"]
  ) %>%
  filter(scenario != "current") %>%
  mutate(
    period = factor(sub(".*_", "", scenario)),
    rcp = factor(
      sub("_.*", "", scenario),
      levels = c("rcp26", "rcp45", "rcp85")
    )
  ) %>%
  left_join(turf_coords, by = "sub_id") %>%
  ungroup()

richness_h2_model <- feols(
  change ~ 0 + period + period:rcp,
  data = richness_regression,
  vcov = conley(100) ~ latitude + longitude
)

summary(richness_h2_model)
confint(richness_h2_model)

# H2 REGRESSION RESULTS TABLE -----------------------------------------------
#AI promted 

tidy_h2 <- function(model, metric) {
  ci <- confint(model)
  ct <- summary(model)$coeftable
  
  data.frame(
    Metric = metric,
    Coefficient = names(coef(model)),
    Estimate = unname(coef(model)),
    SE = ct[, "Std. Error"],
    CI_lower = ci[, 1],
    CI_upper = ci[, 2],
    p_value = ct[, "Pr(>|t|)"],
    N = nobs(model),
    row.names = NULL
  )
}

h2_results <- bind_rows(
  tidy_h2(hsi_h2_model, "Mean HSI"),
  tidy_h2(percent_h2_model, "Percent suitable"),
  tidy_h2(richness_h2_model, "Species richness")
)

View(h2_results)