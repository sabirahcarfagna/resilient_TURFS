# DATA VALIDATION AND EXPLORATORY ANALYSIS ------------------------------------
#
# Goal:
# Inspect and explore the data used in the TURF analysis
# before beginning modeling and regressions.
#
# This script creates reproducible tables and figures for data validation
# and exploratory analysis.


# PACKAGES --------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(sf)
library(terra)


# FILE PATHS ------------------------------------------------------------------

# Cleaned TURF polygons and target species
turf_file <- "data/turfs/mex_turfs_combined.gpkg"

# Processed AquaX rasters
raster_folder <- "data/processed/sdm_rasters"

# Species-specific AquaX HSI cutoffs
cutoff_file <- "data/processed/sdm_cutoffs.csv"

# Previously calculated TURF x species x scenario metrics
metrics_file <- "data/processed/turf_metrics/turf_species_metrics.csv"

# Previously calculated TURF x scenario species richness
richness_file <- "data/processed/turf_metrics/turf_species_richness.csv"


# LOAD DATA -------------------------------------------------------------------

# Load cleaned TURF polygons and target species
turfs <- sf::st_read(turf_file)

# Load species-specific suitability cutoffs
species_cutoffs <- readr::read_csv(
  cutoff_file,
  show_col_types = FALSE
)

# Load TURF x species x scenario metrics
turf_metrics <- readr::read_csv(
  metrics_file,
  show_col_types = FALSE
)

# Load TURF x scenario species richness
species_richness <- readr::read_csv(
  richness_file,
  show_col_types = FALSE
)

# Create a list of all processed raster files
raster_files <- list.files(
  path = raster_folder,
  pattern = "\\.tif$",
  full.names = TRUE
)


# 1. DATA INVENTORY -----------------------------------------------------------

# Summarize the basic structure of the datasets used in the analysis
data_inventory <- tibble::tibble(
  dataset_component = c(
    "Records in cleaned TURF dataset",
    "Unique non-missing TURF IDs",
    "Unique sub-IDs",
    "Aphia IDs in cleaned TURF dataset",
    "Aphia IDs with processed AquaX rasters",
    "Unique scenarios",
    "Processed raster files"
  ),
  n = c(
    nrow(turfs),
    
    # Count unique non-missing TURF IDs
    dplyr::n_distinct(
      turfs$turf_id,
      na.rm = TRUE
    ),
    
    # Count unique sub-IDs
    dplyr::n_distinct(
      turfs$sub_id
    ),
    
    # Count valid Aphia IDs, excluding NA and "NA"
    dplyr::n_distinct(
      turfs$aphia_id[
        !is.na(turfs$aphia_id) &
          turfs$aphia_id != "NA"
      ]
    ),
    
    # Extract Aphia IDs from raster filenames and count unique IDs
    dplyr::n_distinct(
      sub(
        "_.*",
        "",
        basename(raster_files)
      )
    ),
    
    # Extract scenario names from raster filenames and count unique scenarios
    dplyr::n_distinct(
      sub(
        "^[0-9]+_",
        "",
        tools::file_path_sans_ext(
          basename(raster_files)
        )
      )
    ),
    
    # Count total processed raster files
    length(raster_files)
  )
)

data_inventory


# Check how missing species identifiers are stored
missing_species_ids <- turfs |>
  st_drop_geometry() |>
  filter(
    is.na(aphia_id) |
      aphia_id == "NA"
  ) |>
  count(
    aphia_id,
    scientific_name
  )

missing_species_ids


# 2. ANALYTICAL SAMPLE --------------------------------------------------------

# Get the Aphia IDs represented by the processed AquaX rasters
raster_aphia_ids <- raster_files |>
  basename() |>
  sub(
    "_.*",
    "",
    x = _
  ) |>
  unique()


# Classify each sub-ID based on species identification and AquaX availability
sub_id_sample <- turfs |>
  st_drop_geometry() |>
  group_by(sub_id) |>
  summarise(
    has_species = any(
      !is.na(aphia_id) &
        aphia_id != "NA"
    ),
    
    has_aquax = any(
      aphia_id %in% raster_aphia_ids,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) |>
  mutate(
    sample_group = case_when(
      !has_species ~
        "No identified species",
      
      has_aquax ~
        "At least one AquaX species",
      
      TRUE ~
        "Identified species, no AquaX species"
    )
  )


# Count sub-IDs in each group
sub_id_sample_summary <- sub_id_sample |>
  count(
    sample_group,
    name = "n_sub_ids"
  )

sub_id_sample_summary


# Create the final set of modeled sub-ID x species combinations
modeled_combinations <- turfs |>
  st_drop_geometry() |>
  filter(
    !is.na(aphia_id),
    aphia_id != "NA",
    aphia_id %in% raster_aphia_ids
  ) |>
  distinct(
    sub_id,
    turf_id,
    state,
    scientific_name,
    aphia_id
  )


# Count modeled species associated with each sub-ID
sub_id_species_summary <- modeled_combinations |>
  group_by(sub_id) |>
  summarise(
    turf_id = first(turf_id),
    state = first(state),
    n_species = n_distinct(aphia_id),
    
    species = paste(
      sort(
        unique(scientific_name)
      ),
      collapse = "; "
    ),
    
    .groups = "drop"
  ) |>
  arrange(
    desc(n_species)
  )

sub_id_species_summary


# Summarize geographic representation of the modeled sample
aquax_state_summary <- modeled_combinations |>
  group_by(state) |>
  summarise(
    n_sub_ids = n_distinct(sub_id),
    n_turf_ids = n_distinct(
      turf_id,
      na.rm = TRUE
    ),
    n_species = n_distinct(aphia_id),
    .groups = "drop"
  ) |>
  arrange(
    desc(n_sub_ids)
  )

aquax_state_summary

# Map the geographic distribution of the modeled sub-IDs ----------------------

# Keep one spatial geometry for each of the 36 modeled sub-IDs
# and transform to longitude/latitude for mapping
modeled_turfs_map <- turfs |>
  filter(
    sub_id %in% modeled_combinations$sub_id
  ) |>
  distinct(
    sub_id,
    .keep_all = TRUE
  ) |>
  st_transform(
    crs = 4326
  )


# Create a faint outline of Mexico for geographic context
mexico_map <- ggplot2::map_data(
  "world"
) |>
  filter(
    region == "Mexico"
  )


# Plot the 36 modeled sub-IDs
modeled_turfs_plot <- ggplot() +
  
  geom_polygon(
    data = mexico_map,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    fill = "grey98",
    color = "grey90",
    linewidth = 0.25
  ) +
  
  geom_sf(
    data = modeled_turfs_map,
    fill = "#d95f5f",
    color = "black",
    linewidth = 0.15,
    alpha = 0.95,
    inherit.aes = FALSE
  ) +
  
  coord_sf(
    xlim = c(-118, -86),
    ylim = c(14, 33),
    expand = FALSE,
    crs = st_crs(4326)
  ) +
  
  theme_void()


modeled_turfs_plot

# Check that all modeled combinations are represented in every scenario
scenario_summary <- turf_metrics |>
  count(
    scenario,
    name = "n_observations"
  )

scenario_summary


# 3. RASTER CELL SUPPORT ------------------------------------------------------

# Create one row for every sub-ID x species x scenario combination
cell_audit_combinations <- turf_metrics |>
  distinct(
    sub_id,
    aphia_id,
    scenario
  )


# Extract the raster cells supporting each calculated result and record
# how many contain HSI values, are NA, are suitable, or are unsuitable.
# Mean HSI is also recalculated directly from suitable cells so that
# the saved metric can later be checked independently.
get_cell_support <- function(
    sub_id_value,
    aphia_id_value,
    scenario_value
) {
  
  raster_file <- file.path(
    raster_folder,
    paste0(
      aphia_id_value,
      "_",
      scenario_value,
      ".tif"
    )
  )
  
  raster <- terra::rast(
    raster_file
  )
  
  turf <- turfs |>
    filter(
      sub_id == sub_id_value
    ) |>
    slice(1) |>
    st_transform(
      terra::crs(raster)
    )
  
  extracted <- terra::extract(
    raster,
    terra::vect(turf),
    cells = TRUE
  )
  
  hsi <- extracted[[2]]
  
  cutoff_value <- species_cutoffs |>
    filter(
      as.character(aphia_id) ==
        as.character(aphia_id_value)
    ) |>
    pull(cutoff) |>
    first()
  
  suitable_hsi <- hsi[
    !is.na(hsi) &
      hsi > cutoff_value
  ]
  
  tibble::tibble(
    sub_id = sub_id_value,
    aphia_id = as.character(
      aphia_id_value
    ),
    scenario = scenario_value,
    cutoff = cutoff_value,
    
    n_extracted_cells = length(hsi),
    
    n_valid_cells = sum(
      !is.na(hsi)
    ),
    
    n_na_cells = sum(
      is.na(hsi)
    ),
    
    n_zero_hsi_cells = sum(
      hsi == 0,
      na.rm = TRUE
    ),
    
    n_suitable_cells = sum(
      hsi > cutoff_value,
      na.rm = TRUE
    ),
    
    n_unsuitable_cells = sum(
      hsi <= cutoff_value,
      na.rm = TRUE
    ),
    
    mean_hsi_recalculated = if (
      length(suitable_hsi) > 0
    ) {
      mean(suitable_hsi)
    } else {
      NA_real_
    }
  )
}


# Run the raster-cell audit for all modeled combinations
cell_support <- dplyr::bind_rows(
  lapply(
    seq_len(
      nrow(cell_audit_combinations)
    ),
    function(i) {
      
      get_cell_support(
        sub_id_value =
          cell_audit_combinations$sub_id[i],
        
        aphia_id_value =
          cell_audit_combinations$aphia_id[i],
        
        scenario_value =
          cell_audit_combinations$scenario[i]
      )
    }
  )
)


# Check how many HSI-valued raster cells support each calculated result
hsi_cell_summary <- cell_support |>
  mutate(
    hsi_cell_group = case_when(
      n_valid_cells == 0 ~ "0 cells",
      n_valid_cells == 1 ~ "1 cell",
      n_valid_cells <= 3 ~ "2-3 cells",
      n_valid_cells <= 10 ~ "4-10 cells",
      TRUE ~ ">10 cells"
    )
  ) |>
  count(
    hsi_cell_group
  )

hsi_cell_summary

# Plot HSI-valued AquaX cell support -------------------------------------------

# Order cell-support categories from lowest to highest support
hsi_cell_plot_data <- hsi_cell_summary |>
  mutate(
    hsi_cell_group = factor(
      hsi_cell_group,
      levels = c(
        "0 cells",
        "1 cell",
        "2-3 cells",
        "4-10 cells",
        ">10 cells"
      )
    )
  )


# Create categorical cell-support plot
cell_support_plot <- ggplot(
  hsi_cell_plot_data,
  aes(
    x = hsi_cell_group,
    y = n
  )
) +
  
  geom_col(
    fill = "#1769aa",
    width = 0.7
  ) +
  
  geom_text(
    aes(
      label = n
    ),
    vjust = -0.4,
    size = 4
  ) +
  
  labs(
    x = "HSI-valued AquaX cells",
    y = "Number of observations"
  ) +
  
  expand_limits(
    y = max(hsi_cell_plot_data$n) * 1.10
  ) +
  
  theme_minimal() +
  
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(
      size = 12
    ),
    axis.text = element_text(
      size = 11
    )
  )


cell_support_plot

# Summarize AquaX NA coverage by scenario -------------------------------------

na_by_scenario <- cell_support |>
  mutate(
    na_status = case_when(
      n_valid_cells == 0 ~ "All cells NA",
      n_na_cells > 0 ~ "Some cells NA",
      TRUE ~ "No cells NA"
    )
  ) |>
  group_by(scenario) |>
  summarise(
    n_observations = n(),
    
    no_na = sum(na_status == "No cells NA"),
    some_na = sum(na_status == "Some cells NA"),
    all_na = sum(na_status == "All cells NA"),
    
    total_extracted_cells = sum(n_extracted_cells),
    total_na_cells = sum(n_na_cells),
    
    percent_cells_na = 100 *
      total_na_cells / total_extracted_cells,
    
    .groups = "drop"
  )

na_by_scenario

# Identify observations that change from HSI values to all NA -----------------

na_transition_check <- cell_support |>
  select(
    sub_id,
    aphia_id,
    scenario,
    n_extracted_cells,
    n_valid_cells,
    n_na_cells
  ) |>
  filter(
    n_valid_cells == 0
  ) |>
  arrange(
    sub_id,
    aphia_id,
    scenario
  )

na_transition_check

# Classify all-NA patterns across scenarios -----------------------------------

all_na_patterns <- cell_support |>
  group_by(
    sub_id,
    aphia_id
  ) |>
  summarise(
    scenarios_all_na = sum(n_valid_cells == 0),
    scenarios_with_values = sum(n_valid_cells > 0),
    .groups = "drop"
  ) |>
  filter(
    scenarios_all_na > 0
  ) |>
  mutate(
    na_pattern = case_when(
      scenarios_all_na == 7 ~ "All NA in every scenario",
      TRUE ~ "Changes between values and all NA"
    )
  )

print(all_na_patterns, n = Inf, width = Inf)

# Compare partial NA coverage across scenarios --------------------------------

partial_na_by_scenario <- cell_support |>
  filter(
    n_valid_cells > 0,
    n_na_cells > 0
  ) |>
  group_by(
    scenario
  ) |>
  summarise(
    n_observations = n(),
    total_extracted_cells = sum(n_extracted_cells),
    total_valid_cells = sum(n_valid_cells),
    total_na_cells = sum(n_na_cells),
    percent_cells_na = 100 * total_na_cells / total_extracted_cells,
    .groups = "drop"
  )

print(partial_na_by_scenario, n = Inf, width = Inf)
          
# Explain why mean HSI is NA by separating cases with no suitable
# cells from cases where all intersecting raster cells are NA
mean_hsi_na_summary <- turf_metrics |>
  mutate(
    aphia_id =
      as.character(aphia_id)
  ) |>
  left_join(
    cell_support |>
      select(
        sub_id,
        aphia_id,
        scenario,
        n_valid_cells,
        n_suitable_cells
      ),
    by = c(
      "sub_id",
      "aphia_id",
      "scenario"
    )
  ) |>
  filter(
    is.na(mean_hsi)
  ) |>
  mutate(
    mean_hsi_na_type = case_when(
      n_valid_cells == 0 ~
        "All extracted cells are NA",
      
      n_suitable_cells == 0 ~
        "HSI values present, no cells above cutoff",
      
      TRUE ~
        "Other"
    )
  ) |>
  count(
    mean_hsi_na_type
  )

mean_hsi_na_summary


# Calculate the area of each modeled sub-ID to examine whether smaller
# spatial units intersect fewer AquaX raster cells
modeled_subid_areas <- turfs |>
  filter(
    aphia_id %in% raster_aphia_ids
  ) |>
  distinct(
    sub_id,
    .keep_all = TRUE
  ) |>
  mutate(
    area_km2 = as.numeric(
      st_area(geom)
    ) / 1e6
  )


# Add sub-ID area to the raster cell counts
cell_support_with_area <- cell_support |>
  left_join(
    modeled_subid_areas |>
      st_drop_geometry() |>
      select(
        sub_id,
        area_km2
      ),
    by = "sub_id"
  )


# Measure the relationship between sub-ID area and the number of
# raster cells that intersect it
area_cell_correlation <- cor(
  cell_support_with_area$area_km2,
  cell_support_with_area$n_extracted_cells,
  use = "complete.obs"
)

area_cell_correlation


# Identify sub-ID x species combinations that have only 1-3 HSI-valued
# cells in at least one scenario
low_cell_combinations <- cell_support |>
  filter(
    n_valid_cells >= 1,
    n_valid_cells <= 3
  ) |>
  distinct(
    sub_id,
    aphia_id
  )


# Count which sub-IDs contain low-cell combinations
low_cell_subids <- low_cell_combinations |>
  count(
    sub_id,
    name = "n_low_cell_species"
  ) |>
  arrange(
    desc(n_low_cell_species)
  )

low_cell_combinations
low_cell_subids


# 4. METRIC VALIDATION AND PATTERNS -------------------------------------------


# 4A. PERCENT SUITABLE HABITAT ------------------------------------------------

# Independently reconstruct percent suitable habitat from raster-cell
# counts and compare it with the saved metric
percent_suitable_validation <- turf_metrics |>
  mutate(
    aphia_id =
      as.character(aphia_id)
  ) |>
  left_join(
    cell_support |>
      select(
        sub_id,
        aphia_id,
        scenario,
        n_valid_cells,
        n_suitable_cells
      ),
    by = c(
      "sub_id",
      "aphia_id",
      "scenario"
    )
  ) |>
  mutate(
    percent_suitable_recalculated = if_else(
      n_valid_cells > 0,
      100 * n_suitable_cells / n_valid_cells,
      NA_real_
    ),
    
    difference =
      percent_suitable -
      percent_suitable_recalculated
  )

summary(
  percent_suitable_validation$difference
)


# Summarize percent suitable habitat across scenarios
percent_suitable_by_scenario <- turf_metrics |>
  group_by(scenario) |>
  summarise(
    n_with_value = sum(
      !is.na(percent_suitable)
    ),
    
    mean_percent = mean(
      percent_suitable,
      na.rm = TRUE
    ),
    
    median_percent = median(
      percent_suitable,
      na.rm = TRUE
    ),
    
    n_zero = sum(
      percent_suitable == 0,
      na.rm = TRUE
    ),
    
    n_100 = sum(
      percent_suitable == 100,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

percent_suitable_by_scenario


# Show the percentage of results that are 100% suitable in each scenario
percent_100_by_scenario <- cell_support |>
  group_by(scenario) |>
  summarise(
    total_results = n(),
    
    n_100_percent = sum(
      n_valid_cells > 0 &
        n_suitable_cells == n_valid_cells
    ),
    
    percent_100 =
      100 * n_100_percent / total_results,
    
    .groups = "drop"
  )

percent_100_by_scenario

percent_suitable_plot_data <- turf_metrics |>
  mutate(
    suitable_group = case_when(
      is.na(percent_suitable) ~ "NA",
      percent_suitable == 0 ~ "0%",
      percent_suitable == 100 ~ "100%",
      TRUE ~ "Between 0 and 100%"
    )
  ) |>
  count(suitable_group) |>
  mutate(
    suitable_group = factor(
      suitable_group,
      levels = c(
        "0%",
        "Between 0 and 100%",
        "100%",
        "NA"
      )
    )
  )

percent_suitable_plot <- ggplot(
  percent_suitable_plot_data,
  aes(
    x = suitable_group,
    y = n
  )
) +
  geom_col(
    fill = "#1769aa",
    width = 0.7
  ) +
  geom_text(
    aes(label = n),
    vjust = -0.4,
    size = 4
  ) +
  labs(
    x = "% suitable habitat",
    y = "Number of observations"
  ) +
  expand_limits(
    y = max(percent_suitable_plot_data$n) * 1.10
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

percent_suitable_plot

percent_100_plot <- ggplot(
  percent_100_by_scenario,
  aes(
    x = factor(
      scenario,
      levels = c(
        "current",
        "rcp26_2050",
        "rcp26_2100",
        "rcp45_2050",
        "rcp45_2100",
        "rcp85_2050",
        "rcp85_2100"
      )
    ),
    y = percent_100,
    group = 1
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 3
  ) +
  geom_text(
    aes(
      label = paste0(
        round(percent_100, 1),
        "%"
      )
    ),
    vjust = -0.8,
    size = 3.5
  ) +
  labs(
    x = "Scenario",
    y = "Observations with 100% suitable habitat (%)"
  ) +
  scale_x_discrete(
    labels = c(
      "Current",
      "RCP2.6\n2050",
      "RCP2.6\n2100",
      "RCP4.5\n2050",
      "RCP4.5\n2100",
      "RCP8.5\n2050",
      "RCP8.5\n2100"
    )
  ) +
  expand_limits(
    y = c(0, 65)
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

percent_100_plot

# -------------------------------------------------------
# % suitable habitat spatial attention check
# Where is % suitable based on very few HSI-valued cells?
# -------------------------------------------------------

percent_suitable_attention_data <- cell_support |>
  mutate(
    cell_support_group = case_when(
      n_valid_cells == 0 ~ "0 cells",
      n_valid_cells == 1 ~ "1 cell",
      n_valid_cells >= 2 & n_valid_cells <= 3 ~ "2-3 cells",
      TRUE ~ ">3 cells"
    )
  )


# Count observations in each cell-support category

percent_suitable_attention_summary <- percent_suitable_attention_data |>
  count(
    cell_support_group,
    name = "n_observations"
  )

percent_suitable_attention_summary


# See exactly which sub-ID × species combinations
# ever have 0-3 HSI-valued cells

percent_suitable_attention_cases <- percent_suitable_attention_data |>
  filter(
    n_valid_cells <= 3
  ) |>
  select(
    sub_id,
    aphia_id,
    scenario,
    n_valid_cells,
    cell_support_group
  ) |>
  arrange(
    sub_id,
    aphia_id,
    scenario
  )

percent_suitable_attention_cases

# -------------------------------------------------------
# Locate sub-IDs with 0-3 HSI-valued cells
# -------------------------------------------------------

percent_suitable_attention_locations <- percent_suitable_attention_data |>
  filter(
    n_valid_cells <= 3
  ) |>
  group_by(
    sub_id
  ) |>
  summarise(
    has_0_cells = any(n_valid_cells == 0),
    has_1_cell = any(n_valid_cells == 1),
    has_2_3_cells = any(
      n_valid_cells >= 2 &
        n_valid_cells <= 3
    ),
    .groups = "drop"
  ) |>
  left_join(
    turfs |>
      st_drop_geometry() |>
      distinct(
        sub_id,
        state
      ),
    by = "sub_id"
  ) |>
  arrange(
    state,
    sub_id
  )

percent_suitable_attention_locations

# -------------------------------------------------------
# % suitable habitat spatial diagnostic maps
# -------------------------------------------------------

# Assign one mapping category to each affected sub-ID

percent_suitable_map_categories <- percent_suitable_attention_locations |>
  mutate(
    map_category = case_when(
      sub_id == "12_2" ~ "2 cells → 0 cells",
      has_0_cells ~ "0 cells",
      has_1_cell ~ "1 cell",
      has_2_3_cells ~ "2-3 cells",
      TRUE ~ "More than 3 cells"
    )
  )


# Attach categories to all modeled TURF geometries

percent_suitable_attention_map <- turfs |>
  filter(
    sub_id %in% modeled_combinations$sub_id
  ) |>
  distinct(
    sub_id,
    .keep_all = TRUE
  ) |>
  mutate(
    sub_id = as.character(sub_id)
  ) |>
  left_join(
    percent_suitable_map_categories |>
      select(
        sub_id,
        map_category
      ),
    by = "sub_id"
  ) |>
  mutate(
    map_category = if_else(
      is.na(map_category),
      "Other modeled sub-IDs",
      map_category
    )
  ) |>
  st_transform(
    crs = 4326
  )


# -------------------------------------------------------
# Baja California + Baja California Sur
# -------------------------------------------------------

percent_suitable_attention_west <- percent_suitable_attention_map |>
  filter(
    state %in% c("BC", "BCS")
  )

west_problem_bbox <- percent_suitable_attention_west |>
  filter(
    map_category != "Other modeled sub-IDs"
  ) |>
  st_bbox()

percent_suitable_attention_west_plot <- ggplot() +
  geom_sf(
    data = percent_suitable_attention_west,
    aes(fill = map_category),
    color = "black",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = c(
      "0 cells" = "#d95f5f",
      "1 cell" = "#f0ad4e",
      "2-3 cells" = "#f6d55c",
      "2 cells → 0 cells" = "#9b59b6",
      "Other modeled sub-IDs" = "grey80"
    ),
    breaks = c(
      "0 cells",
      "1 cell",
      "2-3 cells",
      "2 cells → 0 cells",
      "Other modeled sub-IDs"
    ),
    drop = FALSE
  ) +
  coord_sf(
    xlim = c(
      west_problem_bbox["xmin"] - 1,
      west_problem_bbox["xmax"] + 1
    ),
    ylim = c(
      west_problem_bbox["ymin"] - 1,
      west_problem_bbox["ymax"] + 1
    ),
    expand = FALSE
  ) +
  labs(
    title = "Baja California + Baja California Sur",
    fill = NULL
  ) +
  theme_void() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "none"
  )


# -------------------------------------------------------
# Veracruz
# -------------------------------------------------------

percent_suitable_attention_ver <- percent_suitable_attention_map |>
  filter(
    state == "VER"
  )

ver_problem_bbox <- percent_suitable_attention_ver |>
  filter(
    map_category != "Other modeled sub-IDs"
  ) |>
  st_bbox()

percent_suitable_attention_ver_plot <- ggplot() +
  geom_sf(
    data = percent_suitable_attention_ver,
    aes(fill = map_category),
    color = "black",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = c(
      "0 cells" = "#d95f5f",
      "1 cell" = "#f0ad4e",
      "2-3 cells" = "#f6d55c",
      "2 cells → 0 cells" = "#9b59b6",
      "Other modeled sub-IDs" = "grey80"
    ),
    breaks = c(
      "0 cells",
      "1 cell",
      "2-3 cells",
      "2 cells → 0 cells",
      "Other modeled sub-IDs"
    ),
    drop = FALSE
  ) +
  coord_sf(
    xlim = c(
      ver_problem_bbox["xmin"] - 1,
      ver_problem_bbox["xmax"] + 1
    ),
    ylim = c(
      ver_problem_bbox["ymin"] - 1,
      ver_problem_bbox["ymax"] + 1
    ),
    expand = FALSE
  ) +
  labs(
    title = "Veracruz",
    fill = NULL
  ) +
  theme_void() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "none"
  )

# Baja California diagnostic map
mean_hsi_attention_bc_plot <- ggplot() +
  geom_polygon(
    data = mexico_map,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    fill = "grey98",
    color = "grey85",
    linewidth = 0.25
  ) +
  geom_sf(
    data = mean_hsi_attention_map |>
      filter(attention_type == "No all-NA issue"),
    fill = "grey80",
    color = "black",
    linewidth = 0.2,
    inherit.aes = FALSE
  ) +
  geom_sf(
    data = mean_hsi_attention_bc,
    aes(fill = attention_type),
    color = "black",
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "All NA in every scenario" = "#d95f5f",
      "Becomes all NA in some scenarios" = "#f0ad4e"
    ),
    name = NULL
  ) +
  coord_sf(
    xlim = bc_xlim,
    ylim = bc_ylim,
    expand = FALSE,
    crs = st_crs(4326)
  ) +
  labs(title = "Baja California") +
  theme_void() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )


# Veracruz diagnostic map
mean_hsi_attention_ver_plot <- ggplot() +
  geom_polygon(
    data = mexico_map,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    fill = "grey98",
    color = "grey85",
    linewidth = 0.25
  ) +
  geom_sf(
    data = mean_hsi_attention_map |>
      filter(attention_type == "No all-NA issue"),
    fill = "grey80",
    color = "black",
    linewidth = 0.2,
    inherit.aes = FALSE
  ) +
  geom_sf(
    data = mean_hsi_attention_ver,
    aes(fill = attention_type),
    color = "black",
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "All NA in every scenario" = "#d95f5f",
      "Becomes all NA in some scenarios" = "#f0ad4e"
    ),
    name = NULL
  ) +
  coord_sf(
    xlim = ver_xlim,
    ylim = ver_ylim,
    expand = FALSE,
    crs = st_crs(4326)
  ) +
  labs(title = "Veracruz") +
  theme_void() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )


# View both diagnostic maps
mean_hsi_attention_bc_plot
mean_hsi_attention_ver_plot

# View both maps

percent_suitable_attention_west_plot

percent_suitable_attention_ver_plot

# ------------------------------------------------------------
# Climate-driven change in % suitable habitat
# ------------------------------------------------------------

percent_suitable_change <- turf_metrics |>
  select(
    sub_id,
    turf_id,
    aphia_id,
    scenario,
    percent_suitable
  ) |>
  tidyr::pivot_wider(
    names_from = scenario,
    values_from = percent_suitable
  ) |>
  tidyr::pivot_longer(
    cols = starts_with("rcp"),
    names_to = "future_scenario",
    values_to = "future_percent_suitable"
  ) |>
  mutate(
    change_percent_suitable =
      future_percent_suitable - current,
    
    change_type = case_when(
      is.na(current) | is.na(future_percent_suitable) ~ "Cannot compare",
      change_percent_suitable > 0 ~ "Gain",
      change_percent_suitable < 0 ~ "Loss",
      TRUE ~ "No change"
    )
  )

percent_suitable_change |>
  count(future_scenario, change_type)

percent_suitable_change |>
  filter(change_type == "Gain") |>
  arrange(future_scenario, desc(change_percent_suitable)) |>
  select(
    sub_id,
    aphia_id,
    future_scenario,
    current,
    future_percent_suitable,
    change_percent_suitable
  ) |>
  as.data.frame()

# ------------------------------------------------------------
# % suitable habitat redistribution under RCP8.5 2100
# ------------------------------------------------------------

percent_suitable_redistribution <- percent_suitable_change |>
  filter(
    future_scenario == "rcp85_2100",
    change_type %in% c("Gain", "Loss")
  ) |>
  mutate(
    combination = paste0(sub_id, " · ", aphia_id)
  ) |>
  arrange(change_percent_suitable) |>
  mutate(
    combination = factor(combination, levels = combination)
  )

percent_suitable_redistribution_plot <-
  ggplot(
    percent_suitable_redistribution,
    aes(
      x = combination,
      y = change_percent_suitable,
      fill = change_type
    )
  ) +
  geom_col(width = 0.8) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  scale_fill_manual(
    values = c(
      "Gain" = "#4C9F70",
      "Loss" = "#D95F5F"
    ),
    name = NULL
  ) +
  scale_y_continuous(
    limits = c(-100, 100),
    breaks = seq(-100, 100, 25)
  ) +
  labs(
    x = NULL,
    y = "Change in suitable habitat (percentage points)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

percent_suitable_redistribution_plot

# ------------------------------------------------------------
# Map direction of % suitable habitat change — RCP8.5 2100
# ------------------------------------------------------------

percent_suitable_direction_85_2100 <- percent_suitable_change |>
  filter(future_scenario == "rcp85_2100") |>
  group_by(sub_id) |>
  summarise(
    has_gain = any(change_type == "Gain"),
    has_loss = any(change_type == "Loss"),
    has_no_change = any(change_type == "No change"),
    all_cannot_compare = all(change_type == "Cannot compare"),
    .groups = "drop"
  ) |>
  mutate(
    change_direction = case_when(
      has_gain & has_loss ~ "Gain + loss",
      has_gain ~ "Gain",
      has_loss ~ "Loss",
      has_no_change ~ "No change",
      all_cannot_compare ~ "Cannot compare",
      TRUE ~ NA_character_
    )
  )


percent_suitable_direction_map_data <- modeled_subid_areas |>
  left_join(
    percent_suitable_direction_85_2100 |>
      select(sub_id, change_direction),
    by = "sub_id"
  )


percent_suitable_direction_map <-
  ggplot() +
  
  # Mexico
  geom_sf(
    data = mex_states,
    fill = "grey95",
    color = "white",
    linewidth = 0.25
  ) +
  
  # TURFs
  geom_sf(
    data = percent_suitable_direction_map_data,
    aes(fill = change_direction),
    color = "black",
    linewidth = 0.35
  ) +
  
  scale_fill_manual(
    values = c(
      "Gain" = "#3f8fc5",
      "Loss" = "#d95f5f",
      "Gain + loss" = "#9b59b6",
      "No change" = "#bdbdbd",
      "Cannot compare" = "white"
    ),
    breaks = c(
      "Gain",
      "Loss",
      "Gain + loss",
      "No change",
      "Cannot compare"
    ),
    name = NULL
  ) +
  
  coord_sf(
    xlim = c(-118, -86),
    ylim = c(14, 33),
    expand = FALSE
  ) +
  
  labs(
    title = "Direction of change in % suitable habitat",
    subtitle = "RCP8.5 · 2100"
  ) +
  
  theme_void() +
  
  theme(
    plot.title = element_text(
      size = 15,
      face = "bold"
    ),
    plot.subtitle = element_text(
      size = 11
    ),
    legend.position = "bottom"
  )


percent_suitable_direction_map

# 4B. MEAN HSI ----------------------------------------------------------------

# Independently reconstruct mean HSI from suitable raster cells and
# compare it with the saved metric
mean_hsi_validation <- turf_metrics |>
  mutate(
    aphia_id =
      as.character(aphia_id)
  ) |>
  left_join(
    cell_support |>
      select(
        sub_id,
        aphia_id,
        scenario,
        mean_hsi_recalculated
      ),
    by = c(
      "sub_id",
      "aphia_id",
      "scenario"
    )
  ) |>
  mutate(
    difference =
      mean_hsi -
      mean_hsi_recalculated
  )

summary(
  mean_hsi_validation$difference
)


# Summarize calculated mean HSI across scenarios
mean_hsi_by_scenario <- turf_metrics |>
  group_by(scenario) |>
  summarise(
    n_with_mean_hsi = sum(
      !is.na(mean_hsi)
    ),
    
    mean_hsi_average = mean(
      mean_hsi,
      na.rm = TRUE
    ),
    
    median_hsi = median(
      mean_hsi,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

mean_hsi_by_scenario

# Mean HSI climate response

# Combine the calculated mean HSI values with the cell-support
# information needed to distinguish structural zeros from true NAs

mean_hsi_climate_data <- turf_metrics |>
  mutate(
    aphia_id = as.character(aphia_id)
  ) |>
  left_join(
    cell_support |>
      mutate(
        aphia_id = as.character(aphia_id)
      ) |>
      select(
        sub_id,
        aphia_id,
        scenario,
        n_valid_cells,
        n_suitable_cells
      ),
    by = c(
      "sub_id",
      "aphia_id",
      "scenario"
    )
  ) |>
  mutate(
    mean_hsi_validated = case_when(
      n_valid_cells == 0 ~ NA_real_,
      n_suitable_cells == 0 ~ 0,
      TRUE ~ mean_hsi
    )
  )


# Summarize validated mean HSI by climate scenario

mean_hsi_climate_summary <- mean_hsi_climate_data |>
  group_by(scenario) |>
  summarise(
    n_with_value = sum(!is.na(mean_hsi_validated)),
    mean_hsi = mean(
      mean_hsi_validated,
      na.rm = TRUE
    ),
    median_hsi = median(
      mean_hsi_validated,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

mean_hsi_climate_summary

# -------------------------------------------------------
# Plot mean HSI climate response
# -------------------------------------------------------

mean_hsi_climate_plot <- ggplot(
  mean_hsi_climate_summary,
  aes(
    x = factor(
      scenario,
      levels = c(
        "current",
        "rcp26_2050",
        "rcp26_2100",
        "rcp45_2050",
        "rcp45_2100",
        "rcp85_2050",
        "rcp85_2100"
      )
    ),
    y = mean_hsi,
    group = 1
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 3
  ) +
  geom_text(
    aes(
      label = round(mean_hsi, 0)
    ),
    vjust = -0.8,
    size = 3.5
  ) +
  labs(
    x = "Scenario",
    y = "Average mean HSI"
  ) +
  scale_x_discrete(
    labels = c(
      "Current",
      "RCP2.6\n2050",
      "RCP2.6\n2100",
      "RCP4.5\n2050",
      "RCP4.5\n2100",
      "RCP8.5\n2050",
      "RCP8.5\n2100"
    )
  ) +
  expand_limits(
    y = c(0, 1000)
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

mean_hsi_climate_plot

# -------------------------------------------------------
# Mean HSI spatial attention check
# Where do all-NA AquaX cases occur?
# -------------------------------------------------------

mean_hsi_attention_data <- cell_support |>
  mutate(
    aphia_id = as.character(aphia_id)
  ) |>
  group_by(
    sub_id,
    aphia_id
  ) |>
  summarise(
    n_scenarios = n(),
    n_all_na = sum(n_valid_cells == 0),
    .groups = "drop"
  ) |>
  mutate(
    attention_type = case_when(
      n_all_na == 7 ~ "All NA in every scenario",
      n_all_na > 0 ~ "Becomes all NA in some scenarios",
      TRUE ~ "No all-NA issue"
    )
  )


# Check the cases that need attention

mean_hsi_attention_data |>
  filter(
    attention_type != "No all-NA issue"
  ) |>
  arrange(
    attention_type,
    sub_id,
    aphia_id
  )

# -------------------------------------------------------
# Mean HSI spatial attention map
# -------------------------------------------------------

# Reduce the attention data to one status per sub-ID
# If any species becomes all NA, that sub-ID is flagged.
# Permanent all-NA takes priority.

mean_hsi_attention_subids <- mean_hsi_attention_data |>
  filter(
    attention_type != "No all-NA issue"
  ) |>
  group_by(
    sub_id
  ) |>
  summarise(
    attention_type = case_when(
      any(attention_type == "All NA in every scenario") ~
        "All NA in every scenario",
      
      any(attention_type == "Becomes all NA in some scenarios") ~
        "Becomes all NA in some scenarios"
    ),
    .groups = "drop"
  )


# Attach the diagnostic status to the TURF geometries

mean_hsi_attention_map <- turfs |>
  filter(
    sub_id %in% modeled_combinations$sub_id
  ) |>
  distinct(
    sub_id,
    .keep_all = TRUE
  ) |>
  mutate(
    sub_id = as.character(sub_id)
  ) |>
  left_join(
    mean_hsi_attention_subids |>
      mutate(
        sub_id = as.character(sub_id)
      ),
    by = "sub_id"
  ) |>
  mutate(
    attention_type = if_else(
      is.na(attention_type),
      "No all-NA issue",
      attention_type
    )
  ) |>
  st_transform(
    crs = 4326
  )


# -------------------------------------------------------
# Mean HSI spatial attention maps
# Two diagnostic zooms: BC and Veracruz
# -------------------------------------------------------

# Separate the two affected clusters

mean_hsi_attention_bc <- mean_hsi_attention_map |>
  filter(
    sub_id %in% c(
      "10_1",
      "12_1",
      "12_2"
    )
  )

mean_hsi_attention_ver <- mean_hsi_attention_map |>
  filter(
    sub_id %in% c(
      "193_1",
      "203_1",
      "205_1",
      "209_1"
    )
  )


# Get bounding boxes around each affected cluster

bc_bbox <- st_bbox(mean_hsi_attention_bc)

ver_bbox <- st_bbox(mean_hsi_attention_ver)


# Add space around each cluster

bc_xlim <- c(
  bc_bbox["xmin"] - 1,
  bc_bbox["xmax"] + 1
)

bc_ylim <- c(
  bc_bbox["ymin"] - 1,
  bc_bbox["ymax"] + 1
)

ver_xlim <- c(
  ver_bbox["xmin"] - 1,
  ver_bbox["xmax"] + 1
)

ver_ylim <- c(
  ver_bbox["ymin"] - 1,
  ver_bbox["ymax"] + 1
)


# -------------------------------------------------------
# Baja California panel
# -------------------------------------------------------

percent_suitable_attention_west_plot <- ggplot() +
  geom_polygon(
    data = mexico_map,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    fill = "grey98",
    color = "grey85",
    linewidth = 0.25
  ) +
  geom_sf(
    data = percent_suitable_attention_map |>
      filter(map_category == "Other modeled sub-IDs"),
    fill = "grey80",
    color = "black",
    linewidth = 0.2,
    inherit.aes = FALSE
  ) +
  geom_sf(
    data = percent_suitable_attention_west |>
      filter(map_category != "Other modeled sub-IDs"),
    aes(
      fill = map_category
    ),
    color = "black",
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "0 cells" = "#d95f5f",
      "1 cell" = "#f0ad4e",
      "2-3 cells" = "#f6d55c",
      "2 cells → 0 cells" = "#9b59b6"
    ),
    name = NULL
  ) +
  coord_sf(
    xlim = c(
      west_problem_bbox["xmin"] - 1,
      west_problem_bbox["xmax"] + 1
    ),
    ylim = c(
      west_problem_bbox["ymin"] - 1,
      west_problem_bbox["ymax"] + 1
    ),
    expand = FALSE,
    crs = st_crs(4326)
  ) +
  labs(
    title = "Baja California + Baja California Sur"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold"
    ),
    legend.position = "none"
  )

# -------------------------------------------------------
# Veracruz panel
# -------------------------------------------------------
percent_suitable_attention_ver_plot <- ggplot() +
  geom_polygon(
    data = mexico_map,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    fill = "grey98",
    color = "grey85",
    linewidth = 0.25
  ) +
  geom_sf(
    data = percent_suitable_attention_map |>
      filter(map_category == "Other modeled sub-IDs"),
    fill = "grey80",
    color = "black",
    linewidth = 0.2,
    inherit.aes = FALSE
  ) +
  geom_sf(
    data = percent_suitable_attention_ver |>
      filter(map_category != "Other modeled sub-IDs"),
    aes(
      fill = map_category
    ),
    color = "black",
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "0 cells" = "#d95f5f",
      "1 cell" = "#f0ad4e",
      "2-3 cells" = "#f6d55c",
      "2 cells → 0 cells" = "#9b59b6"
    ),
    name = NULL
  ) +
  coord_sf(
    xlim = c(
      ver_problem_bbox["xmin"] - 1,
      ver_problem_bbox["xmax"] + 1
    ),
    ylim = c(
      ver_problem_bbox["ymin"] - 1,
      ver_problem_bbox["ymax"] + 1
    ),
    expand = FALSE,
    crs = st_crs(4326)
  ) +
  labs(
    title = "Veracruz"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold"
    ),
    legend.position = "none"
  )

# View both diagnostic maps

mean_hsi_attention_bc_plot
mean_hsi_attention_ver_plot

# MEAN HSI: ZERO SUITABLE-CELL CASES -----------------------------------------
#
# These observations have numeric AquaX HSI values within the TURF,
# but zero cells exceed the species-specific suitability cutoff.
#
# Mean HSI therefore has no suitable cells to average over.
# These are the cases that should be represented as mean HSI = 0.

mean_hsi_zero_cases <- cell_support |>
  dplyr::filter(
    n_valid_cells > 0,
    n_suitable_cells == 0
  )

# Spatial sub-IDs containing at least one of these observations
mean_hsi_zero_subids <- mean_hsi_zero_cases |>
  dplyr::distinct(sub_id)

# Check
mean_hsi_zero_cases |>
  dplyr::summarise(
    observations = dplyr::n(),
    unique_sub_ids = dplyr::n_distinct(sub_id),
    unique_species = dplyr::n_distinct(aphia_id)
  ) |>
  as.data.frame()

# MAP: MEAN HSI ZERO-SUITABLE-CELL CASES -------------------------------------
#
# Highlights the spatial sub-IDs containing observations where
# AquaX has numeric HSI values, but no cells exceed the
# species-specific suitability cutoff.

mean_hsi_zero_map_data <- modeled_turfs_map |>
  dplyr::mutate(
    zero_suitable_cells = dplyr::if_else(
      sub_id %in% mean_hsi_zero_subids$sub_id,
      "No cells above cutoff",
      "No issue"
    )
  )

mean_hsi_zero_map <- ggplot() +
  
  # Mexico
  geom_sf(
    data = mex_states,
    fill = "grey95",
    color = "white",
    linewidth = 0.25
  ) +
  
  # All 36 modeled sub-IDs
  geom_sf(
    data = mean_hsi_zero_map_data,
    aes(fill = zero_suitable_cells),
    color = "black",
    linewidth = 0.25
  ) +
  
  scale_fill_manual(
    values = c(
      "No cells above cutoff" = "#1769aa",
      "No issue" = "grey75"
    ),
    breaks = c(
      "No cells above cutoff",
      "No issue"
    ),
    name = NULL
  ) +
  
  coord_sf(
    xlim = c(-118.5, -86),
    ylim = c(14, 33),
    expand = FALSE
  ) +
  
  theme_void() +
  
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10)
  )

mean_hsi_zero_map

# 4C. SPECIES RICHNESS --------------------------------------------------------

# Independently reconstruct species richness from species presence
# and compare it with the saved richness metric
richness_check <- turf_metrics |>
  group_by(
    sub_id,
    scenario
  ) |>
  summarise(
    richness_recalculated = sum(
      present,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


richness_validation <- species_richness |>
  left_join(
    richness_check,
    by = c(
      "sub_id",
      "scenario"
    )
  ) |>
  mutate(
    difference =
      species_richness -
      richness_recalculated
  )

summary(
  richness_validation$difference
)


# Summarize species richness across scenarios
richness_by_scenario <- species_richness |>
  group_by(scenario) |>
  summarise(
    mean_richness = mean(
      species_richness
    ),
    
    median_richness = median(
      species_richness
    ),
    
    min_richness = min(
      species_richness
    ),
    
    max_richness = max(
      species_richness
    ),
    
    .groups = "drop"
  )

richness_by_scenario


# Store current richness for each sub-ID so every future scenario can
# be compared with the same current baseline
current_richness <- species_richness |>
  filter(
    scenario == "current"
  ) |>
  select(
    sub_id,
    current_richness = species_richness
  )


# Calculate richness change relative to current for every future scenario
richness_changes <- species_richness |>
  filter(
    scenario != "current"
  ) |>
  left_join(
    current_richness,
    by = "sub_id"
  ) |>
  mutate(
    richness_change =
      species_richness -
      current_richness,
    
    change_type = case_when(
      richness_change < 0 ~ "Loss",
      richness_change == 0 ~ "No change",
      richness_change > 0 ~ "Gain"
    )
  )


# Count how many sub-IDs lose, maintain, or gain richness in each
# future scenario
richness_change_summary <- richness_changes |>
  count(
    scenario,
    change_type,
    name = "n_sub_ids"
  )

richness_change_summary


# Identify every sub-ID with a richness change in at least one
# future scenario
richness_changes_nonzero <- richness_changes |>
  filter(
    richness_change != 0
  ) |>
  select(
    sub_id,
    scenario,
    current_richness,
    species_richness,
    richness_change
  ) |>
  arrange(
    sub_id,
    scenario
  )

richness_changes_nonzero


# Examine species-level presence changes between current and RCP8.5 2100
presence_change_85_2100 <- turf_metrics |>
  filter(
    scenario %in% c(
      "current",
      "rcp85_2100"
    )
  ) |>
  select(
    sub_id,
    aphia_id,
    scenario,
    present
  ) |>
  tidyr::pivot_wider(
    names_from = scenario,
    values_from = present
  ) |>
  mutate(
    presence_change =
      rcp85_2100 - current
  ) |>
  filter(
    presence_change != 0
  ) |>
  arrange(
    sub_id,
    presence_change
  )

presence_change_85_2100


# Check whether any sub-ID has species turnover that is hidden by
# zero net change in total richness under RCP8.5 2100
hidden_turnover_85_2100 <- presence_change_85_2100 |>
  count(
    sub_id,
    name = "n_species_changes"
  ) |>
  left_join(
    richness_changes |>
      filter(
        scenario == "rcp85_2100"
      ) |>
      select(
        sub_id,
        current_richness,
        species_richness,
        richness_change
      ),
    by = "sub_id"
  ) |>
  filter(
    richness_change == 0
  )

hidden_turnover_85_2100

# ------------------------------------------------------------
# Species richness climate response
# ------------------------------------------------------------

richness_climate_summary <- species_richness |>
  group_by(scenario) |>
  summarise(
    mean_richness = mean(species_richness, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    scenario = factor(
      scenario,
      levels = c(
        "current",
        "rcp26_2050",
        "rcp26_2100",
        "rcp45_2050",
        "rcp45_2100",
        "rcp85_2050",
        "rcp85_2100"
      ),
      labels = c(
        "Current",
        "RCP2.6\n2050",
        "RCP2.6\n2100",
        "RCP4.5\n2050",
        "RCP4.5\n2100",
        "RCP8.5\n2050",
        "RCP8.5\n2100"
      )
    )
  )


richness_climate_plot <-
  ggplot(
    richness_climate_summary,
    aes(
      x = scenario,
      y = mean_richness,
      group = 1
    )
  ) +
  geom_line(
    linewidth = 0.8
  ) +
  geom_point(
    size = 3
  ) +
  scale_y_continuous(
    limits = c(0, 3),
    breaks = seq(0, 3, 0.5)
  ) +
  labs(
    x = NULL,
    y = "Mean species richness"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      size = 10
    ),
    axis.title.y = element_text(
      size = 11
    )
  )


richness_climate_plot

# TARGET-SPECIES LOSSES: RCP8.5 2100 -----------------------------------------
#
# A target-species loss occurs when:
#   1. the species is currently targeted by the TURF,
#   2. it has >0% suitable habitat currently, and
#   3. it has 0% suitable habitat under RCP8.5 in 2100.
#
# This uses the original target-species metrics, so non-target species
# cannot be counted as losses.

target_species_losses_85_2100 <- turf_metrics |>
  dplyr::filter(
    scenario %in% c("current", "rcp85_2100")
  ) |>
  dplyr::select(
    sub_id,
    turf_id,
    aphia_id,
    scenario,
    percent_suitable
  ) |>
  tidyr::pivot_wider(
    names_from = scenario,
    values_from = percent_suitable
  ) |>
  dplyr::filter(
    !is.na(current),
    !is.na(rcp85_2100),
    current > 0,
    rcp85_2100 == 0
  )

target_species_loss_subids_85_2100 <- target_species_losses_85_2100 |>
  dplyr::count(
    sub_id,
    name = "n_target_species_lost"
  )

# MAP: TARGET-SPECIES LOSSES UNDER RCP8.5 2100 -------------------------------

target_species_loss_map_data <- modeled_turfs_map |>
  dplyr::left_join(
    target_species_loss_subids_85_2100,
    by = "sub_id"
  ) |>
  dplyr::mutate(
    loss_category = dplyr::case_when(
      n_target_species_lost == 3 ~ "3 target species lost",
      n_target_species_lost == 1 ~ "1 target species lost",
      TRUE ~ "No target species lost"
    )
  )

target_species_loss_map <- ggplot() +
  
  # Mexico
  geom_sf(
    data = mex_states,
    fill = "grey95",
    color = "white",
    linewidth = 0.25
  ) +
  
  # All 36 modeled sub-IDs
  geom_sf(
    data = target_species_loss_map_data,
    aes(fill = loss_category),
    color = "black",
    linewidth = 0.25
  ) +
  
  scale_fill_manual(
    values = c(
      "No target species lost" = "grey75",
      "1 target species lost" = "#f0ad4e",
      "3 target species lost" = "#d95f5f"
    ),
    breaks = c(
      "3 target species lost",
      "1 target species lost",
      "No target species lost"
    ),
    name = NULL
  ) +
  
  coord_sf(
    xlim = c(-118.5, -86),
    ylim = c(14, 33),
    expand = FALSE
  ) +
  
  theme_void() +
  
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10)
  )

target_species_loss_map

# MAP: SPECIES RICHNESS RESPONSE UNDER RCP8.5 2100 ---------------------------
#
# Separates interpretable richness losses from changes caused by
# AquaX missing predictions and apparent gains involving species
# already targeted by the TURF.

richness_response_map_data <- modeled_turfs_map |>
  dplyr::left_join(
    richness_85_2100_check |>
      dplyr::select(
        sub_id,
        current,
        rcp85_2100,
        richness_change
      ),
    by = "sub_id"
  ) |>
  dplyr::mutate(
    richness_response = dplyr::case_when(
      
      # These apparent losses occur because S. purpuratus
      # becomes all-NA in future AquaX scenarios
      sub_id %in% c("10_1", "12_2") ~
        "Future AquaX prediction becomes NA",
      
      # These apparent gains involve species already targeted
      # by the TURF that are 0% suitable currently
      sub_id %in% c("22_1", "25_1") ~
        "Existing target: 0% → suitable",
      
      # Interpretable richness losses
      richness_change < 0 ~
        "Target-species richness loss",
      
      # Everything else
      TRUE ~
        "No richness change"
    )
  )


richness_response_map <- ggplot() +
  
  geom_sf(
    data = mex_states,
    fill = "grey95",
    color = "white",
    linewidth = 0.25
  ) +
  
  geom_sf(
    data = richness_response_map_data,
    aes(fill = richness_response),
    color = "black",
    linewidth = 0.25
  ) +
  
  scale_fill_manual(
    values = c(
      "Target-species richness loss" = "#d95f5f",
      "Future AquaX prediction becomes NA" = "#9b59b6",
      "Existing target: 0% → suitable" = "#f0ad4e",
      "No richness change" = "grey75"
    ),
    breaks = c(
      "Target-species richness loss",
      "Future AquaX prediction becomes NA",
      "Existing target: 0% → suitable",
      "No richness change"
    ),
    name = NULL
  ) +
  
  coord_sf(
    xlim = c(-118.5, -86),
    ylim = c(14, 33),
    expand = FALSE
  ) +
  
  theme_void() +
  
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9)
  )

richness_response_map

# SAVE VALIDATION OBJECTS -----------------------------------------------------

saveRDS(
  list(
    turfs = turfs,
    species_cutoffs = species_cutoffs,
    turf_metrics = turf_metrics,
    species_richness = species_richness,
    raster_aphia_ids = raster_aphia_ids,
    modeled_combinations = modeled_combinations,
    sub_id_sample = sub_id_sample,
    sub_id_sample_summary = sub_id_sample_summary,
    sub_id_species_summary = sub_id_species_summary,
    aquax_state_summary = aquax_state_summary,
    modeled_turfs_map = modeled_turfs_map,
    modeled_turfs_plot = modeled_turfs_plot,
    cell_support = cell_support,
    cell_support_plot = cell_support_plot,
    cell_support_with_area = cell_support_with_area,
    na_by_scenario = na_by_scenario,
    all_na_patterns = all_na_patterns,
    partial_na_by_scenario = partial_na_by_scenario,
    modeled_subid_areas = modeled_subid_areas,
    low_cell_combinations = low_cell_combinations,
    low_cell_subids = low_cell_subids,
    mean_hsi_na_summary = mean_hsi_na_summary,
    percent_suitable_plot = percent_suitable_plot,
    percent_100_plot = percent_100_plot,
    percent_suitable_by_scenario = percent_suitable_by_scenario,
    mean_hsi_by_scenario = mean_hsi_by_scenario,
    mean_hsi_climate_data = mean_hsi_climate_data,
    mean_hsi_climate_summary = mean_hsi_climate_summary,
    mean_hsi_climate_plot = mean_hsi_climate_plot,
    mean_hsi_attention_bc_plot = mean_hsi_attention_bc_plot,  
    mean_hsi_attention_ver_plot = mean_hsi_attention_ver_plot, 
    richness_by_scenario = richness_by_scenario,
    mean_hsi_zero_map = mean_hsi_zero_map,
    mean_hsi_zero_cases = mean_hsi_zero_cases,
    mean_hsi_zero_subids = mean_hsi_zero_subids,
    richness_changes = richness_changes,
    percent_suitable_redistribution_plot = percent_suitable_redistribution_plot,
    percent_suitable_direction_map = percent_suitable_direction_map,
    percent_suitable_attention_west_plot = percent_suitable_attention_west_plot,
    percent_suitable_attention_ver_plot = percent_suitable_attention_ver_plot,
    richness_change_summary = richness_change_summary,
    richness_85_2100_check = richness_85_2100_check,
    richness_response_map_data = richness_response_map_data,
    richness_response_map = richness_response_map,
    richness_changes_nonzero = richness_changes_nonzero,
    presence_change_85_2100 = presence_change_85_2100
  ),
  "data/processed/turf_metrics/data_validation_objects.rds"
)