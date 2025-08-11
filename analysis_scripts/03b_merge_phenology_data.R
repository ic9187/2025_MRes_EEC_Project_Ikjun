# ============================================================================
# Script: 03b_merge_phenology_data.R
# Purpose: Create PAD analysis dataset with environmental predictors
# Input: All processed data from scripts 02a-02d
# Output: pad_analysis_dataset with PAD as response variable
# Key approach:
#   - PAD as response variable (continuous)
#   - Same year temporal alignment (no lag)
#   - Host and climate variables as predictors
#   - UPDATED: Added nymph_peak_day, larva_peak_day, and nlcdClass
# ============================================================================

cat("\n==== CREATING PAD ANALYSIS DATASET ====\n")

# ---- Step 1: Load processed datasets ----
cat("Step 1: Loading processed datasets...\n")

phenology_data <- readRDS("data/processed/phenology_data.RDS")
host_data <- readRDS("data/processed/host_data.RDS")
climate_data <- readRDS("data/processed/climate_data.RDS")

cat("  Loaded datasets:\n")
cat("  - Phenology:", nrow(phenology_data), "site_years\n")
cat("  - Host data:", nrow(host_data), "site_years\n")
cat("  - Climate:", nrow(climate_data), "site_years\n")

# ---- Step 1.5: Extract habitat data ----
cat("\nStep 1.5: Extracting habitat data from tick field data...\n")

# Load tick field data for nlcdClass
if(!exists("tick_raw")) {
  stop("Error: tick_raw not found. Please run 00_setup.R first.")
}

# Extract nlcdClass by site (take first non-NA value per site)
habitat_data <- tick_raw %>%
  filter(!is.na(nlcdClass) & nlcdClass != "") %>%
  mutate(siteID = sub("_.*", "", plotID)) %>%
  group_by(siteID) %>%
  summarise(
    nlcdClass = first(nlcdClass),
    .groups = "drop"
  )

cat("  Extracted habitat data for", nrow(habitat_data), "sites\n")
cat("  Habitat types found:", length(unique(habitat_data$nlcdClass)), "\n")
cat("  Habitat types:", paste(unique(habitat_data$nlcdClass), collapse = ", "), "\n")

# ---- Step 2: Start with phenology data (response variable) ----
cat("\nStep 2: Starting with phenology data as base...\n")

# Use phenology data as the foundation since PAD is our response
# UPDATED: Include nymph_peak_day and larva_peak_day
pad_base <- phenology_data %>%
  dplyr::select(site_year, siteID, year, PAD, nymph_peak_day, larva_peak_day, 
         nymph_max, larva_max)

cat("  Base dataset:", nrow(pad_base), "site_years with PAD\n")
cat("  Peak day variables included: nymph_peak_day, larva_peak_day\n")

# ---- Step 3: Merge with host data ----
cat("\nStep 3: Merging with host data...\n")

merge_1 <- pad_base %>%
  left_join(
    host_data %>% 
      dplyr::select(site_year, mean_density_ha, mean_weighted_density_ha, 
                    species_richness, shannon_diversity, simpson_diversity),
    by = "site_year"
  )

# Check merge success
n_with_host <- sum(!is.na(merge_1$mean_density_ha))
cat("  Rows after merge:", nrow(merge_1), "\n")
cat("  Site_years with host data:", n_with_host, "\n")

# ---- Step 4: Merge with climate data ----
cat("\nStep 4: Merging with climate data...\n")

merge_2 <- merge_1 %>%
  left_join(
    climate_data %>% 
      dplyr::select(site_year, lat, long, elevation, mean_annual_temp, 
                    mean_annual_humidity, total_annual_precip, mean_annual_vpdmax, mean_annual_vpdmin),
    by = "site_year"
  )

# Check merge success
n_with_climate <- sum(!is.na(merge_2$mean_annual_temp))
cat("  Rows after merge:", nrow(merge_2), "\n")
cat("  Site_years with climate data:", n_with_climate, "\n")

# ---- Step 5: Merge with habitat data ----
cat("\nStep 5: Merging with habitat data...\n")

pad_analysis_dataset <- merge_2 %>%
  left_join(habitat_data, by = "siteID")

# Check final merge success
n_with_habitat <- sum(!is.na(pad_analysis_dataset$nlcdClass))
cat("  Rows after merge:", nrow(pad_analysis_dataset), "\n")
cat("  Site_years with habitat data:", n_with_habitat, "\n")

# Show habitat distribution
if(n_with_habitat > 0) {
  habitat_counts <- table(pad_analysis_dataset$nlcdClass, useNA = "ifany")
  cat("  Habitat distribution:\n")
  for(i in 1:length(habitat_counts)) {
    cat("    ", names(habitat_counts)[i], ":", habitat_counts[i], "\n")
  }
}

# ---- Step 6: Check data completeness ----
cat("\nStep 6: Checking data completeness...\n")

# Count complete cases for different variable sets
complete_host <- pad_analysis_dataset %>%
  filter(!is.na(PAD), !is.na(mean_density_ha), !is.na(shannon_diversity))

complete_climate <- pad_analysis_dataset %>%
  filter(!is.na(PAD), !is.na(mean_annual_temp), !is.na(mean_annual_humidity))

complete_all <- pad_analysis_dataset %>%
  filter(!is.na(PAD), !is.na(mean_density_ha), !is.na(shannon_diversity),
         !is.na(mean_annual_temp), !is.na(mean_annual_humidity), 
         !is.na(total_annual_precip), !is.na(lat))

complete_peak_days <- pad_analysis_dataset %>%
  filter(!is.na(nymph_peak_day), !is.na(larva_peak_day))

cat("  Complete cases:\n")
cat("  - PAD + Host variables:", nrow(complete_host), "\n")
cat("  - PAD + Climate variables:", nrow(complete_climate), "\n")
cat("  - PAD + All variables:", nrow(complete_all), "\n")
cat("  - Peak day variables:", nrow(complete_peak_days), "\n")
cat("  - Habitat data:", n_with_habitat, "\n")

# ---- Step 7: Summary statistics ----
cat("\nStep 7: Summary statistics...\n")

# PAD summary
pad_summary <- pad_analysis_dataset %>%
  filter(!is.na(PAD)) %>%
  summarise(
    n_observations = n(),
    mean_PAD = mean(PAD),
    sd_PAD = sd(PAD),
    min_PAD = min(PAD),
    max_PAD = max(PAD),
    median_PAD = median(PAD)
  )

cat("  PAD Summary:\n")
print(pad_summary)

# Peak day summaries
if(sum(!is.na(pad_analysis_dataset$nymph_peak_day)) > 0) {
  nymph_summary <- pad_analysis_dataset %>%
    filter(!is.na(nymph_peak_day)) %>%
    summarise(
      mean_nymph_peak = mean(nymph_peak_day),
      sd_nymph_peak = sd(nymph_peak_day),
      min_nymph_peak = min(nymph_peak_day),
      max_nymph_peak = max(nymph_peak_day)
    )
  
  cat("\n  Nymph Peak Day Summary:\n")
  print(nymph_summary)
}

if(sum(!is.na(pad_analysis_dataset$larva_peak_day)) > 0) {
  larva_summary <- pad_analysis_dataset %>%
    filter(!is.na(larva_peak_day)) %>%
    summarise(
      mean_larva_peak = mean(larva_peak_day),
      sd_larva_peak = sd(larva_peak_day),
      min_larva_peak = min(larva_peak_day),
      max_larva_peak = max(larva_peak_day)
    )
  
  cat("\n  Larva Peak Day Summary:\n")
  print(larva_summary)
}

# Variable availability
missing_summary <- pad_analysis_dataset %>%
  summarise(
    total_rows = n(),
    missing_PAD = sum(is.na(PAD)),
    missing_nymph_peak = sum(is.na(nymph_peak_day)),
    missing_larva_peak = sum(is.na(larva_peak_day)),
    missing_host_density = sum(is.na(mean_density_ha)),
    missing_shannon = sum(is.na(shannon_diversity)),
    missing_temp = sum(is.na(mean_annual_temp)),
    missing_humidity = sum(is.na(mean_annual_humidity)),
    missing_precip = sum(is.na(total_annual_precip)),
    missing_habitat = sum(is.na(nlcdClass)),
    missing_lat = sum(is.na(lat))
  )

cat("\n  Missing data summary:\n")
print(missing_summary)

# ---- Step 8: Save final dataset ----
cat("\nStep 8: Saving PAD analysis dataset...\n")

saveRDS(pad_analysis_dataset, "data/processed/pad_analysis_dataset.RDS")
write.csv(pad_analysis_dataset, "data/processed/pad_analysis_dataset.csv", row.names = FALSE)

cat("  Saved to: data/processed/pad_analysis_dataset.RDS\n")

# ---- Step 9: Data preview ----
cat("\nStep 9: Data preview...\n")

cat("  Dataset structure:\n")
cat("  - Response variable: PAD (continuous)\n")
cat("  - Peak day components: nymph_peak_day, larva_peak_day\n")
cat("  - Host predictors: mean_density_ha, mean_weighted_density_ha, species_richness, shannon_diversity, simpson_diversity\n")
cat("  - Climate predictors: mean_annual_temp, mean_annual_humidity, total_annual_precip, mean_annual_vpdmax, mean_annual_vpdmin\n")
cat("  - Spatial predictors: lat, long, elevation\n")
cat("  - Habitat variable: nlcdClass\n")
cat("  - Temporal alignment: Same year (no lag)\n")

# Show first few rows
cat("\n  First 5 rows:\n")
preview <- pad_analysis_dataset %>%
  dplyr::select(site_year, PAD, nymph_peak_day, larva_peak_day, mean_density_ha, shannon_diversity, 
         mean_annual_temp, mean_annual_humidity, nlcdClass, siteID) %>%
  head(5)
print(preview)

cat("\n==== PAD ANALYSIS DATASET COMPLETE ====\n")
cat("UPDATES: Added nymph_peak_day, larva_peak_day, and nlcdClass for additional analyses\n")