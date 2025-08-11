# ============================================================================
# Script: 03c_merge_prevalence_data.R
# Purpose: Create yearly transmission success dataset with temporal alignment
# Input: All processed data from scripts 02a-02d
# Output: yearly_transmission with transmission success variables
# Key approach:
#   - Aggregate Borrelia to site_year level
#   - Calculate transmission success (raw and log difference)
#   - Temporal alignment: predictors year N → transmission year N+1
#   - UPDATED: Added nymph_peak_day, larva_peak_day, and nlcdClass
# ============================================================================

cat("\n==== CREATING YEARLY TRANSMISSION DATASET ====\n")

# ---- Step 1: Load processed datasets ----
cat("Step 1: Loading processed datasets...\n")

borrelia_data <- readRDS("data/processed/borrelia_data.RDS")
phenology_data <- readRDS("data/processed/phenology_data.RDS")
host_data <- readRDS("data/processed/host_data.RDS")
climate_data <- readRDS("data/processed/climate_data.RDS")

cat("  Loaded datasets: borrelia, phenology, host, climate\n")

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

# ---- Step 2: Aggregate Borrelia to site_year level ----
cat("\nStep 2: Aggregating Borrelia to site_year level...\n")

site_year_borrelia <- borrelia_data %>%
  group_by(site_year, siteID, year) %>%
  summarise(
    total_tests = sum(total_tests),
    positive_tests = sum(positive_tests), 
    positivity_rate = positive_tests / total_tests,
    n_subsamples = n(),
    .groups = "drop"
  )

cat("  Aggregated to", nrow(site_year_borrelia), "site_years\n")

# ---- Step 3: Calculate transmission success variables ----
cat("\nStep 3: Calculating transmission success variables...\n")

transmission_borrelia <- site_year_borrelia %>%
  mutate(
    prev_year = year - 1,
    prev_site_year = paste(siteID, prev_year, sep = "_")
  ) %>%
  left_join(
    site_year_borrelia %>% 
      dplyr::select(site_year, positivity_rate) %>%
      rename(prev_positivity_rate = positivity_rate),
    by = c("prev_site_year" = "site_year")
  ) %>%
  mutate(
    # Raw difference (additive effects)
    transmission_raw = positivity_rate - prev_positivity_rate,
    # Ratio (multiplicative effects)
    transmission_ratio = ifelse(prev_positivity_rate == 0, NA, positivity_rate / prev_positivity_rate),
    # Log difference (multiplicative effects, log scale)
    transmission_log = ifelse(is.finite(log(transmission_ratio)), 
                              log(transmission_ratio), 
                              NA)
  ) %>%
  # Keep only rows where we have both years
  filter(!is.na(transmission_raw))

cat("  Transmission success calculated for", nrow(transmission_borrelia), "site_years\n")

# ---- Step 4: Merge with phenology data ----
cat("\nStep 4: Merging with phenology data...\n")
cat("  Alignment: Phenology year N → Transmission year N+1\n")

# UPDATED: Include nymph_peak_day and larva_peak_day
merge_1 <- transmission_borrelia %>%
  left_join(
    phenology_data %>% 
      dplyr::select(site_year, PAD, PAD_category, nymph_peak_day, larva_peak_day, nymph_max, larva_max),
    by = c("prev_site_year" = "site_year")
  ) %>%
  filter(!is.na(PAD))

cat("  Rows after merge:", nrow(merge_1), "\n")
cat("  Peak day variables included: nymph_peak_day, larva_peak_day\n")

# ---- Step 5: Merge with host data ----
cat("\nStep 5: Merging with host data...\n")

merge_2 <- merge_1 %>%
  left_join(
    host_data %>% 
      dplyr::select(site_year, mean_density_ha, mean_weighted_density_ha, 
                    species_richness, shannon_diversity, simpson_diversity),
    by = c("prev_site_year" = "site_year")
  ) %>%
  filter(!is.na(mean_density_ha))

cat("  Rows after merge:", nrow(merge_2), "\n")

# ---- Step 6: Merge with climate data ----
cat("\nStep 6: Merging with climate data...\n")

merge_3 <- merge_2 %>%
  left_join(
    climate_data %>% 
      dplyr::select(site_year, lat, long, elevation, mean_annual_temp, 
                    mean_annual_humidity, total_annual_precip, mean_annual_vpdmax, mean_annual_vpdmin),
    by = c("prev_site_year" = "site_year")
  ) %>%
  filter(!is.na(mean_annual_temp))

cat("  Rows after merge:", nrow(merge_3), "\n")

# ---- Step 7: Merge with habitat data ----
cat("\nStep 7: Merging with habitat data...\n")

yearly_transmission <- merge_3 %>%
  left_join(habitat_data, by = "siteID") %>%
  # Select and order columns
  dplyr::select(
    # Identifiers
    site_year, siteID, year, prev_site_year,
    
    # Response variables
    total_tests, positive_tests, positivity_rate, prev_positivity_rate,
    transmission_raw, transmission_log, transmission_ratio,
    
    # Phenology/Tick measures (UPDATED: added peak days)
    PAD, PAD_category, nymph_peak_day, larva_peak_day, nymph_max, larva_max,
    
    # Host variables  
    mean_density_ha, mean_weighted_density_ha, shannon_diversity, simpson_diversity,
    
    # Climate variables
    mean_annual_temp, mean_annual_humidity, total_annual_precip, mean_annual_vpdmax, mean_annual_vpdmin,
    
    # Geography and habitat
    lat, long, elevation, nlcdClass
  ) %>%
  arrange(siteID, year)

# Check final merge success
n_with_habitat <- sum(!is.na(yearly_transmission$nlcdClass))
cat("  Final rows:", nrow(yearly_transmission), "\n")
cat("  Site_years with habitat data:", n_with_habitat, "\n")

# Show habitat distribution
if(n_with_habitat > 0) {
  habitat_counts <- table(yearly_transmission$nlcdClass, useNA = "ifany")
  cat("  Habitat distribution:\n")
  for(i in 1:length(habitat_counts)) {
    cat("    ", names(habitat_counts)[i], ":", habitat_counts[i], "\n")
  }
}

# ---- Step 8: Summary and save ----
cat("\nStep 8: Summary and saving...\n")

# Quick summary of transmission variables
transmission_summary <- yearly_transmission %>%
  summarise(
    n_site_years = n(),
    mean_transmission_raw = round(mean(transmission_raw), 4),
    sd_transmission_raw = round(sd(transmission_raw), 4),
    mean_transmission_log = round(mean(transmission_log), 4),
    sd_transmission_log = round(sd(transmission_log), 4)
  )

cat("  Transmission success summary:\n")
cat("  - Site-years:", transmission_summary$n_site_years, "\n")
cat("  - Raw transmission (mean ± SD):", transmission_summary$mean_transmission_raw, "±", transmission_summary$sd_transmission_raw, "\n")
cat("  - Log transmission (mean ± SD):", transmission_summary$mean_transmission_log, "±", transmission_summary$sd_transmission_log, "\n")

# Peak day summaries
if(sum(!is.na(yearly_transmission$nymph_peak_day)) > 0) {
  nymph_summary <- yearly_transmission %>%
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

if(sum(!is.na(yearly_transmission$larva_peak_day)) > 0) {
  larva_summary <- yearly_transmission %>%
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

# Data completeness check
complete_cases_summary <- yearly_transmission %>%
  summarise(
    total_rows = n(),
    missing_transmission_log = sum(is.na(transmission_log)),
    missing_PAD = sum(is.na(PAD)),
    missing_nymph_peak = sum(is.na(nymph_peak_day)),
    missing_larva_peak = sum(is.na(larva_peak_day)),
    missing_host_density = sum(is.na(mean_density_ha)),
    missing_shannon = sum(is.na(shannon_diversity)),
    missing_temp = sum(is.na(mean_annual_temp)),
    missing_humidity = sum(is.na(mean_annual_humidity)),
    missing_precip = sum(is.na(total_annual_precip)),
    missing_habitat = sum(is.na(nlcdClass))
  )

cat("\n  Missing data summary:\n")
print(complete_cases_summary)

# Save dataset
saveRDS(yearly_transmission, "data/processed/yearly_transmission.RDS")
write.csv(yearly_transmission, "data/processed/yearly_transmission.csv", row.names = FALSE)

cat("\n  Saved yearly transmission dataset to: data/processed/yearly_transmission.RDS\n")
cat("  Dataset contains year-to-year transmission success with temporal alignment\n")

# ---- Step 9: Data preview ----
cat("\nStep 9: Data preview...\n")

cat("  Dataset structure:\n")
cat("  - Response variable: transmission_log (log ratio)\n")
cat("  - Peak day components: nymph_peak_day, larva_peak_day\n")
cat("  - Phenology predictor: PAD\n")
cat("  - Host predictors: mean_density_ha, shannon_diversity, simpson_diversity\n")
cat("  - Climate predictors: mean_annual_temp, mean_annual_humidity, total_annual_precip, mean_annual_vpdmax, mean_annual_vpdmin\n")
cat("  - Habitat variable: nlcdClass\n")
cat("  - Temporal alignment: Predictors year N → Transmission year N+1\n")

# Show first few rows
cat("\n  First 5 rows:\n")
preview <- yearly_transmission %>%
  dplyr::select(site_year, transmission_log, PAD, nymph_peak_day, larva_peak_day, 
         mean_density_ha, shannon_diversity, nlcdClass, siteID) %>%
  head(5)
print(preview)

cat("\n==== YEARLY TRANSMISSION DATASET COMPLETE ====\n")
cat("UPDATES: Added nymph_peak_day, larva_peak_day, and nlcdClass for additional analyses\n")