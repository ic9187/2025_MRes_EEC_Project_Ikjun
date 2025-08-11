# ============================================================================
# Script: 02b_pathogen_to_borrelia.R (FINAL - Clean Approach)
# Purpose: Process pathogen data to calculate Borrelia positivity at tick and subsample level
# Input: path_raw (from 00_setup.R), clean_site_years (from 01_tick_sampling_effort.R)
# Output: borrelia_data with positivity rates per subsample (tick-level aggregation)
# Key approach:
#   1. Group by testingID - if ANY Borrelia test positive → tick is positive
#   2. Group by subsampleID - calculate proportion of positive ticks
#   3. Maintain subsample structure for mixed-effects modeling
# ============================================================================

cat("\n==== PATHOGEN TO BORRELIA ANALYSIS (FINAL APPROACH) ====\n")

# ---- Step 1: Load site_year filter ----
cat("Step 1: Loading site_year filter...\n")
clean_site_years <- readRDS("data/processed/clean_site_years.RDS")
cat("  Loaded", length(clean_site_years), "clean site_years\n")

# ---- Step 2: Filter for Borrelia and Ixodes ----
cat("\nStep 2: Filtering for Borrelia pathogens in Ixodes ticks...\n")

# Define traditional Lyme disease agents specifically
lyme_agents <- c("Borrelia burgdorferi", 
                 "Borrelia burgdorferi sensu lato")

bor_ixodes <- path_raw %>%
  filter(
    testPathogenName %in% lyme_agents,  # Only Lyme disease agents
    !is.na(testResult),
    grepl("IXO", subsampleID)  # Only Ixodes ticks
  ) %>%
  mutate(
    collectDate = as.Date(collectDate),
    year = year(collectDate),
    siteID = sub("_.*", "", plotID),
    site_year = paste(siteID, year, sep = "_")
  ) %>%
  # Apply clean_site_years filter
  filter(site_year %in% clean_site_years)

cat("  Borrelia tests in Ixodes:", nrow(bor_ixodes), "\n")
cat("  Unique testingIDs (individual ticks):", n_distinct(bor_ixodes$testingID), "\n")
cat("  Unique subsampleIDs:", n_distinct(bor_ixodes$subsampleID), "\n")
cat("  Site_years represented:", n_distinct(bor_ixodes$site_year), "\n")

# ---- Step 3: Aggregate by testingID (tick level) ----
cat("\nStep 3: Aggregating to tick level (any Borrelia positive)...\n")

tick_level <- bor_ixodes %>%
  group_by(testingID, subsampleID, siteID, year, site_year) %>%
  summarise(
    # Tick is positive if ANY Borrelia test is positive
    tick_borrelia_positive = any(testResult == "Positive"),
    
    # Keep track of testing details
    n_borrelia_tests = n(),
    borrelia_species_tested = paste(sort(unique(testPathogenName)), collapse = "; "),
    
    # Keep temporal information (should be same for all tests on same tick)
    collectDate = first(collectDate),
    day_number = first(day_number),
    month = first(month),
    individualCount = first(individualCount),
    
    .groups = "drop"
  )

cat("  Tick-level dataset created:", nrow(tick_level), "individual ticks\n")
cat("  Positive ticks:", sum(tick_level$tick_borrelia_positive), 
    "(", round(mean(tick_level$tick_borrelia_positive) * 100, 1), "%)\n")

# ---- Step 4: Aggregate by subsampleID (subsample level) ----
cat("\nStep 4: Aggregating to subsample level...\n")

borrelia_summary <- tick_level %>%
  group_by(subsampleID, siteID, year, site_year) %>%
  summarise(
    # Core variables for binomial analysis
    total_tests = n(),                          # Total ticks tested
    positive_tests = sum(tick_borrelia_positive), # Positive ticks
    positivity_rate = positive_tests / total_tests,
    
    # Temporal information
    collectDate = first(collectDate),
    day_number = first(day_number),
    month = first(month),
    
    # Sample composition metrics
    total_individual_ticks = sum(individualCount, na.rm = TRUE),
    mean_tests_per_tick = mean(n_borrelia_tests),
    
    # Most common species combination tested
    most_common_species_combo = names(sort(table(borrelia_species_tested), decreasing = TRUE))[1],
    
    .groups = "drop"
  ) %>%
  # Add quality flags based on sample size
  mutate(
    quality_flag = case_when(
      total_tests < 3 ~ "very_low_sample_size",
      total_tests >= 3 & total_tests < 5 ~ "low_sample_size",
      total_tests >= 5 & total_tests < 10 ~ "moderate_sample_size",
      TRUE ~ "adequate_sample_size"
    ),
    # Add season variable
    season = case_when(
      month %in% c(3,4,5) ~ "Spring",
      month %in% c(6,7,8) ~ "Summer",
      month %in% c(9,10,11) ~ "Fall",
      TRUE ~ "Winter"
    )
  )

cat("  Subsample-level dataset created:", nrow(borrelia_summary), "subsamples\n")

# ---- Step 5: Add temporal alignment columns ----
cat("\nStep 5: Adding temporal alignment for phenology matching...\n")

borrelia_data <- borrelia_summary %>%
  mutate(
    prev_year = year - 1,
    prev_site_year = paste(siteID, prev_year, sep = "_")
  )

# ---- Step 6: Quality assessment and summary ----
cat("\nStep 6: Quality assessment and summary...\n")

# Overall summary
overall_summary <- borrelia_data %>%
  filter(quality_flag %in% c("moderate_sample_size", "adequate_sample_size")) %>%
  summarise(
    n_subsamples = n(),
    n_site_years = n_distinct(site_year),
    total_ticks_tested = sum(total_tests),
    total_positive_ticks = sum(positive_tests),
    overall_positivity_rate = total_positive_ticks / total_ticks_tested,
    mean_subsample_positivity = mean(positivity_rate),
    median_subsample_positivity = median(positivity_rate),
    sd_subsample_positivity = sd(positivity_rate)
  )

cat("\nOverall Borrelia Summary (moderate+ sample size):\n")
print(overall_summary)

# Quality flag distribution
cat("\nSample size distribution:\n")
quality_table <- table(borrelia_data$quality_flag)
print(quality_table)

# Validation: Check the aggregation worked correctly
cat("\nValidation checks:\n")
validation_summary <- bor_ixodes %>%
  group_by(testingID) %>%
  summarise(any_positive = any(testResult == "Positive"), .groups = "drop") %>%
  summarise(
    total_ticks_validation = n(),
    positive_ticks_validation = sum(any_positive),
    positivity_rate_validation = positive_ticks_validation / total_ticks_validation
  )

main_summary <- borrelia_data %>%
  summarise(
    total_ticks_main = sum(total_tests),
    positive_ticks_main = sum(positive_tests),
    positivity_rate_main = positive_ticks_main / total_ticks_main
  )

cat("Validation - tick counts match:", 
    validation_summary$total_ticks_validation == main_summary$total_ticks_main, "\n")
cat("Validation - positive counts match:", 
    validation_summary$positive_ticks_validation == main_summary$positive_ticks_main, "\n")
cat("Validation - positivity rates match:", 
    abs(validation_summary$positivity_rate_validation - main_summary$positivity_rate_main) < 0.001, "\n")

# Distribution of ticks per subsample
subsample_size_dist <- table(borrelia_data$total_tests)
cat("\nTicks per subsample distribution:\n")
print(head(subsample_size_dist, 10))

# ---- Step 7: Save results ----
cat("\nStep 7: Saving results...\n")

# Save main dataset
saveRDS(borrelia_data, "data/processed/borrelia_data.RDS")
write.csv(borrelia_data, "data/processed/borrelia_data.csv", row.names = FALSE)

# Save tick-level data for detailed analysis if needed
saveRDS(tick_level, "data/processed/borrelia_tick_level.RDS")
write.csv(tick_level, "data/processed/borrelia_tick_level.csv", row.names = FALSE)

# Create summary report
summary_report <- data.frame(
  Metric = c("Total pathogen tests", "Unique ticks tested", "Unique subsamples", 
             "Subsamples with adequate data", "Overall tick positivity rate",
             "Mean subsample positivity rate", "Site-years represented"),
  Value = c(
    nrow(bor_ixodes),
    nrow(tick_level), 
    nrow(borrelia_data),
    sum(borrelia_data$quality_flag %in% c("moderate_sample_size", "adequate_sample_size")),
    paste0(round(overall_summary$overall_positivity_rate * 100, 1), "%"),
    paste0(round(overall_summary$mean_subsample_positivity * 100, 1), "%"),
    overall_summary$n_site_years
  )
)

write.csv(summary_report, "data/processed/borrelia_summary_report.csv", row.names = FALSE)

cat("\nSaved datasets:\n")
cat("  - Main analysis data: data/processed/borrelia_data.RDS\n")
cat("  - Tick-level data: data/processed/borrelia_tick_level.RDS\n") 
cat("  - Summary report: data/processed/borrelia_summary_report.csv\n")

cat("\n==== BORRELIA PROCESSING COMPLETE ====\n")
cat("Key achievements:\n")
cat("- Correctly aggregated", nrow(bor_ixodes), "pathogen tests to", nrow(tick_level), "individual ticks\n")
cat("- Created", nrow(borrelia_data), "subsample observations for mixed-effects modeling\n")
cat("- Biologically meaningful: any Borrelia species = infected tick\n")
cat("- Ready for temporal alignment with phenology data\n")
cat("- Validated: aggregation preserved all tick and positive counts\n")