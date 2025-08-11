# ============================================================================
# TICK DATA PIPELINE SUMMARY - UPDATED
# ============================================================================
# Comprehensive flow: collected → identified → tested → positive
# Now includes BOTH Ixodes and non-Ixodes Borrelia testing
# Uses filtered data from scripts 01, 02a, and 02b

cat("\n==== UPDATED TICK DATA PIPELINE SUMMARY ====\n")

# ---- Step 1: Load clean site_years filter ----
clean_site_years <- readRDS("data/processed/clean_site_years.RDS")
cat("Using", length(clean_site_years), "clean site_years for consistent sampling effort\n\n")

# ---- Step 2: Total ticks collected (from field counts) ----
cat("=== STEP 1: FIELD COLLECTION ===\n")

# Apply same filtering as script 01
tick_collected <- tick_raw %>%
  filter(
    samplingImpractical == "OK" | is.na(samplingImpractical) | samplingImpractical == "",
    !(is.na(nymphCount) & is.na(larvaCount) & is.na(adultCount)),
    year != 2019,
    sampleCondition == "No known compromise" | is.na(sampleCondition) | sampleCondition == ""
  ) %>%
  mutate(site_year = paste(siteID, year, sep = "_")) %>%
  filter(site_year %in% clean_site_years)

# Calculate total ticks collected
total_collected <- tick_collected %>%
  summarise(
    total_nymphs = sum(nymphCount, na.rm = TRUE),
    total_larvae = sum(larvaCount, na.rm = TRUE), 
    total_adults = sum(adultCount, na.rm = TRUE),
    sampling_events = n(),
    .groups = "drop"
  ) %>%
  mutate(total_ticks_collected = total_nymphs + total_larvae + total_adults)

cat("Total ticks collected in field:\n")
cat("  - Nymphs:", format(total_collected$total_nymphs, big.mark = ","), "\n")
cat("  - Larvae:", format(total_collected$total_larvae, big.mark = ","), "\n")
cat("  - Adults:", format(total_collected$total_adults, big.mark = ","), "\n")
cat("  - TOTAL:", format(total_collected$total_ticks_collected, big.mark = ","), "\n")
cat("  - From", total_collected$sampling_events, "sampling events\n\n")

# ---- Step 3: Ticks sent for taxonomic identification ----
cat("=== STEP 2: TAXONOMIC IDENTIFICATION ===\n")

# Apply same filtering as script 02a
taxo_identified <- taxo_raw %>%
  filter(sampleCondition == "OK") %>%
  mutate(
    collectDate = as.Date(collectDate),
    year = year(collectDate),
    siteID = sub("_.*", "", plotID),
    site_year = paste(siteID, year, sep = "_")
  ) %>%
  filter(site_year %in% clean_site_years)

# Calculate identification numbers
identification_summary <- taxo_identified %>%
  summarise(
    total_identified = sum(individualCount, na.rm = TRUE),
    identification_records = n(),
    unique_samples = n_distinct(sampleID),
    .groups = "drop"
  )

# Ixodes specifically
ixodes_summary <- taxo_identified %>%
  filter(!genus %in% c("Amblyomma", "Haemaphysalis", "Dermacentor")) %>%
  summarise(
    ixodes_identified = sum(individualCount, na.rm = TRUE),
    ixodes_records = n(),
    .groups = "drop"
  )

# Non-Ixodes specifically  
non_ixodes_summary <- taxo_identified %>%
  filter(genus %in% c("Amblyomma", "Haemaphysalis", "Dermacentor")) %>%
  summarise(
    non_ixodes_identified = sum(individualCount, na.rm = TRUE),
    non_ixodes_records = n(),
    .groups = "drop"
  )

cat("Taxonomic identification results:\n")
cat("  - Total ticks identified:", format(identification_summary$total_identified, big.mark = ","), "\n")
cat("  - Ixodes ticks identified:", format(ixodes_summary$ixodes_identified, big.mark = ","), "\n")
cat("  - Non-Ixodes ticks identified:", format(non_ixodes_summary$non_ixodes_identified, big.mark = ","), "\n")
cat("  - Proportion Ixodes:", round(ixodes_summary$ixodes_identified / identification_summary$total_identified * 100, 1), "%\n")
cat("  - Proportion Non-Ixodes:", round(non_ixodes_summary$non_ixodes_identified / identification_summary$total_identified * 100, 1), "%\n")
cat("  - From", identification_summary$unique_samples, "unique samples\n\n")

# ---- Step 4: Ticks tested for Borrelia (ALL species) ----
cat("=== STEP 3: BORRELIA PATHOGEN TESTING (ALL SPECIES) ===\n")

# Apply same filtering as script 02b but INCLUDE non-Ixodes
pathogen_tested_all <- path_raw %>%
  filter(
    grepl("Borrelia", testPathogenName, ignore.case = TRUE),
    !is.na(testResult)
    # REMOVED: grepl("IXO", subsampleID) - now includes non-Ixodes!
  ) %>%
  mutate(
    collectDate = as.Date(collectDate),
    year = year(collectDate),
    siteID = sub("_.*", "", plotID),
    site_year = paste(siteID, year, sep = "_"),
    # Determine tick type from subsampleID
    tick_type = ifelse(grepl("IXO", subsampleID), "Ixodes", "Non-Ixodes")
  ) %>%
  filter(site_year %in% clean_site_years)

# Aggregate by testingID (tick level) for ALL species
tick_level_all <- pathogen_tested_all %>%
  group_by(testingID, tick_type) %>%
  summarise(
    tick_positive = any(testResult == "Positive"),
    n_tests_per_tick = n(),
    .groups = "drop"
  )

# Calculate testing summary for ALL species
testing_summary_all <- tick_level_all %>%
  summarise(
    total_ticks_tested = n(),
    positive_ticks = sum(tick_positive),
    .groups = "drop"
  )

# Calculate testing summary by tick type
testing_by_type <- tick_level_all %>%
  group_by(tick_type) %>%
  summarise(
    ticks_tested = n(),
    positive_ticks = sum(tick_positive),
    positivity_rate = round(mean(tick_positive) * 100, 1),
    .groups = "drop"
  )

cat("Borrelia testing results (ALL species):\n")
cat("  - Total ticks tested:", format(testing_summary_all$total_ticks_tested, big.mark = ","), "\n")
cat("  - Total positive ticks:", format(testing_summary_all$positive_ticks, big.mark = ","), "\n")
cat("  - Overall positivity rate:", round(testing_summary_all$positive_ticks / testing_summary_all$total_ticks_tested * 100, 1), "%\n\n")

cat("Breakdown by tick type:\n")
print(testing_by_type)
cat("\n")

# Extract Ixodes-specific numbers for downstream calculations
ixodes_testing <- testing_by_type %>% filter(tick_type == "Ixodes")
non_ixodes_testing <- testing_by_type %>% filter(tick_type == "Non-Ixodes")

# Handle case where non-Ixodes might be empty
if(nrow(non_ixodes_testing) == 0) {
  non_ixodes_testing <- data.frame(tick_type = "Non-Ixodes", ticks_tested = 0, positive_ticks = 0, positivity_rate = 0)
}

cat("Testing breakdown:\n")
cat("  - Ixodes tested:", format(ixodes_testing$ticks_tested, big.mark = ","), "\n")
cat("  - Non-Ixodes tested:", format(non_ixodes_testing$ticks_tested, big.mark = ","), "\n")
cat("  - Proportion Ixodes tested:", round(ixodes_testing$ticks_tested / testing_summary_all$total_ticks_tested * 100, 1), "%\n")
cat("  - Proportion Non-Ixodes tested:", round(non_ixodes_testing$ticks_tested / testing_summary_all$total_ticks_tested * 100, 1), "%\n\n")

# ---- Step 5: Complete pipeline flow summary ----
cat("=== COMPLETE PIPELINE FLOW ===\n")

# Calculate the requested proportions
prop_identified_collected <- round(identification_summary$total_identified / total_collected$total_ticks_collected * 100, 1)
prop_ixodes_others <- round(ixodes_summary$ixodes_identified / identification_summary$total_identified * 100, 1)
prop_tested_identified <- round(testing_summary_all$total_ticks_tested / identification_summary$total_identified * 100, 1)
prop_ixodes_tested_others <- round(ixodes_testing$ticks_tested / testing_summary_all$total_ticks_tested * 100, 1)
prop_ixodes_positive_tested <- round(ixodes_testing$positive_ticks / ixodes_testing$ticks_tested * 100, 1)

pipeline_flow <- data.frame(
  Metric = c(
    "Total collected",
    "Proportion identified / collected", 
    "Proportion Ixodes identified / others",
    "Proportion ticks tested / identified",
    "Proportion Ixodes tested / others tested",
    "Proportion Ixodes positive / Ixodes tested"
  ),
  Value = c(
    format(total_collected$total_ticks_collected, big.mark = ","),
    paste0(prop_identified_collected, "%"),
    paste0(prop_ixodes_others, "%"),
    paste0(prop_tested_identified, "%"), 
    paste0(prop_ixodes_tested_others, "%"),
    paste0(prop_ixodes_positive_tested, "%")
  ),
  Raw_Numbers = c(
    paste0(format(total_collected$total_ticks_collected, big.mark = ","), " ticks"),
    paste0(format(identification_summary$total_identified, big.mark = ","), " / ", format(total_collected$total_ticks_collected, big.mark = ",")),
    paste0(format(ixodes_summary$ixodes_identified, big.mark = ","), " / ", format(identification_summary$total_identified, big.mark = ",")),
    paste0(format(testing_summary_all$total_ticks_tested, big.mark = ","), " / ", format(identification_summary$total_identified, big.mark = ",")),
    paste0(format(ixodes_testing$ticks_tested, big.mark = ","), " / ", format(testing_summary_all$total_ticks_tested, big.mark = ",")),
    paste0(format(ixodes_testing$positive_ticks, big.mark = ","), " / ", format(ixodes_testing$ticks_tested, big.mark = ","))
  )
)

print(pipeline_flow)

# ---- Step 6: Data attrition visualization ----
cat("\n=== DATA ATTRITION ===\n")

attrition_data <- c(
  "Collected" = total_collected$total_ticks_collected,
  "ID'd" = identification_summary$total_identified,
  "Ixodes" = ixodes_summary$ixodes_identified,
  "Tested" = testing_summary_all$total_ticks_tested,
  "Positive" = testing_summary_all$positive_ticks
)

# Calculate losses between stages
losses <- c(
  total_collected$total_ticks_collected - identification_summary$total_identified,
  identification_summary$total_identified - ixodes_summary$ixodes_identified,
  ixodes_summary$ixodes_identified - ixodes_testing$ticks_tested,  # Note: comparing Ixodes to Ixodes tested
  ixodes_testing$ticks_tested - ixodes_testing$positive_ticks
)

cat("Data attrition between stages:\n")
cat("  - Collection → ID:", format(losses[1], big.mark = ","), "ticks lost\n")
cat("  - All species → Ixodes:", format(losses[2], big.mark = ","), "ticks lost\n") 
cat("  - Ixodes ID'd → Ixodes Tested:", format(losses[3], big.mark = ","), "ticks lost\n")
cat("  - Ixodes Tested → Ixodes Positive:", format(losses[4], big.mark = ","), "ticks lost\n")

# ---- Step 7: Save pipeline data for maps ----
cat("\n=== SAVING PIPELINE DATA ===\n")

pipeline_data <- list(
  total_collected = total_collected,
  identification_summary = identification_summary,
  ixodes_summary = ixodes_summary,
  non_ixodes_summary = non_ixodes_summary,
  testing_summary_all = testing_summary_all,
  testing_by_type = testing_by_type,
  pipeline_flow = pipeline_flow,
  attrition_data = attrition_data
)

# Save for use in maps script
saveRDS(pipeline_data, "data/processed/pipeline_summary_data.RDS")

cat("✓ Pipeline summary data saved for maps script consistency\n")
cat("\n==== UPDATED TICK PIPELINE SUMMARY COMPLETE ====\n")