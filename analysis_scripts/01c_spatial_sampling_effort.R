# ============================================================================
# Script: 01c_spatial_sampling_effort.R
# Purpose: Add spatial sampling effort filtering (≥4 plots) on top of existing 
#          temporal filtering from 01a and 01b
# Input: clean_site_years.RDS, clean_mammal_site_years.RDS (from 01a, 01b)
# Output: Updated clean_site_years.RDS, clean_mammal_site_years.RDS with spatial filtering
# Key decision: ≥4 plots per site (mean ± 2, symmetric filtering)
# ============================================================================

cat("\n==== SPATIAL SAMPLING EFFORT ANALYSIS ====\n")
cat("Adding spatial filtering (≥4 plots) on top of existing temporal filtering\n\n")

library(dplyr)
library(tidyr)

# ---- Step 1: Load existing temporal filters ----
cat("Step 1: Loading existing temporal filters...\n")
clean_site_years <- readRDS("data/processed/clean_site_years.RDS")
clean_mammal_site_years <- readRDS("data/processed/clean_mammal_site_years.RDS")

cat("  Temporally filtered tick site_years:", length(clean_site_years), "\n")
cat("  Temporally filtered mammal site_years:", length(clean_mammal_site_years), "\n")

# ---- Step 2: Create tick plots per site-year table ----
cat("\nStep 2: Creating tick plots per site-year table...\n")

tick_plots_table <- tick_raw %>%
  filter(
    samplingImpractical == "OK" | is.na(samplingImpractical) | samplingImpractical == "",
    !(is.na(nymphCount) & is.na(larvaCount) & is.na(adultCount)),
    year != 2019,
    sampleCondition == "No known compromise" | is.na(sampleCondition) | sampleCondition == ""
  ) %>%
  mutate(
    siteID = sub("_.*", "", plotID),
    site_year = paste(siteID, year, sep = "_")
  ) %>%
  filter(site_year %in% clean_site_years) %>%
  group_by(siteID, year) %>%
  summarise(n_plots = n_distinct(plotID), .groups = "drop")

# Wide format
tick_plots_wide <- tick_plots_table %>%
  pivot_wider(names_from = year, values_from = n_plots, values_fill = NA) %>%
  rowwise() %>%
  mutate(total_plots = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(total_plots)) %>%
  {
    year_cols <- names(.)[!names(.) %in% c("siteID", "total_plots")]
    year_cols_ordered <- sort(as.numeric(year_cols))
    dplyr::select(., siteID, all_of(as.character(year_cols_ordered)), total_plots)
  }

# Long format
tick_plots_long <- tick_plots_table %>%
  mutate(year = as.numeric(year)) %>%
  arrange(siteID, year)

cat("  Tick plots - Wide format:", nrow(tick_plots_wide), "sites\n")
cat("  Tick plots - Long format:", nrow(tick_plots_long), "site-years\n")
cat("  Tick plots mean:", round(mean(tick_plots_long$n_plots), 2), "\n")

# ---- Step 3: Create mammal plots per site-year table ----
cat("\nStep 3: Creating mammal plots per site-year table...\n")

mammal_plots_table <- trap_raw %>%
  filter(samplingImpractical == "OK" | is.na(samplingImpractical) | samplingImpractical == "") %>%
  mutate(
    collectDate = as.Date(collectDate),
    year = year(collectDate),
    site_year = paste(siteID, year, sep = "_")
  ) %>%
  filter(site_year %in% clean_mammal_site_years) %>%
  group_by(siteID, year) %>%
  summarise(n_plots = n_distinct(plotID), .groups = "drop")

# Wide format
mammal_plots_wide <- mammal_plots_table %>%
  pivot_wider(names_from = year, values_from = n_plots, values_fill = NA) %>%
  rowwise() %>%
  mutate(total_plots = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(total_plots)) %>%
  {
    year_cols <- names(.)[!names(.) %in% c("siteID", "total_plots")]
    year_cols_ordered <- sort(as.numeric(year_cols))
    dplyr::select(., siteID, all_of(as.character(year_cols_ordered)), total_plots)
  }

# Long format
mammal_plots_long <- mammal_plots_table %>%
  mutate(year = as.numeric(year)) %>%
  arrange(siteID, year)

cat("  Mammal plots - Wide format:", nrow(mammal_plots_wide), "sites\n")
cat("  Mammal plots - Long format:", nrow(mammal_plots_long), "site-years\n")
cat("  Mammal plots mean:", round(mean(mammal_plots_long$n_plots), 2), "\n")

# ---- Step 4: Apply spatial filtering (≥4 plots) ----
cat("\nStep 4: Applying spatial filtering (≥4 plots)...\n")
cat("Filter logic: Mean ± 2 plots (symmetric filtering)\n")

# Identify site_years with ≥4 plots
tick_spatial_filtered <- tick_plots_long %>%
  filter(n_plots >= 4) %>%
  mutate(site_year = paste(siteID, year, sep = "_")) %>%
  pull(site_year)

mammal_spatial_filtered <- mammal_plots_long %>%
  filter(n_plots >= 4) %>%
  mutate(site_year = paste(siteID, year, sep = "_")) %>%
  pull(site_year)

# Apply spatial filter to existing temporal filters
clean_site_years_final <- clean_site_years[clean_site_years %in% tick_spatial_filtered]
clean_mammal_site_years_final <- clean_mammal_site_years[clean_mammal_site_years %in% mammal_spatial_filtered]

cat("  Original tick site_years:", length(clean_site_years), "\n")
cat("  Spatially filtered tick site_years:", length(clean_site_years_final), "\n")
cat("  Retention rate:", round(length(clean_site_years_final)/length(clean_site_years)*100, 1), "%\n")

cat("  Original mammal site_years:", length(clean_mammal_site_years), "\n")
cat("  Spatially filtered mammal site_years:", length(clean_mammal_site_years_final), "\n")
cat("  Retention rate:", round(length(clean_mammal_site_years_final)/length(clean_mammal_site_years)*100, 1), "%\n")

# ---- Step 5: Update filter files ----
cat("\nStep 5: Updating filter files...\n")

saveRDS(clean_site_years_final, "data/processed/clean_site_years.RDS")
saveRDS(clean_mammal_site_years_final, "data/processed/clean_mammal_site_years.RDS")

cat("  ✓ Updated: data/processed/clean_site_years.RDS\n")
cat("  ✓ Updated: data/processed/clean_mammal_site_years.RDS\n")

# ---- Step 6: Create filtered tables for display ----
cat("\nStep 6: Creating filtered tables...\n")

# Filtered tick table
tick_plots_filtered <- tick_plots_long %>%
  filter(n_plots >= 4) %>%
  group_by(siteID) %>%
  summarise(
    n_years = n(),
    total_plots = sum(n_plots),
    mean_plots = round(mean(n_plots), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(total_plots))

# Filtered mammal table
mammal_plots_filtered <- mammal_plots_long %>%
  filter(n_plots >= 4) %>%
  group_by(siteID) %>%
  summarise(
    n_years = n(),
    total_plots = sum(n_plots),
    mean_plots = round(mean(n_plots), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(total_plots))

# ---- Step 7: Display final results ----
cat("\n", rep("=", 60), "\n")
cat("FINAL SPATIALLY FILTERED TABLES (≥4 plots)\n")
cat(rep("=", 60), "\n")

cat("\nTICK SITES (spatially + temporally filtered):\n")
print(tick_plots_filtered)

cat("\nMAMMAL SITES (spatially + temporally filtered):\n")
print(mammal_plots_filtered)

cat("\nFINAL SUMMARY:\n")
cat("Tick sites:", nrow(tick_plots_filtered), "\n")
cat("Tick site-years:", length(clean_site_years_final), "\n")
cat("Mammal sites:", nrow(mammal_plots_filtered), "\n")
cat("Mammal site-years:", length(clean_mammal_site_years_final), "\n")

cat("\n✓ Spatial sampling effort filtering complete!\n")
cat("✓ Downstream scripts (02a, 02c) will automatically use updated filters\n")

cat("\n==== SPATIAL SAMPLING EFFORT ANALYSIS COMPLETE ====\n")