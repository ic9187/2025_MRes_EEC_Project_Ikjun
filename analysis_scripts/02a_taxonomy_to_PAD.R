# ============================================================================
# Script: 02a_taxonomy_to_PAD.R
# Purpose: Calculate Peak Activity Difference (PAD) from taxonomy data
# Input: taxo_raw (from 00_setup.R), clean_site_years (from 01_tick_sampling_effort.R)
# Output: phenology_data with PAD values per site_year
# Key decisions:
#   - Using maximum counts per day
#   - Taking mean when multiple days have same max count
#   - Filtering for Ixodes only
# ============================================================================

cat("\n==== TAXONOMY TO PAD ANALYSIS ====\n")

# ---- Step 1: Load site_year filter ----
cat("Step 1: Loading site_year filter...\n")
clean_site_years <- readRDS("data/processed/clean_site_years.RDS")
cat("  Loaded", length(clean_site_years), "clean site_years\n")

# ---- Step 2: Filter taxonomy data ----
cat("\nStep 2: Filtering taxonomy data...\n")

# Initial filtering
ixo_filtered <- taxo_raw %>% 
  # Keep only Ixodes genus (exclude other genera)
  filter(!genus %in% c("Amblyomma", "Haemaphysalis", "Dermacentor")) %>%
  filter(sampleCondition == "OK") %>%
  # Add temporal columns
  mutate(
    collectDate = as.Date(collectDate),
    day_number = yday(collectDate),
    year = year(collectDate),
    siteID = sub("_.*", "", plotID),
    site_year = paste(siteID, year, sep = "_")
  ) %>%
  # Apply clean_site_years filter
  filter(site_year %in% clean_site_years)

cat("  Rows after Ixodes filtering:", nrow(ixo_filtered), "\n")
cat("  Site_years represented:", n_distinct(ixo_filtered$site_year), "\n")

# ---- Step 3: Transform to life stage counts ----
cat("\nStep 3: Creating life stage counts...\n")

# Create life stage categories
ixo_lifestage <- ixo_filtered %>%
  filter(!is.na(sexOrAge)) %>%
  mutate(
    life_stage = case_when(
      sexOrAge == "Larva" ~ "larva",
      sexOrAge == "Nymph" ~ "nymph",  
      sexOrAge %in% c("Male", "Female") ~ "adult",
      TRUE ~ NA_character_
    )
  )

# Aggregate to daily counts by site and life stage
daily_counts <- ixo_lifestage %>%
  group_by(siteID, year, site_year, day_number, life_stage) %>%
  summarise(
    count = sum(individualCount, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = life_stage,
    values_from = count,
    values_fill = 0,
    names_prefix = "count_"
  )

cat("  Daily count records created:", nrow(daily_counts), "\n")

# ---- Step 4: Find peak days ----
cat("\nStep 4: Finding peak activity days...\n")

# Calculate peaks for each site_year
peak_data <- daily_counts %>%
  group_by(site_year, siteID, year) %>%
  summarise(
    # Maximum counts
    nymph_max = max(count_nymph, na.rm = TRUE),
    larva_max = max(count_larva, na.rm = TRUE),
    
    # Peak days (handling ties by taking mean)
    nymph_peak_day = ifelse(
      nymph_max > 0,
      mean(day_number[count_nymph == nymph_max]),
      NA
    ),
    larva_peak_day = ifelse(
      larva_max > 0,
      mean(day_number[count_larva == larva_max]),
      NA
    ),
    
    # Number of sampling days
    n_sampling_days = n(),
    
    .groups = "drop"
  ) %>%
  # Filter out site_years with no activity
  filter(nymph_max > 0 & larva_max > 0)

# Calculate PAD
peak_data <- peak_data %>%
  mutate(
    PAD = larva_peak_day - nymph_peak_day,
    # Categorize PAD
    PAD_category = case_when(
      PAD < -10 ~ "Larvae earlier",
      PAD >= -10 & PAD <= 10 ~ "Synchronous",
      PAD > 10 ~ "Nymphs earlier",
      TRUE ~ NA_character_
    )
  )

cat("  Site_years with valid PAD values:", nrow(peak_data), "\n")

# ---- Step 5: Summary and save ----
cat("\nStep 5: Creating summary and saving...\n")

# Summary statistics
summary_stats <- peak_data %>%
  summarise(
    n_site_years = n(),
    mean_PAD = mean(PAD, na.rm = TRUE),
    sd_PAD = sd(PAD, na.rm = TRUE),
    min_PAD = min(PAD, na.rm = TRUE),
    max_PAD = max(PAD, na.rm = TRUE)
  )

cat("\n  PAD Summary:\n")
print(summary_stats)

cat("\n  PAD by category:\n")
print(table(peak_data$PAD_category))

# Save phenology data
phenology_data <- peak_data %>%
  dplyr::select(site_year, siteID, year, PAD, PAD_category, 
         nymph_peak_day, larva_peak_day, nymph_max, larva_max)

saveRDS(phenology_data, "data/processed/phenology_data.RDS")
write.csv(phenology_data, "data/processed/phenology_data.csv", row.names = FALSE)

cat("\n  Saved phenology data to: data/processed/phenology_data.RDS\n")
cat("\n==== PAD CALCULATION COMPLETE ====\n")