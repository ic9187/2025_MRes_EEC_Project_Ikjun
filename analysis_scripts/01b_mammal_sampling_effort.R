# ============================================================================
# Script: 01b_mammal_sampling_effort.R  
# Purpose: Apply bout interval analysis to identify site-years with consistent 
#          mammal sampling effort
# Input: trap_raw (from 00_setup.R)
# Output: clean_mammal_site_years (vector of site_years with consistent sampling)
# ============================================================================

cat("\n==== MAMMAL SAMPLING EFFORT ANALYSIS ====\n")

# ---- Step 1: Initial data cleaning ----
cat("Step 1: Initial data cleaning...\n")

trap_filtered <- trap_raw %>%
  filter(samplingImpractical == "OK" | 
           is.na(samplingImpractical) | 
           samplingImpractical == "") %>%
  mutate(
    collectDate = as.Date(collectDate),
    year = year(collectDate),
    site_year = paste(siteID, year, sep = "_")
  ) %>%
  filter(!is.na(collectDate)) %>%
  dplyr::select(siteID, year, site_year, collectDate) %>%
  distinct() %>%
  arrange(siteID, year, collectDate)

cat("  Rows after filtering:", nrow(trap_filtered), "\n")

# ---- Step 2: Identify bout groupings ----
cat("\nStep 2: Grouping consecutive days into bouts...\n")

bout_identification <- trap_filtered %>%
  group_by(siteID, year) %>%
  mutate(
    days_since_last = as.numeric(collectDate - lag(collectDate)),
    new_bout = is.na(days_since_last) | days_since_last > 7,  # >1 week = new bout
    bout_id = cumsum(new_bout)
  ) %>%
  ungroup()

# ---- Step 3: Create bout-level summary ----
cat("Step 3: Creating bout-level data...\n")

bout_summary <- bout_identification %>%
  group_by(siteID, year, site_year, bout_id) %>%
  summarise(
    bout_date = min(collectDate),  # Use first day as bout date
    bout_duration = as.numeric(max(collectDate) - min(collectDate)) + 1,
    n_trap_nights = n(),
    .groups = "drop"
  )

cat("  Total bouts identified:", nrow(bout_summary), "\n")

# ---- Step 4: Calculate bout intervals ----
cat("Step 4: Calculating bout intervals...\n")

bout_intervals <- bout_summary %>%
  group_by(site_year, siteID, year) %>%
  arrange(bout_date) %>%
  mutate(interval = as.numeric(bout_date - lag(bout_date))) %>%
  filter(!is.na(interval)) %>%
  ungroup()

cat("  Bout intervals calculated:", nrow(bout_intervals), "\n")

# ---- Step 5: Optimized reference level selection ----
cat("Step 5: Testing all possible reference levels...\n")

test_reference_levels <- function(interval_data) {
  site_years <- unique(interval_data$site_year)
  results <- list()
  
  for(ref_level in site_years) {
    cat("  Testing reference:", ref_level, "\r")
    
    # Set this site_year as reference level
    test_data <- interval_data
    test_data$site_year <- factor(test_data$site_year, 
                                  levels = c(ref_level, setdiff(site_years, ref_level)))
    
    # Run model
    tryCatch({
      model <- lm(interval ~ site_year, data = test_data)
      coef_table <- summary(model)$coefficients
      p_values <- coef_table[,4]
      term_names <- rownames(coef_table)
      
      # Find non-significant terms (p >= 0.05)
      non_sig_terms <- term_names[p_values >= 0.05 & term_names != "(Intercept)"]
      non_sig_site_years <- gsub("site_year", "", non_sig_terms)
      
      # Include reference level
      clean_site_years <- c(ref_level, non_sig_site_years)
      
      results[[ref_level]] <- list(
        reference = ref_level,
        n_retained = length(clean_site_years),
        retention_rate = length(clean_site_years) / length(site_years),
        clean_site_years = clean_site_years
      )
    }, error = function(e) {
      results[[ref_level]] <- list(
        reference = ref_level,
        n_retained = 0,
        retention_rate = 0,
        clean_site_years = character(0)
      )
    })
  }
  
  return(results)
}

# Test all reference levels
reference_tests <- test_reference_levels(bout_intervals)

# Find the best reference level (maximum retention)
retention_summary <- sapply(reference_tests, function(x) x$n_retained)
best_reference <- names(retention_summary)[which.max(retention_summary)]
best_result <- reference_tests[[best_reference]]

cat("\n")
cat("  Best reference level:", best_reference, "\n")
cat("  Site-years retained:", best_result$n_retained, 
    "out of", length(unique(bout_intervals$site_year)), "\n")
cat("  Retention rate:", round(best_result$retention_rate * 100, 1), "%\n")

# Use the optimal result
clean_mammal_site_years <- best_result$clean_site_years

# Show retention comparison
cat("\nRetention rate comparison (top 5):\n")
top_5 <- sort(retention_summary, decreasing = TRUE)[1:5]
for(i in 1:length(top_5)) {
  cat("  ", names(top_5)[i], ":", top_5[i], "site-years\n")
}

# ---- Step 6: Save results ----
cat("\nStep 6: Saving results...\n")

# Save the clean mammal site_years list
saveRDS(clean_mammal_site_years, "data/processed/clean_mammal_site_years.RDS")
write.csv(data.frame(site_year = clean_mammal_site_years), 
          "data/processed/clean_mammal_site_years.csv", 
          row.names = FALSE)

# Save bout summary for inspection
write.csv(bout_summary, "data/processed/bout_summary.csv", row.names = FALSE)

cat("\n  Saved clean_mammal_site_years to: data/processed/clean_mammal_site_years.RDS\n")
cat("\n==== MAMMAL SAMPLING EFFORT ANALYSIS COMPLETE ====\n")