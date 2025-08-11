# ============================================================================
# Script: 01a_tick_sampling_effort.R
# Purpose: Apply interval analysis to identify site-years with consistent 
#          sampling effort in tick count data
# Input: tick_raw (from 00_setup.R)
# Output: clean_site_years (vector of site_years with consistent sampling)
# Key decisions: 
#   - Using linear model to compare sampling intervals
#   - ABBY_2016 as reference site_year
#   - p >= 0.05 indicates similar sampling pattern
# ============================================================================

cat("\n==== TICK SAMPLING EFFORT ANALYSIS ====\n")

# ---- Step 1: Initial data cleaning ----
cat("Step 1: Initial data cleaning...\n")

# Filter by samplingImpractical
tick_filtered <- tick_raw %>%
  filter(samplingImpractical == "OK" | 
           is.na(samplingImpractical) | 
           samplingImpractical == "")

cat("  Rows after filtering samplingImpractical:", nrow(tick_filtered), "\n")

# Remove rows with all NA counts
tick_filtered <- tick_filtered %>%
  filter(!(is.na(nymphCount) & is.na(larvaCount) & is.na(adultCount)))

cat("  Rows after removing all-NA counts:", nrow(tick_filtered), "\n")

# Remove 2019 data (protocol changes)
tick_filtered <- tick_filtered %>% 
  filter(year != 2019)

cat("  Rows after removing 2019:", nrow(tick_filtered), "\n")

# Keep only samples with no compromise in condition
tick_filtered <- tick_filtered %>%
  filter(sampleCondition == "No known compromise" | 
           is.na(sampleCondition) | 
           sampleCondition == "")

cat("  Final filtered rows:", nrow(tick_filtered), "\n")

# ---- Step 2: Interval analysis ----
cat("\nStep 2: Calculating sampling intervals...\n")

# Calculate intervals between consecutive sampling days
intervals <- tick_filtered %>%
  dplyr::select(siteID, year, day_number) %>% 
  distinct() %>%
  mutate(site_year = paste(siteID, year, sep="_")) %>%
  arrange(site_year, day_number) %>%
  group_by(site_year, siteID, year) %>%
  mutate(interval = day_number - lag(day_number)) %>%
  filter(!is.na(interval))

cat("  Total intervals calculated:", nrow(intervals), "\n")

# ---- Step 3: Optimized reference level selection for interval comparison ----
cat("\nStep 3: Testing all reference levels for maximum retention...\n")

# Function to test all possible reference levels
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
reference_tests <- test_reference_levels(intervals)

# Find the best reference level (maximum retention)
retention_summary <- sapply(reference_tests, function(x) x$n_retained)
best_reference <- names(retention_summary)[which.max(retention_summary)]
best_result <- reference_tests[[best_reference]]

cat("\n")
cat("  Best reference level:", best_reference, "\n")
cat("  Site-years retained:", best_result$n_retained, 
    "out of", length(unique(intervals$site_year)), "\n")
cat("  Retention rate:", round(best_result$retention_rate * 100, 1), "%\n")

# Use the optimal result
clean_site_years <- best_result$clean_site_years

# Show retention comparison
cat("\nRetention rate comparison (top 5):\n")
top_5 <- sort(retention_summary, decreasing = TRUE)[1:5]
for(i in 1:length(top_5)) {
  cat("  ", names(top_5)[i], ":", top_5[i], "site-years\n")
}

# ---- Step 3b: Re-run best model for detailed examination ----
cat("\nRe-running best model for detailed summary...\n")

# Re-run the model with the optimal reference level
intervals_best <- intervals
intervals_best$site_year <- factor(intervals_best$site_year, 
                                   levels = c(best_reference, 
                                              setdiff(unique(intervals$site_year), 
                                                      best_reference)))

# Fit the best model
best_model <- lm(interval ~ site_year, data = intervals_best)

# Display key results
cat("  Model R-squared:", round(summary(best_model)$r.squared, 4), "\n")
cat("  Reference site_year:", best_reference, "\n")

# ---- Step 4: Quality check and save (ENHANCED) ----
cat("\nStep 4: Saving results...\n")

# Original summary
all_site_years <- unique(intervals$site_year)
excluded_site_years <- setdiff(all_site_years, clean_site_years)

summary_df <- data.frame(
  total_site_years = length(all_site_years),
  retained_site_years = length(clean_site_years),
  excluded_site_years = length(excluded_site_years),
  retention_rate = length(clean_site_years) / length(all_site_years) * 100,
  best_reference = best_reference  # NEW: Record which reference was optimal
)

cat("  Summary:\n")
print(summary_df)

# Save the clean site_years list (same as before)
saveRDS(clean_site_years, "data/processed/clean_site_years.RDS")
write.csv(data.frame(site_year = clean_site_years), 
          "data/processed/clean_site_years.csv", 
          row.names = FALSE)

# NEW: Save optimization details
optimization_results <- data.frame(
  reference_tested = names(retention_summary),
  site_years_retained = retention_summary,
  retention_rate = retention_summary / length(all_site_years) * 100
) %>%
  arrange(desc(site_years_retained))

write.csv(optimization_results, "data/processed/tick_reference_optimization.csv", row.names = FALSE)

cat("\n  Saved clean_site_years to: data/processed/clean_site_years.RDS\n")
cat("  Saved optimization details to: data/processed/tick_reference_optimization.csv\n")