# ============================================================================
# Peak Day Models Sensitivity Analysis - UPDATED
# Testing decomposed models robustness by comparing original vs outlier-free models
# Goal: Validate decomposition analysis by checking which life stage drives VPD effects
# Updated: Focus on nymph and larva peak day models with/without nymph outliers
# ============================================================================

cat("=== PEAK DAY DECOMPOSITION SENSITIVITY ANALYSIS ===\n")
cat("Comparing original vs outlier-free peak day models to validate decomposition\n")
cat("Focus: How nymph outliers affect both nymph and larva peak day models\n\n")

# Load required libraries
library(dplyr)
library(lme4)
library(lmerTest)
library(performance)

# ============================================================================
# SECTION 1: IDENTIFY NYMPH OUTLIERS
# ============================================================================

cat("=== SECTION 1: IDENTIFYING NYMPH OUTLIERS ===\n")

# Extract outlier row numbers from nymph model only
nymph_outliers <- outlier_results[["nymph_peak"]]$observation[outlier_results[["nymph_peak"]]$is_any_outlier]

cat("Nymph model outliers:", length(nymph_outliers), "observations\n")
if(length(nymph_outliers) > 0) {
  cat("  Row numbers:", paste(nymph_outliers, collapse = ", "), "\n")
  
  # Show site_year information for nymph outliers
  nymph_outlier_site_years <- outlier_results[["nymph_peak"]]$site_year[outlier_results[["nymph_peak"]]$is_any_outlier]
  cat("  Site_years:", paste(nymph_outlier_site_years, collapse = ", "), "\n")
} else {
  cat("  No nymph outliers detected\n")
}

# ============================================================================
# SECTION 2: CREATE CLEAN DATASET
# ============================================================================

cat("\n=== SECTION 2: CREATING CLEAN DATASET ===\n")

# Create clean dataset by removing nymph outliers
if(length(nymph_outliers) > 0) {
  analysis_data_clean <- analysis_data[-nymph_outliers, ]
} else {
  analysis_data_clean <- analysis_data
  cat("No outliers to remove - using original dataset\n")
}

cat("Original dataset size:", nrow(analysis_data), "observations\n")
cat("Clean dataset size:", nrow(analysis_data_clean), "observations\n")
cat("Observations removed:", nrow(analysis_data) - nrow(analysis_data_clean), "\n")

# Verify required variables exist
required_vars <- c("PAD", "mean_annual_vpdmin_std", "siteID")
missing_vars <- required_vars[!required_vars %in% names(analysis_data_clean)]
if(length(missing_vars) > 0) {
  stop("Error: Missing required variables in clean dataset: ", paste(missing_vars, collapse = ", "))
}

# ============================================================================
# SECTION 3: BUILD CLEAN PEAK DAY MODELS
# ============================================================================

cat("\n=== SECTION 3: BUILDING CLEAN PEAK DAY MODELS ===\n")

# Verify required variables exist for peak day models
required_peak_vars <- c("nymph_peak_day", "larva_peak_day", "mean_annual_vpdmin_std", "siteID")
missing_peak_vars <- required_peak_vars[!required_peak_vars %in% names(analysis_data_clean)]
if(length(missing_peak_vars) > 0) {
  stop("Error: Missing required peak day variables in clean dataset: ", paste(missing_peak_vars, collapse = ", "))
}

cat("Building clean nymph peak model without outliers...\n")
nymph_model_clean <- lmer(nymph_peak_day ~ mean_annual_vpdmin_std + (1|siteID), 
                          data = analysis_data_clean)

cat("Building clean larva peak model without outliers...\n")
larva_model_clean <- lmer(larva_peak_day ~ mean_annual_vpdmin_std + (1|siteID), 
                          data = analysis_data_clean)

cat("✓ Both clean peak day models built successfully\n")

# ============================================================================
# SECTION 4: DIAGNOSTIC FUNCTION
# ============================================================================

run_quick_diagnostics <- function(model, model_name) {
  cat("Diagnosing:", model_name, "\n")
  
  n_obs <- nobs(model)
  converged <- model@optinfo$conv$opt == 0
  singular <- isSingular(model)
  
  # Normality test
  residuals_vals <- residuals(model)
  normality_p <- NA
  if (n_obs >= 3 && n_obs <= 5000) {
    tryCatch({
      shapiro_result <- shapiro.test(residuals_vals)
      normality_p <- shapiro_result$p.value
    }, error = function(e) {
      normality_p <- NA
    })
  }
  
  # R-squared
  r2_marginal <- NA
  r2_conditional <- NA
  tryCatch({
    r2_results <- performance::r2(model)
    r2_marginal <- r2_results$R2_marginal
    r2_conditional <- r2_results$R2_conditional
  }, error = function(e) {
    # Skip if fails
  })
  
  reliable <- converged && !singular && (is.na(normality_p) || normality_p >= 0.05)
  
  cat("  Convergence:", ifelse(converged, "✓", "✗"), "\n")
  cat("  Singularity:", ifelse(singular, "✗", "✓"), "\n") 
  cat("  Normality p-value:", round(normality_p, 4), ifelse(is.na(normality_p) || normality_p >= 0.05, "✓", "✗"), "\n")
  cat("  R² marginal:", round(r2_marginal, 4), "\n")
  cat("  R² conditional:", round(r2_conditional, 4), "\n")
  cat("  Overall reliable:", ifelse(reliable, "✓", "✗"), "\n\n")
  
  return(list(
    n_obs = n_obs,
    converged = converged,
    singular = singular,
    normality_p = normality_p,
    r2_marginal = r2_marginal,
    r2_conditional = r2_conditional,
    reliable = reliable
  ))
}

# ============================================================================
# SECTION 5: COMPARE DIAGNOSTICS
# ============================================================================

cat("=== SECTION 5: DIAGNOSTIC COMPARISON ===\n")

cat("ORIGINAL PEAK DAY MODELS (with nymph outliers):\n")
nymph_orig_diag <- run_quick_diagnostics(peak_models[["nymph_peak"]], "Nymph Peak (Original)")
larva_orig_diag <- run_quick_diagnostics(peak_models[["larva_peak"]], "Larva Peak (Original)")

cat("CLEAN PEAK DAY MODELS (without nymph outliers):\n")
nymph_clean_diag <- run_quick_diagnostics(nymph_model_clean, "Nymph Peak (Clean)")
larva_clean_diag <- run_quick_diagnostics(larva_model_clean, "Larva Peak (Clean)")

# ============================================================================
# SECTION 6: COMPARE MODEL SUMMARIES
# ============================================================================

cat("=== SECTION 6: MODEL SUMMARY COMPARISON ===\n")

cat("NYMPH PEAK DAY MODEL COMPARISON:\n")
cat(rep("=", 80), "\n")
cat("ORIGINAL NYMPH MODEL (with outliers):\n")
print(summary(peak_models[["nymph_peak"]]))

cat("\nCLEAN NYMPH MODEL (without outliers):\n")
print(summary(nymph_model_clean))

cat("\n", rep("=", 80), "\n")
cat("LARVA PEAK DAY MODEL COMPARISON:\n")
cat(rep("=", 80), "\n")
cat("ORIGINAL LARVA MODEL (with outliers):\n")
print(summary(peak_models[["larva_peak"]]))

cat("\nCLEAN LARVA MODEL (without outliers):\n")
print(summary(larva_model_clean))

# ============================================================================
# SECTION 7: KEY EFFECTS COMPARISON
# ============================================================================

cat("\n=== SECTION 7: KEY EFFECTS COMPARISON ===\n")

# Extract key coefficients for comparison
extract_key_effects <- function(model, model_name) {
  coef_table <- summary(model)$coefficients
  
  effects <- list(
    model_name = model_name,
    n_obs = nobs(model)
  )
  
  if("mean_annual_vpdmin_std" %in% rownames(coef_table)) {
    effects$vpdmin_coef <- coef_table["mean_annual_vpdmin_std", "Estimate"]
    effects$vpdmin_p <- coef_table["mean_annual_vpdmin_std", "Pr(>|t|)"]
  }
  
  return(effects)
}

# ============================================================================
# SECTION 7: KEY EFFECTS COMPARISON
# ============================================================================

cat("\n=== SECTION 7: KEY EFFECTS COMPARISON ===\n")

# Extract key coefficients for comparison
extract_key_effects <- function(model, model_name) {
  coef_table <- summary(model)$coefficients
  
  effects <- list(
    model_name = model_name,
    n_obs = nobs(model)
  )
  
  if("mean_annual_vpdmin_std" %in% rownames(coef_table)) {
    effects$vpdmin_coef <- coef_table["mean_annual_vpdmin_std", "Estimate"]
    effects$vpdmin_p <- coef_table["mean_annual_vpdmin_std", "Pr(>|t|)"]
  }
  
  return(effects)
}

# Extract effects for all models
nymph_orig_effects <- extract_key_effects(peak_models[["nymph_peak"]], "Nymph Original")
nymph_clean_effects <- extract_key_effects(nymph_model_clean, "Nymph Clean")
larva_orig_effects <- extract_key_effects(peak_models[["larva_peak"]], "Larva Original") 
larva_clean_effects <- extract_key_effects(larva_model_clean, "Larva Clean")

# Create comparison table
cat("COEFFICIENT COMPARISON TABLE:\n")
cat("Model                    | N   | VPD Min Coef | VPD Min P\n")
cat(rep("-", 55), "\n")

print_effects_row <- function(effects) {
  cat(sprintf("%-24s | %-3d | %-12.3f | %-9.4f\n",
              effects$model_name,
              effects$n_obs,
              ifelse(is.null(effects$vpdmin_coef), NA, effects$vpdmin_coef),
              ifelse(is.null(effects$vpdmin_p), NA, effects$vpdmin_p)))
}

print_effects_row(nymph_orig_effects)
print_effects_row(nymph_clean_effects)
print_effects_row(larva_orig_effects)
print_effects_row(larva_clean_effects)

# ============================================================================
# SECTION 8: DECOMPOSITION VALIDATION SUMMARY
# ============================================================================

cat("\n=== SECTION 8: DECOMPOSITION VALIDATION SUMMARY ===\n")

# Compare R² marginal between models
nymph_orig_r2 <- nymph_orig_diag$r2_marginal
nymph_clean_r2 <- nymph_clean_diag$r2_marginal
larva_orig_r2 <- larva_orig_diag$r2_marginal
larva_clean_r2 <- larva_clean_diag$r2_marginal

cat("R² MARGINAL COMPARISON:\n")
cat("  Nymph peak - Original:", round(nymph_orig_r2, 4), "| Clean:", round(nymph_clean_r2, 4), "\n")
cat("  Larva peak - Original:", round(larva_orig_r2, 4), "| Clean:", round(larva_clean_r2, 4), "\n")

# Determine which life stage is more VPD-sensitive
orig_stronger <- ifelse(nymph_orig_r2 > larva_orig_r2, "nymph", "larva")
clean_stronger <- ifelse(nymph_clean_r2 > larva_clean_r2, "nymph", "larva")

cat("\nVPD SENSITIVITY BY LIFE STAGE:\n")
cat("  Original models: ", orig_stronger, " more VPD-sensitive (R² = ", 
    round(ifelse(orig_stronger == "nymph", nymph_orig_r2, larva_orig_r2), 4), ")\n", sep = "")
cat("  Clean models: ", clean_stronger, " more VPD-sensitive (R² = ", 
    round(ifelse(clean_stronger == "nymph", nymph_clean_r2, larva_clean_r2), 4), ")\n", sep = "")

# Check consistency of decomposition pattern
pattern_consistent <- orig_stronger == clean_stronger
cat("  Decomposition pattern consistent:", ifelse(pattern_consistent, "YES ✓", "NO ✗"), "\n")

# Check coefficient direction consistency
nymph_direction_consistent <- sign(nymph_orig_effects$vpdmin_coef) == sign(nymph_clean_effects$vpdmin_coef)
larva_direction_consistent <- sign(larva_orig_effects$vpdmin_coef) == sign(larva_clean_effects$vpdmin_coef)

cat("  Nymph VPD effect direction consistent:", ifelse(nymph_direction_consistent, "YES ✓", "NO ✗"), "\n")
cat("  Larva VPD effect direction consistent:", ifelse(larva_direction_consistent, "YES ✓", "NO ✗"), "\n")

# Diagnostic improvement
nymph_improved <- nymph_clean_diag$reliable && !nymph_orig_diag$reliable
larva_improved <- larva_clean_diag$reliable && !larva_orig_diag$reliable

cat("\nDIAGNOSTIC IMPROVEMENT:\n")
cat("  Nymph model reliability: Original =", ifelse(nymph_orig_diag$reliable, "✓", "✗"), 
    "| Clean =", ifelse(nymph_clean_diag$reliable, "✓", "✗"), 
    ifelse(nymph_improved, " (IMPROVED)", ""), "\n")
cat("  Larva model reliability: Original =", ifelse(larva_orig_diag$reliable, "✓", "✗"), 
    "| Clean =", ifelse(larva_clean_diag$reliable, "✓", "✗"), 
    ifelse(larva_improved, " (IMPROVED)", ""), "\n")

# Final recommendation for decomposition analysis
cat("\nDECOMPOSITION ANALYSIS VALIDATION:\n")
if(pattern_consistent && nymph_direction_consistent && larva_direction_consistent) {
  cat("✓ DECOMPOSITION ANALYSIS VALIDATED\n")
  cat("  → Consistent life stage VPD sensitivity patterns with and without outliers\n")
  cat("  → VPD effects on", clean_stronger, "peak timing drive PAD patterns\n")
  cat("  → Outliers do not change biological interpretation\n")
  cat("  → Safe to proceed with decomposition conclusions\n")
} else if(pattern_consistent) {
  cat("✓ DECOMPOSITION PATTERN VALIDATED\n") 
  cat("  → Same life stage drives VPD sensitivity in both datasets\n")
  cat("  → Some coefficient magnitude changes but pattern holds\n")
  cat("  → Decomposition interpretation remains valid\n")
} else {
  cat("✗ DECOMPOSITION PATTERN INCONSISTENT\n")
  cat("  → Different life stages show stronger VPD sensitivity with/without outliers\n")
  cat("  → Decomposition interpretation may be outlier-driven\n")
  cat("  → Consider robust modeling approaches\n")
}

# ============================================================================
# SECTION 9: DECOMPOSITION CONCLUSIONS
# ============================================================================

cat("\n=== SECTION 9: DECOMPOSITION CONCLUSIONS ===\n")

# Calculate effect magnitude changes
nymph_effect_change <- abs(nymph_clean_effects$vpdmin_coef - nymph_orig_effects$vpdmin_coef)
larva_effect_change <- abs(larva_clean_effects$vpdmin_coef - larva_orig_effects$vpdmin_coef)

cat("EFFECT MAGNITUDE CHANGES:\n")
cat("  Nymph VPD coefficient change:", round(nymph_effect_change, 4), "\n")
cat("  Larva VPD coefficient change:", round(larva_effect_change, 4), "\n")

# R² changes
nymph_r2_change <- abs(nymph_clean_r2 - nymph_orig_r2)
larva_r2_change <- abs(larva_clean_r2 - larva_orig_r2)

cat("  Nymph R² change:", round(nymph_r2_change, 4), "\n")
cat("  Larva R² change:", round(larva_r2_change, 4), "\n")

# Which life stage shows more sensitivity to outlier removal?
more_sensitive_to_outliers <- ifelse(nymph_effect_change > larva_effect_change, "nymph", "larva")
cat("  Life stage more sensitive to outliers:", more_sensitive_to_outliers, "\n")

cat("\nFINAL DECOMPOSITION INTERPRETATION:\n")
if(pattern_consistent) {
  cat("✓ ROBUST DECOMPOSITION FINDINGS:\n")
  cat("  → VPD min effects on PAD are consistently driven by", clean_stronger, "timing changes\n")
  cat("  → This biological pattern holds with and without outliers\n")
  cat("  → Decomposition analysis supports VPD-only model findings\n")
  
  if(clean_stronger == "nymph") {
    cat("  → VPD min primarily affects nymph questing timing\n")
    cat("  → Larva timing is less responsive to VPD min\n")
  } else {
    cat("  → VPD min primarily affects larva questing timing\n") 
    cat("  → Nymph timing is less responsive to VPD min\n")
  }
  
} else {
  cat("⚠ INCONSISTENT DECOMPOSITION FINDINGS:\n")
  cat("  → Different life stages drive VPD sensitivity with/without outliers\n")
  cat("  → Decomposition interpretation may be outlier-dependent\n")
  cat("  → Focus on VPD-only model for PAD rather than decomposition\n")
}

cat("\nCONCLUSION:\n")
cat("The decomposition analysis is", 
    ifelse(pattern_consistent, "VALIDATED", "QUESTIONABLE"), 
    "- VPD effects on PAD are", 
    ifelse(pattern_consistent, paste("consistently driven by", clean_stronger, "timing"), "inconsistently driven by different life stages"), ".\n")

cat("\n=== SENSITIVITY ANALYSIS COMPLETE ===\n")