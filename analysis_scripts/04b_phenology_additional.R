# ============================================================================
# 04b_peak_day_decomposition.R
# Peak day decomposition analysis - what drives PAD?
# Models: VPD effects on nymph and larva peak days
# Goal: Understand which life stage drives VPD sensitivity in PAD
# ============================================================================

cat("=== PEAK DAY DECOMPOSITION ANALYSIS ===\n")
cat("Decomposing PAD into nymph and larva peak day components\n")
cat("Using VPD-only model: VPD → each peak day\n\n")

# Load required libraries
library(dplyr)
library(lme4)
library(lmerTest)
library(car)
library(performance)

# ============================================================================
# SECTION 1: VERIFY DATA IS LOADED
# ============================================================================

cat("=== SECTION 1: DATA VERIFICATION ===\n")

# Check that analysis_data is already loaded from main phenology script
if (!exists("analysis_data")) {
  stop("Error: analysis_data not found. Please run 04_simplified_phenology_analysis.R first.")
}

# Check required variables (including peak day variables)
required_vars <- c("nymph_peak_day", "larva_peak_day", "mean_annual_vpdmin_std", "siteID")
missing_vars <- required_vars[!required_vars %in% names(analysis_data)]
if(length(missing_vars) > 0) {
  stop("Error: Missing required variables in analysis_data: ", paste(missing_vars, collapse = ", "))
}

cat("✓ Using same analysis_data from main phenology script\n")
cat("Data specifications:\n")
cat("  Observations:", nrow(analysis_data), "\n")
cat("  Sites:", length(unique(analysis_data$siteID)), "\n")
cat("  Years:", length(unique(analysis_data$year)), "\n")

# Check peak day variable availability
n_with_peak_days <- sum(!is.na(analysis_data$nymph_peak_day) & !is.na(analysis_data$larva_peak_day))
cat("  Observations with both peak days:", n_with_peak_days, "\n")

if(n_with_peak_days == 0) {
  stop("Error: No observations have both nymph_peak_day and larva_peak_day data")
}

# Check peak day distributions
cat("\nPeak day variable summaries:\n")
cat("  Nymph peak day: mean =", round(mean(analysis_data$nymph_peak_day, na.rm = TRUE), 1), 
    ", range =", min(analysis_data$nymph_peak_day, na.rm = TRUE), "-", max(analysis_data$nymph_peak_day, na.rm = TRUE), "\n")
cat("  Larva peak day: mean =", round(mean(analysis_data$larva_peak_day, na.rm = TRUE), 1),
    ", range =", min(analysis_data$larva_peak_day, na.rm = TRUE), "-", max(analysis_data$larva_peak_day, na.rm = TRUE), "\n")

# Calculate PAD from components to verify
calculated_PAD <- analysis_data$larva_peak_day - analysis_data$nymph_peak_day
cat("  PAD (calculated): mean =", round(mean(calculated_PAD, na.rm = TRUE), 1), "\n")
if("PAD" %in% names(analysis_data)) {
  cat("  PAD (original): mean =", round(mean(analysis_data$PAD, na.rm = TRUE), 1), "\n")
}
cat("\n")

# ============================================================================
# SECTION 2: BUILD PEAK DAY MODELS
# ============================================================================

cat("=== SECTION 2: BUILDING PEAK DAY DECOMPOSITION MODELS ===\n")

# Store all models
peak_models <- list()

# Model 1: Nymph peak day ~ VPD only
cat("Building Model 1: Nymph peak day ~ VPD only...\n")
peak_models[["nymph_peak"]] <- lmer(nymph_peak_day ~ mean_annual_vpdmin_std + (1|siteID), 
                                   data = analysis_data)

# Model 2: Larva peak day ~ VPD only  
cat("Building Model 2: Larva peak day ~ VPD only...\n")
peak_models[["larva_peak"]] <- lmer(larva_peak_day ~ mean_annual_vpdmin_std + (1|siteID), 
                                   data = analysis_data)

cat("✓ Built 2 peak day models successfully\n\n")

# ============================================================================
# SECTION 3: MIXED MODEL DIAGNOSTICS
# ============================================================================

cat("=== SECTION 3: MIXED MODEL DIAGNOSTICS ===\n")

# Function for mixed model diagnostics (same as main script)
run_mixed_diagnostics <- function(model, model_name) {
  cat("Diagnosing:", model_name, "\n")
  
  # Basic model info
  n_obs <- nobs(model)
  
  # 1. Convergence check
  converged <- model@optinfo$conv$opt == 0
  
  # 2. Singularity check  
  singular <- isSingular(model)
  
  # 3. Normality of residuals (Shapiro-Wilk)
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
  
  # 4. Model fit statistics
  model_summary <- summary(model)
  r2_marginal <- NA
  r2_conditional <- NA
  
  tryCatch({
    r2_results <- performance::r2(model)
    r2_marginal <- r2_results$R2_marginal
    r2_conditional <- r2_results$R2_conditional
  }, error = function(e) {
    # If performance package fails, skip R²
  })
  
  # Extract p-values for fixed effects
  coef_table <- coef(model_summary)
  fixed_effects_p <- if(nrow(coef_table) > 1) coef_table[-1, "Pr(>|t|)"] else NA
  
  # Determine reliability
  reliable <- converged && !singular && 
             (is.na(normality_p) || normality_p >= 0.05)
  
  cat("  Convergence:", ifelse(converged, "✓", "✗"), "\n")
  cat("  Singularity:", ifelse(singular, "✗", "✓"), "\n")
  cat("  Normality p-value:", round(normality_p, 4), ifelse(is.na(normality_p) || normality_p >= 0.05, "✓", "✗"), "\n")
  cat("  R² marginal:", round(r2_marginal, 4), "\n")
  cat("  R² conditional:", round(r2_conditional, 4), "\n")
  cat("  Reliable:", ifelse(reliable, "✓", "✗"), "\n\n")
  
  return(data.frame(
    model_name = model_name,
    n_obs = n_obs,
    converged = converged,
    singular = singular,
    normality_p = normality_p,
    r2_marginal = r2_marginal,
    r2_conditional = r2_conditional,
    reliable = reliable,
    stringsAsFactors = FALSE
  ))
}

# Run diagnostics on all models
diagnostic_results <- data.frame()

for(i in 1:length(peak_models)) {
  model_name <- names(peak_models)[i]
  model <- peak_models[[i]]
  
  result <- run_mixed_diagnostics(model, model_name)
  diagnostic_results <- rbind(diagnostic_results, result)
}

# ============================================================================
# SECTION 4: MODEL COMPARISON AND INTERPRETATION
# ============================================================================

cat("=== SECTION 4: MODEL COMPARISON AND INTERPRETATION ===\n")

# Calculate AIC for model comparison
model_comparison <- diagnostic_results %>%
  mutate(
    AIC = sapply(names(peak_models)[match(model_name, names(peak_models))], 
                 function(x) AIC(peak_models[[x]])),
    delta_AIC = AIC - min(AIC, na.rm = TRUE)
  ) %>%
  arrange(AIC)

cat("Peak day models comparison:\n")
print(model_comparison[c("model_name", "AIC", "delta_AIC", "r2_marginal", "r2_conditional", "reliable")])

# Detailed coefficient comparison
cat("\n=== COEFFICIENT COMPARISON ===\n")

for(i in 1:length(peak_models)) {
  model_name <- names(peak_models)[i]
  model <- peak_models[[i]]
  
  cat("\n", toupper(model_name), "MODEL:\n")
  cat("Formula:", deparse(formula(model)), "\n")
  
  # Extract and display key coefficients
  coef_table <- summary(model)$coefficients
  
  if("mean_annual_vpdmin_std" %in% rownames(coef_table)) {
    vpd_coef <- coef_table["mean_annual_vpdmin_std", "Estimate"]
    vpd_p <- coef_table["mean_annual_vpdmin_std", "Pr(>|t|)"]
    cat("  VPD effect:", round(vpd_coef, 3), "(p =", round(vpd_p, 4), ")\n")
  }
}

# ============================================================================
# SECTION 5: BIOLOGICAL INTERPRETATION
# ============================================================================

cat("\n=== SECTION 5: BIOLOGICAL INTERPRETATION ===\n")

# Compare model fit
nymph_r2 <- diagnostic_results$r2_marginal[diagnostic_results$model_name == "nymph_peak"]
larva_r2 <- diagnostic_results$r2_marginal[diagnostic_results$model_name == "larva_peak"]

cat("EFFECT SIZE COMPARISON:\n")
cat("  Nymph peak day model R² marginal:", round(nymph_r2, 4), "\n")
cat("  Larva peak day model R² marginal:", round(larva_r2, 4), "\n")

stronger_effect <- ifelse(nymph_r2 > larva_r2, "nymph", "larva")
cat("  Stronger VPD effect on:", stronger_effect, "peak timing\n\n")

# Interpretation of PAD drivers
cat("PAD DECOMPOSITION INSIGHTS:\n")
cat("→ PAD = larva_peak_day - nymph_peak_day\n")
cat("→ VPD effects on PAD are driven by changes in:", stronger_effect, "timing\n")

if(stronger_effect == "nymph") {
  cat("→ VPD primarily affects nymph questing timing\n")
  cat("→ Larva timing is less responsive to VPD\n")
} else {
  cat("→ VPD primarily affects larva questing timing\n") 
  cat("→ Nymph timing is less responsive to VPD\n")
}

# Check if both models are reliable
both_reliable <- all(diagnostic_results$reliable)
cat("→ Both models reliable:", ifelse(both_reliable, "YES", "NO"), "\n")

if(!both_reliable) {
  unreliable_models <- diagnostic_results$model_name[!diagnostic_results$reliable]
  cat("→ Unreliable models:", paste(unreliable_models, collapse = ", "), "\n")
  cat("→ Interpret results with caution\n")
}

# ============================================================================
# SECTION 6: FINAL MODEL SUMMARIES
# ============================================================================

cat("\n=== SECTION 6: DETAILED MODEL SUMMARIES ===\n")

for(i in 1:length(peak_models)) {
  model_name <- names(peak_models)[i]
  model <- peak_models[[i]]
  
  cat("\n", rep("=", 50), "\n")
  cat(toupper(model_name), "MODEL SUMMARY\n")
  cat(rep("=", 50), "\n")
  print(summary(model))
}

cat("\n")
cat("DIAGNOSTIC SUMMARY:\n")
print(diagnostic_results)

# ============================================================================
# SECTION 7: SAVE RESULTS
# ============================================================================

cat("\n=== SECTION 7: SAVING RESULTS ===\n")

# Create output directory
if (!dir.exists("output")) {
  dir.create("output", recursive = TRUE)
}

# Save all models
saveRDS(peak_models, "output/04b_peak_day_models.RDS")

# Save diagnostic results
write.csv(diagnostic_results, "output/04b_peak_day_diagnostics.csv", row.names = FALSE)

# Save model comparison
write.csv(model_comparison, "output/04b_peak_day_comparison.csv", row.names = FALSE)

# Save analysis info
peak_analysis_info <- list(
  n_observations = nrow(analysis_data),
  n_sites = length(unique(analysis_data$siteID)),
  n_with_peak_days = n_with_peak_days,
  models_tested = names(peak_models),
  approach = "Peak day decomposition analysis - VPD only models",
  variables = c("mean_annual_vpdmin_std"),
  responses = c("nymph_peak_day", "larva_peak_day"),
  random_effects = "siteID",
  stronger_effect_on = stronger_effect,
  both_models_reliable = both_reliable,
  biological_insight = paste("VPD effects on PAD primarily driven by", stronger_effect, "timing changes"),
  note = "Uses VPD-only models to match best model from main phenology analysis"
)

saveRDS(peak_analysis_info, "output/04b_peak_day_analysis_info.RDS")

cat("✓ Results saved to output/ directory\n")
cat("  - Models: 04b_peak_day_models.RDS\n")
cat("  - Diagnostics: 04b_peak_day_diagnostics.csv\n")
cat("  - Comparison: 04b_peak_day_comparison.csv\n")
cat("  - Analysis info: 04b_peak_day_analysis_info.RDS\n\n")

cat("=== PEAK DAY DECOMPOSITION ANALYSIS COMPLETE ===\n")
cat("\nKEY FINDINGS:\n")
cat("→ VPD effects on PAD are primarily driven by:", stronger_effect, "timing\n")
cat("→ This reveals the biological mechanism behind VPD-only model results\n")
cat("→ Use these insights to interpret main PAD analysis results\n")