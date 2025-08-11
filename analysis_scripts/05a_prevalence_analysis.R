# ============================================================================
# 05_simplified_prevalence_analysis.R
# Simplified prevalence analysis with exactly 5 individual models
# Models: PAD, log_nymph_count, log_larva_count, mean_host_density, vpd_min
# Response: transmission_log (using lm, not GLM)
# Time lag: year n prevalence predicted by year n-1 predictors
# ============================================================================

cat("=== SIMPLIFIED PREVALENCE ANALYSIS ===\n")
cat("Exactly 5 individual models: PAD, log_nymph, log_larva, host_density, vpd_min\n")
cat("Response: transmission_log (linear regression with log-transformed response)\n")
cat("Time lag: prevalence(n) ~ predictors(n-1)\n\n")

# Load required libraries
library(dplyr)
library(lmtest)
library(car)

# ============================================================================
# SECTION 1: LOAD AND PREPARE DATA
# ============================================================================

cat("=== SECTION 1: DATA PREPARATION ===\n")

# Load prevalence data
if (file.exists("data/processed/yearly_transmission.RDS")) {
  prevalence_data <- readRDS("data/processed/yearly_transmission.RDS")
  cat("✓ Loaded yearly_transmission.RDS\n")
} else {
  stop("Error: yearly_transmission.RDS not found. Please run 03c_merge_prevalence_data.R first.")
}

cat("Initial dataset:", nrow(prevalence_data), "observations\n")

# Check required variables
required_vars <- c("transmission_log", "PAD", "nymph_max", "larva_max", 
                   "mean_density_ha", "mean_annual_vpdmin", "nlcdClass", "lat", "long")
missing_vars <- required_vars[!required_vars %in% names(prevalence_data)]
if(length(missing_vars) > 0) {
  warning("Missing variables (will be excluded from model_data): ", paste(missing_vars, collapse = ", "))
  cat("  Note: Missing variables are acceptable for core modeling but may affect plotting\n")
}

# Prepare modeling dataset
model_data <- prevalence_data %>%
  # Select required variables (including habitat and coordinates for plotting)
  dplyr::select(
    site_year, siteID, year,
    transmission_log,
    PAD,
    nymph_max,
    larva_max,
    mean_density_ha,
    mean_annual_vpdmin,
    nlcdClass,  # For habitat analysis
    lat, long   # For coordinate-based analysis and latitude ordering
  ) %>%
  # Remove rows with missing values in CORE modeling variables only
  # (nlcdClass, lat, long can be missing for some observations)
  filter(!is.na(transmission_log), !is.na(PAD), !is.na(nymph_max), !is.na(larva_max),
         !is.na(mean_density_ha), !is.na(mean_annual_vpdmin)) %>%
  # Create log-transformed tick count variables
  mutate(
    log_nymph_count = log(nymph_max + 1),  # Adding 1 to handle zeros
    log_larva_count = log(larva_max + 1)   # Adding 1 to handle zeros
  ) %>%
  # Standardize all predictors (z-scores)
  mutate(
    PAD_std = as.numeric(scale(PAD)),
    log_nymph_count_std = as.numeric(scale(log_nymph_count)),
    log_larva_count_std = as.numeric(scale(log_larva_count)),
    mean_density_ha_std = as.numeric(scale(mean_density_ha)),
    mean_annual_vpdmin_std = as.numeric(scale(mean_annual_vpdmin))
  )

cat("Complete cases for modeling:", nrow(model_data), "\n")
cat("Sites represented:", length(unique(model_data$siteID)), "\n")
cat("Years represented:", paste(range(model_data$year), collapse = "-"), "\n")
cat("✓ All variables standardized\n")

# Check availability of additional variables for plotting
n_with_habitat <- sum(!is.na(model_data$nlcdClass))
n_with_coords <- sum(!is.na(model_data$lat) & !is.na(model_data$long))
cat("Additional variables for plotting:\n")
cat("  Observations with habitat data (nlcdClass):", n_with_habitat, "\n")
cat("  Observations with coordinates (lat, long):", n_with_coords, "\n\n")

# ============================================================================
# SECTION 2: BUILD 5 INDIVIDUAL MODELS
# ============================================================================

cat("=== SECTION 2: BUILDING 5 INDIVIDUAL MODELS ===\n")

# Initialize model storage
prevalence_models <- list()

# Build exactly 5 individual predictor models
cat("Building 5 individual predictor models...\n")
prevalence_models[["pad_only"]] <- lm(transmission_log ~ PAD_std, data = model_data)
prevalence_models[["nymph_only"]] <- lm(transmission_log ~ log_nymph_count_std, data = model_data)
prevalence_models[["larva_only"]] <- lm(transmission_log ~ log_larva_count_std, data = model_data)
prevalence_models[["host_only"]] <- lm(transmission_log ~ mean_density_ha_std, data = model_data)
prevalence_models[["vpd_only"]] <- lm(transmission_log ~ mean_annual_vpdmin_std, data = model_data)

cat("✓ Built exactly", length(prevalence_models), "individual models\n\n")

# ============================================================================
# SECTION 3: LINEAR MODEL DIAGNOSTICS
# ============================================================================

cat("=== SECTION 3: MODEL DIAGNOSTICS ===\n")

# Function for linear model diagnostics
run_lm_diagnostics <- function(model, model_name) {
  cat("Diagnosing:", model_name, "\n")
  
  # Basic model info
  n_obs <- nobs(model)
  n_predictors <- length(coef(model)) - 1  # Exclude intercept
  
  # 1. Normality test (Shapiro-Wilk on residuals)
  normality_p <- NA
  if (n_obs >= 3 && n_obs <= 5000) {
    tryCatch({
      shapiro_result <- shapiro.test(residuals(model))
      normality_p <- shapiro_result$p.value
    }, error = function(e) {
      normality_p <- NA
    })
  }
  
  # 2. Homoscedasticity test (Breusch-Pagan)
  homoscedasticity_p <- NA
  tryCatch({
    bp_result <- bptest(model)
    homoscedasticity_p <- bp_result$p.value
  }, error = function(e) {
    homoscedasticity_p <- NA
  })
  
  # 3. Autocorrelation test (Durbin-Watson)
  autocorrelation_p <- NA
  tryCatch({
    dw_result <- durbinWatsonTest(model)
    autocorrelation_p <- dw_result$p
  }, error = function(e) {
    autocorrelation_p <- NA
  })
  
  # 4. Model fit statistics
  model_summary <- summary(model)
  r_squared <- model_summary$r.squared
  adj_r_squared <- model_summary$adj.r.squared
  f_statistic <- model_summary$fstatistic[1]
  f_p_value <- pf(f_statistic, n_predictors, n_obs - n_predictors - 1, lower.tail = FALSE)
  
  # 5. Extract coefficient p-values (excluding intercept)
  coef_table <- summary(model)$coefficients
  if(nrow(coef_table) > 1) {
    min_coef_p <- min(coef_table[-1, "Pr(>|t|)"], na.rm = TRUE)
  } else {
    min_coef_p <- NA
  }
  
  # Determine reliability (pass all diagnostic tests)
  reliable <- (is.na(normality_p) || normality_p >= 0.05) &&
             (is.na(homoscedasticity_p) || homoscedasticity_p >= 0.05) &&
             (is.na(autocorrelation_p) || autocorrelation_p >= 0.05)
  
  cat("  Normality p-value:", round(normality_p, 4), ifelse(is.na(normality_p) || normality_p >= 0.05, "✓", "✗"), "\n")
  cat("  Homoscedasticity p-value:", round(homoscedasticity_p, 4), ifelse(is.na(homoscedasticity_p) || homoscedasticity_p >= 0.05, "✓", "✗"), "\n")
  cat("  Autocorrelation p-value:", round(autocorrelation_p, 4), ifelse(is.na(autocorrelation_p) || autocorrelation_p >= 0.05, "✓", "✗"), "\n")
  cat("  R²:", round(r_squared, 4), "\n")
  cat("  F p-value:", round(f_p_value, 4), "\n")
  cat("  Reliable:", ifelse(reliable, "✓", "✗"), "\n\n")
  
  return(data.frame(
    model_name = model_name,
    n_obs = n_obs,
    n_predictors = n_predictors,
    normality_p = normality_p,
    homoscedasticity_p = homoscedasticity_p,
    autocorrelation_p = autocorrelation_p,
    r_squared = r_squared,
    adj_r_squared = adj_r_squared,
    f_p_value = f_p_value,
    min_coef_p = min_coef_p,
    reliable = reliable,
    stringsAsFactors = FALSE
  ))
}

# Run diagnostics on all models
diagnostic_results <- data.frame()

for(i in 1:length(prevalence_models)) {
  model_name <- names(prevalence_models)[i]
  model <- prevalence_models[[i]]
  
  result <- run_lm_diagnostics(model, model_name)
  diagnostic_results <- rbind(diagnostic_results, result)
}

# ============================================================================
# SECTION 4: MODEL RANKING AND SELECTION
# ============================================================================

cat("=== SECTION 4: MODEL RANKING ===\n")

# Calculate AIC for all models and rank
model_comparison <- diagnostic_results %>%
  mutate(
    AIC = sapply(names(prevalence_models)[match(model_name, names(prevalence_models))], 
                 function(x) AIC(prevalence_models[[x]])),
    delta_AIC = AIC - min(AIC, na.rm = TRUE)
  ) %>%
  arrange(AIC)

cat("All models ranked by AIC:\n")
print(model_comparison[c("model_name", "AIC", "delta_AIC", "r_squared", "reliable")])

# Identify reliable models
reliable_models <- model_comparison %>%
  filter(reliable == TRUE)

if(nrow(reliable_models) > 0) {
  cat("\nReliable models (passed all diagnostic tests):\n")
  print(reliable_models[c("model_name", "AIC", "delta_AIC", "r_squared")])
  
  best_model_name <- reliable_models$model_name[1]
  cat("\n✓ Best reliable model:", best_model_name, "\n")
} else {
  cat("\n⚠ No models passed all diagnostic tests\n")
  cat("Selecting best model by AIC despite diagnostic failures\n")
  best_model_name <- model_comparison$model_name[1]
  cat("✓ Best model (by AIC):", best_model_name, "\n")
}

best_model <- prevalence_models[[best_model_name]]

cat("  AIC:", round(model_comparison$AIC[model_comparison$model_name == best_model_name], 2), "\n")
cat("  R²:", round(model_comparison$r_squared[model_comparison$model_name == best_model_name], 4), "\n\n")

# ============================================================================
# SECTION 5: FINAL MODEL SUMMARY
# ============================================================================

cat("=== SECTION 5: FINAL MODEL SUMMARY ===\n")

cat("SIMPLIFIED PREVALENCE ANALYSIS RESULTS\n")
cat("======================================\n\n")

cat("Models tested:", length(prevalence_models), "(5 individual predictors only)\n")
cat("Individual predictors: PAD, log_nymph_count, log_larva_count, host_density, vpd_min\n")
cat("No combination models - clean, focused approach\n\n")

cat("Selected model:", best_model_name, "\n")
cat("Model formula:", deparse(formula(best_model)), "\n\n")

# Print model summary
cat("MODEL SUMMARY:\n")
print(summary(best_model))

cat("\n")
cat("DIAGNOSTIC SUMMARY:\n")
print(diagnostic_results)

# ============================================================================
# SECTION 6: SAVE RESULTS
# ============================================================================

cat("\n=== SECTION 6: SAVING RESULTS ===\n")

# Create output directory
if (!dir.exists("output")) {
  dir.create("output", recursive = TRUE)
}

# Save all models
saveRDS(prevalence_models, "output/05_simplified_prevalence_models.RDS")

# Save diagnostic results
write.csv(diagnostic_results, "output/05_simplified_prevalence_diagnostics.csv", row.names = FALSE)

# Save model comparison
write.csv(model_comparison, "output/05_simplified_prevalence_ranking.csv", row.names = FALSE)

# Save best model
saveRDS(best_model, "output/05_simplified_best_model.RDS")

# Save model data for further analysis
saveRDS(model_data, "output/05_simplified_model_data.RDS")

# Save analysis info
analysis_info <- list(
  n_observations = nrow(model_data),
  n_sites = length(unique(model_data$siteID)),
  n_with_habitat = n_with_habitat,
  n_with_coords = n_with_coords,
  models_tested = names(prevalence_models),
  best_model = best_model_name,
  approach = "Simplified prevalence analysis - 5 individual predictors only",
  predictors = c("PAD_std", "log_nymph_count_std", "log_larva_count_std", 
                "mean_density_ha_std", "mean_annual_vpdmin_std"),
  response = "transmission_log",
  time_lag = "year n prevalence predicted by year n-1 predictors",
  reliable_models = sum(diagnostic_results$reliable),
  additional_variables = c("nlcdClass", "lat", "long"),
  notes = "Includes habitat and coordinate data for downstream plotting"
)

saveRDS(analysis_info, "output/05_simplified_analysis_info.RDS")

cat("✓ Results saved to output/ directory\n")
cat("  - Models: 05_simplified_prevalence_models.RDS\n")
cat("  - Diagnostics: 05_simplified_prevalence_diagnostics.csv\n")
cat("  - Rankings: 05_simplified_prevalence_ranking.csv\n")
cat("  - Best model: 05_simplified_best_model.RDS\n")
cat("  - Model data: 05_simplified_model_data.RDS\n")
cat("  - Analysis info: 05_simplified_analysis_info.RDS\n\n")

cat("=== SIMPLIFIED PREVALENCE ANALYSIS COMPLETE ===\n")