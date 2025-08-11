# ============================================================================
# 04_simplified_phenology_analysis.R
# Simplified phenology analysis with only 3 models
# Models: mean host density, VPD min, and their interaction
# All with siteID random effects, all variables standardized
# ============================================================================

cat("=== SIMPLIFIED PHENOLOGY ANALYSIS ===\n")
cat("Building 3 models only: host density, VPD min, and interaction\n")
cat("All models with siteID random effects\n\n")

# Load required libraries
library(dplyr)
library(lme4)
library(lmerTest)
library(car)
library(performance)

# ============================================================================
# SECTION 1: LOAD AND PREPARE DATA
# ============================================================================

cat("=== SECTION 1: DATA PREPARATION ===\n")

# Load analysis data
if (!exists("analysis_data")) {
  if (file.exists("data/processed/pad_analysis_dataset.RDS")) {
    analysis_data <- readRDS("data/processed/pad_analysis_dataset.RDS")
    cat("✓ Loaded analysis_data from RDS file\n")
  } else {
    stop("Error: pad_analysis_dataset.RDS not found. Please run 03b_merge_phenology_data.R first.")
  }
}

# Check required variables
required_vars <- c("PAD", "mean_density_ha", "mean_annual_vpdmin", "siteID")
missing_vars <- required_vars[!required_vars %in% names(analysis_data)]
if(length(missing_vars) > 0) {
  stop("Error: Missing required variables: ", paste(missing_vars, collapse = ", "))
}

# Filter and prepare data
analysis_data <- analysis_data %>%
  filter(!is.na(PAD) & !is.na(siteID)) %>%
#  filter(PAD >= 0) %>%  # Remove negative PAD values
  filter(!is.na(mean_density_ha) & !is.na(mean_annual_vpdmin)) %>%
  mutate(siteID = as.factor(siteID))

cat("Data after filtering:\n")
cat("  Observations:", nrow(analysis_data), "\n")
cat("  Sites:", length(unique(analysis_data$siteID)), "\n")
cat("  Years:", length(unique(analysis_data$year)), "\n\n")

# Standardize variables (z-scores)
analysis_data <- analysis_data %>%
  mutate(
    mean_density_ha_std = as.numeric(scale(mean_density_ha)),
    mean_annual_vpdmin_std = as.numeric(scale(mean_annual_vpdmin))
  )

cat("✓ Variables standardized: mean_density_ha_std, mean_annual_vpdmin_std\n\n")

# ============================================================================
# SECTION 2: BUILD THE 3 MODELS
# ============================================================================

cat("=== SECTION 2: BUILDING 3 SIMPLIFIED MODELS ===\n")

# Store all models
phenology_models <- list()

# Model 1: Host density only
cat("Building Model 1: Host density only...\n")
phenology_models[["host_only"]] <- lmer(PAD ~ mean_density_ha_std + (1|siteID), 
                                      data = analysis_data)

# Model 2: VPD min only  
cat("Building Model 2: VPD min only...\n")
phenology_models[["vpd_only"]] <- lmer(PAD ~ mean_annual_vpdmin_std + (1|siteID), 
                                     data = analysis_data)

# Model 3: Host density * VPD min interaction
cat("Building Model 3: Host density * VPD min interaction...\n")
phenology_models[["host_vpd_interaction"]] <- lmer(PAD ~ mean_density_ha_std * mean_annual_vpdmin_std + (1|siteID), 
                                                 data = analysis_data)

cat("✓ Built 3 models successfully\n\n")

# ============================================================================
# SECTION 3: MIXED MODEL DIAGNOSTICS
# ============================================================================

cat("=== SECTION 3: MIXED MODEL DIAGNOSTICS ===\n")

# Function for mixed model diagnostics
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

for(i in 1:length(phenology_models)) {
  model_name <- names(phenology_models)[i]
  model <- phenology_models[[i]]
  
  result <- run_mixed_diagnostics(model, model_name)
  diagnostic_results <- rbind(diagnostic_results, result)
}

# ============================================================================
# SECTION 4: MODEL RANKING AND SELECTION
# ============================================================================

cat("=== SECTION 4: MODEL RANKING ===\n")

# Calculate AIC for reliable models
model_comparison <- diagnostic_results %>%
  filter(reliable == TRUE) %>%
  mutate(
    AIC = sapply(names(phenology_models)[match(model_name, names(phenology_models))], 
                 function(x) AIC(phenology_models[[x]])),
    delta_AIC = AIC - min(AIC, na.rm = TRUE)
  ) %>%
  arrange(AIC)

if(nrow(model_comparison) > 0) {
  cat("Reliable models ranked by AIC:\n")
  print(model_comparison[c("model_name", "AIC", "delta_AIC", "r2_marginal", "r2_conditional")])
  
  # Select best model
  best_model_name <- model_comparison$model_name[1]
  best_model <- phenology_models[[best_model_name]]
  
  cat("\n✓ Best model:", best_model_name, "\n")
  cat("  AIC:", round(model_comparison$AIC[1], 2), "\n")
  cat("  R² marginal:", round(model_comparison$r2_marginal[1], 4), "\n")
  cat("  R² conditional:", round(model_comparison$r2_conditional[1], 4), "\n\n")
  
} else {
  cat("⚠ No reliable models found - all models failed diagnostic criteria\n")
  cat("Proceeding with model comparison despite diagnostic failures\n\n")
  
  # If no reliable models, rank all models by AIC
  model_comparison <- diagnostic_results %>%
    mutate(
      AIC = sapply(names(phenology_models)[match(model_name, names(phenology_models))], 
                   function(x) AIC(phenology_models[[x]])),
      delta_AIC = AIC - min(AIC, na.rm = TRUE)
    ) %>%
    arrange(AIC)
  
  best_model_name <- model_comparison$model_name[1]
  best_model <- phenology_models[[best_model_name]]
  
  cat("All models ranked by AIC (diagnostic failures noted):\n")
  print(model_comparison[c("model_name", "AIC", "delta_AIC", "reliable")])
}

# ============================================================================
# SECTION 5: FINAL MODEL SUMMARY
# ============================================================================

cat("\n=== SECTION 5: FINAL MODEL SUMMARY ===\n")

cat("SIMPLIFIED PHENOLOGY ANALYSIS RESULTS\n")
cat("=====================================\n\n")

cat("Models tested:\n")
cat("1. PAD ~ mean_density_ha_std + (1|siteID)\n")
cat("2. PAD ~ mean_annual_vpdmin_std + (1|siteID)\n") 
cat("3. PAD ~ mean_density_ha_std * mean_annual_vpdmin_std + (1|siteID)\n\n")

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
saveRDS(phenology_models, "output/04_simplified_phenology_models.RDS")

# Save diagnostic results
write.csv(diagnostic_results, "output/04_simplified_phenology_diagnostics.csv", row.names = FALSE)

# Save model comparison
write.csv(model_comparison, "output/04_simplified_phenology_ranking.csv", row.names = FALSE)

# Save best model
saveRDS(best_model, "output/04_simplified_best_model.RDS")

# Save analysis info
analysis_info <- list(
  n_observations = nrow(analysis_data),
  n_sites = length(unique(analysis_data$siteID)),
  models_tested = names(phenology_models),
  best_model = best_model_name,
  approach = "Simplified phenology analysis with 3 mixed effects models",
  variables = c("mean_density_ha_std", "mean_annual_vpdmin_std"),
  response = "PAD",
  random_effects = "siteID"
)

saveRDS(analysis_info, "output/04_simplified_analysis_info.RDS")

cat("✓ Results saved to output/ directory\n")
cat("  - Models: 04_simplified_phenology_models.RDS\n")
cat("  - Diagnostics: 04_simplified_phenology_diagnostics.csv\n")
cat("  - Rankings: 04_simplified_phenology_ranking.csv\n")
cat("  - Best model: 04_simplified_best_model.RDS\n")
cat("  - Analysis info: 04_simplified_analysis_info.RDS\n\n")

cat("=== SIMPLIFIED PHENOLOGY ANALYSIS COMPLETE ===\n")