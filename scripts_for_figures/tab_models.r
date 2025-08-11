# ============================================================================
# Making Tables and Plots for the Results Section
# ============================================================================
library(sjPlot)

# phenology analysis models
summary(phenology_models[["host_only"]])           # Host density only
summary(phenology_models[["vpd_only"]])            # VPD max only  
summary(phenology_models[["host_vpd_interaction"]])  # Full interaction model

# Or see all model names first
names(phenology_models)

# phenology final model
summary(phenology_models[["vpd_only"]])
performance::r2(phenology_models[["vpd_only"]])

tab_model(phenology_models[["vpd_only"]], 
          dv.labels = "Peak Activity Difference (PAD)")



# prevalence analysis models
summary(prevalence_models[["pad_only"]])     # PAD predictor
summary(prevalence_models[["nymph_only"]])   # Nymph count predictor
summary(prevalence_models[["host_only"]])    # Host density predictor
summary(prevalence_models[["vpd_only"]])     # VPD max predictor

# Or see all model names
names(prevalence_models)

# tab model
tab_model(prevalence_models[["pad_only"]],
          prevalence_models[["host_only"]],
          prevalence_models[["vpd_only"]]
          )


# ============================================================================
# Effect size tabling for VPD-only PAD analysis
# Extract effect size metrics for the reliable vpd_only model
# ============================================================================

# Get the vpd_only model (the reliable model)
vpd_model <- phenology_models[["vpd_only"]]
model_summary <- summary(vpd_model)
coef_table <- model_summary$coefficients

# Get SD of PAD to standardize coefficients
pad_sd <- sd(analysis_data$PAD, na.rm = TRUE)

# Calculate R² for the model
model_r2 <- performance::r2(vpd_model)$R2_marginal

# Calculate Cohen's f² for the model
cohens_f2 <- model_r2 / (1 - model_r2)

# Calculate standardized coefficient (SD change in PAD per 1 SD predictor change)
std_beta_vpd <- coef_table["mean_annual_vpdmin_std", "Estimate"] / pad_sd

# Calculate standardized CI bounds
vpd_ci_lower <- (coef_table["mean_annual_vpdmin_std", "Estimate"] - 1.96 * coef_table["mean_annual_vpdmin_std", "Std. Error"]) / pad_sd
vpd_ci_upper <- (coef_table["mean_annual_vpdmin_std", "Estimate"] + 1.96 * coef_table["mean_annual_vpdmin_std", "Std. Error"]) / pad_sd

# Create table with standardized coefficients for VPD min only
effect_table_dash <- data.frame(
  Measures = c("β (standardized)", "95% CI", "R²", "Cohen's f²"),
  `mean annual vpdmin std` = c(
    round(std_beta_vpd, 3),
    paste0(round(vpd_ci_lower, 3), " – ", round(vpd_ci_upper, 3)),
    round(model_r2, 3),
    round(cohens_f2, 3)
  ),
  check.names = FALSE
)

# Create table
tab_df(effect_table_dash, 
       title = NULL,
       col.header = c("Measures", "mean annual vpdmin std"))










# ----- old one. so don't use!!!! --------------
# effect size tabling for PAD analysis
# Extract data for the two significant terms
vpd_data <- phenology_decomposition[phenology_decomposition$Effect == "VPD Main", ]
host_data <- phenology_decomposition[phenology_decomposition$Effect == "Density Main", ]

# Get coefficients from model summary
model_summary <- summary(phenology_models[["host_vpd_interaction"]])
coef_table <- model_summary$coefficients

# Get SD of PAD to standardize coefficients
pad_sd <- sd(analysis_data$PAD, na.rm = TRUE)

# Calculate standardized coefficients (SD change in PAD per 1 SD predictor change)
std_beta_vpd <- vpd_data$Std_Beta / pad_sd
std_beta_host <- host_data$Std_Beta / pad_sd

# Calculate standardized CI bounds
vpd_ci_lower <- (coef_table["mean_annual_vpdmax_std", "Estimate"] - 1.96 * coef_table["mean_annual_vpdmax_std", "Std. Error"]) / pad_sd
vpd_ci_upper <- (coef_table["mean_annual_vpdmax_std", "Estimate"] + 1.96 * coef_table["mean_annual_vpdmax_std", "Std. Error"]) / pad_sd
host_ci_lower <- (coef_table["mean_density_ha_std", "Estimate"] - 1.96 * coef_table["mean_density_ha_std", "Std. Error"]) / pad_sd
host_ci_upper <- (coef_table["mean_density_ha_std", "Estimate"] + 1.96 * coef_table["mean_density_ha_std", "Std. Error"]) / pad_sd

# Create table with standardized coefficients
effect_table_dash <- data.frame(
  Measures = c("β (standardized)", "95% CI", "ΔR²", "Cohen's f²"),
  `mean annual vpdmax std` = c(
    round(std_beta_vpd, 3),
    paste0(round(vpd_ci_lower, 3), " – ", round(vpd_ci_upper, 3)),
    round(vpd_data$Delta_R2, 3),
    round(vpd_data$Cohens_f2, 3)
  ),
  `mean density ha std` = c(
    round(std_beta_host, 3),
    paste0(round(host_ci_lower, 3), " – ", round(host_ci_upper, 3)),
    round(host_data$Delta_R2, 3),
    round(host_data$Cohens_f2, 3)
  ),
  check.names = FALSE
)

# Create table
tab_df(effect_table_dash, 
       title = NULL,
       col.header = c("Measures", "mean annual vpdmax std", "mean density ha std"))



