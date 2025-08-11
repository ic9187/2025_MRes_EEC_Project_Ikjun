# ============================================================================
# PAD VPD-Only Effect Visualization - UPDATED VERSION
# Using the VPD-ONLY MODEL (the reliable model)
# Shows VPD min effect on PAD with sites colored
# ============================================================================

library(ggplot2)
library(dplyr)
library(lme4)
library(ggeffects)

# Ensure data is loaded (assuming analysis_data and phenology_models exist)
if (!exists("analysis_data") || !exists("phenology_models")) {
  stop("Please run phenology analysis scripts first to load analysis_data and phenology_models")
}

# Use the VPD-ONLY MODEL (the reliable model!)
vpd_model <- phenology_models[["vpd_only"]]

cat("=== USING VPD-ONLY MODEL FOR MAIN EFFECT PLOT ===\n")
cat("Model formula:", deparse(formula(vpd_model)), "\n")

# Generate predictions using ggeffects
cat("\nGenerating predictions from vpd-only model...\n")

# VPD main effect (standardized)
vpd_predictions_std <- ggpredict(vpd_model,
                               terms = "mean_annual_vpdmin_std [all]")

vpd_pred_data_std <- data.frame(
  mean_annual_vpdmin_std = vpd_predictions_std$x,
  predicted = vpd_predictions_std$predicted,
  ci_lower = vpd_predictions_std$conf.low,
  ci_upper = vpd_predictions_std$conf.high
)

# For non-standardized plot, we need to convert predictions back to original scale
# Create non-standardized VPD sequence
vpd_min_range <- range(analysis_data$mean_annual_vpdmin, na.rm = TRUE)
vpd_min_seq <- seq(vpd_min_range[1], vpd_min_range[2], length.out = 100)

# Standardize the sequence for prediction
vpd_min_mean <- mean(analysis_data$mean_annual_vpdmin, na.rm = TRUE)
vpd_min_sd <- sd(analysis_data$mean_annual_vpdmin, na.rm = TRUE)
vpd_min_seq_std <- (vpd_min_seq - vpd_min_mean) / vpd_min_sd

# Generate predictions with confidence intervals for raw scale
new_data_raw <- data.frame(mean_annual_vpdmin_std = vpd_min_seq_std)
predictions_raw <- predict(vpd_model, newdata = new_data_raw, re.form = NA)

# Calculate confidence intervals manually
# Get model matrix and variance-covariance matrix
X <- model.matrix(~ mean_annual_vpdmin_std, data = new_data_raw)
vcov_matrix <- vcov(vpd_model)
se_pred <- sqrt(diag(X %*% vcov_matrix %*% t(X)))

vpd_pred_data_raw <- data.frame(
  mean_annual_vpdmin = vpd_min_seq,
  predicted = predictions_raw,
  ci_lower = predictions_raw - 1.96 * se_pred,
  ci_upper = predictions_raw + 1.96 * se_pred
)

cat("✓ Predictions generated for both standardized and raw VPD min\n")

# Plot 1: PAD vs VPD Min Standardized (from vpd-only model)
p_vpd_std <- ggplot(analysis_data, aes(x = mean_annual_vpdmin_std, y = PAD)) +
  geom_ribbon(data = vpd_pred_data_std,
              aes(x = mean_annual_vpdmin_std, ymin = ci_lower, ymax = ci_upper),
              alpha = 0.3, fill = "grey40", inherit.aes = FALSE) +
  geom_point(aes(color = siteID), alpha = 0.7, size = 2) +
  geom_line(data = vpd_pred_data_std, 
            aes(x = mean_annual_vpdmin_std, y = predicted),
            color = "black", size = 1, linetype = "solid", inherit.aes = FALSE) +
  labs(
    x = "Vapor Pressure Deficit Min (standardized)",
    y = "Peak Activity Difference (PAD)",
    title = "Climate Stress Effect on Tick Phenology", 
    subtitle = "From vpd-only model: reliable single-predictor model (standardized)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray50"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  scale_color_viridis_d(alpha = 0.8)

# Plot 2: PAD vs VPD Min Raw (from vpd-only model)
p_vpd_raw <- ggplot(analysis_data, aes(x = mean_annual_vpdmin, y = PAD)) +
  geom_ribbon(data = vpd_pred_data_raw,
              aes(x = mean_annual_vpdmin, ymin = ci_lower, ymax = ci_upper),
              alpha = 0.3, fill = "grey40", inherit.aes = FALSE) +
  geom_point(aes(color = siteID), alpha = 0.7, size = 2) +
  geom_line(data = vpd_pred_data_raw, 
            aes(x = mean_annual_vpdmin, y = predicted),
            color = "black", size = 1, linetype = "solid", inherit.aes = FALSE) +
  labs(
    x = "Vapor Pressure Deficit Min (kPa)",
    y = "Peak Activity Difference (PAD)",
    title = "Climate Stress Effect on Tick Phenology", 
    subtitle = "From vpd-only model: reliable single-predictor model (raw values)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray50"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  scale_color_viridis_d(alpha = 0.8)

# Display plots
cat("\n=== VPD-ONLY EFFECT PLOTS ===\n")
print("Plot 1: VPD Min Effect - Standardized (from vpd-only model)")
print(p_vpd_std)

print("Plot 2: VPD Min Effect - Raw Values (from vpd-only model)")
print(p_vpd_raw)

# Print model information for context
cat("\n=== MODEL SUMMARY ===\n")
model_summary <- summary(vpd_model)
cat("VPD min effect: β =", round(model_summary$coefficients["mean_annual_vpdmin_std", "Estimate"], 2),
    ", p =", round(model_summary$coefficients["mean_annual_vpdmin_std", "Pr(>|t|)"], 3), "\n")

# Get R² for context
r2_results <- performance::r2(vpd_model)
cat("R² marginal:", round(r2_results$R2_marginal, 3), "\n")
cat("R² conditional:", round(r2_results$R2_conditional, 3), "\n")

# Save plots
if (!dir.exists("figures")) {
  dir.create("figures", recursive = TRUE)
}

ggsave("figures/PAD_vpd_min_effect_standardized.pdf", p_vpd_std, width = 6, height = 5, dpi = 300)
ggsave("figures/PAD_vpd_min_effect_standardized.png", p_vpd_std, width = 6, height = 5, dpi = 300)

ggsave("figures/PAD_vpd_min_effect_raw.pdf", p_vpd_raw, width = 6, height = 5, dpi = 300)
ggsave("figures/PAD_vpd_min_effect_raw.png", p_vpd_raw, width = 6, height = 5, dpi = 300)

cat("\n=== PLOTS SAVED ===\n")
cat("VPD min effect plots:\n")
cat("  - Standardized: figures/PAD_vpd_min_effect_standardized.*\n")
cat("  - Raw values: figures/PAD_vpd_min_effect_raw.*\n")

cat("\n=== VPD-ONLY VISUALIZATION COMPLETE ===\n")
cat("✓ Used reliable vpd-only model\n")
cat("✓ Shows VPD min effect on PAD\n")
cat("✓ Sites colored for spatial context\n")