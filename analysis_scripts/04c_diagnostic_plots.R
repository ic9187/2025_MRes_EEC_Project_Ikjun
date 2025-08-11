# ============================================================================
# Simple Diagnostic Plots for Phenology Models
# Creates diagnostic plots for all 3 phenology models: host_only, vpd_only, interaction
# No saving, just visual diagnostics
# ============================================================================

cat("=== CREATING DIAGNOSTIC PLOTS FOR 3 PHENOLOGY MODELS ===\n")

# Load required libraries
library(ggplot2)
library(gridExtra)

# ============================================================================
# VERIFY MODELS EXIST
# ============================================================================

if (!exists("phenology_models")) {
  stop("Error: phenology_models not found. Please run 04a_phenology_analysis.R first.")
}

required_models <- c("host_only", "vpd_only", "host_vpd_interaction")
missing_models <- required_models[!required_models %in% names(phenology_models)]
if(length(missing_models) > 0) {
  stop("Error: Missing models: ", paste(missing_models, collapse = ", "))
}

cat("✓ All 3 models found: host_only, vpd_only, host_vpd_interaction\n\n")

# ============================================================================
# FUNCTION TO CREATE DIAGNOSTIC PLOTS
# ============================================================================

create_model_diagnostics <- function(model, model_name) {
  
  # Extract model components
  fitted_vals <- fitted(model)
  residuals_vals <- residuals(model)
  random_effects <- ranef(model)$siteID[,1]
  
  # Plot 1: Residuals vs Fitted
  p1 <- ggplot(data.frame(fitted = fitted_vals, residuals = residuals_vals), 
               aes(x = fitted, y = residuals)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "loess", se = TRUE, color = "blue", linewidth = 0.8) +
    labs(title = paste(model_name, "- Residuals vs Fitted"),
         x = "Fitted Values", y = "Residuals") +
    theme_minimal()
  
  # Plot 2: Q-Q plot for residuals
  p2 <- ggplot(data.frame(residuals = residuals_vals), aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line(color = "red", linewidth = 0.8) +
    labs(title = paste(model_name, "- Q-Q Plot (Residuals)"),
         x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal()
  
  # Plot 3: Scale-Location plot
  p3 <- ggplot(data.frame(fitted = fitted_vals, sqrt_abs_resid = sqrt(abs(residuals_vals))), 
               aes(x = fitted, y = sqrt_abs_resid)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = TRUE, color = "blue", linewidth = 0.8) +
    labs(title = paste(model_name, "- Scale-Location"),
         x = "Fitted Values", y = "√|Residuals|") +
    theme_minimal()
  
  # Plot 4: Q-Q plot for random effects
  p4 <- ggplot(data.frame(random_effects = random_effects), aes(sample = random_effects)) +
    stat_qq() +
    stat_qq_line(color = "red", linewidth = 0.8) +
    labs(title = paste(model_name, "- Q-Q Plot (Random Effects)"),
         x = "Theoretical Quantiles", y = "Random Effects") +
    theme_minimal()
  
  return(list(p1, p2, p3, p4))
}

# ============================================================================
# CREATE AND DISPLAY PLOTS FOR ALL MODELS
# ============================================================================

cat("=== CREATING DIAGNOSTIC PLOTS ===\n")

# Create plots for each model
for(model_name in names(phenology_models)) {
  cat("Creating plots for:", model_name, "\n")
  
  model <- phenology_models[[model_name]]
  plots <- create_model_diagnostics(model, toupper(model_name))
  
  # Display plots in 2x2 grid
  grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], 
               ncol = 2, 
               top = paste("Diagnostic Plots -", toupper(model_name), "Model"))
}

cat("\n=== DIAGNOSTIC PLOTS COMPLETE ===\n")
cat("Review plots for:\n")
cat("  - Residuals vs Fitted: Should show random scatter\n") 
cat("  - Q-Q Plot (Residuals): Points should follow line\n")
cat("  - Scale-Location: Should show horizontal trend\n")
cat("  - Q-Q Plot (Random Effects): Points should follow line\n")