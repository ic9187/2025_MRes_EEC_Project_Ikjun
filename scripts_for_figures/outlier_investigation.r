# ============================================================================
# 04c_outlier_investigation.R - UPDATED
# Outlier investigation for VPD-only phenology model and nymph peak model
# Goal: Identify and understand outliers in the reliable vpd_only model
# Updated: Focus on vpd_only model (not interaction) and nymph peak outliers only
# ============================================================================

cat("=== OUTLIER INVESTIGATION FOR VPD-ONLY MODEL ===\n")
cat("Investigating outliers in vpd_only phenology model and nymph peak model\n\n")

# Load required libraries
library(dplyr)
library(ggplot2)
library(lme4)
library(influence.ME)  # For mixed model influence measures
library(gridExtra)

# ============================================================================
# SECTION 1: VERIFY MODELS AND DATA
# ============================================================================

cat("=== SECTION 1: VERIFYING MODELS AND DATA ===\n")

# Check that required models exist
if (!exists("phenology_models")) {
  stop("Error: phenology_models not found. Please run 04a_phenology_analysis.R first.")
}

if (!exists("peak_models")) {
  stop("Error: peak_models not found. Please run 04b_phenology_additional.R first.")
}

if (!exists("analysis_data")) {
  stop("Error: analysis_data not found. Please run your phenology analysis scripts first.")
}

# Check that vpd_only model exists (the reliable model)
if (!"vpd_only" %in% names(phenology_models)) {
  stop("Error: vpd_only model not found in phenology_models.")
}

# Check that nymph peak model exists
if (!"nymph_peak" %in% names(peak_models)) {
  stop("Error: nymph_peak model not found in peak_models.")
}

# Focus on only these two models since larva model was fine
models_to_analyze <- list(
  "vpd_only" = phenology_models[["vpd_only"]],
  "nymph_peak" = peak_models[["nymph_peak"]]
)

cat("✓ VPD-only model loaded from phenology_models\n")
cat("✓ Nymph peak model loaded from peak_models\n")
cat("✓ Analysis focuses on 2 models (larva model was fine)\n")
cat("✓ Data available: n =", nrow(analysis_data), "observations\n\n")

# ============================================================================
# SECTION 2: OUTLIER DETECTION FUNCTIONS
# ============================================================================

cat("=== SECTION 2: OUTLIER DETECTION FUNCTIONS ===\n")

# Function to identify outliers using multiple criteria
identify_outliers <- function(model, model_name, data) {
  cat("Analyzing outliers for:", model_name, "\n")
  
  # Get basic model info
  n_obs <- nobs(model)
  fitted_vals <- fitted(model)
  residuals_vals <- residuals(model)
  
  # 1. Standardized residuals
  std_residuals <- residuals_vals / sd(residuals_vals)
  std_outliers <- which(abs(std_residuals) > 2.5)
  
  # 2. Cook's distance equivalent for mixed models
  # Using influence.ME package for proper mixed model influence
  influence_data <- NULL
  cooks_outliers <- c()
  
  tryCatch({
    # This might fail for some model structures
    influence_data <- influence(model, obs = TRUE)
    cooks_d <- cooks.distance(influence_data)
    cooks_threshold <- 4 / n_obs
    cooks_outliers <- which(cooks_d > cooks_threshold)
  }, error = function(e) {
    cat("  Note: Cook's distance calculation failed for", model_name, "\n")
  })
  
  # Combine outlier indices
  all_outlier_indices <- unique(c(std_outliers, cooks_outliers))
  
  # Create results dataframe
  outlier_results <- data.frame(
    observation = 1:n_obs,
    fitted = fitted_vals,
    residual = residuals_vals,
    std_residual = std_residuals,
    abs_std_residual = abs(std_residuals),
    is_std_outlier = 1:n_obs %in% std_outliers,
    is_cooks_outlier = 1:n_obs %in% cooks_outliers,
    is_any_outlier = 1:n_obs %in% all_outlier_indices,
    stringsAsFactors = FALSE
  )
  
  # Add Cook's distance if available
  if (!is.null(influence_data)) {
    outlier_results$cooks_distance <- cooks.distance(influence_data)
  }
  
  # Add site information
  if ("siteID" %in% names(data)) {
    outlier_results$siteID <- data$siteID[1:n_obs]
  }
  
  # Add year information if available
  if ("year" %in% names(data)) {
    outlier_results$year <- data$year[1:n_obs]
  }
  
  # Create site_year identifier for clear identification
  if ("siteID" %in% names(data) && "year" %in% names(data)) {
    outlier_results$site_year <- paste(data$siteID[1:n_obs], data$year[1:n_obs], sep = "_")
  } else if ("siteID" %in% names(data)) {
    outlier_results$site_year <- as.character(data$siteID[1:n_obs])
  } else {
    outlier_results$site_year <- paste("obs", 1:n_obs, sep = "_")
  }
  
  # Summary statistics
  n_std_outliers <- length(std_outliers)
  n_cooks_outliers <- length(cooks_outliers)
  n_total_outliers <- length(all_outlier_indices)
  
  cat("  Standardized residual outliers (|z| > 2.5):", n_std_outliers, "\n")
  if (n_std_outliers > 0) {
    cat("    Row numbers:", paste(std_outliers, collapse = ", "), "\n")
  }
  cat("  Cook's distance outliers:", n_cooks_outliers, "\n")
  if (n_cooks_outliers > 0) {
    cat("    Row numbers:", paste(cooks_outliers, collapse = ", "), "\n")
  }
  cat("  Total unique outliers:", n_total_outliers, "\n")
  if (n_total_outliers > 0) {
    cat("    Row numbers:", paste(all_outlier_indices, collapse = ", "), "\n")
  }
  cat("  Percentage of outliers:", round(100 * n_total_outliers / n_obs, 1), "%\n\n")
  
  return(outlier_results)
}

# Function to create outlier diagnostic plots
create_outlier_plots <- function(outlier_data, model_name) {
  
  # Plot 1: Residuals vs Fitted
  p1 <- ggplot(outlier_data, aes(x = fitted, y = residual)) +
    geom_point(aes(color = is_any_outlier, size = is_any_outlier)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = c(-2.5, 2.5) * sd(outlier_data$residual), 
               linetype = "dotted", color = "red") +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
    scale_size_manual(values = c("FALSE" = 1, "TRUE" = 2)) +
    labs(title = paste(model_name, "- Residuals vs Fitted"),
         x = "Fitted Values", y = "Residuals") +
    theme_minimal() +
    guides(color = "none", size = "none")
  
  # Plot 2: Q-Q plot of residuals
  p2 <- ggplot(outlier_data, aes(sample = std_residual)) +
    stat_qq(aes(color = is_any_outlier, size = is_any_outlier)) +
    stat_qq_line() +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
    scale_size_manual(values = c("FALSE" = 1, "TRUE" = 2)) +
    labs(title = paste(model_name, "- Q-Q Plot"),
         x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal() +
    guides(color = "none", size = "none")
  
  # Plot 3: Standardized residuals
  p3 <- ggplot(outlier_data, aes(x = observation, y = abs_std_residual)) +
    geom_point(aes(color = is_any_outlier, size = is_any_outlier)) +
    geom_hline(yintercept = 2.5, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
    scale_size_manual(values = c("FALSE" = 1, "TRUE" = 2)) +
    labs(title = paste(model_name, "- Standardized Residuals"),
         x = "Observation", y = "|Standardized Residual|") +
    theme_minimal() +
    guides(color = "none", size = "none")
  
  # Plot 4: Cook's distance (if available)
  p4 <- NULL
  if ("cooks_distance" %in% names(outlier_data)) {
    p4 <- ggplot(outlier_data, aes(x = observation, y = cooks_distance)) +
      geom_point(aes(color = is_cooks_outlier, size = is_cooks_outlier)) +
      geom_hline(yintercept = 4/nrow(outlier_data), linetype = "dashed", color = "red") +
      scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
      scale_size_manual(values = c("FALSE" = 1, "TRUE" = 2)) +
      labs(title = paste(model_name, "- Cook's Distance"),
           x = "Observation", y = "Cook's Distance") +
      theme_minimal() +
      guides(color = "none", size = "none")
  }
  
  return(list(residuals_fitted = p1, qq_plot = p2, std_residuals = p3, cooks_distance = p4))
}

# ============================================================================
# SECTION 3: RUN OUTLIER ANALYSIS
# ============================================================================

cat("=== SECTION 3: OUTLIER ANALYSIS ===\n")

# Store outlier results for models
outlier_results <- list()

# Analyze each model
for (model_name in names(models_to_analyze)) {
  model <- models_to_analyze[[model_name]]
  
  # Run outlier detection
  outlier_data <- identify_outliers(model, model_name, analysis_data)
  outlier_results[[model_name]] <- outlier_data
  
  # Create diagnostic plots
  plots <- create_outlier_plots(outlier_data, model_name)
  
  # Display plots
  if (!is.null(plots$cooks_distance)) {
    grid.arrange(plots$residuals_fitted, plots$qq_plot, 
                 plots$std_residuals, plots$cooks_distance, 
                 ncol = 2, top = paste("Diagnostic Plots -", model_name))
  } else {
    grid.arrange(plots$residuals_fitted, plots$qq_plot, plots$std_residuals, 
                 ncol = 2, top = paste("Diagnostic Plots -", model_name))
  }
}

# ============================================================================
# SECTION 4: OUTLIER CHARACTERIZATION
# ============================================================================

cat("\n=== SECTION 4: OUTLIER CHARACTERIZATION ===\n")

# Function to characterize outliers
characterize_outliers <- function(outlier_data, model_name, original_data) {
  
  outliers <- outlier_data[outlier_data$is_any_outlier, ]
  
  if (nrow(outliers) == 0) {
    cat(model_name, "- No outliers detected\n\n")
    return(NULL)
  }
  
  cat(model_name, "- Detailed outlier analysis:\n")
  cat("Number of outliers:", nrow(outliers), "\n")
  
  # Show specific site_year combinations that are outliers
  if ("site_year" %in% names(outliers)) {
    cat("OUTLIER SITE_YEAR COMBINATIONS:\n")
    outlier_site_years <- outliers$site_year
    for(i in 1:length(outlier_site_years)) {
      cat("  Row", outliers$observation[i], ":", outlier_site_years[i], 
          "(std residual =", round(outliers$std_residual[i], 2), ")\n")
    }
    cat("\n")
  }
  
  # Site patterns
  if ("siteID" %in% names(outliers)) {
    outlier_sites <- table(outliers$siteID)
    cat("Outliers by site:\n")
    print(outlier_sites)
    
    # Check if any sites have multiple outliers
    multi_outlier_sites <- names(outlier_sites)[outlier_sites > 1]
    if (length(multi_outlier_sites) > 0) {
      cat("Sites with multiple outliers:", paste(multi_outlier_sites, collapse = ", "), "\n")
    }
  }
  
  # Year patterns (if available)
  if ("year" %in% names(outliers)) {
    outlier_years <- table(outliers$year)
    cat("Outliers by year:\n")
    print(outlier_years)
    
    # Check if any years have multiple outliers
    multi_outlier_years <- names(outlier_years)[outlier_years > 1]
    if (length(multi_outlier_years) > 0) {
      cat("Years with multiple outliers:", paste(multi_outlier_years, collapse = ", "), "\n")
    }
  }
  
  # Residual characteristics
  cat("Outlier residual statistics:\n")
  cat("  Mean residual:", round(mean(outliers$residual), 3), "\n")
  cat("  Max |residual|:", round(max(abs(outliers$residual)), 3), "\n")
  cat("  Range std residuals:", round(range(outliers$std_residual), 3), "\n")
  
  # Get original variable values for outliers
  outlier_indices <- outliers$observation
  if (nrow(original_data) >= max(outlier_indices)) {
    outlier_original <- original_data[outlier_indices, ]
    
    cat("Environmental characteristics of outliers:\n")
    if ("mean_density_ha_std" %in% names(outlier_original)) {
      cat("  Host density (std): mean =", round(mean(outlier_original$mean_density_ha_std, na.rm = TRUE), 3),
          ", range =", paste(round(range(outlier_original$mean_density_ha_std, na.rm = TRUE), 3), collapse = " to "), "\n")
    }
    if ("mean_annual_vpdmin_std" %in% names(outlier_original)) {
      cat("  VPD min (std): mean =", round(mean(outlier_original$mean_annual_vpdmin_std, na.rm = TRUE), 3),
          ", range =", paste(round(range(outlier_original$mean_annual_vpdmin_std, na.rm = TRUE), 3), collapse = " to "), "\n")
    }
    
    # Response variable values for outliers
    response_var <- ifelse(model_name == "nymph_peak", "nymph_peak_day", "PAD")
    if (response_var %in% names(outlier_original)) {
      cat("  ", response_var, ": mean =", round(mean(outlier_original[[response_var]], na.rm = TRUE), 1),
          ", range =", paste(round(range(outlier_original[[response_var]], na.rm = TRUE), 1), collapse = " to "), "\n")
    }
  }
  
  cat("\n")
  
  return(outliers)
}

# Characterize outliers for each model
outlier_summaries <- list()

for (model_name in names(outlier_results)) {
  outlier_summary <- characterize_outliers(outlier_results[[model_name]], model_name, analysis_data)
  outlier_summaries[[model_name]] <- outlier_summary
}

# ============================================================================
# SECTION 5: CROSS-MODEL OUTLIER COMPARISON
# ============================================================================

cat("=== SECTION 5: CROSS-MODEL OUTLIER COMPARISON ===\n")

# Find observations that are outliers across both models
if (length(outlier_results) == 2) {
  model_names <- names(outlier_results)
  
  # Get outlier observation numbers for each model
  outlier_obs <- lapply(outlier_results, function(x) x$observation[x$is_any_outlier])
  
  # Find overlaps between models
  cat("OUTLIER OVERLAPS BETWEEN MODELS:\n")
  
  vpd_outliers <- outlier_obs[["vpd_only"]]
  nymph_outliers <- outlier_obs[["nymph_peak"]]
  common_outliers <- intersect(vpd_outliers, nymph_outliers)
  
  cat("→ VPD-only model outliers:", length(vpd_outliers), "observations")
  if (length(vpd_outliers) > 0) {
    cat(" (rows:", paste(vpd_outliers, collapse = ", "), ")")
  }
  cat("\n")
  
  cat("→ Nymph peak model outliers:", length(nymph_outliers), "observations")
  if (length(nymph_outliers) > 0) {
    cat(" (rows:", paste(nymph_outliers, collapse = ", "), ")")
  }
  cat("\n")
  
  cat("→ Common outliers:", length(common_outliers), "observations")
  if (length(common_outliers) > 0) {
    cat(" (rows:", paste(common_outliers, collapse = ", "), ")")
    
    # Get site_year information for common outliers
    if ("site_year" %in% names(analysis_data)) {
      common_site_years <- analysis_data$site_year[common_outliers]
      cat(" | site_years:", paste(unique(common_site_years), collapse = ", "))
    }
  }
  cat("\n")
  
  # Model-specific outliers
  cat("\nMODEL-SPECIFIC OUTLIERS:\n")
  vpd_specific <- setdiff(vpd_outliers, nymph_outliers)
  nymph_specific <- setdiff(nymph_outliers, vpd_outliers)
  
  cat("→ VPD-only model specific:", length(vpd_specific), "outliers")
  if (length(vpd_specific) > 0) {
    cat(" (rows:", paste(vpd_specific, collapse = ", "), ")")
  }
  cat("\n")
  
  cat("→ Nymph peak model specific:", length(nymph_specific), "outliers")
  if (length(nymph_specific) > 0) {
    cat(" (rows:", paste(nymph_specific, collapse = ", "), ")")
  }
  cat("\n")
}

# ============================================================================
# SECTION 6: DETAILED OUTLIER SUMMARY
# ============================================================================

cat("\n=== SECTION 6: DETAILED OUTLIER SUMMARY ===\n")

# Create detailed outlier summary for user review
cat("DETAILED OUTLIER IDENTIFICATION:\n\n")

for (model_name in names(outlier_results)) {
  outlier_data <- outlier_results[[model_name]]
  outliers <- outlier_data[outlier_data$is_any_outlier, ]
  
  if (nrow(outliers) > 0) {
    cat("", toupper(model_name), "MODEL OUTLIERS:\n")
    
    # Show specific site_year and row information prominently
    cat("OUTLIER SITE_YEAR COMBINATIONS:\n")
    if ("site_year" %in% names(outliers)) {
      for(i in 1:nrow(outliers)) {
        cat("  • Row", outliers$observation[i], ":", outliers$site_year[i], "\n")
        cat("    Standardized residual =", round(outliers$std_residual[i], 3), "\n")
        if ("cooks_distance" %in% names(outliers)) {
          cat("    Cook's distance =", round(outliers$cooks_distance[i], 4), "\n")
        }
        cat("\n")
      }
    }
    
    # Show detailed outlier info in table format
    cat("DETAILED OUTLIER TABLE:\n")
    outlier_summary <- outliers[, c("observation", "site_year", "fitted", "residual", "std_residual")]
    if ("cooks_distance" %in% names(outliers)) {
      outlier_summary$cooks_distance <- round(outliers$cooks_distance, 4)
    }
    
    print(outlier_summary)
    cat("\n")
  }
}

cat("=== SUMMARY: ALL OUTLIER SITE_YEARS ===\n")

# Collect all outlier site_years across models
all_outlier_site_years <- c()
outlier_summary_table <- data.frame()

for (model_name in names(outlier_results)) {
  outlier_data <- outlier_results[[model_name]]
  outliers <- outlier_data[outlier_data$is_any_outlier, ]
  
  if (nrow(outliers) > 0 && "site_year" %in% names(outliers)) {
    model_outlier_site_years <- outliers$site_year
    all_outlier_site_years <- c(all_outlier_site_years, model_outlier_site_years)
    
    # Create summary table
    temp_summary <- data.frame(
      model = model_name,
      row_number = outliers$observation,
      site_year = outliers$site_year,
      std_residual = round(outliers$std_residual, 3),
      stringsAsFactors = FALSE
    )
    
    if ("cooks_distance" %in% names(outliers)) {
      temp_summary$cooks_distance <- round(outliers$cooks_distance, 4)
    }
    
    outlier_summary_table <- rbind(outlier_summary_table, temp_summary)
  }
}

# Show unique outlier site_years
unique_outlier_site_years <- unique(all_outlier_site_years)
cat("UNIQUE OUTLIER SITE_YEARS ACROSS BOTH MODELS:\n")
if (length(unique_outlier_site_years) > 0) {
  for(i in 1:length(unique_outlier_site_years)) {
    site_year <- unique_outlier_site_years[i]
    models_affected <- unique(outlier_summary_table$model[outlier_summary_table$site_year == site_year])
    cat("  •", site_year, "- affects:", paste(models_affected, collapse = ", "), "\n")
  }
} else {
  cat("  No outliers detected\n")
}

cat("\nCOMPLETE OUTLIER SUMMARY TABLE:\n")
if (nrow(outlier_summary_table) > 0) {
  print(outlier_summary_table)
} else {
  cat("No outliers detected in any model\n")
}

# ============================================================================
# SECTION 7: SENSITIVITY ANALYSIS CODE
# ============================================================================

cat("\n=== SECTION 7: SENSITIVITY ANALYSIS CODE ===\n")

# Example code for user to run sensitivity analysis
cat("EXAMPLE CODE FOR SENSITIVITY ANALYSIS:\n")
cat("# Investigate these outlying observations before deciding on any action\n\n")

# Show which site_years are outliers for each model
for (model_name in names(outlier_results)) {
  outlier_indices <- which(outlier_results[[model_name]]$is_any_outlier)
  
  if (length(outlier_indices) > 0) {
    # Get site_year info for outliers
    outlier_site_years <- outlier_results[[model_name]]$site_year[outlier_indices]
    cat("# ", toupper(model_name), " model outliers detected:\n")
    cat("# Site_years:", paste(outlier_site_years, collapse = ", "), "\n")
    cat("# Row numbers:", paste(outlier_indices, collapse = ", "), "\n")
    
    # Generate appropriate model code based on model name
    var_name <- paste0(tolower(model_name), "_outliers")
    data_name <- paste0("clean_data_", tolower(model_name))
    
    cat("# Option 1: Examine outlying observations\n")
    cat("outlying_observations_", tolower(model_name), " <- analysis_data[c(", paste(outlier_indices, collapse = ", "), "), ]\n")
    cat("View(outlying_observations_", tolower(model_name), ")\n\n")
    
    cat("# Option 2: Sensitivity analysis (if outliers deemed problematic)\n")
    if (model_name == "vpd_only") {
      cat(var_name, " <- c(", paste(outlier_indices, collapse = ", "), ")  # Rows to exclude\n")
      cat(data_name, " <- analysis_data[-", var_name, ", ]\n")
      cat("vpd_model_clean <- lmer(PAD ~ mean_annual_vpdmin_std + (1|siteID), data = ", data_name, ")\n")
      cat("summary(vpd_model_clean)\n")
      cat("# Compare: summary(phenology_models[['vpd_only']])  # Original model\n\n")
    } else if (model_name == "nymph_peak") {
      cat(var_name, " <- c(", paste(outlier_indices, collapse = ", "), ")  # Rows to exclude\n")
      cat(data_name, " <- analysis_data[-", var_name, ", ]\n")
      cat("nymph_model_clean <- lmer(nymph_peak_day ~ mean_annual_vpdmin_std + (1|siteID), data = ", data_name, ")\n")
      cat("summary(nymph_model_clean)\n")
      cat("# Compare: summary(peak_models[['nymph_peak']])  # Original model\n\n")
    }
  } else {
    cat("# No outliers detected in", toupper(model_name), "model\n\n")
  }
}

cat("\nRECOMMENDED APPROACH:\n")
cat("1. INVESTIGATE FIRST: Examine outlying observations for data quality issues or extreme conditions\n")
cat("2. SENSITIVITY ANALYSIS: Test model robustness by comparing results with/without outliers\n")
cat("3. ROBUST METHODS: Consider alternatives to outlier removal:\n")
cat("   → Robust regression: Use rlmer() from robustlmm package\n")
cat("   → Winsorizing: Cap extreme values at 95th/5th percentiles\n")
cat("   → Transformation: Log or square root transform response variables\n")
cat("   → Bootstrap: Use bootstrap confidence intervals\n")
cat("4. TRANSPARENT REPORTING: Document outlier handling decisions in methods section\n\n")

# ============================================================================
# SECTION 8: SAVE RESULTS
# ============================================================================

cat("=== SECTION 8: SAVING RESULTS ===\n")

# Create output directory
if (!dir.exists("output")) {
  dir.create("output", recursive = TRUE)
}

# Save outlier analysis results
saveRDS(outlier_results, "output/04c_vpd_nymph_outlier_analysis.RDS")

# Create summary report
outlier_summary_report <- list(
  analysis_date = Sys.Date(),
  models_analyzed = names(models_to_analyze),
  approach = "VPD-only and nymph peak models (larva model was fine)",
  total_observations = nrow(analysis_data),
  outlier_counts = sapply(outlier_results, function(x) sum(x$is_any_outlier, na.rm = TRUE)),
  outlier_percentages = sapply(outlier_results, function(x) 100 * sum(x$is_any_outlier, na.rm = TRUE) / nrow(x)),
  recommendations = "Investigate outliers before removal - see detailed output for specific recommendations",
  method = "Standardized residuals > 2.5 and Cook's distance > 4/n",
  focus = "Updated analysis focusing on reliable vpd_only model and nymph peak outliers only"
)

saveRDS(outlier_summary_report, "output/04c_vpd_nymph_outlier_summary.RDS")

# Save detailed outlier data
for (model_name in names(outlier_results)) {
  filename <- paste0("output/04c_outliers_", model_name, ".csv")
  write.csv(outlier_results[[model_name]], filename, row.names = FALSE)
}

cat("✓ Outlier analysis results saved:\n")
cat("  - Complete analysis: 04c_vpd_nymph_outlier_analysis.RDS\n")
cat("  - Summary report: 04c_vpd_nymph_outlier_summary.RDS\n")
cat("  - Individual model outliers: 04c_outliers_[model_name].csv\n\n")

cat("=== OUTLIER INVESTIGATION COMPLETE ===\n")
cat("Analyzed 2 models: vpd_only, nymph_peak (larva model was fine)\n")
cat("Review the diagnostic plots and outlier characteristics above.\n")
cat("NEXT STEPS:\n")
cat("→ Investigate ecological context of outlying site_years\n") 
cat("→ Check for data quality issues in flagged observations\n")
cat("→ Run sensitivity analysis to test model robustness\n")
cat("→ Consider robust modeling approaches if needed\n")
cat("→ Document outlier handling decisions transparently\n")