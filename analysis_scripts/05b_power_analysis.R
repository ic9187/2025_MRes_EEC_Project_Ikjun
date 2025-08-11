# ============================================================================
# 06_prevalence_power_analysis.R
# Power analysis for prevalence models only
# Models: PAD, host density, VPD min
# Clean summary table with unified power curve
# ============================================================================

cat("=== PREVALENCE POWER AND EFFECT SIZE ANALYSIS ===\n")
cat("Analysis: PAD, host density, and VPD min effects on transmission\n\n")

# Load required libraries
library(WebPower)
library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================================================
# SECTION 1: VERIFY REQUIRED OBJECTS ARE LOADED
# ============================================================================

cat("=== SECTION 1: VERIFYING REQUIRED OBJECTS ===\n")

required_objects <- c("prevalence_models", "model_data")
missing_objects <- required_objects[!sapply(required_objects, exists)]

if(length(missing_objects) > 0) {
  stop("Missing required objects. Please run analyses first:\n",
       "  Missing: ", paste(missing_objects, collapse = ", "), "\n",
       "  Run: 05_prevalence_analysis.R")
}

cat("✓ Required objects found:\n")
cat("  Prevalence models:", length(prevalence_models), "\n")
cat("  Prevalence sample size:", nrow(model_data), "\n\n")

# ============================================================================
# SECTION 2: PREVALENCE SUMMARY ANALYSIS  
# ============================================================================

cat("=== SECTION 2: PREVALENCE SUMMARY ANALYSIS ===\n")

# Define which models to analyze (excluding larva_only and nymph_only)
prevalence_model_names <- c("pad_only", "host_only", "vpd_only")
missing_models <- prevalence_model_names[!prevalence_model_names %in% names(prevalence_models)]

if(length(missing_models) > 0) {
  cat("Warning: Missing prevalence models:", paste(missing_models, collapse = ", "), "\n")
  prevalence_model_names <- prevalence_model_names[prevalence_model_names %in% names(prevalence_models)]
}

cat("Analyzing", length(prevalence_model_names), "prevalence models\n")

# Function to extract model metrics
extract_model_metrics <- function(model, model_name) {
  # Get R²
  model_summary <- summary(model)
  r2 <- model_summary$r.squared
  
  # Calculate Cohen's f²
  cohens_f2 <- r2 / (1 - r2)
  
  # Effect size category
  effect_category <- case_when(
    cohens_f2 < 0.02 ~ "Negligible",
    cohens_f2 < 0.15 ~ "Small", 
    cohens_f2 < 0.35 ~ "Medium",
    TRUE ~ "Large"
  )
  
  # Current power
  n_obs <- nrow(model_data)
  current_power <- wp.regression(n = n_obs, p1 = 1, p2 = 0, f2 = cohens_f2)$power
  
  # Sample size for 80% power
  n_for_80 <- wp.regression(p1 = 1, p2 = 0, f2 = cohens_f2, power = 0.8)$n
  
  # Variable category
  variable_category <- case_when(
    model_name == "pad_only" ~ "Tick Measure",
    model_name == "host_only" ~ "Host",
    model_name == "vpd_only" ~ "Climate",
    TRUE ~ "Other"
  )
  
  # Clean model name
  clean_name <- case_when(
    model_name == "pad_only" ~ "PAD",
    model_name == "host_only" ~ "Host Density", 
    model_name == "vpd_only" ~ "VPD Min",
    TRUE ~ model_name
  )
  
  return(data.frame(
    Model = clean_name,
    Variable_Category = variable_category,
    R_squared = r2,
    Cohens_f2 = cohens_f2,
    Effect_Category = effect_category,
    Current_Power = current_power,
    N_for_80_Power = ceiling(n_for_80),
    stringsAsFactors = FALSE
  ))
}

# Extract metrics for all prevalence models
prevalence_results <- data.frame()

for(model_name in prevalence_model_names) {
  model <- prevalence_models[[model_name]]
  metrics <- extract_model_metrics(model, model_name)
  prevalence_results <- rbind(prevalence_results, metrics)
}

# Sort by effect size
prevalence_results <- prevalence_results %>%
  arrange(desc(Cohens_f2)) %>%
  mutate(Rank = 1:n())

cat("PREVALENCE ANALYSIS SUMMARY:\n")
prevalence_display <- prevalence_results %>%
  mutate(across(c(R_squared, Cohens_f2, Current_Power), ~ round(.x, 4))) %>%
  dplyr::select(Rank, Model, Variable_Category, R_squared, Cohens_f2, Effect_Category, Current_Power, N_for_80_Power)
print(prevalence_display)

# ============================================================================
# SECTION 3: PREVALENCE POWER CURVES
# ============================================================================

cat("\n=== SECTION 3: PREVALENCE POWER CURVES ===\n")

# Sample size range for power curves
sample_sizes <- seq(10, 200, by = 5)

# Generate prevalence power curves
prevalence_curves <- data.frame()

for(i in 1:nrow(prevalence_results)) {
  model_name <- prevalence_results$Model[i]
  cohens_f2 <- prevalence_results$Cohens_f2[i]
  
  power_values <- sapply(sample_sizes, function(n) {
    wp.regression(n = n, p1 = 1, p2 = 0, f2 = cohens_f2)$power
  })
  
  curve_data <- data.frame(
    model = model_name,
    category = prevalence_results$Variable_Category[i],
    sample_size = sample_sizes,
    power = power_values
  )
  
  prevalence_curves <- rbind(prevalence_curves, curve_data)
}

cat("Power curves calculated for", length(prevalence_model_names), "models\n")

# ============================================================================
# SECTION 4: VISUALIZATION
# ============================================================================

cat("\n=== SECTION 4: CREATING POWER CURVE VISUALIZATION ===\n")

# Prevalence unified power curve
cat("Creating prevalence power curve...\n")

p_prevalence <- ggplot(prevalence_curves, aes(x = sample_size, y = power, color = model)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_vline(xintercept = nrow(model_data), linetype = "dotted", color = "blue", alpha = 0.7) +
  
  scale_color_manual(
    name = "Variable",
    values = c(
      "PAD" = "#1F78B4",
      "Host Density" = "#33A02C",
      "VPD Min" = "#FF7F00"
    )
  ) +
  
  labs(
    title = "Prevalence Analysis: Power Curves by Variable",
    subtitle = paste("Individual predictor comparison | Current n =", nrow(model_data)),
    x = "Sample Size",
    y = "Statistical Power", 
    caption = "Red line = 80% power threshold | Blue line = current sample size"
  ) +
  
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::percent) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
  )

# Display plot
print(p_prevalence)

# ============================================================================
# SECTION 5: COMPREHENSIVE SUMMARY
# ============================================================================

cat("\n=== SECTION 5: COMPREHENSIVE SUMMARY ===\n")

# Prevalence insights
cat("PREVALENCE ANALYSIS INSIGHTS:\n")
strongest_prev_effect <- prevalence_results$Model[1]  # Already sorted by effect size
cat("  Strongest effect:", strongest_prev_effect,
    "(f² =", round(prevalence_results$Cohens_f2[1], 4), ")\n")

adequate_power_prev <- sum(prevalence_results$Current_Power >= 0.8)
cat("  Models with adequate power (≥80%):", adequate_power_prev, "of", nrow(prevalence_results), "\n")

# Sample size recommendations
max_n_prev <- max(prevalence_results$N_for_80_Power)
cat("  Sample size needed for 80% power:", max_n_prev, "(currently n =", nrow(model_data), ")\n")

# Variable category insights
cat("\nBY VARIABLE CATEGORY:\n")
category_summary <- prevalence_results %>%
  group_by(Variable_Category) %>%
  summarise(
    max_f2 = max(Cohens_f2),
    max_power = max(Current_Power),
    .groups = "drop"
  )
print(category_summary)

# ============================================================================
# SECTION 6: SAVE RESULTS
# ============================================================================

cat("\n=== SECTION 6: SAVING RESULTS ===\n")

# Create output directory
if (!dir.exists("output")) {
  dir.create("output", recursive = TRUE)
}

# Save plot
ggsave("output/prevalence_power_curves.png", p_prevalence, width = 12, height = 8, dpi = 300)

# Save summary tables
write.csv(prevalence_results, "output/prevalence_power_summary.csv", row.names = FALSE)
write.csv(prevalence_curves, "output/prevalence_power_curves_data.csv", row.names = FALSE)

# Save comprehensive results
prevalence_power_results <- list(
  prevalence_summary = prevalence_results,
  prevalence_curves = prevalence_curves,
  sample_size = nrow(model_data),
  insights = list(
    strongest_effect = strongest_prev_effect,
    adequate_power_count = adequate_power_prev,
    max_n_needed = max_n_prev
  ),
  plot = p_prevalence
)

saveRDS(prevalence_power_results, "output/06_prevalence_power_results.RDS")

cat("✓ Results saved:\n")
cat("  - Power curve plot: output/prevalence_power_curves.png\n")
cat("  - Summary table: output/prevalence_power_summary.csv\n")
cat("  - Power curve data: output/prevalence_power_curves_data.csv\n")
cat("  - Complete results: output/06_prevalence_power_results.RDS\n")

cat("\n=== PREVALENCE POWER ANALYSIS COMPLETE ===\n")
cat("\nKEY INSIGHTS:\n")
cat("→ Clean comparison of PAD, host density, and VPD min effects\n")
cat("→ Effect sizes quantified with Cohen's f²\n") 
cat("→ Power curves show sample size requirements\n")
cat("→ Strongest predictor:", strongest_prev_effect, "\n")