# ============================================================================
# habitat_plots_phenology_direct.R
# Purpose: Create habitat and site boxplots for phenology analysis
# Input: pad_analysis_dataset.RDS (loaded directly for complete PAD data)
# Output: Boxplots showing PAD, nymph_peak_day, larva_peak_day patterns
# ENHANCED: Show all sites ordered by latitude, include ALL PAD values
# ============================================================================

library(dplyr)
library(ggplot2)
library(gridExtra)

cat("=== HABITAT AND SITE PLOTS: PHENOLOGY ANALYSIS ===\n")
cat("Creating boxplots for PAD and peak day patterns\n")
cat("ENHANCED: All sites shown, ordered by latitude, ALL PAD values included\n\n")

# ============================================================================
# SECTION 1: LOAD DATA DIRECTLY
# ============================================================================

cat("=== SECTION 1: LOADING DATA DIRECTLY ===\n")

# Load pad_analysis_dataset directly for complete visualization
if(file.exists("data/processed/pad_analysis_dataset.RDS")) {
  cat("Loading data directly from pad_analysis_dataset.RDS...\n")
  viz_data <- readRDS("data/processed/pad_analysis_dataset.RDS")
  
  # Check if nlcdClass column exists
  if(!"nlcdClass" %in% names(viz_data)) {
    stop("ERROR: nlcdClass column not found in pad_analysis_dataset.RDS\n",
         "SOLUTION: Please run the updated 03b_merge_phenology_data.R script first\n",
         "This will add habitat data (nlcdClass) to your dataset.")
  }
  
  cat("✓ nlcdClass column found in data\n")
  
  # Apply MINIMAL filtering - keep ALL PAD values including negatives
  viz_data <- viz_data %>%
    filter(!is.na(PAD) & !is.na(siteID)) %>%
    # NO PAD >= 0 filter - keep ALL PAD values for complete visualization
    filter(!is.na(mean_density_ha)) %>%
    mutate(siteID = as.factor(siteID))
  
  cat("✓ Loaded visualization data (ALL PAD values included)\n")
} else {
  stop("ERROR: pad_analysis_dataset.RDS not found\n",
       "SOLUTION: Please run scripts 03b_merge_phenology_data.R first")
}

cat("  Observations:", nrow(viz_data), "\n")
cat("  Sites:", length(unique(viz_data$siteID)), "\n")

# Check habitat data availability
n_with_habitat <- sum(!is.na(viz_data$nlcdClass))
cat("  Observations with habitat data:", n_with_habitat, "\n")

# Check peak day data availability
n_with_peaks <- sum(!is.na(viz_data$nymph_peak_day) & !is.na(viz_data$larva_peak_day))
cat("  Observations with peak day data:", n_with_peaks, "\n")

# Check coordinate data for latitude ordering
n_with_coords <- sum(!is.na(viz_data$lat) & !is.na(viz_data$long))
cat("  Observations with coordinates:", n_with_coords, "\n")

# More detailed habitat data check
if("nlcdClass" %in% names(viz_data)) {
  habitat_values <- unique(viz_data$nlcdClass)
  habitat_values <- habitat_values[!is.na(habitat_values)]
  cat("  Unique habitat types found:", length(habitat_values), "\n")
  if(length(habitat_values) > 0) {
    cat("  Habitat types:", paste(habitat_values, collapse = ", "), "\n")
  }
} else {
  cat("  nlcdClass column: NOT FOUND\n")
}

if(n_with_habitat == 0) {
  stop("ERROR: No habitat data available (all nlcdClass values are NA)\n",
       "POSSIBLE CAUSES:\n",
       "1. The updated 03b_merge_phenology_data.R script hasn't been run\n",
       "2. No habitat data exists for your sites in tick_raw\n",
       "3. Data filtering removed all records with habitat data\n",
       "\nSOLUTION:\n",
       "1. Run updated 03b_merge_phenology_data.R first\n",
       "2. Check that tick_raw contains nlcdClass data\n",
       "3. Verify your site IDs match between datasets")
}

# ============================================================================
# SECTION 2: DATA SUMMARY
# ============================================================================

cat("\n=== SECTION 2: DATA SUMMARY ===\n")

# Response variables summary (including negative PADs)
response_summary <- viz_data %>%
  filter(!is.na(PAD), !is.na(nymph_peak_day), !is.na(larva_peak_day)) %>%
  summarise(
    PAD_mean = mean(PAD),
    PAD_sd = sd(PAD),
    PAD_min = min(PAD),
    PAD_max = max(PAD),
    PAD_negative_count = sum(PAD < 0),
    PAD_zero_count = sum(PAD == 0),
    PAD_positive_count = sum(PAD > 0),
    Nymph_peak_mean = mean(nymph_peak_day),
    Nymph_peak_sd = sd(nymph_peak_day),
    Larva_peak_mean = mean(larva_peak_day),
    Larva_peak_sd = sd(larva_peak_day)
  )

cat("Response Variables Summary (ALL PAD values):\n")
print(response_summary)

# PAD distribution breakdown
total_pad_obs <- sum(!is.na(viz_data$PAD))
cat("\nPAD Distribution:\n")
cat("  Negative PAD (larvae earlier):", response_summary$PAD_negative_count, 
    paste0("(", round(response_summary$PAD_negative_count/total_pad_obs*100, 1), "%)"), "\n")
cat("  Zero PAD (perfect synchrony):", response_summary$PAD_zero_count, 
    paste0("(", round(response_summary$PAD_zero_count/total_pad_obs*100, 1), "%)"), "\n")
cat("  Positive PAD (larvae later):", response_summary$PAD_positive_count, 
    paste0("(", round(response_summary$PAD_positive_count/total_pad_obs*100, 1), "%)"), "\n")

# Habitat distribution
habitat_summary <- viz_data %>%
  filter(!is.na(nlcdClass)) %>%
  group_by(nlcdClass) %>%
  summarise(
    N_observations = n(),
    N_sites = n_distinct(siteID),
    Mean_PAD = mean(PAD, na.rm = TRUE),
    SD_PAD = sd(PAD, na.rm = TRUE),
    Mean_nymph_peak = mean(nymph_peak_day, na.rm = TRUE),
    Mean_larva_peak = mean(larva_peak_day, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(N_observations))

cat("\nHabitat Distribution:\n")
print(habitat_summary)

# Site distribution
site_summary <- viz_data %>%
  group_by(siteID) %>%
  summarise(
    N_observations = n(),
    Mean_PAD = mean(PAD, na.rm = TRUE),
    SD_PAD = sd(PAD, na.rm = TRUE),
    Mean_nymph_peak = mean(nymph_peak_day, na.rm = TRUE),
    Mean_larva_peak = mean(larva_peak_day, na.rm = TRUE),
    Habitat = first(nlcdClass),
    Latitude = first(lat),
    Longitude = first(long),
    .groups = "drop"
  ) %>%
  arrange(desc(N_observations))

cat("\nSite Summary (top 10):\n")
print(head(site_summary, 10))
cat("  Total sites:", nrow(site_summary), "\n\n")

# ============================================================================
# SECTION 3: HABITAT BOXPLOTS
# ============================================================================

cat("=== SECTION 3: CREATING HABITAT BOXPLOTS ===\n")

# Filter data with habitat information
habitat_data <- viz_data %>%
  filter(!is.na(nlcdClass))

cat("Creating habitat boxplots with", nrow(habitat_data), "observations\n")

# PAD by habitat
pad_habitat_plot <- ggplot(habitat_data, aes(x = nlcdClass, y = PAD)) +
  geom_boxplot(aes(fill = nlcdClass), alpha = 0.7, outlier.alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  
  # Styling
  scale_fill_viridis_d(name = "Habitat Type") +
  
  # Labels
  labs(
    title = "Peak Activity Difference (PAD) by Habitat Type",
    subtitle = "Temporal synchronicity of nymph and larva questing across land cover classes",
    x = "Habitat Type (NLCD Class)",
    y = "PAD (days)",
    caption = paste("n =", nrow(habitat_data), "observations from", 
                   length(unique(habitat_data$siteID)), "sites")
  ) +
  
  # Theme
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12, color = "gray40"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  
  # Add horizontal line at zero (perfect synchrony)
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7)

# Add sample sizes
habitat_counts <- habitat_data %>%
  count(nlcdClass) %>%
  mutate(label = paste0("n=", n))

pad_habitat_plot <- pad_habitat_plot +
  geom_text(data = habitat_counts, 
            aes(x = nlcdClass, y = max(habitat_data$PAD, na.rm = TRUE) * 0.9, 
                label = label),
            size = 3, color = "black", fontface = "bold")

print(pad_habitat_plot)

# Nymph peak day by habitat (if data available)
if(sum(!is.na(habitat_data$nymph_peak_day)) > 0) {
  
  nymph_habitat_plot <- ggplot(habitat_data, aes(x = nlcdClass, y = nymph_peak_day)) +
    geom_boxplot(aes(fill = nlcdClass), alpha = 0.7, outlier.alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
    
    # Styling
    scale_fill_viridis_d(name = "Habitat Type") +
    
    # Labels
    labs(
      title = "Nymph Peak Day by Habitat Type",
      subtitle = "Day of year when maximum nymph questing occurs",
      x = "Habitat Type (NLCD Class)",
      y = "Day of Year",
      caption = paste("n =", sum(!is.na(habitat_data$nymph_peak_day)), "observations")
    ) +
    
    # Theme
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray40"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
  
  print(nymph_habitat_plot)
} else {
  nymph_habitat_plot <- NULL
  cat("No nymph peak day data available - skipping nymph habitat plot\n")
}

# Larva peak day by habitat (if data available)
if(sum(!is.na(habitat_data$larva_peak_day)) > 0) {
  
  larva_habitat_plot <- ggplot(habitat_data, aes(x = nlcdClass, y = larva_peak_day)) +
    geom_boxplot(aes(fill = nlcdClass), alpha = 0.7, outlier.alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
    
    # Styling
    scale_fill_viridis_d(name = "Habitat Type") +
    
    # Labels
    labs(
      title = "Larva Peak Day by Habitat Type",
      subtitle = "Day of year when maximum larva questing occurs",
      x = "Habitat Type (NLCD Class)",
      y = "Day of Year",
      caption = paste("n =", sum(!is.na(habitat_data$larva_peak_day)), "observations")
    ) +
    
    # Theme
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray40"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
  
  print(larva_habitat_plot)
} else {
  larva_habitat_plot <- NULL
  cat("No larva peak day data available - skipping larva habitat plot\n")
}

# ============================================================================
# SECTION 4: PAD DISTRIBUTION HISTOGRAM
# ============================================================================

cat("\n=== SECTION 4: CREATING PAD DISTRIBUTION HISTOGRAM ===\n")

# Create PAD histogram
pad_data_for_hist <- viz_data %>%
  filter(!is.na(PAD))

cat("Creating PAD histogram with", nrow(pad_data_for_hist), "observations\n")

# Calculate summary stats for annotations
pad_mean <- mean(pad_data_for_hist$PAD)
pad_median <- median(pad_data_for_hist$PAD)
pad_sd <- sd(pad_data_for_hist$PAD)
n_negative <- sum(pad_data_for_hist$PAD < 0)
n_positive <- sum(pad_data_for_hist$PAD > 0)

# Create histogram
pad_histogram <- ggplot(pad_data_for_hist, aes(x = PAD)) +
  geom_histogram(bins = 25, fill = "skyblue", color = "white", alpha = 0.7) +
  
  # Add reference lines
  geom_vline(xintercept = 0, color = "red", linetype = "solid", size = 1, alpha = 0.8) +
  geom_vline(xintercept = pad_mean, color = "darkblue", linetype = "dashed", size = 0.8) +
  
  # Labels and styling
  labs(
    title = "Distribution of Peak Activity Difference (PAD)",
    subtitle = paste0("n = ", nrow(pad_data_for_hist), " | Mean = ", round(pad_mean, 1), 
                     " days | ", n_negative, " negative, ", n_positive, " positive"),
    x = "PAD (days)",
    y = "Frequency",
    caption = "Red line = perfect synchrony (PAD = 0) | Blue dashed line = mean PAD\nNegative = larvae earlier | Positive = larvae later"
  ) +
  
  # Theme
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12, color = "gray40"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.caption = element_text(size = 9, color = "gray50"),
    panel.grid.minor = element_blank()
  )

print(pad_histogram)

# ============================================================================
# SECTION 5: SITE BOXPLOTS (ALL SITES, ORDERED BY LATITUDE)
# ============================================================================

cat("\n=== SECTION 5: CREATING SITE BOXPLOTS ===\n")
cat("ENHANCED: Showing ALL sites, ordered by latitude (south to north)\n")

# Create site coordinate summary for ordering
site_coords <- viz_data %>%
  group_by(siteID) %>%
  summarise(
    latitude = first(lat),
    longitude = first(long),
    habitat = first(nlcdClass),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(latitude)) %>%  # Only sites with coordinates
  arrange(latitude)  # Order by latitude (south to north)

cat("Sites with coordinates:", nrow(site_coords), "\n")
if(nrow(site_coords) > 0) {
  cat("Latitude range:", round(min(site_coords$latitude, na.rm = TRUE), 2), "to", 
      round(max(site_coords$latitude, na.rm = TRUE), 2), "\n")
}

# Filter viz_data to sites with coordinates and create ordered factor
sites_with_coords <- viz_data %>%
  filter(siteID %in% site_coords$siteID) %>%
  mutate(siteID = factor(siteID, levels = site_coords$siteID))

cat("Creating site boxplot with", nrow(sites_with_coords), "observations from", 
    length(unique(sites_with_coords$siteID)), "sites\n")

if(nrow(sites_with_coords) > 0) {
  
  # PAD by site
  pad_site_plot <- ggplot(sites_with_coords, aes(x = siteID, y = PAD)) +
    geom_boxplot(aes(fill = nlcdClass), alpha = 0.7, outlier.alpha = 0.6) +
    geom_jitter(width = 0.3, alpha = 0.5, size = 0.8) +
    
    # Styling
    scale_fill_viridis_d(name = "Habitat Type") +
    
    # Labels
    labs(
      title = "Peak Activity Difference (PAD) by Site",
      subtitle = "ordered by latitude (mean = 40.0 (SD = 5.10), 29.7 to 46.2)",
      x = "Site ID",
      y = "PAD (days)",
      caption = paste("n =", nrow(sites_with_coords), "observations from", 
                     length(unique(sites_with_coords$siteID)), "sites")
    ) +
    
    # Theme
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray40"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    
    # Add horizontal line at zero
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7)
  
  print(pad_site_plot)
  
  # Show latitude information for first and last few sites
  cat("\nLatitude ordering verification:\n")
  lat_check <- site_coords %>%
    dplyr::select(siteID, latitude, habitat) %>%
    mutate(position = case_when(
      row_number() <= 3 ~ "Southernmost",
      row_number() >= n() - 2 ~ "Northernmost", 
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(position))
  
  print(lat_check)
  
} else {
  cat("No sites with coordinate data found - skipping site boxplots\n")
  pad_site_plot <- NULL
}

# ============================================================================
# SECTION 6: STATISTICAL SUMMARIES
# ============================================================================

cat("\n=== SECTION 6: STATISTICAL SUMMARIES ===\n")

# ANOVA for habitat differences in PAD
if(length(unique(habitat_data$nlcdClass)) > 1) {
  cat("Testing for PAD habitat differences (ANOVA):\n")
  
  pad_anova <- aov(PAD ~ nlcdClass, data = habitat_data)
  pad_anova_summary <- summary(pad_anova)
  print(pad_anova_summary)
  
  pad_anova_p <- pad_anova_summary[[1]][["Pr(>F)"]][1]
  cat("  PAD ANOVA p-value:", format.pval(pad_anova_p, digits = 3), "\n")
  
  if(pad_anova_p < 0.05) {
    cat("  Significant PAD habitat differences detected\n")
    
    # Post-hoc test if significant
    if(length(unique(habitat_data$nlcdClass)) > 2) {
      cat("\nPAD post-hoc pairwise comparisons:\n")
      pad_posthoc <- TukeyHSD(pad_anova)
      print(pad_posthoc)
    }
  } else {
    cat("  No significant PAD habitat differences\n")
  }
  
  # ANOVA for nymph peak day (if available)
  if(sum(!is.na(habitat_data$nymph_peak_day)) > 10) {
    cat("\nTesting for nymph peak day habitat differences (ANOVA):\n")
    
    nymph_data <- habitat_data %>% filter(!is.na(nymph_peak_day))
    nymph_anova <- aov(nymph_peak_day ~ nlcdClass, data = nymph_data)
    nymph_anova_summary <- summary(nymph_anova)
    print(nymph_anova_summary)
    
    nymph_anova_p <- nymph_anova_summary[[1]][["Pr(>F)"]][1]
    cat("  Nymph peak ANOVA p-value:", format.pval(nymph_anova_p, digits = 3), "\n")
  }
  
  # ANOVA for larva peak day (if available)
  if(sum(!is.na(habitat_data$larva_peak_day)) > 10) {
    cat("\nTesting for larva peak day habitat differences (ANOVA):\n")
    
    larva_data <- habitat_data %>% filter(!is.na(larva_peak_day))
    larva_anova <- aov(larva_peak_day ~ nlcdClass, data = larva_data)
    larva_anova_summary <- summary(larva_anova)
    print(larva_anova_summary)
    
    larva_anova_p <- larva_anova_summary[[1]][["Pr(>F)"]][1]
    cat("  Larva peak ANOVA p-value:", format.pval(larva_anova_p, digits = 3), "\n")
  }
  
} else {
  cat("Only one habitat type present - no statistical tests possible\n")
}

# Summary statistics by habitat
cat("\nPhenology statistics by habitat:\n")
phenology_stats <- habitat_data %>%
  group_by(nlcdClass) %>%
  summarise(
    N = n(),
    PAD_mean = round(mean(PAD, na.rm = TRUE), 2),
    PAD_sd = round(sd(PAD, na.rm = TRUE), 2),
    Nymph_peak_mean = round(mean(nymph_peak_day, na.rm = TRUE), 1),
    Larva_peak_mean = round(mean(larva_peak_day, na.rm = TRUE), 1),
    .groups = "drop"
  )
print(phenology_stats)

# ============================================================================
# SECTION 7: SAVE OUTPUTS
# ============================================================================

cat("\n=== SECTION 7: SAVING OUTPUTS ===\n")

# Create output directory
if (!dir.exists("output")) {
  dir.create("output", recursive = TRUE)
}

# Save PAD habitat plot
ggsave("output/habitat_PAD_phenology.png", pad_habitat_plot, 
       width = 10, height = 6, dpi = 300)
cat("✓ PAD habitat plot saved: output/habitat_PAD_phenology.png\n")

# Save PAD histogram
ggsave("output/PAD_histogram_phenology.png", pad_histogram, 
       width = 8, height = 6, dpi = 300)
cat("✓ PAD histogram saved: output/PAD_histogram_phenology.png\n")

# Save peak day habitat plots (if available)
if(!is.null(nymph_habitat_plot)) {
  ggsave("output/habitat_nymph_peak_phenology.png", nymph_habitat_plot, 
         width = 10, height = 6, dpi = 300)
  cat("✓ Nymph peak habitat plot saved: output/habitat_nymph_peak_phenology.png\n")
}

if(!is.null(larva_habitat_plot)) {
  ggsave("output/habitat_larva_peak_phenology.png", larva_habitat_plot, 
         width = 10, height = 6, dpi = 300)
  cat("✓ Larva peak habitat plot saved: output/habitat_larva_peak_phenology.png\n")
}

# Save site plot (if available)
if(!is.null(pad_site_plot)) {
  n_sites <- length(unique(sites_with_coords$siteID))
  plot_width <- max(12, n_sites * 0.4)  # Scale with number of sites
  
  ggsave("output/site_PAD_phenology.png", pad_site_plot, 
         width = plot_width, height = 6, dpi = 300)
  cat("✓ PAD site plot saved: output/site_PAD_phenology.png\n")
}

# Save data summaries
write.csv(habitat_summary, "output/habitat_summary_phenology.csv", row.names = FALSE)
write.csv(site_summary, "output/site_summary_phenology.csv", row.names = FALSE)
write.csv(phenology_stats, "output/phenology_stats_by_habitat.csv", row.names = FALSE)

# Save site coordinates for reference
if(exists("site_coords")) {
  write.csv(site_coords, "output/site_coordinates_phenology.csv", row.names = FALSE)
  cat("✓ Site coordinates saved: output/site_coordinates_phenology.csv\n")
}

cat("✓ Data summaries saved: output/*_summary_phenology.csv\n")

# Save statistical results
stat_results <- list(
  pad_anova = if(exists("pad_anova")) pad_anova else NULL,
  pad_anova_summary = if(exists("pad_anova_summary")) pad_anova_summary else NULL,
  pad_anova_p = if(exists("pad_anova_p")) pad_anova_p else NULL,
  nymph_anova = if(exists("nymph_anova")) nymph_anova else NULL,
  larva_anova = if(exists("larva_anova")) larva_anova else NULL,
  pad_posthoc = if(exists("pad_posthoc")) pad_posthoc else NULL,
  phenology_stats = phenology_stats,
  site_coords = if(exists("site_coords")) site_coords else NULL,
  response_summary = response_summary
)

saveRDS(stat_results, "output/habitat_statistics_phenology.RDS")
cat("✓ Statistical results saved: output/habitat_statistics_phenology.RDS\n")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n", rep("=", 60), "\n")
cat("HABITAT AND SITE PLOTS: PHENOLOGY ANALYSIS COMPLETE\n")
cat(rep("=", 60), "\n")

cat("PLOTS CREATED:\n")
cat("✓ PAD habitat boxplot: PAD by nlcdClass\n")
cat("✓ PAD distribution histogram: Overall PAD patterns\n")
if(!is.null(nymph_habitat_plot)) {
  cat("✓ Nymph peak habitat boxplot: nymph_peak_day by nlcdClass\n")
} else {
  cat("✗ Nymph peak habitat boxplot: No peak day data available\n")
}
if(!is.null(larva_habitat_plot)) {
  cat("✓ Larva peak habitat boxplot: larva_peak_day by nlcdClass\n")
} else {
  cat("✗ Larva peak habitat boxplot: No peak day data available\n")
}
if(!is.null(pad_site_plot)) {
  cat("✓ PAD site boxplot: PAD by siteID (latitude-ordered)\n")
} else {
  cat("✗ PAD site boxplot: No coordinate data available\n")
}

cat("\nENHANCEMENTS:\n")
cat("✓ Direct load from pad_analysis_dataset.RDS\n")
cat("✓ All sites included (no observation threshold)\n")
cat("✓ All PAD values included (including negative PADs)\n")
if(!is.null(pad_site_plot)) {
  cat("✓ Sites ordered by latitude (south → north)\n")
  cat("✓ Latitude range:", round(min(site_coords$latitude, na.rm = TRUE), 2), 
      "to", round(max(site_coords$latitude, na.rm = TRUE), 2), "\n")
}

cat("\nDATA SUMMARY:\n")
cat("✓ Total observations:", nrow(viz_data), "(including negative PADs)\n")
cat("✓ Observations with habitat:", n_with_habitat, "\n")
cat("✓ Observations with peak days:", n_with_peaks, "\n")
cat("✓ Observations with coordinates:", n_with_coords, "\n")
cat("✓ Unique habitats:", length(unique(habitat_data$nlcdClass)), "\n")
cat("✓ Unique sites:", length(unique(viz_data$siteID)), "\n")

cat("\nPAD BREAKDOWN:\n")
if(exists("response_summary")) {
  cat("✓ Negative PADs:", response_summary$PAD_negative_count, 
      paste0("(", round(response_summary$PAD_negative_count/total_pad_obs*100, 1), "%)"), "\n")
  cat("✓ Positive PADs:", response_summary$PAD_positive_count, 
      paste0("(", round(response_summary$PAD_positive_count/total_pad_obs*100, 1), "%)"), "\n")
}

if(exists("pad_anova_p")) {
  cat("\nSTATISTICAL RESULTS:\n")
  cat("✓ PAD habitat ANOVA p-value:", format.pval(pad_anova_p, digits = 3), "\n")
  cat("✓ Significant PAD habitat effect:", ifelse(pad_anova_p < 0.05, "YES", "NO"), "\n")
}

cat("\nOUTPUTS SAVED:\n")
cat("✓ Plots: output/habitat_*_phenology.png, output/site_*_phenology.png\n")
cat("         output/PAD_histogram_phenology.png\n")
cat("✓ Summaries: output/*_summary_phenology.csv\n")
cat("✓ Coordinates: output/site_coordinates_phenology.csv\n")
cat("✓ Statistics: output/habitat_statistics_phenology.RDS\n")

cat("\nINTERPRETATION:\n")
cat("→ PAD = 0: Perfect synchrony between nymph and larva peaks\n")
cat("→ PAD > 0: Larva peak later than nymph peak\n")
cat("→ PAD < 0: Larva peak earlier than nymph peak\n")
cat("→ Sites ordered south to north by latitude\n")

cat("\n", rep("=", 60), "\n")