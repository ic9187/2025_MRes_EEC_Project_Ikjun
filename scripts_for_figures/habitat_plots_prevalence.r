# ============================================================================
# habitat_site_plots_prevalence.R
# Purpose: Create habitat and site boxplots for prevalence analysis
# Input: model_data from 05 prevalence analysis pipeline
# Output: Boxplots showing transmission_log patterns by habitat and site
# ENHANCED: Show all sites ordered by latitude
# ============================================================================

library(dplyr)
library(ggplot2)
library(gridExtra)

cat("=== HABITAT AND SITE PLOTS: PREVALENCE ANALYSIS ===\n")
cat("Creating boxplots for transmission_log patterns\n")
cat("ENHANCED: All sites shown, ordered by latitude\n\n")

# ============================================================================
# SECTION 1: LOAD DATA
# ============================================================================

cat("=== SECTION 1: LOADING DATA ===\n")

# Load model_data from prevalence analysis
if(!exists("model_data")) {
  # Try to load from 05a output first (updated with habitat data)
  if(file.exists("output/05a_model_data.RDS")) {
    cat("Loading data from 05a_model_data.RDS...\n")
    model_data <- readRDS("output/05a_model_data.RDS")
    
    if(!"nlcdClass" %in% names(model_data)) {
      stop("ERROR: nlcdClass column not found in 05a_model_data.RDS\n",
           "SOLUTION: Please run the updated 05a_setup_log.R script first")
    }
    
    cat("✓ Loaded model_data from 05a output\n")
    
  } else {
    # Fallback to yearly_transmission
    if(file.exists("data/processed/yearly_transmission.RDS")) {
      cat("Loading data from yearly_transmission.RDS...\n")
      full_data <- readRDS("data/processed/yearly_transmission.RDS")
      
      # Check if nlcdClass column exists
      if(!"nlcdClass" %in% names(full_data)) {
        stop("ERROR: nlcdClass column not found in yearly_transmission.RDS\n",
             "SOLUTION: Please run the updated 03c_merge_prevalence_data.R script first\n",
             "This will add habitat data (nlcdClass) to your dataset.")
      }
      
      cat("✓ nlcdClass column found in data\n")
      
      # Recreate model_data filtering (matching 05a setup)
      model_data <- full_data %>%
        dplyr::select(
          site_year, siteID, year,
          transmission_log,
          PAD,
          nymph_max,
          larva_max,
          mean_density_ha,
          simpson_diversity, 
          mean_annual_temp,
          mean_annual_humidity,
          total_annual_precip,
          nlcdClass,  # Include habitat data
          lat, long   # Include coordinates for latitude ordering
        ) %>%
        # Remove any rows with missing values in core variables
        filter(!is.na(transmission_log), !is.na(PAD), !is.na(mean_density_ha), 
               !is.na(simpson_diversity), !is.na(mean_annual_temp),
               !is.na(mean_annual_humidity), !is.na(total_annual_precip))
      
      cat("✓ Recreated model_data from yearly_transmission\n")
    } else {
      stop("ERROR: Neither 05a_model_data.RDS nor yearly_transmission.RDS found\n",
           "SOLUTION: Please run scripts 05a_setup_log.R first")
    }
  }
} else {
  cat("✓ Using existing model_data\n")
  
  # Check if nlcdClass exists in existing model_data
  if(!"nlcdClass" %in% names(model_data)) {
    stop("ERROR: nlcdClass column not found in existing model_data\n",
         "SOLUTION: Please run the updated 05a_setup_log.R script first")
  }
}

cat("  Observations:", nrow(model_data), "\n")
cat("  Sites:", length(unique(model_data$siteID)), "\n")

# Check habitat data availability
n_with_habitat <- sum(!is.na(model_data$nlcdClass))
cat("  Observations with habitat data:", n_with_habitat, "\n")

# Check coordinate data for latitude ordering
n_with_coords <- sum(!is.na(model_data$lat) & !is.na(model_data$long))
cat("  Observations with coordinates:", n_with_coords, "\n")

# More detailed habitat data check
if("nlcdClass" %in% names(model_data)) {
  habitat_values <- unique(model_data$nlcdClass)
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
       "1. The updated 05a_setup_log.R script hasn't been run\n",
       "2. No habitat data exists for your sites in tick_raw\n",
       "3. Data filtering removed all records with habitat data\n",
       "\nSOLUTION:\n",
       "1. Run updated 05a_setup_log.R first\n",
       "2. Check that tick_raw contains nlcdClass data\n",
       "3. Verify your site IDs match between datasets")
}

# ============================================================================
# SECTION 2: DATA SUMMARY
# ============================================================================

cat("\n=== SECTION 2: DATA SUMMARY ===\n")

# Response variable summary
response_summary <- model_data %>%
  summarise(
    Variable = "transmission_log",
    N = n(),
    Mean = mean(transmission_log),
    SD = sd(transmission_log),
    Min = min(transmission_log),
    Max = max(transmission_log),
    Q25 = quantile(transmission_log, 0.25),
    Median = median(transmission_log),
    Q75 = quantile(transmission_log, 0.75)
  )

cat("Response Variable Summary:\n")
print(response_summary)

# Habitat distribution
habitat_summary <- model_data %>%
  filter(!is.na(nlcdClass)) %>%
  group_by(nlcdClass) %>%
  summarise(
    N_observations = n(),
    N_sites = n_distinct(siteID),
    Mean_transmission = mean(transmission_log),
    SD_transmission = sd(transmission_log),
    .groups = "drop"
  ) %>%
  arrange(desc(N_observations))

cat("\nHabitat Distribution:\n")
print(habitat_summary)

# Site distribution
site_summary <- model_data %>%
  group_by(siteID) %>%
  summarise(
    N_observations = n(),
    Mean_transmission = mean(transmission_log),
    SD_transmission = sd(transmission_log),
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
habitat_data <- model_data %>%
  filter(!is.na(nlcdClass))

cat("Creating habitat boxplot with", nrow(habitat_data), "observations\n")

# Create habitat boxplot
habitat_plot <- ggplot(habitat_data, aes(x = nlcdClass, y = transmission_log)) +
  geom_boxplot(aes(fill = nlcdClass), alpha = 0.7, outlier.alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  
  # Styling
  scale_fill_viridis_d(name = "Habitat Type") +
  
  # Labels
  labs(
    title = "Pathogen Transmission by Habitat Type",
    subtitle = "Log(transmission ratio) across different land cover classes",
    x = "Habitat Type (NLCD Class)",
    y = "Log(Transmission Ratio)",
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
    legend.position = "none",  # Remove legend since x-axis labels are clear
    panel.grid.minor = element_blank()
  ) +
  
  # Add horizontal line at zero (no change in transmission)
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7)

# Add sample sizes to plot
habitat_counts <- habitat_data %>%
  count(nlcdClass) %>%
  mutate(label = paste0("n=", n))

habitat_plot <- habitat_plot +
  geom_text(data = habitat_counts, 
            aes(x = nlcdClass, y = max(habitat_data$transmission_log) * 0.9, 
                label = label),
            size = 3, color = "black", fontface = "bold")

print(habitat_plot)

# ============================================================================
# SECTION 4: SITE BOXPLOTS (ALL SITES, ORDERED BY LATITUDE)
# ============================================================================

cat("\n=== SECTION 4: CREATING SITE BOXPLOTS ===\n")
cat("ENHANCED: Showing ALL sites, ordered by latitude (south to north)\n")

# Create site coordinate summary for ordering
site_coords <- model_data %>%
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
cat("Latitude range:", round(min(site_coords$latitude, na.rm = TRUE), 2), "to", 
    round(max(site_coords$latitude, na.rm = TRUE), 2), "\n")

# Filter model_data to sites with coordinates and create ordered factor
sites_with_coords <- model_data %>%
  filter(siteID %in% site_coords$siteID) %>%
  mutate(siteID = factor(siteID, levels = site_coords$siteID))

cat("Creating site boxplot with", nrow(sites_with_coords), "observations from", 
    length(unique(sites_with_coords$siteID)), "sites\n")

if(nrow(sites_with_coords) > 0) {
  
  # Create site boxplot
  site_plot <- ggplot(sites_with_coords, aes(x = siteID, y = transmission_log)) +
    geom_boxplot(aes(fill = nlcdClass), alpha = 0.7, outlier.alpha = 0.6) +
    geom_jitter(width = 0.3, alpha = 0.5, size = 0.8) +
    
    # Styling
    scale_fill_viridis_d(name = "Habitat Type") +
    
    # Labels
    labs(
      title = "Pathogen Transmission by Site",
      subtitle = "Log(transmission ratio) across NEON sites (ordered by latitude: south → north)",
      x = "Site ID (Ordered by Latitude)",
      y = "Log(Transmission Ratio)",
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
  
  print(site_plot)
  
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
  cat("No sites with coordinate data found - skipping site boxplot\n")
  site_plot <- NULL
}

# ============================================================================
# SECTION 5: STATISTICAL SUMMARIES
# ============================================================================

cat("\n=== SECTION 5: STATISTICAL SUMMARIES ===\n")

# ANOVA for habitat differences
if(length(unique(habitat_data$nlcdClass)) > 1) {
  cat("Testing for habitat differences (ANOVA):\n")
  
  habitat_anova <- aov(transmission_log ~ nlcdClass, data = habitat_data)
  anova_summary <- summary(habitat_anova)
  print(anova_summary)
  
  anova_p <- anova_summary[[1]][["Pr(>F)"]][1]
  cat("  ANOVA p-value:", format.pval(anova_p, digits = 3), "\n")
  
  if(anova_p < 0.05) {
    cat("  Significant habitat differences detected\n")
    
    # Post-hoc test if significant
    if(length(unique(habitat_data$nlcdClass)) > 2) {
      cat("\nPost-hoc pairwise comparisons:\n")
      posthoc <- TukeyHSD(habitat_anova)
      print(posthoc)
    }
  } else {
    cat("  No significant habitat differences\n")
  }
} else {
  cat("Only one habitat type present - no statistical test possible\n")
}

# Summary statistics by habitat
cat("\nTransmission statistics by habitat:\n")
habitat_stats <- habitat_data %>%
  group_by(nlcdClass) %>%
  summarise(
    N = n(),
    Mean = round(mean(transmission_log), 3),
    SD = round(sd(transmission_log), 3),
    Median = round(median(transmission_log), 3),
    Min = round(min(transmission_log), 3),
    Max = round(max(transmission_log), 3),
    .groups = "drop"
  )
print(habitat_stats)

# ============================================================================
# SECTION 6: SAVE OUTPUTS
# ============================================================================

cat("\n=== SECTION 6: SAVING OUTPUTS ===\n")

# Create output directory
if (!dir.exists("output")) {
  dir.create("output", recursive = TRUE)
}

# Save plots
ggsave("output/habitat_transmission_prevalence.png", habitat_plot, 
       width = 10, height = 6, dpi = 300)
cat("✓ Habitat plot saved: output/habitat_transmission_prevalence.png\n")

if(!is.null(site_plot)) {
  # Adjust width based on number of sites
  n_sites <- length(unique(sites_with_coords$siteID))
  plot_width <- max(12, n_sites * 0.4)  # Scale with number of sites
  
  ggsave("output/site_transmission_prevalence.png", site_plot, 
         width = plot_width, height = 6, dpi = 300)
  cat("✓ Site plot saved: output/site_transmission_prevalence.png\n")
}

# Save data summaries
write.csv(habitat_summary, "output/habitat_summary_prevalence.csv", row.names = FALSE)
write.csv(site_summary, "output/site_summary_prevalence.csv", row.names = FALSE)

# Save site coordinates for reference
if(exists("site_coords")) {
  write.csv(site_coords, "output/site_coordinates_prevalence.csv", row.names = FALSE)
  cat("✓ Site coordinates saved: output/site_coordinates_prevalence.csv\n")
}

cat("✓ Data summaries saved: output/*_summary_prevalence.csv\n")

# Save statistical results
if(exists("anova_summary")) {
  stat_results <- list(
    habitat_anova = habitat_anova,
    anova_summary = anova_summary,
    anova_p_value = anova_p,
    habitat_stats = habitat_stats,
    posthoc = if(exists("posthoc")) posthoc else NULL,
    site_coords = if(exists("site_coords")) site_coords else NULL
  )
  
  saveRDS(stat_results, "output/habitat_statistics_prevalence.RDS")
  cat("✓ Statistical results saved: output/habitat_statistics_prevalence.RDS\n")
}

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n", rep("=", 60), "\n")
cat("HABITAT AND SITE PLOTS: PREVALENCE ANALYSIS COMPLETE\n")
cat(rep("=", 60), "\n")

cat("PLOTS CREATED:\n")
cat("✓ Habitat boxplot: transmission_log by nlcdClass\n")
if(!is.null(site_plot)) {
  cat("✓ Site boxplot: transmission_log by siteID (latitude-ordered)\n")
} else {
  cat("✗ Site boxplot: No coordinate data available\n")
}

cat("\nENHANCEMENTS:\n")
cat("✓ All sites included (no observation threshold)\n")
if(!is.null(site_plot)) {
  cat("✓ Sites ordered by latitude (south → north)\n")
  cat("✓ Latitude range:", round(min(site_coords$latitude, na.rm = TRUE), 2), 
      "to", round(max(site_coords$latitude, na.rm = TRUE), 2), "\n")
}

cat("\nDATA SUMMARY:\n")
cat("✓ Total observations:", nrow(model_data), "\n")
cat("✓ Observations with habitat:", n_with_habitat, "\n")
cat("✓ Observations with coordinates:", n_with_coords, "\n")
cat("✓ Unique habitats:", length(unique(habitat_data$nlcdClass)), "\n")
cat("✓ Unique sites:", length(unique(model_data$siteID)), "\n")

if(exists("anova_p")) {
  cat("\nSTATISTICAL RESULTS:\n")
  cat("✓ Habitat ANOVA p-value:", format.pval(anova_p, digits = 3), "\n")
  cat("✓ Significant habitat effect:", ifelse(anova_p < 0.05, "YES", "NO"), "\n")
}

cat("\nOUTPUTS SAVED:\n")
cat("✓ Plots: output/habitat_transmission_prevalence.png\n")
if(!is.null(site_plot)) {
  cat("         output/site_transmission_prevalence.png\n")
}
cat("✓ Summaries: output/*_summary_prevalence.csv\n")
cat("✓ Coordinates: output/site_coordinates_prevalence.csv\n")
if(exists("stat_results")) {
  cat("✓ Statistics: output/habitat_statistics_prevalence.RDS\n")
}

cat("\nINTERPRETATION:\n")
cat("→ Positive values indicate increased transmission (year-to-year)\n")
cat("→ Negative values indicate decreased transmission\n")
cat("→ Zero line represents no change in transmission\n")
cat("→ Sites ordered south to north by latitude\n")

cat("\n", rep("=", 60), "\n")