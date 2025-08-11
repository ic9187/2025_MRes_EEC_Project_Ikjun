# ============================================================================
# UKFS TICK INDIVIDUAL SAMPLING EVENTS PLOT
# ============================================================================
# Plot individual sampling events (not aggregated by day) for UKFS
# Shows log(tick count + 1) vs day of year for each individual sampling event
# Uses field count data (all ticks, not genus-specific)
# Follows same filtering pipeline as main analysis scripts
# ============================================================================

library(dplyr)
library(ggplot2)

cat("Creating UKFS tick sampling events plot...\n")

# ---- Step 1: Load clean site_years filter ----
clean_site_years <- readRDS("data/processed/clean_site_years.RDS")
cat("Loaded", length(clean_site_years), "clean site_years\n")

# ---- Step 2: Apply same filtering as pipeline ----
tick_filtered <- tick_raw %>%
  # Filter by samplingImpractical
  filter(samplingImpractical == "OK" | 
           is.na(samplingImpractical) | 
           samplingImpractical == "") %>%
  # Remove rows with all NA counts
  filter(!(is.na(nymphCount) & is.na(larvaCount) & is.na(adultCount))) %>%
  # Remove 2019 data (protocol changes)
  filter(year != 2019) %>%
  # Keep only samples with no compromise in condition
  filter(sampleCondition == "No known compromise" | 
           is.na(sampleCondition) | 
           sampleCondition == "") %>%
  # Add site_year variable
  mutate(site_year = paste(siteID, year, sep = "_")) %>%
  # Apply clean_site_years filter
  filter(site_year %in% clean_site_years)

cat("Filtered to", nrow(tick_filtered), "sampling events\n")

# ---- Step 3: Filter for UKFS only ----
ukfs_events <- tick_filtered %>%
  filter(siteID == "UKFS") %>%
  # Create log transformed counts for individual events
  mutate(
    log_nymph = log(nymphCount + 1),
    log_larva = log(larvaCount + 1),
    log_adult = log(adultCount + 1)
  ) %>%
  # Keep only events with at least some ticks
  filter(nymphCount > 0 | larvaCount > 0 | adultCount > 0)

cat("Found", nrow(ukfs_events), "UKFS sampling events with ticks\n")
cat("Years:", paste(sort(unique(ukfs_events$year)), collapse = ", "), "\n")
cat("Total nymphs:", sum(ukfs_events$nymphCount, na.rm = TRUE), "\n")
cat("Total larvae:", sum(ukfs_events$larvaCount, na.rm = TRUE), "\n")
cat("Total adults:", sum(ukfs_events$adultCount, na.rm = TRUE), "\n")

# ---- Step 4: Find peak days and values for annotations ----
peak_nymph_idx <- which.max(ukfs_events$nymphCount)
peak_larva_idx <- which.max(ukfs_events$larvaCount)

peak_nymph_day <- ukfs_events$day_number[peak_nymph_idx]
peak_larva_day <- ukfs_events$day_number[peak_larva_idx]
peak_nymph_log_value <- ukfs_events$log_nymph[peak_nymph_idx]
peak_larva_log_value <- ukfs_events$log_larva[peak_larva_idx]

peak_difference <- abs(peak_larva_day - peak_nymph_day)

cat("Peak nymph day:", peak_nymph_day, "| Log value:", round(peak_nymph_log_value, 2), "\n")
cat("Peak larva day:", peak_larva_day, "| Log value:", round(peak_larva_log_value, 2), "\n")
cat("Peak difference:", peak_difference, "days\n")

# ---- Step 5: Create scatter plot of individual events ----
ukfs_scatter <- ggplot(ukfs_events, aes(x = day_number)) +
  # Add points for nymphs and larvae (each individual sampling event)
  geom_point(aes(y = log_nymph, color = "Nymph"), alpha = 0.7, size = 1.5) +
  geom_point(aes(y = log_larva, color = "Larva"), alpha = 0.7, size = 1.5) +
  
  # Add solid vertical lines at peak days
  geom_vline(xintercept = peak_nymph_day, linetype = "solid", color = "#0f62fe", size = 1) +
  geom_vline(xintercept = peak_larva_day, linetype = "solid", color = "#ff832b", size = 1) +
  
  # Add horizontal dotted lines from peak points to y-axis
  geom_hline(yintercept = peak_nymph_log_value, linetype = "dotted", color = "#0f62fe", size = 0.8) +
  geom_hline(yintercept = peak_larva_log_value, linetype = "dotted", color = "#ff832b", size = 0.8) +
  
  # Add double-headed arrow between peak days (narrower)
  geom_segment(
    x = min(peak_nymph_day, peak_larva_day), 
    xend = max(peak_nymph_day, peak_larva_day),
    y = max(ukfs_events$log_nymph, ukfs_events$log_larva, na.rm = TRUE) * 0.9,
    yend = max(ukfs_events$log_nymph, ukfs_events$log_larva, na.rm = TRUE) * 0.9,
    arrow = arrow(ends = "both", length = unit(0.3, "cm")),
    color = "black", size = 0.5
  ) +
  
  # Add annotation for time difference
  annotate(
    "text",
    x = (peak_nymph_day + peak_larva_day) / 2,
    y = max(ukfs_events$log_nymph, ukfs_events$log_larva, na.rm = TRUE) * 0.95,
    label = paste0("Δt"),
    # insert , peak_difference, " days" above to paste the days also
    size = 4, fontface = "bold"
  ) +
  
  # Color scheme
  scale_color_manual(values = c("Nymph" = "#0f62fe", "Larva" = "#ff832b")) +
  
  # Labels and theme
  labs(
    title = "UKFS Tick Individual Sampling Events",
    subtitle = paste("Each point = one sampling event •", nrow(ukfs_events), "events total"),
    x = "Day of Year",
    y = "Log(Count + 1)",
    color = "Life Stage"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray60"),
    legend.position = "right"
  ) +
  
  # Fixed x-axis from 0 to 366
  xlim(0, 366)

# Display the plot
print(ukfs_scatter)

# ---- Step 6: Summary by year ----
cat("\n=== UKFS TICK EVENTS BY YEAR ===\n")

yearly_summary <- ukfs_events %>%
  group_by(year) %>%
  summarise(
    n_events = n(),
    total_nymphs = sum(nymphCount, na.rm = TRUE),
    total_larvae = sum(larvaCount, na.rm = TRUE),
    first_day = min(day_number),
    last_day = max(day_number),
    sampling_span = last_day - first_day,
    .groups = "drop"
  )

print(yearly_summary)

# ---- Step 7: Tick phenology peaks from individual events ----
cat("\n=== TICK PHENOLOGY FROM INDIVIDUAL EVENTS ===\n")

# Find events with maximum tick counts
max_nymph_event <- ukfs_events[which.max(ukfs_events$nymphCount), ]
max_larva_event <- ukfs_events[which.max(ukfs_events$larvaCount), ]

cat("Peak nymph tick event:\n")
cat("  Day:", max_nymph_event$day_number, "| Count:", max_nymph_event$nymphCount, 
    "| Date:", as.character(max_nymph_event$collectDate), "\n")

cat("Peak larva tick event:\n")
cat("  Day:", max_larva_event$day_number, "| Count:", max_larva_event$larvaCount, 
    "| Date:", as.character(max_larva_event$collectDate), "\n")

cat("Peak difference:", max_larva_event$day_number - max_nymph_event$day_number, "days\n")

cat("\n✓ Individual tick sampling events plot complete!\n")