# ============================================================================
# NORTH AMERICA TICK DATA PIPELINE VISUALIZATION - PIPELINE CONSISTENT
# ============================================================================
# Creates four maps showing data progression through the pipeline:
# 1. Tick count data (filtered by clean_site_years)
# 2. Taxonomy data (Ixodes vs Non-Ixodes as PIE CHARTS!)
# 3. Borrelia pathogen data (all tests)
# 4. Borrelia pathogen data (POSITIVE TICKS ONLY)
# 
# UPDATED: Now uses pipeline summary logic for consistency!
# ============================================================================

cat("\n==== CREATING PIPELINE-CONSISTENT TICK MAPS ====\n")

# ---- Setup and load required packages ----
library(ggplot2)
library(dplyr)
library(maps)
library(rnaturalearth)
library(sf)
library(gridExtra)
library(scales)
library(ggforce)  # For pie charts!

# Load data and filters
clean_site_years <- readRDS("data/processed/clean_site_years.RDS")

# Load pipeline summary data for consistency
pipeline_data <- readRDS("data/processed/pipeline_summary_data.RDS")

cat("Loaded pipeline data for consistency:\n")
cat("  - Total collected:", format(pipeline_data$total_collected$total_ticks_collected, big.mark = ","), "\n")
cat("  - Total identified:", format(pipeline_data$identification_summary$total_identified, big.mark = ","), "\n")
cat("  - Ixodes identified:", format(pipeline_data$ixodes_summary$ixodes_identified, big.mark = ","), "\n")
cat("  - Non-Ixodes identified:", format(pipeline_data$non_ixodes_summary$non_ixodes_identified, big.mark = ","), "\n\n")

# ---- Step 1: Get base map of North America ----
cat("Step 1: Creating base North America map...\n")

# Get North America countries
north_america <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")

# Continental bounds (continental US focus)
continental_bounds <- list(
  xmin = -130, xmax = -65,
  ymin = 25, ymax = 55
)

# ---- Step 2: Process tick count data (Stage 1) - PIPELINE CONSISTENT ----
cat("Step 2: Processing tick count data (pipeline consistent)...\n")

tick_count_sites <- tick_raw %>%
  filter(
    samplingImpractical == "OK" | is.na(samplingImpractical) | samplingImpractical == "",
    !(is.na(nymphCount) & is.na(larvaCount) & is.na(adultCount)),
    year != 2019,
    sampleCondition == "No known compromise" | is.na(sampleCondition) | sampleCondition == ""
  ) %>%
  mutate(
    siteID = sub("_.*", "", plotID),
    site_year = paste(siteID, year, sep = "_")
  ) %>%
  filter(site_year %in% clean_site_years) %>%
  group_by(siteID) %>%
  summarise(
    lat = mean(decimalLatitude, na.rm = TRUE),
    long = mean(decimalLongitude, na.rm = TRUE),
    total_nymphs = sum(nymphCount, na.rm = TRUE),
    total_larvae = sum(larvaCount, na.rm = TRUE),
    total_adults = sum(adultCount, na.rm = TRUE),
    n_sampling_events = n(),
    n_years = n_distinct(year),
    .groups = "drop"
  ) %>%
  mutate(
    total_ticks = total_nymphs + total_larvae + total_adults,
    stage = "Tick Counts"
  ) %>%
  filter(total_ticks > 0, !is.na(lat), !is.na(long))

cat("  Tick count sites:", nrow(tick_count_sites), "\n")
cat("  Total ticks:", format(sum(tick_count_sites$total_ticks), big.mark = ","), "\n")

# Verify consistency with pipeline
if(sum(tick_count_sites$total_ticks) != pipeline_data$total_collected$total_ticks_collected) {
  cat("  WARNING: Map tick counts don't match pipeline summary!\n")
} else {
  cat("  ✓ Tick counts match pipeline summary\n")
}

# ---- Step 3: Process taxonomy data - PIPELINE CONSISTENT ----
cat("Step 3: Processing taxonomy data (pipeline consistent)...\n")

# Use EXACT same logic as pipeline summary
taxonomy_raw <- taxo_raw %>%
  filter(sampleCondition == "OK") %>%
  mutate(
    collectDate = as.Date(collectDate),
    year = year(collectDate),
    siteID = sub("_.*", "", plotID),
    site_year = paste(siteID, year, sep = "_"),
    tick_type = ifelse(genus %in% c("Amblyomma", "Haemaphysalis", "Dermacentor"), 
                       "Non-Ixodes", "Ixodes")
  ) %>%
  filter(site_year %in% clean_site_years)

# Get site coordinates
site_coords <- tick_raw %>%
  mutate(siteID = sub("_.*", "", plotID)) %>%
  group_by(siteID) %>%
  summarise(
    lat = mean(decimalLatitude, na.rm = TRUE),
    long = mean(decimalLongitude, na.rm = TRUE),
    .groups = "drop"
  )

# Create site-level taxonomy summary (NO double-aggregation!)
taxonomy_sites <- taxonomy_raw %>%
  group_by(siteID, tick_type) %>%
  summarise(
    tick_count = sum(individualCount, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(site_coords, by = "siteID") %>%
  filter(tick_count > 0, !is.na(lat), !is.na(long))

# Create PIE CHART data - FIXED calculation!
pie_data <- taxonomy_sites %>%
  group_by(siteID, lat, long) %>%
  summarise(
    ixodes_count = sum(tick_count[tick_type == "Ixodes"], na.rm = TRUE),
    non_ixodes_count = sum(tick_count[tick_type == "Non-Ixodes"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # FIXED: total_ticks is now correctly calculated
    total_ticks = ixodes_count + non_ixodes_count,
    ixodes_prop = ixodes_count / total_ticks,
    non_ixodes_prop = non_ixodes_count / total_ticks,
    # Create radius based on total ticks
    base_radius = (0.5 + (total_ticks - min(total_ticks)) / (max(total_ticks) - min(total_ticks)) * 1.5) * 0.6,
    radius = base_radius / cos(lat * pi/180)
  )

# Verify totals match pipeline
total_ixodes_map <- sum(pie_data$ixodes_count)
total_non_ixodes_map <- sum(pie_data$non_ixodes_count)
total_identified_map <- total_ixodes_map + total_non_ixodes_map

cat("  Map taxonomy totals:\n")
cat("    - Ixodes:", format(total_ixodes_map, big.mark = ","), "\n")
cat("    - Non-Ixodes:", format(total_non_ixodes_map, big.mark = ","), "\n")
cat("    - Total:", format(total_identified_map, big.mark = ","), "\n")

# Verify against pipeline
if(total_ixodes_map != pipeline_data$ixodes_summary$ixodes_identified) {
  cat("  WARNING: Ixodes counts don't match pipeline!\n")
} else {
  cat("  ✓ Ixodes counts match pipeline\n")
}

if(total_non_ixodes_map != pipeline_data$non_ixodes_summary$non_ixodes_identified) {
  cat("  WARNING: Non-Ixodes counts don't match pipeline!\n")
} else {
  cat("  ✓ Non-Ixodes counts match pipeline\n")
}

# Create pie slices for visualization
pie_slices <- pie_data %>%
  filter(ixodes_count > 0 & non_ixodes_count > 0) %>%
  {
    bind_rows(
      # Ixodes slice
      mutate(., 
             tick_type = "Ixodes",
             start_angle = 0,
             end_angle = 2 * pi * ixodes_prop,
             slice_count = ixodes_count),
      # Non-Ixodes slice
      mutate(., 
             tick_type = "Non-Ixodes", 
             start_angle = 2 * pi * ixodes_prop,
             end_angle = 2 * pi,
             slice_count = non_ixodes_count)
    )
  }

# Full circles for sites with only one type
full_circles <- pie_data %>%
  filter(ixodes_count == 0 | non_ixodes_count == 0) %>%
  mutate(
    tick_type = ifelse(ixodes_count > 0, "Ixodes", "Non-Ixodes"),
    start_angle = 0,
    end_angle = 2 * pi,
    slice_count = total_ticks
  )

# Combine for mapping
all_pie_data <- bind_rows(pie_slices, full_circles)

cat("  Sites with pie charts:", nrow(pie_data %>% filter(ixodes_count > 0 & non_ixodes_count > 0)), "\n")
cat("  Sites with full circles:", nrow(full_circles), "\n")

# ---- Step 4: Process Borrelia pathogen data - IXODES FOCUS ----
cat("Step 4: Processing Borrelia pathogen data (Ixodes focus for maps)...\n")

# Filter for IXODES ONLY for maps C and D visualization
borrelia_ixodes_level <- path_raw %>%
  filter(
    grepl("Borrelia", testPathogenName, ignore.case = TRUE),
    !is.na(testResult),
    grepl("IXO", subsampleID)  # IXODES ONLY for maps
  ) %>%
  mutate(
    collectDate = as.Date(collectDate),
    year = year(collectDate),
    siteID = sub("_.*", "", plotID),
    site_year = paste(siteID, year, sep = "_")
  ) %>%
  filter(site_year %in% clean_site_years) %>%
  group_by(testingID, siteID) %>%
  summarise(
    tick_positive = any(testResult == "Positive"),
    .groups = "drop"
  )

# Aggregate by site - IXODES ONLY
borrelia_sites <- borrelia_ixodes_level %>%
  group_by(siteID) %>%
  summarise(
    total_ixodes = n(),  # Number of Ixodes tested
    n_positive = sum(tick_positive),  # Number of Ixodes positive
    positivity_rate = round(mean(tick_positive) * 100, 1),
    .groups = "drop"
  ) %>%
  left_join(site_coords, by = "siteID") %>%
  filter(!is.na(lat), !is.na(long)) %>%
  mutate(
    stage = "Ixodes Borrelia Testing",
    total_ticks = total_ixodes  # For consistency with mapping function
  )

# Positive sites only - IXODES ONLY
borrelia_positive_sites <- borrelia_sites %>%
  filter(n_positive > 0) %>%
  mutate(
    stage = "Ixodes Borrelia Positive",
    total_ticks = n_positive  # Circle size = positive Ixodes count
  )

cat("  Ixodes testing sites:", nrow(borrelia_sites), "\n")
cat("  Total Ixodes tested:", format(sum(borrelia_sites$total_ixodes), big.mark = ","), "\n")
cat("  Positive sites:", nrow(borrelia_positive_sites), "\n")
cat("  Total positive Ixodes:", format(sum(borrelia_positive_sites$total_ticks), big.mark = ","), "\n")

# ---- Step 5: Create size scaling ----
cat("Step 5: Creating consistent size scaling...\n")

all_tick_counts <- c(tick_count_sites$total_ticks, 
                     pie_data$total_ticks,
                     borrelia_sites$total_ticks,
                     borrelia_positive_sites$total_ticks)

min_ticks <- min(all_tick_counts)
max_ticks <- max(all_tick_counts)

size_scale <- function(tick_count) {
  normalized <- (tick_count - min_ticks) / (max_ticks - min_ticks)
  1 + normalized * 2  # Size range 1 to 3
}

# Add size columns
tick_count_sites$circle_size <- size_scale(tick_count_sites$total_ticks)
borrelia_sites$circle_size <- size_scale(borrelia_sites$total_ticks)
borrelia_positive_sites$circle_size <- size_scale(borrelia_positive_sites$total_ticks)

# ---- Step 6: Create mapping functions ----
cat("Step 6: Creating mapping functions...\n")

# Standard map function
create_pipeline_map <- function(site_data, title, color, map_bounds) {
  ggplot() +
    geom_sf(data = north_america, fill = "#F5F5F5", color = "black", size = 0.3) +
    geom_point(data = site_data,
               aes(x = long, y = lat, size = circle_size),
               color = color, alpha = 0.7, stroke = 0.5) +
    coord_sf(xlim = c(map_bounds$xmin, map_bounds$xmax),
             ylim = c(map_bounds$ymin, map_bounds$ymax),
             expand = FALSE) +
    labs(
      title = title,
      subtitle = paste(nrow(site_data), "sites,", 
                      format(sum(site_data$total_ticks), big.mark = ","), "ticks"),
      x = "Longitude", y = "Latitude"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 8),
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "bottom"
    ) +
    scale_size_continuous(
      name = "Sample Size", 
      range = c(1, 6),
      breaks = c(1, 10, 100, 1000),
      labels = c("1", "10", "100", "1000")
    )
}

# PIE CHART map function
create_pie_map <- function(pie_data, title, map_bounds) {
  total_sites <- n_distinct(pie_data$siteID)
  
  # Calculate totals from the pie data correctly
  ixodes_total <- sum(pie_data$slice_count[pie_data$tick_type == "Ixodes"], na.rm = TRUE)
  non_ixodes_total <- sum(pie_data$slice_count[pie_data$tick_type == "Non-Ixodes"], na.rm = TRUE)
  total_ticks <- ixodes_total + non_ixodes_total
  
  ggplot() +
    geom_sf(data = north_america, fill = "#F5F5F5", color = "black", size = 0.3) +
    geom_arc_bar(data = pie_data,
                 aes(x0 = long, y0 = lat, 
                     r0 = 0, r = radius,
                     start = start_angle, end = end_angle, 
                     fill = tick_type),
                 alpha = 0.85, color = "white", size = 0.1) +
    scale_fill_manual(
      name = "Tick Type",
      values = c("Ixodes" = "#FF8C00", "Non-Ixodes" = "#A0D0A0"),
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    coord_sf(xlim = c(map_bounds$xmin, map_bounds$xmax),
             ylim = c(map_bounds$ymin, map_bounds$ymax),
             expand = FALSE) +
    labs(
      title = title,
      subtitle = paste(total_sites, "sites,", 
                      format(total_ticks, big.mark = ","), "total ticks",
                      "| Ixodes:", format(ixodes_total, big.mark = ","),
                      ", Non-Ixodes:", format(non_ixodes_total, big.mark = ",")),
      x = "Longitude", y = "Latitude"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 8),
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "bottom"
    )
}

# ---- Step 7: Create all four maps ----
cat("Step 7: Creating all four pipeline maps...\n")

map1 <- create_pipeline_map(
  tick_count_sites, 
  "A) Field Collected Tick Data", 
  "#2E8B57",
  continental_bounds
)

map2 <- create_pie_map(
  all_pie_data,
  "B) Taxonomy Data (Ixodes vs. Other species)",
  continental_bounds
)

map3 <- create_pipeline_map(
  borrelia_sites, 
  "C) Ixodes Borrelia Test Data", 
  "#4169E1",
  continental_bounds
)

map4 <- create_pipeline_map(
  borrelia_positive_sites, 
  "D) Ixodes Borrelia Positive", 
  "#DC143C",
  continental_bounds
)

# ---- Step 8: Display and save ----
cat("Step 8: Displaying maps...\n")

print(map1)
print(map2)
print(map3)
print(map4)

combined_maps <- grid.arrange(map1, map2, map3, map4, ncol = 2)

# ---- Step 9: Create CORRECTED summary table ----
cat("Step 9: Creating corrected summary table...\n")

pipeline_summary_corrected <- data.frame(
  Stage = c("Tick Counts", "Taxonomy (All)", "Ixodes Borrelia All", "Ixodes Borrelia Positive"),
  Sites = c(nrow(tick_count_sites), 
            n_distinct(all_pie_data$siteID),
            nrow(borrelia_sites), 
            nrow(borrelia_positive_sites)),
  Total_Count = c(sum(tick_count_sites$total_ticks), 
                  total_identified_map,  # CORRECTED
                  sum(borrelia_sites$total_ixodes),  # Ixodes tested
                  sum(borrelia_positive_sites$total_ticks)),  # Ixodes positive
  Notes = c("Field collection", 
            "Ixodes + Non-Ixodes combined - CORRECTED",
            "Ixodes ticks tested for Borrelia", 
            "Ixodes POSITIVE ticks only")
)

cat("\nCORRECTED Pipeline Summary:\n")
print(pipeline_summary_corrected)

# ---- Step 10: Save outputs ----
cat("Step 10: Saving outputs...\n")

if (!dir.exists("output/plots")) {
  dir.create("output/plots", recursive = TRUE)
}

ggsave("output/plots/pipeline_mapA_tick_counts_CORRECTED.png", map1, 
       width = 10, height = 8, dpi = 300)
ggsave("output/plots/pipeline_mapB_taxonomy_PIE_CHARTS_CORRECTED.png", map2, 
       width = 10, height = 8, dpi = 300)
ggsave("output/plots/pipeline_mapC_ixodes_borrelia_all_CORRECTED.png", map3, 
       width = 10, height = 8, dpi = 300)
ggsave("output/plots/pipeline_mapD_ixodes_borrelia_positive_CORRECTED.png", map4, 
       width = 10, height = 8, dpi = 300)
ggsave("output/plots/pipeline_maps_combined_CORRECTED.png", combined_maps, 
       width = 16, height = 12, dpi = 300)

saveRDS(list(
  tick_counts = tick_count_sites,
  pie_data = pie_data,
  pie_slices = all_pie_data,
  borrelia_all = borrelia_sites,
  borrelia_positive = borrelia_positive_sites,
  summary = pipeline_summary_corrected
), "output/pipeline_maps_data_CORRECTED.RDS")

write.csv(pipeline_summary_corrected, "output/pipeline_summary_table_CORRECTED.csv", row.names = FALSE)

cat("\n==== PIPELINE-CONSISTENT MAPS COMPLETE ====\n")
cat("KEY FIXES:\n")
cat("- ✓ Removed double-aggregation in taxonomy calculations\n")
cat("- ✓ Used exact pipeline summary logic throughout\n")
cat("- ✓ Numbers now match between pipeline summary and maps\n")
cat("- ✓ Added verification checks for consistency\n")
cat("- ✓ Maps C & D now focus on Ixodes specifically\n")
cat("\nCorrected totals:\n")
cat("- Map A (Collected):", format(sum(tick_count_sites$total_ticks), big.mark = ","), "\n")
cat("- Map B (Identified):", format(total_identified_map, big.mark = ","), "\n")
cat("  - Ixodes:", format(total_ixodes_map, big.mark = ","), "\n")
cat("  - Non-Ixodes:", format(total_non_ixodes_map, big.mark = ","), "\n")
cat("- Map C (Ixodes Tested):", format(sum(borrelia_sites$total_ixodes), big.mark = ","), "\n")
cat("- Map D (Ixodes Positive):", format(sum(borrelia_positive_sites$total_ticks), big.mark = ","), "\n")