# ============================================================================
# SITE LOCATION MAPS - TICK AND MAMMAL SAMPLING SITES
# ============================================================================
# Creates three simple location maps showing:
# 1. Tick collection sites only
# 2. Mammal trapping sites only  
# 3. Both tick and mammal sites combined
# ============================================================================

cat("\n==== CREATING SITE LOCATION MAPS ====\n")

# ---- Setup and load required packages ----
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(sf)

# ---- Load data ----
# Check if raw data objects exist (from 00_setup.R)
if (!exists("tick_raw") || !exists("trap_raw")) {
  cat("Loading raw data objects...\n")
  
  # Check if data files exist
  if (file.exists("data/raw/tick_data.csv") && file.exists("data/raw/trap_data.csv")) {
    tick_raw <- read.csv("data/raw/tick_data.csv", stringsAsFactors = FALSE)
    trap_raw <- read.csv("data/raw/trap_data.csv", stringsAsFactors = FALSE)
    cat("  ✓ Raw data loaded successfully\n")
  } else {
    stop("Error: Raw data files not found. Please run 00_setup.R first or check file paths.")
  }
} else {
  cat("  ✓ Raw data objects already loaded\n")
}

# Load filters
clean_site_years <- readRDS("data/processed/clean_site_years.RDS")
clean_mammal_site_years <- readRDS("data/processed/clean_mammal_site_years.RDS")

# ---- Get base map of North America ----
cat("Step 1: Creating base North America map...\n")

# Get North America countries
north_america <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")

# Continental bounds (same as existing maps)
continental_bounds <- list(
  xmin = -130, xmax = -65,
  ymin = 25, ymax = 55
)

# ---- Process tick site locations ----
cat("Step 2: Processing tick site locations...\n")

tick_sites <- tick_raw %>%
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
    n_years = n_distinct(year),
    n_sampling_events = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(lat), !is.na(long)) %>%
  mutate(data_type = "Tick")

cat("  Tick sites:", nrow(tick_sites), "\n")

# ---- Process mammal site locations ----
cat("Step 3: Processing mammal site locations...\n")

mammal_sites <- trap_raw %>%
  filter(samplingImpractical == "OK" | is.na(samplingImpractical) | samplingImpractical == "") %>%
  mutate(
    collectDate = as.Date(collectDate),
    year = year(collectDate),
    site_year = paste(siteID, year, sep = "_")
  ) %>%
  filter(site_year %in% clean_mammal_site_years) %>%
  group_by(siteID) %>%
  summarise(
    n_years = n_distinct(year),
    n_sampling_events = n(),
    .groups = "drop"
  ) %>%
  mutate(data_type = "Mammal")

# Get coordinates from tick data (all NEON sites have same coordinates)
site_coords <- tick_raw %>%
  mutate(siteID = sub("_.*", "", plotID)) %>%
  group_by(siteID) %>%
  summarise(
    lat = mean(decimalLatitude, na.rm = TRUE),
    long = mean(decimalLongitude, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(lat), !is.na(long))

# Join coordinates to mammal sites
mammal_sites <- mammal_sites %>%
  left_join(site_coords, by = "siteID") %>%
  filter(!is.na(lat), !is.na(long))

cat("  Mammal sites:", nrow(mammal_sites), "\n")

# ---- Create base map function ----
create_site_map <- function(site_data, title, colors, shapes, map_bounds) {
  ggplot() +
    geom_sf(data = north_america, fill = "#F5F5F5", color = "black", size = 0.3) +
    geom_point(data = site_data,
               aes(x = long, y = lat, color = data_type, shape = data_type),
               size = 2.5,
               alpha = 0.7) +
    scale_color_manual(
      name = "Data Type",
      values = colors,
      guide = guide_legend(override.aes = list(alpha = 1, size = 3))
    ) +
    scale_shape_manual(
      name = "Data Type",
      values = shapes,
      guide = guide_legend(override.aes = list(alpha = 1, size = 3))
    ) +
    coord_sf(xlim = c(map_bounds$xmin, map_bounds$xmax),
             ylim = c(map_bounds$ymin, map_bounds$ymax),
             expand = FALSE) +
    labs(
      title = title,
      subtitle = paste(nrow(site_data), "sites"),
      x = "Longitude", 
      y = "Latitude"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
}

# ---- Create the three maps ----
cat("Step 4: Creating site location maps...\n")

# Map 1: Tick sites only
map_tick <- create_site_map(
  tick_sites,
  "Tick Collection Sites",
  c("Tick" = "orange"),
  c("Tick" = 15),  # Square
  continental_bounds
)

# Map 2: Mammal sites only  
map_mammal <- create_site_map(
  mammal_sites,
  "Mammal Trapping Sites", 
  c("Mammal" = "#8B4513"),  # Brown color
  c("Mammal" = 17),  # Triangle
  continental_bounds
)

# Map 3: Combined sites
combined_sites <- bind_rows(tick_sites, mammal_sites)

map_combined <- create_site_map(
  combined_sites,
  "Tick and Mammal Sampling Sites",
  c("Tick" = "orange", "Mammal" = "#8B4513"),
  c("Tick" = 15, "Mammal" = 17),  # Square and Triangle
  continental_bounds
) +
  labs(subtitle = paste("Tick sites:", nrow(tick_sites), 
                       "| Mammal sites:", nrow(mammal_sites)
                       ))

# ---- Create summary statistics ----
cat("Step 5: Creating summary statistics...\n")

# Check for overlapping sites
tick_site_ids <- tick_sites$siteID
mammal_site_ids <- mammal_sites$siteID
overlapping_sites <- intersect(tick_site_ids, mammal_site_ids)

site_summary <- data.frame(
  Category = c("Tick sites only", "Mammal sites only", "Both tick and mammal", "Total unique sites"),
  Count = c(
    length(setdiff(tick_site_ids, mammal_site_ids)),
    length(setdiff(mammal_site_ids, tick_site_ids)), 
    length(overlapping_sites),
    length(union(tick_site_ids, mammal_site_ids))
  )
)

cat("\nSite Summary:\n")
print(site_summary)

cat("\nOverlapping sites (both tick and mammal data):", length(overlapping_sites), "\n")
cat("Overlapping site IDs:", paste(overlapping_sites, collapse = ", "), "\n")

# ---- Save outputs ----
cat("Step 6: Saving outputs...\n")

# Create output directory
if (!dir.exists("output/plots")) {
  dir.create("output/plots", recursive = TRUE)
}

# Save maps
ggsave("output/plots/site_locations_tick.png", map_tick, 
       width = 10, height = 8, dpi = 300)
ggsave("output/plots/site_locations_mammal.png", map_mammal, 
       width = 10, height = 8, dpi = 300)
ggsave("output/plots/site_locations_combined.png", map_combined, 
       width = 10, height = 8, dpi = 300)

# Save data
saveRDS(list(
  tick_sites = tick_sites,
  mammal_sites = mammal_sites,
  combined_sites = combined_sites,
  overlapping_sites = overlapping_sites,
  site_summary = site_summary
), "output/site_locations_data.RDS")

write.csv(site_summary, "output/site_summary_table.csv", row.names = FALSE)

# Save overlapping sites info
overlapping_details <- tick_sites %>%
  filter(siteID %in% overlapping_sites) %>%
  dplyr::select(siteID, lat, long) %>%
  mutate(has_tick_data = TRUE) %>%
  full_join(
    mammal_sites %>%
      filter(siteID %in% overlapping_sites) %>%
      dplyr::select(siteID, lat, long) %>%
      mutate(has_mammal_data = TRUE),
    by = c("siteID", "lat", "long")
  )

write.csv(overlapping_details, "output/overlapping_sites_details.csv", row.names = FALSE)

cat("Maps and data saved!\n")

# ---- Display maps in R plot window ----
cat("Step 7: Displaying maps in R plot window...\n")

# Display each map in sequence
cat("  Displaying tick sites map...\n")
print(map_tick)

cat("  Displaying mammal sites map...\n") 
print(map_mammal)

cat("  Displaying combined sites map...\n")
print(map_combined)

cat("\n==== SITE LOCATION MAPS COMPLETE ====\n")
cat("Created three maps:\n")
cat("- Tick sites:", nrow(tick_sites), "locations\n")
cat("- Mammal sites:", nrow(mammal_sites), "locations\n") 
cat("- Combined view with", length(overlapping_sites), "overlapping sites\n")
cat("- Black squares = tick sites, brown squares = mammal sites\n")
cat("- Transparency allows visualization of overlapping sampling efforts\n")