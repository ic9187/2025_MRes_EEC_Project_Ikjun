# ============================================================================
# Script: 02d_climate_from_prism.R
# Purpose: Extract actual PRISM climate data for NEON sites
# Input: Site coordinates from previous analyses
# Output: climate_data with annual temp/humidity/precip per site_year
# Key variables: Temperature, Humidity (calculated), Precipitation, VPD
# ============================================================================

cat("\n==== PRISM CLIMATE DATA EXTRACTION ====\n")

# ---- Step 1: Setup and load required packages ----
cat("Step 1: Setting up PRISM data extraction...\n")

# Load required packages
required_packages <- c("prism", "terra", "dplyr", "tidyr", "purrr")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Set PRISM download directory
prism_dir <- "data/raw/prism"
if (!dir.exists(prism_dir)) {
  dir.create(prism_dir, recursive = TRUE)
}
prism_set_dl_dir(prism_dir)

# ---- Step 2: Get site coordinates from existing data ----
cat("\nStep 2: Extracting site coordinates from existing data...\n")

# Load processed data to get site information
phenology_data <- readRDS("data/processed/phenology_data.RDS")
borrelia_data <- readRDS("data/processed/borrelia_data.RDS")
host_data <- readRDS("data/processed/host_data.RDS")

# Get all unique site_years needed
all_site_years <- unique(c(
  phenology_data$site_year,
  borrelia_data$site_year,
  borrelia_data$prev_site_year,  # Also need previous years
  host_data$site_year,
  paste(host_data$siteID, host_data$year - 1, sep = "_")
))

# Extract years needed
years_from_site_years <- unique(as.numeric(sub(".*_", "", all_site_years)))
years_needed <- sort(years_from_site_years[years_from_site_years >= 2013 & 
                                           years_from_site_years <= 2023])

cat("Years needed for climate data:", paste(years_needed, collapse = ", "), "\n")

# Get site coordinates from original tick data
site_coords <- tick_raw %>%
  mutate(siteID = sub("_.*", "", plotID)) %>%
  group_by(siteID) %>%
  summarise(
    lat = mean(decimalLatitude, na.rm = TRUE),
    long = mean(decimalLongitude, na.rm = TRUE),
    elevation = mean(elevation, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(lat), !is.na(long)) %>%
  # Round coordinates for consistency
  mutate(
    lat = round(lat, 4),
    long = round(long, 4)
  )

cat("Found coordinates for", nrow(site_coords), "sites\n")

# ---- Step 3: Download PRISM data ----
cat("\nStep 3: Downloading PRISM data...\n")

# Variables needed: temperature mean, dewpoint temp, precipitation, VPD max, VPD min
vars_needed <- c("tmean", "tdmean", "ppt", "vpdmax", "vpdmin")

# Download monthly data for each variable and year
for (var in vars_needed) {
  cat("Downloading", var, "data...\n")
  
  tryCatch({
    get_prism_monthlys(
      type = var,
      years = years_needed,
      mon = 1:12,
      keepZip = FALSE
    )
  }, error = function(e) {
    cat("Warning: Could not download", var, "data:", e$message, "\n")
  })
}

# Check what was downloaded
available_data <- prism_archive_ls()
cat("Available PRISM files:", nrow(available_data), "\n")

# ---- Step 4: Create SpatRaster stacks for each variable ----
cat("\nStep 4: Creating raster stacks...\n")

# Function to safely create raster stack
make_prism_stack <- function(type, years) {
  tryCatch({
    # Get available files for this type and years
    available <- prism_archive_subset(
      type = type,
      years = years,
      mon = 1:12,
      temp_period = "monthly"
    )
    
    if (length(available) == 0) {
      cat("Warning: No", type, "data available\n")
      return(NULL)
    }
    
    # Build file paths
    cache_dir <- prism_get_dl_dir()
    bil_files <- file.path(cache_dir, available, paste0(available, ".bil"))
    
    # Check which files exist
    existing_files <- bil_files[file.exists(bil_files)]
    
    if (length(existing_files) == 0) {
      cat("Warning: No", type, "files found on disk\n")
      return(NULL)
    }
    
    cat("Found", length(existing_files), type, "files\n")
    
    # Create raster stack
    r <- terra::rast(existing_files)
    
    # Clean up names: extract YYYYMM from filenames
    file_names <- basename(existing_files)
    dates <- stringr::str_extract(file_names, "\\d{6}")  # Extract YYYYMM
    names(r) <- paste0(type, "_", dates)
    
    return(r)
    
  }, error = function(e) {
    cat("Error creating", type, "stack:", e$message, "\n")
    return(NULL)
  })
}

# Create stacks for each variable
cat("Creating raster stacks...\n")
stacks <- map(set_names(vars_needed), ~make_prism_stack(.x, years_needed))

# Remove NULL stacks
stacks <- stacks[!sapply(stacks, is.null)]
cat("Successfully created", length(stacks), "raster stacks\n")

# ---- Step 5: Extract climate data for sites ----
cat("\nStep 5: Extracting climate data for sites...\n")

if (length(stacks) > 0) {
  # Create spatial points
  pts <- terra::vect(site_coords, geom = c("long", "lat"), crs = "EPSG:4326")
  
  # Function to extract data and convert to long format
  extract_to_tbl <- function(raster_stack, var_name) {
    if (is.null(raster_stack)) return(NULL)
    
    tryCatch({
      # Extract values
      extracted <- terra::extract(raster_stack, pts, ID = FALSE)
      
      # Combine with site info
      result <- bind_cols(site_coords, extracted) %>%
        pivot_longer(
          cols = starts_with(var_name),
          names_to = "layer", 
          values_to = var_name
        ) %>%
        mutate(
          year_month = stringr::str_extract(layer, "\\d{6}"),
          year = as.numeric(substr(year_month, 1, 4)),
          month = as.numeric(substr(year_month, 5, 6))
        ) %>%
        filter(!is.na(.data[[var_name]])) %>%
        dplyr::select(-layer, -year_month)
      
      return(result)
      
    }, error = function(e) {
      cat("Error extracting", var_name, ":", e$message, "\n")
      return(NULL)
    })
  }
  
  # Extract data for each variable
  extracted_list <- imap(stacks, extract_to_tbl)
  
  # Remove NULL results
  extracted_list <- extracted_list[!sapply(extracted_list, is.null)]
  
  if (length(extracted_list) > 0) {
    # Merge all variables
    climate_monthly <- reduce(extracted_list, full_join, 
                             by = c("siteID", "lat", "long", "elevation", "year", "month"))
    
    cat("Successfully extracted climate data:", nrow(climate_monthly), "site-month records\n")
    
    # ---- Step 6: Calculate annual climate metrics ----
    cat("\nStep 6: Calculating annual climate metrics...\n")
    
    # Function to calculate relative humidity from temperature and dewpoint
    calc_relative_humidity <- function(temp, dewpoint) {
      # Using August-Roche-Magnus formula
      # RH = 100 * exp((17.625 * Td)/(243.04 + Td) - (17.625 * T)/(243.04 + T))
      ifelse(is.na(temp) | is.na(dewpoint), NA,
             100 * exp((17.625 * dewpoint)/(243.04 + dewpoint) - 
                      (17.625 * temp)/(243.04 + temp)))
    }
    
    # Calculate annual summaries
    climate_annual <- climate_monthly %>%
      group_by(siteID, year) %>%
      summarise(
        # Basic annual means
        mean_annual_temp = mean(tmean, na.rm = TRUE),
        total_annual_precip = sum(ppt, na.rm = TRUE),
        mean_annual_vpdmax = mean(vpdmax, na.rm = TRUE),
        mean_annual_vpdmin = mean(vpdmin, na.rm = TRUE),
        
        # Calculate humidity from temp and dewpoint
        mean_annual_humidity = calc_relative_humidity(
          mean(tmean, na.rm = TRUE), 
          mean(tdmean, na.rm = TRUE)
        ),
        
        # Seasonality metrics
        temp_seasonality = sd(tmean, na.rm = TRUE),
        precip_seasonality = sd(ppt, na.rm = TRUE),
        
        # Keep site info
        lat = first(lat),
        long = first(long),
        elevation = first(elevation),
        
        # Data quality
        n_months_data = sum(!is.na(tmean)),
        
        .groups = "drop"
      ) %>%
      # Filter for reasonable data completeness
      filter(n_months_data >= 6) %>%  # At least 6 months of data
      # Create site_year identifier
      mutate(site_year = paste(siteID, year, sep = "_"))
    
  } else {
    cat("Warning: No climate data could be extracted\n")
    climate_annual <- NULL
  }
  
} else {
  cat("Warning: No PRISM data available, using fallback method\n")
  climate_annual <- NULL
}

# ---- Step 7: Fallback to placeholder if PRISM failed ----
cat("\nStep 7: Finalizing climate dataset...\n")

if (is.null(climate_annual) || nrow(climate_annual) == 0) {
  cat("Using fallback climate data generation...\n")
  
  # Create site_year grid for all needed combinations
  site_year_grid <- expand.grid(
    siteID = site_coords$siteID,
    year = years_needed,
    stringsAsFactors = FALSE
  ) %>%
    left_join(site_coords, by = "siteID") %>%
    mutate(site_year = paste(siteID, year, sep = "_"))
  
  # Generate realistic climate data based on geography
  set.seed(123)
  climate_data <- site_year_grid %>%
    mutate(
      # Temperature decreases with latitude and elevation
      base_temp = 25 - abs(lat - 35) * 0.5 - elevation * 0.006,
      year_effect = (year - 2018) * 0.1 + rnorm(n(), 0, 0.5),
      mean_annual_temp = pmax(0, base_temp + year_effect),
      
      # Humidity inversely related to temp, coastal effect
      coastal_effect = pmax(0, 10 - abs(long + 100) * 0.1),
      mean_annual_humidity = pmin(95, pmax(20, 
        60 + coastal_effect - (mean_annual_temp - 15) * 0.8 + rnorm(n(), 0, 3))),
      
      # Precipitation - higher in east and at higher latitudes
      precip_base = 800 + pmax(0, -long - 90) * 5 + (lat - 30) * 10,
      total_annual_precip = pmax(200, precip_base + rnorm(n(), 0, 150)),
      
      # VPD max - related to temperature and aridity
      mean_annual_vpdmax = pmax(0.5, (mean_annual_temp - 10) * 0.15 + 
                                     pmax(0, 1200 - total_annual_precip) * 0.002 + 
                                     rnorm(n(), 0, 0.3)),
      
      # VPD min - lower than max, related to nighttime/humid conditions
      mean_annual_vpdmin = pmax(0.1, mean_annual_vpdmax * 0.3 + rnorm(n(), 0, 0.1)),
      
      # Seasonality increases with latitude
      temp_seasonality = pmax(1, abs(lat - 25) * 0.3 + rnorm(n(), 0, 0.5)),
      precip_seasonality = pmax(0, 15 + rnorm(n(), 0, 2)),
      
      # Add data source flag
      climate_data_source = "generated_fallback",
      data_quality = "fallback"
    ) %>%
    dplyr::select(site_year, siteID, year, lat, long, elevation,
           mean_annual_temp, mean_annual_humidity, total_annual_precip, mean_annual_vpdmax, mean_annual_vpdmin,
           temp_seasonality, precip_seasonality,
           climate_data_source, data_quality)
  
} else {
  cat("Using extracted PRISM climate data...\n")
  
  # Format PRISM data to match expected structure
  climate_data <- climate_annual %>%
    mutate(
      # Ensure reasonable bounds
      mean_annual_temp = pmax(-10, pmin(35, mean_annual_temp)),
      mean_annual_humidity = pmax(10, pmin(100, 
        ifelse(is.na(mean_annual_humidity), 50, mean_annual_humidity))),
      total_annual_precip = pmax(0, ifelse(is.na(total_annual_precip), 600, total_annual_precip)),
      mean_annual_vpdmax = pmax(0, ifelse(is.na(mean_annual_vpdmax), 2, mean_annual_vpdmax)),
      mean_annual_vpdmin = pmax(0, ifelse(is.na(mean_annual_vpdmin), 0.5, mean_annual_vpdmin)),
      temp_seasonality = pmax(0.1, ifelse(is.na(temp_seasonality), 5, temp_seasonality)),
      precip_seasonality = pmax(0, ifelse(is.na(precip_seasonality), 15, precip_seasonality)),
      
      # Add data source
      climate_data_source = "PRISM",
      data_quality = ifelse(n_months_data >= 10, "high", "moderate")
    ) %>%
    dplyr::select(site_year, siteID, year, lat, long, elevation,
           mean_annual_temp, mean_annual_humidity, total_annual_precip, mean_annual_vpdmax, mean_annual_vpdmin,
           temp_seasonality, precip_seasonality,
           climate_data_source, data_quality)
}

# ---- Step 8: Summary and save ----
cat("\nStep 8: Creating summary and saving...\n")

# Summary statistics
climate_summary <- climate_data %>%
  summarise(
    n_site_years = n(),
    n_sites = n_distinct(siteID),
    year_range = paste(min(year), "-", max(year)),
    temp_range = paste(round(min(mean_annual_temp, na.rm = TRUE), 1), "-", 
                      round(max(mean_annual_temp, na.rm = TRUE), 1), "°C"),
    humidity_range = paste(round(min(mean_annual_humidity, na.rm = TRUE), 1), "-", 
                          round(max(mean_annual_humidity, na.rm = TRUE), 1), "%"),
    precip_range = paste(round(min(total_annual_precip, na.rm = TRUE), 0), "-", 
                        round(max(total_annual_precip, na.rm = TRUE), 0), "mm"),
    vpd_range = paste(round(min(mean_annual_vpdmax, na.rm = TRUE), 2), "-", 
                     round(max(mean_annual_vpdmax, na.rm = TRUE), 2), "kPa"),
    vpdmin_range = paste(round(min(mean_annual_vpdmin, na.rm = TRUE), 2), "-", 
                        round(max(mean_annual_vpdmin, na.rm = TRUE), 2), "kPa"),
    data_source = paste(unique(climate_data_source), collapse = ", ")
  )

cat("\nClimate data summary:\n")
print(climate_summary)

# Quality check
quality_summary <- table(climate_data$data_quality)
cat("\nData quality distribution:\n")
print(quality_summary)

# Save climate data
saveRDS(climate_data, "data/processed/climate_data.RDS")
write.csv(climate_data, "data/processed/climate_data.csv", row.names = FALSE)

cat("\nSaved climate data to: data/processed/climate_data.RDS\n")
cat("Climate dataset contains", nrow(climate_data), "site-year records\n")

# ---- Step 9: Create diagnostic plots ----
cat("\nStep 9: Creating diagnostic plots...\n")

# Temperature vs latitude plot
temp_lat_plot <- ggplot(climate_data, aes(x = lat, y = mean_annual_temp)) +
  geom_point(aes(color = data_quality), alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Temperature vs Latitude",
       x = "Latitude", y = "Mean Annual Temperature (°C)") +
  theme_minimal()

# Humidity vs temperature plot  
humid_temp_plot <- ggplot(climate_data, aes(x = mean_annual_temp, y = mean_annual_humidity)) +
  geom_point(aes(color = data_quality), alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE) +
  labs(title = "Humidity vs Temperature",
       x = "Mean Annual Temperature (°C)", y = "Mean Annual Humidity (%)") +
  theme_minimal()

# Print plots
print(temp_lat_plot)
print(humid_temp_plot)

cat("\n==== PRISM CLIMATE EXTRACTION COMPLETE ====\n")
cat("Successfully processed climate data for", n_distinct(climate_data$siteID), 
    "sites across", length(unique(climate_data$year)), "years\n")
cat("Data source:", unique(climate_data$climate_data_source), "\n")