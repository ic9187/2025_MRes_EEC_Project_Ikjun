# ============================================================================
# Script: 02c_mam_to_host.R
# Purpose: Calculate mammal host variables with phylogenetic competence weighting
# Input: mammal_raw, trap_raw (from 00_setup.R)
# Output: host_data with density and diversity metrics per site_year
# Key decisions:
#   - Hierarchical phylogenetic competence weights (LoGiudice → family → order means)
#   - Minimum 50 traps set for reliable density estimate
#   - Including recapture rate summary
#   - Host diversity metrics
# ============================================================================

cat("\n==== MAMMAL TO HOST DENSITY & DIVERSITY ANALYSIS ====\n")

# ============================================================================
# CONFIGURATION SECTION - Easy to modify
# ============================================================================

# Analysis parameters
trap_threshold_high <- 80  # High-quality data threshold
trap_threshold_low <- 50   # Minimum acceptable threshold

# Host competence weights (LoGiudice 2003 hierarchical approach)
host_competence_weights <- list(
  
  # ========== EXACT LOGIUDICE EMPIRICAL VALUES ==========
  
  # Genus Peromyscus (LoGiudice empirical: 92.1/92.1 = 1.000)
  peromyscus = list(
    weight = 1.000,
    source = "LoGiudice_empirical_92.1/92.1",
    species = c(
      "Peromyscus leucopus",                # White-footed mouse
      "Peromyscus maniculatus",             # Deer mouse  
      "Peromyscus gossypinus",              # Cotton mouse
      "Peromyscus sp.",                     # Peromyscus in general
      "Peromyscus gossypinus/leucopus",     # Mixed ID
      "Peromyscus leucopus/maniculatus",    # Mixed ID
      "Peromyscus polionotus"               # Oldfield mouse
    )
  ),
  
  # Genus Tamias (LoGiudice empirical: 55.0/92.1 = 0.597)
  tamias = list(
    weight = 0.597,
    source = "LoGiudice_empirical_55.0/92.1",
    species = c(
      "Tamias striatus",                    # Eastern chipmunk
      "Tamias minimus",                     # Least chipmunk
      "Tamias sp."                          # Tamias in general
    )
  ),
  
  # Genus Sorex (LoGiudice empirical: 51.2/92.1 = 0.556)
  sorex = list(
    weight = 0.556,
    source = "LoGiudice_empirical_51.2/92.1",
    species = c(
      "Sorex cinereus",                     # Masked shrew
      "Sorex arcticus",                     # Arctic shrew
      "Sorex fumeus",                       # Smoky shrew
      "Sorex hoyi",                         # Pygmy shrew
      "Sorex palustris",                    # Water shrew
      "Sorex longirostris",                 # Southeastern shrew
      "Sorex sp."                           # Sorex in general
    )
  ),
  
  # Genus Blarina (LoGiudice empirical: 41.8/92.1 = 0.454)
  blarina = list(
    weight = 0.454,
    source = "LoGiudice_empirical_41.8/92.1",
    species = c(
      "Blarina brevicauda",                 # Northern short-tailed shrew
      "Blarina carolinensis",               # Southern short-tailed shrew
      "Blarina hylophaga",                  # Elliot's short-tailed shrew
      "Blarina sp."                         # Blarina in general
    )
  ),
  
  # Genus Sciurus (LoGiudice empirical: 14.7/92.1 = 0.160)
  sciurus = list(
    weight = 0.160,
    source = "LoGiudice_empirical_14.7/92.1",
    species = c(
      "Sciurus carolinensis"                # Eastern gray squirrel
    )
  ),
  
  # Genus Didelphis (LoGiudice empirical: 2.6/92.1 = 0.028)
  didelphis = list(
    weight = 0.028,
    source = "LoGiudice_empirical_2.6/92.1",
    species = c(
      "Didelphis virginiana"                # Virginia opossum
    )
  ),
  
  # ========== FAMILY LEVEL MEANS ==========
  
  # Family Cricetidae uncovered genera (family mean of Peromyscus only = 1.000)
  cricetidae_uncovered = list(
    weight = 1.000,
    source = "Family_Cricetidae_mean_Peromyscus_1.000",
    species = c(
      "Sigmodon hispidus",                  # Cotton rat
      "Sigmodon hispidus eremicus",         # Cotton rat subspecies
      "Myodes gapperi",                     # Red-backed vole
      "Microtus ochrogaster",               # Prairie vole
      "Microtus pennsylvanicus",            # Meadow vole
      "Microtus pinetorum",                 # Woodland vole
      "Neotoma floridana",                  # Eastern woodrat
      "Onychomys leucogaster",              # Northern grasshopper mouse
      "Perognathus flavus",                 # Silky pocket mouse
      "Perognathus flavescens",             # Plains pocket mouse
      "Reithrodontomys megalotis",          # Western harvest mouse
      "Reithrodontomys fulvescens",         # Fulvous harvest mouse
      "Reithrodontomys montanus",           # Plains harvest mouse
      "Reithrodontomys humulis",            # Eastern harvest mouse
      "Ochrotomys nuttalli",                # Golden mouse
      "Oryzomys palustris",                 # Marsh rice rat
      "Synaptomys cooperi",                 # Southern bog lemming
      "Ictidomys tridecemlineatus",         # Thirteen-lined ground squirrel
      "Cricetidae sp.",                     # Cricetidae in general
      "Rodentia sp."                        # Rodentia in general
    )
  ),
  
  # Family Sciuridae uncovered genera (family mean: (0.597+0.160)/2 = 0.379)
  sciuridae_uncovered = list(
    weight = 0.379,
    source = "Family_Sciuridae_mean_(0.597+0.160)/2",
    species = c(
      "Glaucomys volans",                   # Southern flying squirrel
      "Glaucomys sabrinus",                 # Northern flying squirrel
      "Glaucomys sp.",                      # Flying squirrels in general
      "Tamiasciurus hudsonicus"             # Red squirrel
    )
  ),
  
  # Family Soricidae uncovered genera (family mean: (0.556+0.454)/2 = 0.505)
  soricidae_uncovered = list(
    weight = 0.505,
    source = "Family_Soricidae_mean_(0.556+0.454)/2",
    species = c(
      "Cryptotis parva",                    # Least shrew
      "Soricidae sp."                       # Soricidae in general
    )
  ),
  
  # ========== ORDER LEVEL MEANS ==========
  
  # Family Dipodidae (Order Rodentia mean: (1.000+0.597+0.160)/3 = 0.586)
  dipodidae = list(
    weight = 0.586,
    source = "Order_Rodentia_mean_(1.000+0.597+0.160)/3",
    species = c(
      "Zapus hudsonius",                    # Meadow jumping mouse
      "Napaeozapus insignis",               # Woodland jumping mouse
      "Dipodidae sp."                       # Jumping mice in general
    )
  ),
  
  # Family Muridae (Order Rodentia mean: 0.586)
  muridae = list(
    weight = 0.586,
    source = "Order_Rodentia_mean_(1.000+0.597+0.160)/3",
    species = c(
      "Mus musculus",                       # House mouse
      "Rattus rattus"                       # Black rat
    )
  ),
  
  # ========== HIGHER TAXONOMIC LEVELS ==========
  
  # Remaining species (overall vertebrate small mammal mean)
  remaining_species = list(
    weight = 0.568,
    source = "All_LoGiudice_species_mean_(1.000+0.597+0.556+0.454+0.160+0.028)/6",
    species = c(
      "Podomys floridanus",                 # Florida mouse
      "Chaetodipus hispidus",               # Hispid pocket mouse
      "Sylvilagus floridanus",              # Eastern cottontail
      "Mustela frenata",                    # Long-tailed weasel
      "Mustela erminea",                    # Ermine
      "Mustela sp.",                        # Weasels in general
      "Mammalia sp."                        # Mammalia in general
    )
  )
)

# Default weight for any species not covered above
default_weight <- 0.568  # Overall mean of all LoGiudice empirical values

# Function to assign competence weight
assign_competence_weight <- function(species_name) {
  if (is.na(species_name) || species_name == "") {
    return(default_weight)
  }
  
  for (group in host_competence_weights) {
    if (species_name %in% group$species) {
      return(group$weight)
    }
  }
  
  # If not found in any group, use default
  return(default_weight)
}

# ============================================================================
# DATA PROCESSING
# ============================================================================

# ---- Step 1: Load and prepare mammal data ----
cat("\nStep 1: Loading and preparing mammal data...\n")

# Check if mammal_raw exists
if (!exists("mammal_raw")) {
  stop("Error: mammal_raw not found. Please run 00_setup.R first.")
}

# Prepare mammal data with competence weights
mam_data <- mammal_raw %>%
  # Create site_year identifier
  mutate(
    year = year(collectDate),
    site_year = paste(siteID, year, sep = "_"),
    trap_date = as.Date(collectDate),
    
    # Clean trap status (remove text, keep numbers)
    trapStatus_num = as.numeric(gsub("[^0-9]", "", trapStatus)),
    # Binary capture (4,5 = capture)
    capture_success = ifelse(trapStatus_num %in% c(4, 5), 1, 0)
  ) %>%
  # Assign competence weights
  mutate(
    competence_weight = sapply(scientificName, assign_competence_weight)
  )

# ---- Step 1.5: Apply mammal sampling effort filter ----
cat("\n=== APPLYING MAMMAL SAMPLING EFFORT FILTER ===\n")

# Load mammal sampling effort filter
clean_mammal_site_years <- readRDS("data/processed/clean_mammal_site_years.RDS")
cat("Loaded", length(clean_mammal_site_years), "clean mammal site_years\n")

# Filter mammal data before processing
mam_data_filtered <- mam_data %>%
  filter(site_year %in% clean_mammal_site_years)

cat("Mammal data before filter:", nrow(mam_data), "rows\n")
cat("Mammal data after filter:", nrow(mam_data_filtered), "rows\n")
cat("Retention rate:", round(nrow(mam_data_filtered)/nrow(mam_data)*100, 1), "%\n")

# Replace the original mam_data with filtered version
mam_data <- mam_data_filtered

# Check site-year coverage
filtered_site_years <- unique(mam_data$site_year)
cat("Site-years retained:", length(filtered_site_years), "\n")
cat("Sites represented:", n_distinct(mam_data$siteID), "\n\n")

# ---- Step 2: Enhanced data quality diagnostics ----
cat("\nStep 2: Enhanced data quality diagnostics...\n")

# Diagnostic Check 1: Species identification coverage
cat("\n--- DIAGNOSTIC: Species Identification Coverage ---\n")
species_coverage <- mam_data %>%
  filter(capture_success == 1) %>%
  group_by(site_year) %>%
  summarise(
    total_captures = n(),
    identified_to_species = sum(!is.na(scientificName)),
    coverage_rate = identified_to_species / total_captures,
    .groups = "drop"
  ) %>%
  arrange(coverage_rate)

cat("Overall species ID rate:", 
    round(sum(species_coverage$identified_to_species) / sum(species_coverage$total_captures) * 100, 1), 
    "%\n")
cat("Site-years with <80% species ID:", 
    sum(species_coverage$coverage_rate < 0.8), 
    "out of", nrow(species_coverage), "\n")

# Diagnostic Check 2: Competence weight coverage
cat("\n--- DIAGNOSTIC: Competence Weight Coverage ---\n")
weight_coverage <- mam_data %>%
  filter(capture_success == 1) %>%
  mutate(
    weight_source = case_when(
      competence_weight == 1.000 ~ "LoGiudice_Peromyscus",
      competence_weight == 0.597 ~ "LoGiudice_Tamias", 
      competence_weight == 0.556 ~ "LoGiudice_Sorex",
      competence_weight == 0.454 ~ "LoGiudice_Blarina",
      competence_weight == 0.160 ~ "LoGiudice_Sciurus",
      competence_weight == 0.028 ~ "LoGiudice_Didelphis",
      competence_weight == 0.379 ~ "Family_Sciuridae",
      competence_weight == 0.505 ~ "Family_Soricidae",
      competence_weight == 0.586 ~ "Order_Rodentia",
      competence_weight == 0.568 ~ "Overall_mean",
      TRUE ~ "Other"
    )
  ) %>%
  count(weight_source, sort = TRUE) %>%
  mutate(percentage = round(n / sum(n) * 100, 1))

print(weight_coverage)

# ---- Step 3: Temporal coverage ----
cat("\n--- DIAGNOSTIC: Temporal Coverage ---\n")
mammal_years <- mam_data %>%
  group_by(siteID, year) %>%
  summarise(n_records = n(), .groups = "drop") %>%
  group_by(siteID) %>%
  summarise(
    first_year = min(year),
    last_year = max(year),
    n_years = n_distinct(year),
    .groups = "drop"
  ) %>%
  arrange(first_year)

cat("Temporal coverage summary:\n")
cat("Years represented:", min(mammal_years$first_year), "to", max(mammal_years$last_year), "\n")
cat("Mean years per site:", round(mean(mammal_years$n_years), 1), "\n")

# ---- Step 4: Calculate plot-night densities with enhanced quality assessment ----
cat("\nStep 4: Calculating plot-night densities with enhanced quality assessment...\n")

# Plot-night summaries with enhanced quality tiers
plot_night_summary <- mam_data %>%
  group_by(plotID, trap_date, site_year, siteID, year) %>%
  summarise(
    # Trap effort
    traps_set = n_distinct(trapCoordinate[trapStatus_num != 1]),
    
    # Captures
    total_captures = sum(capture_success, na.rm = TRUE),
    unique_individuals = n_distinct(tagID[capture_success == 1 & !is.na(tagID)]),
    
    # Weighted captures
    weighted_captures = sum(competence_weight[capture_success == 1], na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  # Calculate densities (assuming 100m² per trap)
  mutate(
    area_ha = (traps_set * 100) / 10000,
    density_per_ha = total_captures / area_ha,
    individual_density_per_ha = unique_individuals / area_ha,
    weighted_density_per_ha = weighted_captures / area_ha,
    
    # Enhanced quality flag with three tiers
    quality_tier = case_when(
      traps_set >= trap_threshold_high ~ "high",
      traps_set >= trap_threshold_low ~ "acceptable",
      TRUE ~ "low"
    )
  )

# Report trap effort distribution
trap_summary <- plot_night_summary %>%
  group_by(quality_tier) %>%
  summarise(
    n_nights = n(),
    mean_traps = round(mean(traps_set), 1),
    .groups = "drop"
  )

cat("Trap effort distribution:\n")
print(trap_summary)

# ---- Step 5: Aggregate to site_year level with additional metrics ----
cat("\nStep 5: Aggregating to site_year level with enhanced metrics...\n")

# Site-year summaries with additional metrics
site_year_density <- plot_night_summary %>%
  group_by(site_year, siteID, year) %>%
  summarise(
    # Sampling effort
    n_trap_nights = n(),
    total_traps_set = sum(traps_set),
    
    # Enhanced density metrics
    mean_density_ha = mean(density_per_ha, na.rm = TRUE),
    median_density_ha = median(density_per_ha, na.rm = TRUE),
    max_density_ha = max(density_per_ha, na.rm = TRUE),
    
    # Individual-based (accounts for recaptures)
    mean_individual_density_ha = mean(individual_density_per_ha, na.rm = TRUE),
    
    # Weighted by phylogenetic competence
    mean_weighted_density_ha = mean(weighted_density_per_ha, na.rm = TRUE),
    
    # Basic capture metrics
    total_captures = sum(total_captures),
    
    # Quality metrics
    prop_high_quality = mean(quality_tier == "high"),
    
    .groups = "drop"
  ) %>%
  # Enhanced quality flag
  mutate(
    quality_flag = case_when(
      n_trap_nights < 3 ~ "insufficient_nights",
      total_traps_set < 150 ~ "low_effort",
      total_traps_set >= 150 & total_traps_set < 300 ~ "moderate_effort",
      TRUE ~ "high_effort"
    ),
    # Add temporal alignment columns
    next_year = year + 1,
    next_site_year = paste(siteID, next_year, sep = "_")
  )

# ---- Step 6: Calculate host diversity metrics ----
cat("\nStep 6: Calculating host diversity metrics...\n")

# Load vegan for diversity calculations
if (!require("vegan", quietly = TRUE)) {
  install.packages("vegan")
  library(vegan)
}

# Species captures for diversity calculations
species_captures <- mam_data %>%
  filter(capture_success == 1, !is.na(scientificName)) %>%
  group_by(site_year, siteID, year, scientificName) %>%
  summarise(
    total_captures = n(),
    .groups = "drop"
  )

# Calculate diversity metrics at site-year level
host_diversity_data <- species_captures %>%
  group_by(site_year, siteID, year) %>%
  summarise(
    # Basic diversity metrics
    species_richness = n_distinct(scientificName),
    shannon_diversity = ifelse(n_distinct(scientificName) > 1,
                               diversity(table(scientificName), index = "shannon"),
                               0),
    simpson_diversity = ifelse(n_distinct(scientificName) > 1,
                               diversity(table(scientificName), index = "simpson"), 
                               0),
    
    .groups = "drop"
  )

# ---- Step 7: Calculate recapture rates ----
cat("\nStep 7: Calculating recapture rates...\n")

# Recapture analysis
individual_summary <- mam_data %>%
  filter(capture_success == 1, !is.na(tagID)) %>%
  group_by(site_year, tagID) %>%
  summarise(
    n_captures = n(),
    .groups = "drop"
  ) %>%
  mutate(is_recapture = n_captures > 1)

recapture_rates <- individual_summary %>%
  group_by(site_year) %>%
  summarise(
    total_individuals = n(),
    recaptured_individuals = sum(is_recapture),
    recapture_rate = ifelse(total_individuals > 0, recaptured_individuals / total_individuals, 0),
    .groups = "drop"
  )

# ---- Step 8: Combine all host data ----
cat("\nStep 8: Combining all host data...\n")

# Merge all components
host_data <- site_year_density %>%
  left_join(host_diversity_data, by = c("site_year", "siteID", "year")) %>%
  left_join(recapture_rates, by = "site_year") %>%
  # Handle missing diversity data
  mutate(
    species_richness = ifelse(is.na(species_richness), 0, species_richness),
    shannon_diversity = ifelse(is.na(shannon_diversity), 0, shannon_diversity),
    simpson_diversity = ifelse(is.na(simpson_diversity), 0, simpson_diversity),
    total_individuals = ifelse(is.na(total_individuals), 0, total_individuals),
    recapture_rate = ifelse(is.na(recapture_rate), 0, recapture_rate)
  )

# ---- Step 9: Final summary and diagnostic output ----
cat("\n=== FINAL HOST DATA SUMMARY ===\n")

# Summary statistics
summary_stats <- host_data %>%
  filter(quality_flag %in% c("moderate_effort", "high_effort")) %>%
  summarise(
    n_site_years = n(),
    mean_density = mean(mean_density_ha, na.rm = TRUE),
    mean_weighted_density = mean(mean_weighted_density_ha, na.rm = TRUE),
    mean_individual_density = mean(mean_individual_density_ha, na.rm = TRUE),
    mean_recapture_rate = mean(recapture_rate, na.rm = TRUE),
    mean_species_richness = mean(species_richness, na.rm = TRUE),
    mean_shannon_diversity = mean(shannon_diversity, na.rm = TRUE)
  )

cat("\n  Enhanced Host Summary:\n")
print(summary_stats)

cat("\n  Effort distribution:\n")
print(table(host_data$quality_flag))

cat("\n  Species richness distribution:\n")
print(table(host_data$species_richness, useNA = "ifany"))

# Competence weight summary
cat("\n  Competence weighting summary:\n")
species_summary <- species_captures %>%
  left_join(
    mam_data %>% 
      filter(capture_success == 1) %>%
      dplyr::select(scientificName, competence_weight) %>%
      distinct(),
    by = "scientificName"
  ) %>%
  group_by(scientificName) %>%
  summarise(
    total_captures = sum(total_captures),
    weight = first(competence_weight),
    .groups = "drop"
  ) %>%
  arrange(desc(total_captures))

cat("Top 10 species by captures:\n")
print(head(species_summary, 10))

# Save host data
saveRDS(host_data, "data/processed/host_data.RDS")
write.csv(host_data, "data/processed/host_data.csv", row.names = FALSE)

cat("\n  Saved host data to: data/processed/host_data.RDS\n")
cat("\n==== HOST DENSITY & DIVERSITY CALCULATION COMPLETE ====\n")

# Additional diagnostic output
cat("\n=== FINAL DIAGNOSTICS ===\n")
cat("Total site-years processed:", nrow(host_data), "\n")
cat("Site-years with diversity data:", sum(!is.na(host_data$species_richness)), "\n")
cat("High-quality site-years:", sum(host_data$quality_flag == "high_effort", na.rm = TRUE), "\n")
cat("Species in competence database:", 
    length(unlist(lapply(host_competence_weights, function(x) x$species))), "\n")
cat("Mean weighted density (all sites):", 
    round(mean(host_data$mean_weighted_density_ha, na.rm = TRUE), 3), "per ha\n")
cat("Mean species richness (all sites):", 
    round(mean(host_data$species_richness, na.rm = TRUE), 1), "\n")