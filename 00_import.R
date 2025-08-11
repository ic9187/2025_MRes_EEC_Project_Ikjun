# =============================================================================
# Script: 00_import.R
# Purpose: Import raw NEON data and perform initial processing
# Input: Raw NEON CSV files from stackedFiles directories
# Output: Processed CSV files with temporal variables added
# Key approach:
#   - Load raw NEON data using neonUtilities::readTableNEON
#   - Add day_number and year variables for temporal analysis
#   - Save processed files for analysis pipeline
# NEON packages installing
# https://www.neonscience.org/resources/learning-hub/tutorials/download-explore-neon-data?check_logged_in=1
# =============================================================================

cat("\n==== NEON DATA IMPORTING AND INITIAL PROCESSING ====\n")

# ---- Step 1: Load required packages ----
cat("Step 1: Loading NEON packages...\n")
library(neonUtilities)
library(neonOS)

# ---- Step 2: Import Tick Field Data (DP1.10093) ----
cat("\nStep 2: Processing tick field data...\n")

tick_data <- readTableNEON(
  dataFile="data/NEON_count-ticks/stackedFiles/tck_fielddata.csv", 
  varFile="data/NEON_count-ticks/stackedFiles/variables_10093.csv")

tick_data$day_number <- yday(as.Date(tick_data$collectDate))
tick_data$year <- year(as.Date(tick_data$collectDate))
write.csv(tick_data, "tick_data.csv", row.names = F)
cat("  ✓ Saved tick_data.csv (", nrow(tick_data), "rows )\n")

# ---- Step 3: Import Tick Taxonomy Data (DP1.10093) ----
cat("\nStep 3: Processing tick taxonomy data...\n")

tick_taxonomy <- readTableNEON(
  dataFile="data/NEON_count-ticks/stackedFiles/tck_taxonomyProcessed.csv", 
  varFile="data/NEON_count-ticks/stackedFiles/variables_10093.csv")

write.csv(tick_taxonomy, "tick_taxonomy.csv", row.names = F)
cat("  ✓ Saved tick_taxonomy.csv (", nrow(tick_taxonomy), "rows )\n")

# ---- Step 4: Import Pathogen Status Data (DP1.10092) ----
cat("\nStep 4: Processing pathogen data...\n")

pathogen_data <- readTableNEON(
  dataFile="data/NEON_pathogens-tick/stackedFiles/tck_pathogen.csv", 
  varFile="data/NEON_pathogens-tick/stackedFiles/variables_10092.csv")

pathogen_data$day_number <- yday(as.Date(pathogen_data$collectDate))
pathogen_data$year <- year(as.Date(pathogen_data$collectDate))
write.csv(pathogen_data, "pathogen_data.csv", row.names = F)
cat("  ✓ Saved pathogen_data.csv (", nrow(pathogen_data), "rows )\n")

# ---- Step 5: Import Mammal Trap Data (DP1.10072) ----
cat("\nStep 5: Processing mammal trap data...\n")

mammal_data <- readTableNEON(
  dataFile="data/NEON_count-small-mammals/stackedFiles/mam_pertrapnight.csv", 
  varFile="data/NEON_count-small-mammals/stackedFiles/variables_10072.csv")

mammal_data$day_number <- yday(as.Date(mammal_data$collectDate))
mammal_data$year <- year(as.Date(mammal_data$collectDate))
write.csv(mammal_data, "mammal_data.csv", row.names = F)
cat("  ✓ Saved mammal_data.csv (", nrow(mammal_data), "rows )\n")

# ---- Step 6: Import Mammal Plot Data (DP1.10072) ----
cat("\nStep 6: Processing mammal plot data...\n")

trap_data <- readTableNEON(
  dataFile="data/NEON_count-small-mammals/stackedFiles/mam_perplotnight.csv", 
  varFile="data/NEON_count-small-mammals/stackedFiles/variables_10072.csv")

write.csv(trap_data, "trap_data.csv", row.names = F)
cat("  ✓ Saved trap_data.csv (", nrow(trap_data), "rows )\n")

cat("\n==== NEON data import complete ====\n")