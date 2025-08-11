# ==== NEON Tick Data Analysis - Setup Script ====
# This script sets up the working environment and loads all raw data
# Source this at the beginning of all analysis scripts

# ==== Working Directory ====
# Set the working directory to the main project folder
project_dir <- "~/Desktop/MRes_EEC_24/70067_PROJECT"
setwd(project_dir)
cat("Working directory set to:", getwd(), "\n")

# ==== Load Packages ====
# Function to load required packages (and install if missing)
load_packages <- function(pkg_list) {
  new_packages <- pkg_list[!(pkg_list %in% installed.packages()[,"Package"])]
  if(length(new_packages) > 0) {
    cat("Installing missing packages:", paste(new_packages, collapse=", "), "\n")
    install.packages(new_packages)
  }
  
  # Load all packages with messages suppressed
  invisible(lapply(pkg_list, function(pkg) {
    suppressMessages(library(pkg, character.only = TRUE))
    cat("Loaded package:", pkg, "\n")
  }))
}

# Package lists by category
basic_pkgs <- c("ggplot2", "lattice", "tidyr", "dplyr", "lubridate")
gis_pkgs <- c("terra", "sf", "sp", "raster", "geodata", "lwgeom", "openxlsx",
              "gridExtra", "dismo", "ncf", "fastmatrix", "SpatialPack", "spData",
              "spdep", "spatialreg", "nlme", "spgwr", "landscapemetrics", "vegan",
              "units", "rnaturalearth", "prism", "purrr")
stats_pkgs <- c("lme4", "WebPower", "MASS", "ggpubr", "usdm", "psych",
                "lmerTest", "lmtest", "sjPlot", "merTools", "DHARMa", "ggeffects",
                "viridis", "patchwork", "ggnewscale", "mgcv", "MuMIn", "glmmTMB")
neon_pkgs <- c("neonUtilities", "neonOS")

# Load packages by category
cat("\nLoading required packages...\n")
load_packages(basic_pkgs)
load_packages(stats_pkgs)
load_packages(neon_pkgs)

# ==== Load Raw Data ====
cat("\n==== LOADING RAW DATA ====\n")

# 1. Tick count data (sampling effort)
cat("Loading tick count data...")
tick_raw <- read.csv("data/raw/tick_data.csv", stringsAsFactors = FALSE)
cat(" Done. (", nrow(tick_raw), "rows )\n")
# not that the day_number and year columns were added before this right after
# importing the original data from NEON! Importing from NEON is on NEON website

# 2. Tick taxonomy data
cat("Loading tick taxonomy data...")
taxo_raw <- read.csv("data/raw/tick_taxonomy.csv", stringsAsFactors = FALSE)
cat(" Done. (", nrow(taxo_raw), "rows )\n")

# 3. Pathogen data
cat("Loading pathogen data...")
path_raw <- read.csv("data/raw/pathogen_data.csv", stringsAsFactors = FALSE)
cat(" Done. (", nrow(path_raw), "rows )\n")

# 4. Small mammal data
cat("Loading small mammal data...")
mammal_raw <- read.csv("data/raw/mammal_data.csv", stringsAsFactors = FALSE)
cat(" Done. (", nrow(mammal_raw), "rows )\n")

# 5. Trap data (for mammal density calculations)
cat("Loading trap data...")
trap_raw <- read.csv("data/raw/trap_data.csv", stringsAsFactors = FALSE)
cat(" Done. (", nrow(trap_raw), "rows )\n")

cat("\n==== All raw data loaded successfully! ====\n")
cat("Setup complete. Ready for analysis.\n")