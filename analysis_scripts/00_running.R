# ============================================================================
# Analysis overview: running the scripts (and any temporary codes)
# ============================================================================

# setting up the environment (working directory, packages, raw data)
source("~/Desktop/MRes_EEC_24/70067_PROJECT/R_scripts/00_setup.R")

# handling sampling effort for tick field data (for PAD, prevalence and mammal)
source("R_scripts/01a_tick_sampling_effort.R")
source("R_scripts/01b_mammal_sampling_effort.R")
source("R_scripts/01c_spatial_sampling_effort.R")

# variable construction from raw data
source("R_scripts/02a_taxonomy_to_PAD.R")
source("R_scripts/02b_pathogen_to_borrelia.R")
source("R_scripts/02c_mam_to_host.R")
source("R_scripts/02d_climate_from_prism.R")

# checking and merging the variable data
source("R_scripts/03a_tick_pipeline_summary.R")
source("R_scripts/03b_merge_phenology_data.R")
source("R_scripts/03c_merge_prevalence_data.R")

# phenology analysis
source("R_scripts/04a_phenology_analysis.R")
source("R_scripts/phenology_additional_vpd_only.R")
summary(phenology_models[["vpd_only"]])
performance::r2(phenology_models[["vpd_only"]])
names(phenology_models)
source("R_scripts/04c_diagnostic_plots.R")
source("R_scripts/for_thesis/outlier_investigation.r")
source("R_scripts/for_thesis/peak_sensitivity.r")

# prevalence analysis
source("R_scripts/05a_prevalence_analysis.R")
summary(prevalence_models[["pad_only"]])
names(prevalence_models)
source("R_scripts/05b_power_analysis.R")

# ----------------------------------------------------------------------------

# visualization (figures and tables for the thesis) # [if needed -> par(mfrow = c(1, 1))]
# for Introduction and Methods
source("R_scripts/for_thesis/ukfs_density_plot.r")
source("R_scripts/for_thesis/site_location_maps.r")
source("R_scripts/for_thesis/processed_maps.r")

# for Results in general
source("R_scripts/for_thesis/descriptive_stats.r")
source("R_scripts/for_thesis/tab_models.r")

# phenology analysis
source("R_scripts/for_thesis/habitat_plots_phenology.r")
source("R_scripts/for_thesis/phenology_latitude_gam.r")
source("R_scripts/for_thesis/phenology_plots.r")

# prevalence analysis
source("R_scripts/for_thesis/habitat_plots_prevalence.r")
source("R_scripts/for_thesis/prevalence_plots.r")
source("R_scripts/05b_power_analysis.R")
