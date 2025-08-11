# ============================================================================
# DESCRIPTIVE STATISTICS - SIMPLIFIED ANALYSIS
# ============================================================================

# ============================================================================
# PHENOLOGY ANALYSIS DATA
# ============================================================================

str(analysis_data)
unique(analysis_data$siteID)
unique(analysis_data$year)

# PAD
mean(analysis_data$PAD)
sd(analysis_data$PAD)
min(analysis_data$PAD)
max(analysis_data$PAD)

# Host density
mean(analysis_data$mean_density_ha)
sd(analysis_data$mean_density_ha)
min(analysis_data$mean_density_ha)
max(analysis_data$mean_density_ha)

# VPD min
mean(analysis_data$mean_annual_vpdmin)
sd(analysis_data$mean_annual_vpdmin)
min(analysis_data$mean_annual_vpdmin)
max(analysis_data$mean_annual_vpdmin)

# Histograms
par(mfrow = c(2, 2))
hist(analysis_data$PAD, main = "PAD")
hist(analysis_data$mean_density_ha, main = "Host Density")
hist(analysis_data$mean_annual_vpdmin, main = "VPD Min")
par(mfrow = c(1, 1))

# ============================================================================
# PREVALENCE ANALYSIS DATA
# ============================================================================

str(model_data)
unique(model_data$siteID)
unique(model_data$year)

# Response variable
mean(model_data$transmission_log)
sd(model_data$transmission_log)
min(model_data$transmission_log)
max(model_data$transmission_log)

# PAD
mean(model_data$PAD)
sd(model_data$PAD)
min(model_data$PAD)
max(model_data$PAD)

# Nymph count (raw)
mean(model_data$nymph_max)
sd(model_data$nymph_max)
min(model_data$nymph_max)
max(model_data$nymph_max)

# Nymph count (log)
mean(model_data$log_nymph_count)
sd(model_data$log_nymph_count)
min(model_data$log_nymph_count)
max(model_data$log_nymph_count)

# Larva count (raw)
mean(model_data$larva_max)
sd(model_data$larva_max)
min(model_data$larva_max)
max(model_data$larva_max)

# Larva count (log)
mean(model_data$log_larva_count)
sd(model_data$log_larva_count)
min(model_data$log_larva_count)
max(model_data$log_larva_count)

# Host density
mean(model_data$mean_density_ha)
sd(model_data$mean_density_ha)
min(model_data$mean_density_ha)
max(model_data$mean_density_ha)

# VPD max
mean(model_data$mean_annual_vpdmin)
sd(model_data$mean_annual_vpdmin)
min(model_data$mean_annual_vpdmin)
max(model_data$mean_annual_vpdmin)

# Histograms
par(mfrow = c(2, 3))
hist(model_data$transmission_log, main = "Transmission Log")
hist(model_data$PAD, main = "PAD")
hist(model_data$nymph_max, main = "Nymph Count (Raw)")
hist(model_data$log_nymph_count, main = "Nymph Count (Log)")
hist(model_data$larva_max, main = "Larva Count (Raw)")
hist(model_data$log_larva_count, main = "Larva Count (Log)")
par(mfrow = c(1, 1))

par(mfrow = c(1, 2))
hist(model_data$mean_density_ha, main = "Host Density")
hist(model_data$mean_annual_vpdmin, main = "VPD Min")
par(mfrow = c(1, 1))
