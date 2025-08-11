library(dplyr)
library(usdm)
library(psych)

# Collinearity check: host density & VPD
vif_data <- data.frame(
  mean_density_ha = as.numeric(analysis_data$mean_density_ha),
  mean_annual_vpdmin = as.numeric(analysis_data$mean_annual_vpdmin)
) %>%
  filter(!is.na(mean_density_ha), !is.na(mean_annual_vpdmin))

# VIF test
vif_results <- usdm::vif(vif_data)
print(vif_results)

# Pairs panel
psych::pairs.panels(vif_data, main = "Host Density & VPD")
