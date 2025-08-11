library(ggplot2)

# Load the model data (assuming you've run the script)
model_data <- readRDS("output/05_simplified_model_data.RDS")

# Scatter plot: transmission_log vs PAD (colored by site)
p1 <- ggplot(model_data, aes(x = PAD, y = transmission_log, color = siteID)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_point(alpha = 0.7, size = 2) +
  labs(x = "PAD", y = "log(prevalence change)", title = "Prevalence vs PAD") +
  theme_minimal() +
  scale_color_viridis_d(alpha = 0.8) +
  theme(legend.position = "right")

# Scatter plot: transmission_log vs Host Density (colored by site)
p2 <- ggplot(model_data, aes(x = mean_density_ha, y = transmission_log, color = siteID)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_point(alpha = 0.7, size = 2) +
  labs(x = "Mean Host Density (per ha)", y = "log(prevalence change)", title = "Prevalence vs Host Density") +
  theme_minimal() +
  scale_color_viridis_d(alpha = 0.8) +
  theme(legend.position = "right")

# Scatter plot: transmission_log vs VPD Min (colored by site)
p3 <- ggplot(model_data, aes(x = mean_annual_vpdmin, y = transmission_log, color = siteID)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_point(alpha = 0.7, size = 2) +
  labs(x = "Vapor Pressure Deficit Min (VPD)", y = "log(prevalence change)", title = "Prevalence vs VPD Min") +
  theme_minimal() +
  scale_color_viridis_d(alpha = 0.8) +
  theme(legend.position = "right")

# Display plots
p1
p2
p3
