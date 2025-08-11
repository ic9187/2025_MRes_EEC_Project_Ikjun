# ============================================================================
# SECTION 8: GENERALIZED ADDITIVE MODEL (GAM) ANALYSIS
# ============================================================================

cat("\n=== SECTION 8: GAM ANALYSIS: PAD ~ LATITUDE ===\n")

# Load mgcv for GAM
library(mgcv)

# Prepare data for GAM (use sites_with_coords from Section 5)
gam_data <- sites_with_coords %>%
  filter(!is.na(PAD), !is.na(lat), !is.na(siteID)) %>%
  # Get latitude from the lat column
  rename(latitude = lat) %>%
  mutate(siteID = as.factor(siteID))

cat("GAM analysis data:\n")
cat("  Observations:", nrow(gam_data), "\n")
cat("  Sites:", length(unique(gam_data$siteID)), "\n")
cat("  Latitude range:", round(min(gam_data$latitude), 2), "to", round(max(gam_data$latitude), 2), "\n")
cat("  PAD range:", round(min(gam_data$PAD), 2), "to", round(max(gam_data$PAD), 2), "\n")

# ============================================================================
# FIT GAM MODEL
# ============================================================================

cat("\nFitting GAM: PAD ~ s(latitude) + s(siteID, bs='re')\n")

# Fit the GAM model
gam_model <- gam(PAD ~ s(latitude, k = 8) + s(siteID, bs = "re"), 
                 data = gam_data,
                 method = "REML")

cat("✓ GAM model fitted successfully\n")

# Extract model statistics
gam_summary <- summary(gam_model)
latitude_p <- gam_summary$s.table["s(latitude)", "p-value"]
site_re_p <- gam_summary$s.table["s(siteID)", "p-value"]
r_squared <- gam_summary$r.sq
dev_explained <- gam_summary$dev.expl

cat("\n--- KEY RESULTS ---\n")
cat("R-squared:", round(r_squared, 3), "\n")
cat("Deviance explained:", round(dev_explained * 100, 1), "%\n")
cat("Latitude smooth p-value:", format.pval(latitude_p, digits = 3), "\n")
cat("Site random effect p-value:", format.pval(site_re_p, digits = 3), "\n")

print(summary(gam_model))

# ============================================================================
# CREATE PREDICTIONS
# ============================================================================

cat("\n=== CREATING PREDICTIONS ===\n")

# Generate prediction data
lat_range <- range(gam_data$latitude)
pred_lat <- seq(lat_range[1], lat_range[2], length.out = 30)

# Create prediction dataframe
pred_data <- data.frame(
  latitude = pred_lat,
  siteID = factor(rep(levels(gam_data$siteID)[1], 30), levels = levels(gam_data$siteID))
)

# Get predictions
pred_result <- predict(gam_model, newdata = pred_data, se.fit = TRUE)
pred_data$fit <- pred_result$fit
pred_data$se <- pred_result$se.fit
pred_data$lower <- pred_data$fit - 1.96 * pred_data$se
pred_data$upper <- pred_data$fit + 1.96 * pred_data$se

cat("Predictions created successfully\n")

# ============================================================================
# FINAL GAM PLOT
# ============================================================================

cat("\n=== CREATING FINAL GAM PLOT ===\n")

# Create the plot with YOUR REQUESTED COLORS
gam_plot <- ggplot(gam_data, aes(x = latitude, y = PAD)) +
  # Add grey confidence band first
  geom_ribbon(data = pred_data, 
              aes(x = latitude, ymin = lower, ymax = upper),
              fill = "grey40", alpha = 0.3, inherit.aes = FALSE) +
  geom_point(aes(color = nlcdClass), size = 2, alpha = 0.7) +
  geom_line(data = pred_data, aes(x = latitude, y = fit), 
            color = "black", size = 1, inherit.aes = FALSE) +
  scale_color_viridis_d(name = "Habitat Type") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "GAM: Peak Activity Difference vs Latitude",
    subtitle = paste0("PAD ~ s(latitude) + s(siteID, bs='re') | R² = ", 
                     round(r_squared, 3), " | p = ", format.pval(latitude_p, digits = 3)),
    x = "Latitude (°N)",
    y = "PAD (days)",
    caption = paste("n =", nrow(gam_data), "observations |", 
                   length(unique(gam_data$siteID)), "sites")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
  )

print(gam_plot)

# ============================================================================
# GAM DIAGNOSTICS
# ============================================================================

cat("\n=== GAM MODEL DIAGNOSTICS ===\n")
gam.check(gam_model)

cat("\n--- GAM ANALYSIS COMPLETE ---\n")
cat("✓ Latitude effect p-value:", format.pval(latitude_p, digits = 3), "\n")
cat("✓ Model explains", round(dev_explained * 100, 1), "% of deviance\n")

# Find latitude with maximum predicted PAD
max_pad_idx <- which.max(pred_data$fit)
max_pad_latitude <- pred_data$latitude[max_pad_idx]
max_pad_value <- pred_data$fit[max_pad_idx]

cat("\n--- PEAK PAD LOCATION ---\n")
cat("✓ Maximum predicted PAD occurs at latitude:", round(max_pad_latitude, 2), "°N\n")
cat("✓ Maximum predicted PAD value:", round(max_pad_value, 1), "days\n")