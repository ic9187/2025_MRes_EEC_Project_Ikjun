# Individual models
summary(phenology_models[["host_only"]])           # Host density only
summary(phenology_models[["vpd_only"]])            # VPD max only  
summary(phenology_models[["host_vpd_interaction"]])  # Full interaction model

# Or see all model names first
names(phenology_models)



# Individual models
summary(prevalence_models[["pad_only"]])     # PAD predictor
summary(prevalence_models[["nymph_only"]])   # Nymph count predictor
summary(prevalence_models[["host_only"]])    # Host density predictor
summary(prevalence_models[["vpd_only"]])     # VPD max predictor

# Or see all model names
names(prevalence_models)



# For mixed models (phenology)
performance::r2(phenology_models[["vpd_only"]])  # R² marginal/conditional
anova(phenology_models[["vpd_only"]])           # ANOVA table

# For linear models (prevalence)  
summary(prevalence_models[["pad_only"]])$r.squared          # R²
anova(prevalence_models[["pad_only"]])                      # ANOVA table

# Coefficients only
coef(summary(phenology_models[["vpd_only"]]))   # Mixed model coefficients
coef(summary(prevalence_models[["pad_only"]]))              # Linear model coefficients

