# ==============================================================================
# Weed Parameters Analysis: Density, Biomass, and Emergence
# ==============================================================================
#
# DESCRIPTION:
# This script analyzes weed response parameters as a function of metribuzin dose:
#   - Figure 5: Weed density GAM response
#   - Figure 6: Weed biomass GAM response  
#   - Figure 7: Weed emergence timing GAM response
#
# Also includes comparison of metribuzin vs sulfentrazone and S-metolachlor
# using linear mixed models with Tukey HSD tests.
#
# ASSOCIATED PUBLICATION:
# Singh et al. (2025). Optimizing metribuzin rates for herbicide-resistant 
# Amaranthus weed control in soybean. Weed Technology, 39, e97.
# DOI: https://doi.org/10.1017/wet.2025.10047
#
# STATISTICAL METHODS:
# 1. GAM (Generalized Additive Model) for dose-response visualization
# 2. Linear Mixed Model (LMER) for treatment comparisons
#    - Fixed effect: Treatment
#    - Random effects: Location, Weed species
# 3. Tukey HSD for pairwise comparisons
#
# DATA FILE REQUIRED:
# - data/weed_parameters.xlsx
#
# VARIABLES:
# - density: Weed density (plants per m²) at 28 DAA
# - biomass: Dry aboveground biomass (g per m²) at 28 DAA
# - emergence: Days after application to first weed emergence
#
# OUTPUT:
# - Figure 5: Weed density vs metribuzin dose
# - Figure 6: Weed biomass vs metribuzin dose
# - Figure 7: Weed emergence timing vs metribuzin dose
# - Boxplots with Tukey letters for treatment comparisons
#
# AUTHOR: Rishabh Singh
# DATE: 2025
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. SETUP: Load Required Packages
# ------------------------------------------------------------------------------

# Clear workspace
rm(list = ls())

# Load required packages
library(readxl)       # Read Excel files
library(tidyverse)    # Data manipulation and ggplot2
library(dplyr)        # Data wrangling
library(ggplot2)      # Visualization
library(lme4)         # Linear mixed models
library(lmerTest)     # P-values for lmer
library(emmeans)      # Estimated marginal means
library(multcomp)     # Multiple comparisons (cld function)
library(broom)        # Tidy model outputs

# Set seed for reproducibility
set.seed(123)

# ------------------------------------------------------------------------------
# 2. DATA IMPORT
# ------------------------------------------------------------------------------

# Import weed parameters data
weed_data <- read_excel("data/weed_parameters.xlsx", sheet = 1)

# Preview data structure
cat("========== WEED PARAMETERS DATA ==========\n")
head(weed_data)
str(weed_data)
names(weed_data)

# ------------------------------------------------------------------------------
# 3. DATA PREPROCESSING
# ------------------------------------------------------------------------------

# Convert variables to appropriate types
weed_data$trt <- as.factor(weed_data$trt)
weed_data$rate <- as.numeric(weed_data$rate)
weed_data$weed <- as.factor(weed_data$weed)
weed_data$location <- as.factor(weed_data$location)
weed_data$reps <- as.factor(weed_data$reps)
weed_data$density <- as.numeric(weed_data$density)
weed_data$biomass <- as.numeric(weed_data$biomass)
weed_data$emergence <- as.numeric(weed_data$emergence)

# Calculate dose in g ai ha⁻¹
weed_data$dose <- weed_data$rate * 453.6 * 2.47

# Data summary
cat("\n========== DATA SUMMARY ==========\n")
cat("Total observations:", nrow(weed_data), "\n")
cat("Locations:", paste(unique(weed_data$location), collapse = ", "), "\n")
cat("Treatments:", paste(levels(weed_data$trt), collapse = ", "), "\n")
str(weed_data)

# ------------------------------------------------------------------------------
# 4. CHECK DATA DISTRIBUTION
# ------------------------------------------------------------------------------

# Histograms to assess data distribution
par(mfrow = c(1, 3))
hist(weed_data$density, main = "Weed Density Distribution", 
     xlab = "Density (plants/m²)", col = "lightblue", breaks = 30)
hist(weed_data$biomass, main = "Weed Biomass Distribution", 
     xlab = "Biomass (g/m²)", col = "lightgreen", breaks = 30)
hist(weed_data$emergence, main = "Weed Emergence Distribution", 
     xlab = "Days to Emergence", col = "lightyellow", breaks = 20)
par(mfrow = c(1, 1))

# Note: Density and biomass are typically right-skewed
# GAM handles non-normal distributions well

# ------------------------------------------------------------------------------
# 5. EXPLORATORY BOXPLOTS BY LOCATION
# ------------------------------------------------------------------------------

# Weed density by location (diagnostic)
p_density_loc <- ggplot(weed_data, aes(x = trt, y = density)) +
  geom_boxplot() +
  facet_wrap(~ location, scales = "free_y") +
  labs(title = "Weed Density by Treatment and Location",
       x = "Treatment", y = "Density (plants/m²)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

# Weed biomass by location
p_biomass_loc <- ggplot(weed_data, aes(x = trt, y = biomass)) +
  geom_boxplot() +
  facet_wrap(~ location, scales = "free_y") +
  labs(title = "Weed Biomass by Treatment and Location",
       x = "Treatment", y = "Biomass (g/m²)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

# Weed emergence by location
p_emergence_loc <- ggplot(weed_data, aes(x = trt, y = emergence)) +
  geom_boxplot() +
  facet_wrap(~ location, scales = "free_y") +
  labs(title = "Weed Emergence by Treatment and Location",
       x = "Treatment", y = "Days to Emergence") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

# Display diagnostic plots (uncomment to view)
# print(p_density_loc)
# print(p_biomass_loc)
# print(p_emergence_loc)

# ------------------------------------------------------------------------------
# 6. PREPARE DATA FOR GAM ANALYSIS
# ------------------------------------------------------------------------------

# Filter data for metribuzin treatments only (exclude controls and comparators)
# trt 1 = Non-treated, trt 2 = Weed-free (excluded from this file)
# trt 3-15 = Metribuzin doses
# trt 16 = Sulfentrazone, trt 17 = S-metolachlor

# For GAM plots: metribuzin treatments only
sp_mtz <- weed_data %>% 
  filter(!trt %in% c("1", "2", "16", "17"))

cat("\n========== GAM DATA ==========\n")
cat("Metribuzin-only observations:", nrow(sp_mtz), "\n")

# ------------------------------------------------------------------------------
# 7. FIGURE 5: WEED DENSITY GAM
# ------------------------------------------------------------------------------

cat("\n========== GENERATING FIGURE 5: WEED DENSITY ==========\n")

# Calculate reference lines (mean values for comparison herbicides)
# This requires fitting a model to get treatment means
sp_comparison <- weed_data %>% filter(trt %in% c("1", "16", "17"))

if (nrow(sp_comparison) > 0) {
  density_means <- sp_comparison %>%
    group_by(trt) %>%
    summarize(mean_density = mean(density, na.rm = TRUE), .groups = "drop")
  cat("\nDensity means for comparison:\n")
  print(density_means)
}

fig5_density <- sp_mtz %>%
  ggplot(aes(x = dose, y = density)) +
  # GAM smooth curve
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              color = "blue", fill = "lightblue", alpha = 0.3, linewidth = 1.2) +
  # Individual data points
  geom_jitter(color = "brown", width = 5, height = 0, alpha = 0.4, size = 2) +
  
  # Reference lines for sulfentrazone and S-metolachlor (if available)
  # Uncomment and adjust yintercept values based on your data
  # geom_hline(yintercept = 21.4, alpha = 0.8, color = "darkgreen", linetype = "dashed") +
  # geom_hline(yintercept = 66.2, alpha = 0.8, color = "purple", linetype = "dashed") +
  
  # Labels and theme
  labs(
    title = "Weed Density Response to Metribuzin Dose",
    subtitle = "Amaranthus density at 28 days after application",
    x = expression("Metribuzin dose (g ai ha"^-1*")"),
    y = expression("Weed density (plants m"^-2*")")
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 11),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.6)
  ) +
  
  scale_y_continuous(limits = c(0, 1500)) +
  scale_x_continuous(
    breaks = c(210, 263, 315, 368, 420, 473, 525, 578, 630, 683, 736, 788, 841)
  )

print(fig5_density)

# ------------------------------------------------------------------------------
# 8. FIGURE 6: WEED BIOMASS GAM
# ------------------------------------------------------------------------------

cat("\n========== GENERATING FIGURE 6: WEED BIOMASS ==========\n")

fig6_biomass <- sp_mtz %>%
  ggplot(aes(x = dose, y = biomass)) +
  # GAM smooth curve
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              color = "blue", fill = "lightblue", alpha = 0.3, linewidth = 1.2) +
  # Individual data points
  geom_jitter(color = "brown", width = 5, height = 0, alpha = 0.3, size = 2) +
  
  # Labels and theme
  labs(
    title = "Weed Biomass Response to Metribuzin Dose",
    subtitle = "Amaranthus dry biomass at 28 days after application",
    x = expression("Metribuzin dose (g ai ha"^-1*")"),
    y = expression("Dry biomass (g m"^-2*")")
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 11),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.6)
  ) +
  
  scale_y_continuous(limits = c(0, 150)) +
  scale_x_continuous(
    breaks = c(210, 263, 315, 368, 420, 473, 525, 578, 630, 683, 736, 788, 841)
  )

print(fig6_biomass)

# ------------------------------------------------------------------------------
# 9. FIGURE 7: WEED EMERGENCE GAM
# ------------------------------------------------------------------------------

cat("\n========== GENERATING FIGURE 7: WEED EMERGENCE ==========\n")

fig7_emergence <- sp_mtz %>%
  ggplot(aes(x = dose, y = emergence)) +
  # GAM smooth curve
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              color = "blue", fill = "lightblue", alpha = 0.3, linewidth = 1.2) +
  # Individual data points
  geom_jitter(color = "brown", width = 5, height = 0, alpha = 0.4, size = 2) +
  
  # Labels and theme
  labs(
    title = "Weed Emergence Delay Response to Metribuzin Dose",
    subtitle = "Days to first Amaranthus emergence after herbicide application",
    x = expression("Metribuzin dose (g ai ha"^-1*")"),
    y = "Days to emergence"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 11),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.6)
  ) +
  
  scale_y_continuous(limits = c(0, 60)) +
  scale_x_continuous(
    breaks = c(210, 263, 315, 368, 420, 473, 525, 578, 630, 683, 736, 788, 841)
  )

print(fig7_emergence)

# ------------------------------------------------------------------------------
# 10. TREATMENT COMPARISON: LMER + TUKEY HSD
# ------------------------------------------------------------------------------

cat("\n========== TREATMENT COMPARISONS (LMER + TUKEY) ==========\n")

# Prepare data for all treatments including controls
# Include trt 1 (non-treated) and metribuzin treatments
wp_all <- weed_data %>% 
  filter(trt %in% c("1", as.character(3:17)))

# --- WEED DENSITY ---
cat("\n--- Weed Density: LMER Analysis ---\n")

model_density <- lmer(density ~ trt + (1|location) + (1|weed), data = wp_all)
summary(model_density)

# Estimated marginal means and Tukey groups
res_density <- emmeans(model_density, ~ trt)
res_density_cld <- cld(res_density, Letters = letters, reversed = TRUE) %>%
  as.data.frame() %>%
  mutate(letters = trimws(.group))

cat("\nWeed Density - Tukey HSD Results:\n")
print(res_density_cld)

# Create boxplot with Tukey letters
wp_density_plot <- wp_all %>%
  left_join(res_density_cld, by = "trt")

fig_density_boxplot <- ggplot(wp_density_plot, aes(x = factor(trt), y = density)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue", alpha = 0.7) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.4, size = 1.5) +
  geom_text(data = res_density_cld, 
            aes(x = trt, y = max(wp_all$density, na.rm = TRUE) * 0.95, label = letters),
            vjust = -0.5, size = 4, fontface = "bold") +
  labs(
    title = "Weed Density by Treatment",
    subtitle = "Letters indicate Tukey HSD groupings (α = 0.05)",
    x = "Treatment",
    y = expression("Weed density (plants m"^-2*")")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(linewidth = 0.8)
  )

print(fig_density_boxplot)

# --- WEED BIOMASS ---
cat("\n--- Weed Biomass: LMER Analysis ---\n")

model_biomass <- lmer(biomass ~ trt + (1|location) + (1|weed), data = wp_all)
summary(model_biomass)

res_biomass <- emmeans(model_biomass, ~ trt)
res_biomass_cld <- cld(res_biomass, Letters = letters, reversed = TRUE) %>%
  as.data.frame() %>%
  mutate(letters = trimws(.group))

cat("\nWeed Biomass - Tukey HSD Results:\n")
print(res_biomass_cld)

# Create boxplot with Tukey letters
wp_biomass_plot <- wp_all %>%
  left_join(res_biomass_cld, by = "trt")

fig_biomass_boxplot <- ggplot(wp_biomass_plot, aes(x = factor(trt), y = biomass)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgreen", alpha = 0.7) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.4, size = 1.5) +
  geom_text(data = res_biomass_cld, 
            aes(x = trt, y = max(wp_all$biomass, na.rm = TRUE) * 0.95, label = letters),
            vjust = -0.5, size = 4, fontface = "bold") +
  labs(
    title = "Weed Biomass by Treatment",
    subtitle = "Letters indicate Tukey HSD groupings (α = 0.05)",
    x = "Treatment",
    y = expression("Dry biomass (g m"^-2*")")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(linewidth = 0.8)
  )

print(fig_biomass_boxplot)

# --- WEED EMERGENCE ---
cat("\n--- Weed Emergence: LMER Analysis ---\n")

model_emergence <- lmer(emergence ~ trt + (1|location) + (1|weed) + (1|reps), data = wp_all)
summary(model_emergence)

res_emergence <- emmeans(model_emergence, ~ trt)
res_emergence_cld <- cld(res_emergence, Letters = letters, reversed = TRUE) %>%
  as.data.frame() %>%
  mutate(letters = trimws(.group))

cat("\nWeed Emergence - Tukey HSD Results:\n")
print(res_emergence_cld)

# Create boxplot with Tukey letters
wp_emergence_plot <- wp_all %>%
  left_join(res_emergence_cld, by = "trt")

fig_emergence_boxplot <- ggplot(wp_emergence_plot, aes(x = factor(trt), y = emergence)) +
  geom_boxplot(outlier.shape = NA, fill = "lightyellow", alpha = 0.7) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.4, size = 1.5) +
  geom_text(data = res_emergence_cld, 
            aes(x = trt, y = max(wp_all$emergence, na.rm = TRUE) * 0.95, label = letters),
            vjust = -0.5, size = 4, fontface = "bold") +
  labs(
    title = "Weed Emergence Timing by Treatment",
    subtitle = "Letters indicate Tukey HSD groupings (α = 0.05)",
    x = "Treatment",
    y = "Days to emergence"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(linewidth = 0.8)
  )

print(fig_emergence_boxplot)

# ------------------------------------------------------------------------------
# 11. COMPARISON: METRIBUZIN vs SULFENTRAZONE vs S-METOLACHLOR
# ------------------------------------------------------------------------------

cat("\n========== HERBICIDE COMPARISON ==========\n")

# Compare non-treated (1), sulfentrazone (16), and S-metolachlor (17)
sp_herb_comparison <- weed_data %>% 
  filter(trt %in% c("1", "16", "17"))

cat("Observations for comparison:", nrow(sp_herb_comparison), "\n")

# Density comparison
model_herb_density <- lm(density ~ trt, data = sp_herb_comparison)
cat("\n--- Density: Herbicide Comparison ---\n")
summary(model_herb_density)

res_herb_density <- emmeans(model_herb_density, ~ trt)
res_herb_density_cld <- cld(res_herb_density, Letters = letters, reversed = TRUE) %>%
  as.data.frame() %>%
  mutate(letters = trimws(.group))
print(res_herb_density_cld)

# Emergence comparison
model_herb_emergence <- lm(emergence ~ trt, data = sp_herb_comparison)
cat("\n--- Emergence: Herbicide Comparison ---\n")
summary(model_herb_emergence)

res_herb_emergence <- emmeans(model_herb_emergence, ~ trt)
res_herb_emergence_cld <- cld(res_herb_emergence, Letters = letters, reversed = TRUE) %>%
  as.data.frame() %>%
  mutate(letters = trimws(.group))
print(res_herb_emergence_cld)

# ------------------------------------------------------------------------------
# 12. SAVE OUTPUTS
# ------------------------------------------------------------------------------

# Create output directories
dir.create("figures", showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)

# Save GAM figures
ggsave("figures/Fig5_Weed_Density_GAM.png", plot = fig5_density,
       width = 10, height = 7, dpi = 300)
ggsave("figures/Fig5_Weed_Density_GAM.pdf", plot = fig5_density,
       width = 10, height = 7)

ggsave("figures/Fig6_Weed_Biomass_GAM.png", plot = fig6_biomass,
       width = 10, height = 7, dpi = 300)
ggsave("figures/Fig6_Weed_Biomass_GAM.pdf", plot = fig6_biomass,
       width = 10, height = 7)

ggsave("figures/Fig7_Weed_Emergence_GAM.png", plot = fig7_emergence,
       width = 10, height = 7, dpi = 300)
ggsave("figures/Fig7_Weed_Emergence_GAM.pdf", plot = fig7_emergence,
       width = 10, height = 7)

# Save boxplots
ggsave("figures/Weed_Density_Boxplot_Tukey.png", plot = fig_density_boxplot,
       width = 12, height = 7, dpi = 300)
ggsave("figures/Weed_Biomass_Boxplot_Tukey.png", plot = fig_biomass_boxplot,
       width = 12, height = 7, dpi = 300)
ggsave("figures/Weed_Emergence_Boxplot_Tukey.png", plot = fig_emergence_boxplot,
       width = 12, height = 7, dpi = 300)

# Save Tukey results
write.csv(res_density_cld, "output/Weed_Density_Tukey_Results.csv", row.names = FALSE)
write.csv(res_biomass_cld, "output/Weed_Biomass_Tukey_Results.csv", row.names = FALSE)
write.csv(res_emergence_cld, "output/Weed_Emergence_Tukey_Results.csv", row.names = FALSE)

cat("\n========== FILES SAVED ==========\n")
cat("Figures saved to: figures/\n")
cat("Tables saved to: output/\n")

# ------------------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------------------

cat("\n========== SESSION INFO ==========\n")
sessionInfo()

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
