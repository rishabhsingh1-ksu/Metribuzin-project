# ==============================================================================
# Generalized Additive Model (GAM) Analysis: Weed Control & Crop Injury
# ==============================================================================
#
# DESCRIPTION:
# This script performs GAM analysis to evaluate:
#   1. Amaranthus weed control response to metribuzin dose (Fig 3, Table 5)
#   2. Soybean crop injury response to metribuzin dose (Fig 4, Table 6)
#
# Figure 3 & Table 5: All site-years (excluding PCA outliers IL'23, MI'23)
# Figure 4 & Table 6: PCA outlier locations only (IL'23, MI'23)
#
# ASSOCIATED PUBLICATION:
# Singh et al. (2025). Optimizing metribuzin rates for herbicide-resistant 
# Amaranthus weed control in soybean. Weed Technology, 39, e97.
# DOI: https://doi.org/10.1017/wet.2025.10047
#
# STATISTICAL METHOD:
# Generalized Additive Model (GAM) with cubic spline smoothers
# - Response variable: Weed control OR Crop injury (% of non-treated)
# - Predictor: Metribuzin dose (g ai ha⁻¹)
# - Grouped by: Days after application (14, 28, 42 DAA)
#
# DATA FILES REQUIRED:
# - data/weed_control.xlsx (weed control data)
# - data/soy_para.xlsx (crop injury data)
#
# OUTPUT:
# - Figure 3: Weed control GAM (all locations except IL'23, MI'23)
# - Table 5: Predicted weed control values
# - Figure 4: Weed control GAM (IL'23, MI'23 only - PCA outliers)
# - Table 6: Predicted weed control for outlier locations
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
# Install if needed: install.packages(c("readxl", "tidyverse", "mgcv", "RColorBrewer", "lubridate"))

library(readxl)       # Read Excel files
library(tidyverse)    # Data manipulation and ggplot2
library(dplyr)        # Data wrangling
library(ggplot2)      # Visualization
library(RColorBrewer) # Color palettes
library(mgcv)         # GAM modeling (for predictions)
library(lubridate)    # Date handling

# Set seed for reproducibility
set.seed(123)

# ------------------------------------------------------------------------------
# 2. DATA IMPORT
# ------------------------------------------------------------------------------

# Import weed control data
# Sheet 1 ("weed") contains weed control observations across all site-years
raw_data <- read_excel("data/weed.control.xlsx", sheet = 1)

# Import crop injury data
soy_data <- read_excel("data/soy_para.xlsx", sheet = 1)

# Preview data structure
cat("========== WEED CONTROL DATA ==========\n")
head(raw_data)
str(raw_data)

cat("\n========== CROP INJURY DATA ==========\n")
head(soy_data)
str(soy_data)

# ------------------------------------------------------------------------------
# 3. DATA PREPROCESSING - WEED CONTROL
# ------------------------------------------------------------------------------

# Filter and prepare data for GAM analysis
# - Keep only metribuzin treatments (trt 3-15)
# - trt 1 = Non-treated control
# - trt 2 = Weed-free control
# - trt 3-15 = Metribuzin doses (210-841 g ai/ha)
# - trt 16 = Sulfentrazone (420 g ai/ha)
# - trt 17 = S-metolachlor (1790 g ai/ha)

dat_hrbcd <- raw_data[complete.cases(raw_data), ] %>%
  filter(trt %in% c(3:15)) %>%
  mutate(
    rate = as.numeric(if_else(rate == "NT", "0.00", as.character(rate))),
    location = as.factor(location),
    moisture = as.factor(moisture),
    julian_day = yday(ymd(sowing))
  )

# Convert variables to appropriate types
dat_hrbcd$trt <- as.factor(dat_hrbcd$trt)
dat_hrbcd$rate <- as.numeric(dat_hrbcd$rate)
dat_hrbcd$weed <- as.factor(dat_hrbcd$weed)
dat_hrbcd$location <- as.factor(dat_hrbcd$location)
dat_hrbcd$control <- as.numeric(dat_hrbcd$control)

# Calculate metribuzin dose in g ai ha⁻¹
dat_hrbcd$dose <- dat_hrbcd$rate * 453.6 * 2.47

# Data summary
cat("\n========== WEED CONTROL DATA SUMMARY ==========\n")
cat("Total observations:", nrow(dat_hrbcd), "\n")
cat("Locations:", paste(unique(dat_hrbcd$location), collapse = ", "), "\n")
cat("Dose range:", min(dat_hrbcd$dose), "-", max(dat_hrbcd$dose), "g ai/ha\n")

# ------------------------------------------------------------------------------
# 4. DATA CLEANING FOR GAM
# ------------------------------------------------------------------------------

# Remove infinite/NA values and replace zeros with small values
dat_glm_all <- dat_hrbcd %>%
  filter(is.finite(rate) & !is.na(rate) & !is.infinite(rate)) %>%
  filter(if_any(where(is.numeric), ~ is.finite(.))) %>%
  mutate(
    rate = case_when(rate == 0 ~ 0.0001, TRUE ~ rate),
    control = case_when(control == 0 ~ 0.0001, TRUE ~ control)
  )

# ------------------------------------------------------------------------------
# 5. SPLIT DATA: MAIN LOCATIONS vs PCA OUTLIERS
# ------------------------------------------------------------------------------

# PCA analysis identified IL'23 and MI'23 as outlier locations
# These locations had distinct soil/environmental conditions
# See PCA_analysis.R for details

# Main dataset: Exclude PCA outlier locations (for Fig 3, Table 5)
dat_glm <- dat_glm_all %>% 
  filter(!location %in% c("IL'23", "MI'23"))

# Outlier dataset: Only PCA outlier locations (for Fig 4, Table 6)
dat_glm_outliers <- dat_glm_all %>% 
  filter(location %in% c("IL'23", "MI'23"))

cat("\n========== DATA SPLIT ==========\n")
cat("Main locations (Fig 3):", nrow(dat_glm), "observations\n")
cat("Outlier locations (Fig 4):", nrow(dat_glm_outliers), "observations\n")

# ------------------------------------------------------------------------------
# 6. FIGURE 3: WEED CONTROL GAM - MAIN LOCATIONS
# ------------------------------------------------------------------------------

cat("\n========== GENERATING FIGURE 3 ==========\n")

fig3_weed_control <- ggplot(dat_glm, aes(x = dose, y = control)) +
  # Individual data points (semi-transparent)
  geom_point(aes(color = as.factor(days)), alpha = 0.1, size = 1.5) +
  
  # GAM smooth curves with confidence intervals
  geom_smooth(
    aes(y = control, color = as.factor(days), fill = as.factor(days)),
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    linewidth = 1,
    alpha = 0.2
  ) +
  
  # Labels and theme
  labs(
    title = "Amaranthus Weed Control Response to Metribuzin Dose",
    subtitle = "All site-years excluding PCA outliers (IL'23, MI'23)",
    x = expression("Metribuzin dose (g ai ha"^-1*")"),
    y = "Weed control (% of non-treated)",
    color = "Days After\nApplication",
    fill = "Days After\nApplication"
  ) +
  
  # Theme customization
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 11),
    legend.position = "right"
  ) +
  
  # Color palette
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  
  # Axis scales
  scale_y_continuous(
    breaks = seq(0, 100, by = 20),
    limits = c(0, 110)
  ) +
  scale_x_continuous(
    breaks = c(210, 263, 315, 368, 420, 473, 525, 578, 630, 683, 736, 788, 841)
  )

print(fig3_weed_control)

# ------------------------------------------------------------------------------
# 7. TABLE 5: PREDICTED WEED CONTROL - MAIN LOCATIONS
# ------------------------------------------------------------------------------

cat("\n========== GENERATING TABLE 5 ==========\n")

# Subset data by evaluation timing
dat_14 <- dat_glm %>% filter(days == 14)
dat_28 <- dat_glm %>% filter(days == 28)
dat_42 <- dat_glm %>% filter(days == 42)

# Fit GAM models for each timing
gam_14 <- gam(control ~ s(dose, bs = "cs"), data = dat_14)
gam_28 <- gam(control ~ s(dose, bs = "cs"), data = dat_28)
gam_42 <- gam(control ~ s(dose, bs = "cs"), data = dat_42)

# Model summaries
cat("\n--- GAM Model Summary: 14 DAA ---\n")
summary(gam_14)

cat("\n--- GAM Model Summary: 28 DAA ---\n")
summary(gam_28)

cat("\n--- GAM Model Summary: 42 DAA ---\n")
summary(gam_42)

# Create prediction dataframe
dose_levels <- c(210, 263, 315, 368, 420, 473, 525, 578, 630, 683, 736, 788, 841)
pred_df <- data.frame(dose = dose_levels)

# Generate predictions with standard errors
pred_14 <- predict(gam_14, newdata = pred_df, se.fit = TRUE)
pred_28 <- predict(gam_28, newdata = pred_df, se.fit = TRUE)
pred_42 <- predict(gam_42, newdata = pred_df, se.fit = TRUE)

# Compile Table 5
table5 <- data.frame(
  Dose_g_ai_ha = dose_levels,
  Control_14DAA = round(pred_14$fit, 1),
  SE_14DAA = round(pred_14$se.fit, 2),
  Control_28DAA = round(pred_28$fit, 1),
  SE_28DAA = round(pred_28$se.fit, 2),
  Control_42DAA = round(pred_42$fit, 1),
  SE_42DAA = round(pred_42$se.fit, 2)
)

cat("\n========== TABLE 5: Predicted Weed Control (%) - Main Locations ==========\n")
print(table5)

# Find doses achieving target control levels
find_dose_for_control <- function(model, target, dose_range = c(200, 850)) {
  doses <- seq(dose_range[1], dose_range[2], by = 1)
  preds <- predict(model, newdata = data.frame(dose = doses))
  idx <- which.min(abs(preds - target))
  return(doses[idx])
}

cat("\n--- Doses (g ai/ha) for Target Control Levels ---\n")
targets <- c(80, 90, 95)
target_doses <- data.frame(
  Target = paste0(targets, "%"),
  Dose_14DAA = sapply(targets, function(t) find_dose_for_control(gam_14, t)),
  Dose_28DAA = sapply(targets, function(t) find_dose_for_control(gam_28, t)),
  Dose_42DAA = sapply(targets, function(t) find_dose_for_control(gam_42, t))
)
print(target_doses)

# ------------------------------------------------------------------------------
# 8. FIGURE 4: WEED CONTROL GAM - PCA OUTLIER LOCATIONS (IL'23, MI'23)
# ------------------------------------------------------------------------------

cat("\n========== GENERATING FIGURE 4 ==========\n")

fig4_weed_control_outliers <- ggplot(dat_glm_outliers, aes(x = dose, y = control)) +
  # Individual data points
  geom_point(aes(color = as.factor(days)), alpha = 0.15, size = 1.5) +
  
  # GAM smooth curves
  geom_smooth(
    aes(y = control, color = as.factor(days), fill = as.factor(days)),
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    linewidth = 1,
    alpha = 0.2
  ) +
  
  # Labels and theme
  labs(
    title = "Amaranthus Weed Control - PCA Outlier Locations",
    subtitle = "Illinois 2023 and Michigan 2023 site-years",
    x = expression("Metribuzin dose (g ai ha"^-1*")"),
    y = "Weed control (% of non-treated)",
    color = "Days After\nApplication",
    fill = "Days After\nApplication"
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 11),
    legend.position = "right"
  ) +
  
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  
  scale_y_continuous(
    breaks = seq(0, 100, by = 20),
    limits = c(0, 110)
  ) +
  scale_x_continuous(
    breaks = c(210, 263, 315, 368, 420, 473, 525, 578, 630, 683, 736, 788, 841)
  )

print(fig4_weed_control_outliers)

# ------------------------------------------------------------------------------
# 9. TABLE 6: PREDICTED WEED CONTROL - OUTLIER LOCATIONS
# ------------------------------------------------------------------------------

cat("\n========== GENERATING TABLE 6 ==========\n")

# Subset outlier data by timing
dat_out_14 <- dat_glm_outliers %>% filter(days == 14)
dat_out_28 <- dat_glm_outliers %>% filter(days == 28)
dat_out_42 <- dat_glm_outliers %>% filter(days == 42)

# Check if sufficient data exists for GAM fitting
cat("Outlier data points: 14 DAA =", nrow(dat_out_14), 
    ", 28 DAA =", nrow(dat_out_28), 
    ", 42 DAA =", nrow(dat_out_42), "\n")

# Fit GAM models for outlier locations
# Using try() to handle cases with insufficient data
gam_out_14 <- tryCatch(
  gam(control ~ s(dose, bs = "cs", k = 5), data = dat_out_14),
  error = function(e) lm(control ~ dose, data = dat_out_14)
)

gam_out_28 <- tryCatch(
  gam(control ~ s(dose, bs = "cs", k = 5), data = dat_out_28),
  error = function(e) lm(control ~ dose, data = dat_out_28)
)

gam_out_42 <- tryCatch(
  gam(control ~ s(dose, bs = "cs", k = 5), data = dat_out_42),
  error = function(e) lm(control ~ dose, data = dat_out_42)
)

# Generate predictions
pred_out_14 <- predict(gam_out_14, newdata = pred_df, se.fit = TRUE)
pred_out_28 <- predict(gam_out_28, newdata = pred_df, se.fit = TRUE)
pred_out_42 <- predict(gam_out_42, newdata = pred_df, se.fit = TRUE)

# Compile Table 6
table6 <- data.frame(
  Dose_g_ai_ha = dose_levels,
  Control_14DAA = round(pred_out_14$fit, 1),
  SE_14DAA = round(pred_out_14$se.fit, 2),
  Control_28DAA = round(pred_out_28$fit, 1),
  SE_28DAA = round(pred_out_28$se.fit, 2),
  Control_42DAA = round(pred_out_42$fit, 1),
  SE_42DAA = round(pred_out_42$se.fit, 2)
)

cat("\n========== TABLE 6: Predicted Weed Control (%) - Outlier Locations ==========\n")
print(table6)

# ------------------------------------------------------------------------------
# 10. CROP INJURY GAM ANALYSIS
# ------------------------------------------------------------------------------

cat("\n========== CROP INJURY ANALYSIS ==========\n")

# Prepare crop injury data
soy_data$reps <- as.factor(soy_data$reps)
soy_data$trt <- as.factor(soy_data$trt)
soy_data$rate <- as.numeric(soy_data$rate)
soy_data$weed <- as.factor(soy_data$weed)
soy_data$location <- as.factor(soy_data$location)

# Filter for metribuzin treatments only
sp_mtz <- soy_data %>% 
  filter(trt %in% c(3:15))

# Calculate dose
sp_mtz$dose <- sp_mtz$rate * 453.6 * 2.47

# Create crop injury plot (similar to Fig 4 in publication)
fig_crop_injury <- sp_mtz %>%
  ggplot(aes(x = dose)) +
  # 14 DAA injury curve
  geom_smooth(aes(y = ci.14dat, color = "14 DAA"), 
              method = "gam", formula = y ~ s(x, bs = "cs"),
              linewidth = 1, se = TRUE, alpha = 0.2) +
  # 28 DAA injury curve
  geom_smooth(aes(y = ci.28dat, color = "28 DAA"), 
              method = "gam", formula = y ~ s(x, bs = "cs"),
              linewidth = 1, se = TRUE, alpha = 0.2) +
  # 42 DAA injury curve
  geom_smooth(aes(y = ci.42dat, color = "42 DAA"), 
              method = "gam", formula = y ~ s(x, bs = "cs"),
              linewidth = 1, se = TRUE, alpha = 0.2) +
  
  # 5% injury threshold line
  geom_hline(yintercept = 5, alpha = 1, color = "black", linetype = "dotted") +
  annotate("text", x = 250, y = 7, label = "5% injury threshold", 
           size = 3, hjust = 0) +
  
  # Labels and theme
  labs(
    title = "Soybean Crop Injury Response to Metribuzin Dose",
    subtitle = "Injury remained ≤5% even at highest dose (841 g ai/ha)",
    x = expression("Metribuzin dose (g ai ha"^-1*")"),
    y = "Crop injury (% of non-treated)",
    color = "Days After\nApplication"
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 11),
    legend.position = "right"
  ) +
  
  scale_color_manual(values = c("14 DAA" = "#E41A1C", 
                                 "28 DAA" = "#377EB8", 
                                 "42 DAA" = "#4DAF4A")) +
  
  scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 100)) +
  scale_x_continuous(
    breaks = c(210, 263, 315, 368, 420, 473, 525, 578, 630, 683, 736, 788, 841)
  )

print(fig_crop_injury)

# ------------------------------------------------------------------------------
# 11. COMPARISON WITH SULFENTRAZONE AND S-METOLACHLOR
# ------------------------------------------------------------------------------

cat("\n========== HERBICIDE COMPARISON ANALYSIS ==========\n")

# Load required packages for mean comparisons
library(emmeans)
library(multcomp)

# Prepare comparison data (trt 16 = Sulfentrazone, trt 17 = S-metolachlor)
dat_comparison <- raw_data[complete.cases(raw_data), ] %>%
  filter(trt %in% c(9, 16, 17)) %>%  # trt 9 = 525 g/ha metribuzin for comparison
  mutate(
    rate = as.numeric(if_else(rate == "NT", "0.00", as.character(rate))),
    trt = as.factor(trt)
  )

# Compare weed control at 14 DAA
sp_14 <- dat_comparison %>% filter(days == 14)
model_14 <- lm(control ~ trt, data = sp_14)

cat("\n--- Weed Control Comparison at 14 DAA ---\n")
summary(model_14)

res_14 <- emmeans(model_14, ~ trt)
res_14_cld <- cld(res_14, Letters = letters, reversed = TRUE) %>%
  as.data.frame() %>%
  mutate(letters = trimws(.group))
print(res_14_cld)

# Compare weed control at 28 DAA
sp_28 <- dat_comparison %>% filter(days == 28)
model_28 <- lm(control ~ trt, data = sp_28)

cat("\n--- Weed Control Comparison at 28 DAA ---\n")
res_28 <- emmeans(model_28, ~ trt)
res_28_cld <- cld(res_28, Letters = letters, reversed = TRUE) %>%
  as.data.frame() %>%
  mutate(letters = trimws(.group))
print(res_28_cld)

# Compare weed control at 42 DAA
sp_42 <- dat_comparison %>% filter(days == 42)
model_42 <- lm(control ~ trt, data = sp_42)

cat("\n--- Weed Control Comparison at 42 DAA ---\n")
res_42 <- emmeans(model_42, ~ trt)
res_42_cld <- cld(res_42, Letters = letters, reversed = TRUE) %>%
  as.data.frame() %>%
  mutate(letters = trimws(.group))
print(res_42_cld)

# ------------------------------------------------------------------------------
# 12. SAVE OUTPUTS
# ------------------------------------------------------------------------------

# Create output directories if they don't exist
dir.create("figures", showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)

# Save Figure 3
ggsave("figures/Fig3_Weed_Control_GAM_Main.png", 
       plot = fig3_weed_control,
       width = 10, height = 7, dpi = 300)
ggsave("figures/Fig3_Weed_Control_GAM_Main.pdf", 
       plot = fig3_weed_control,
       width = 10, height = 7)

# Save Figure 4
ggsave("figures/Fig4_Weed_Control_GAM_Outliers.png", 
       plot = fig4_weed_control_outliers,
       width = 10, height = 7, dpi = 300)
ggsave("figures/Fig4_Weed_Control_GAM_Outliers.pdf", 
       plot = fig4_weed_control_outliers,
       width = 10, height = 7)

# Save Crop Injury Figure
ggsave("figures/Crop_Injury_GAM.png", 
       plot = fig_crop_injury,
       width = 10, height = 7, dpi = 300)
ggsave("figures/Crop_Injury_GAM.pdf", 
       plot = fig_crop_injury,
       width = 10, height = 7)

# Save Tables
write.csv(table5, "output/Table5_Weed_Control_Predictions_Main.csv", row.names = FALSE)
write.csv(table6, "output/Table6_Weed_Control_Predictions_Outliers.csv", row.names = FALSE)

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
