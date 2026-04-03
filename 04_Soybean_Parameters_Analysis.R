# ==============================================================================
# Soybean Parameters Analysis: Height and Yield
# ==============================================================================
#
# DESCRIPTION:
# This script analyzes soybean crop response parameters:
#   - Figure 8: Soybean height by treatment (boxplot with Tukey HSD)
#   - Figure 9: Soybean yield by treatment (boxplot with Tukey HSD)
#   - Crop injury GAM analysis (supplementary)
#
# ASSOCIATED PUBLICATION:
# Singh et al. (2025). Optimizing metribuzin rates for herbicide-resistant 
# Amaranthus weed control in soybean. Weed Technology, 39, e97.
# DOI: https://doi.org/10.1017/wet.2025.10047
#
# STATISTICAL METHODS:
# 1. Linear Model (LM) for soybean height analysis
# 2. Linear Mixed Model (LMER) for soybean yield analysis
#    - Fixed effect: Treatment
#    - Random effects: Location, Weed species, Replication
# 3. Tukey HSD for pairwise treatment comparisons
# 4. GAM for crop injury visualization
#
# DATA FILE REQUIRED:
# - data/soy_para.xlsx
#
# VARIABLES:
# - ci.14dat, ci.28dat, ci.42dat: Crop injury (%) at 14, 28, 42 DAA
# - soy.ht.cm: Soybean height (cm) at ~28 DAA
# - soy.yield.bu/a: Soybean yield (bushels per acre)
#
# KEY FINDINGS:
# - Soybean height and yield were not significantly affected by metribuzin
# - Crop injury remained ≤5% even at highest doses
#
# OUTPUT:
# - Figure 8: Soybean height boxplot with Tukey letters
# - Figure 9: Soybean yield boxplot with Tukey letters
# - Crop injury GAM figure
# - Tukey HSD results tables
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
library(car)          # Anova function for lmer
library(broom)        # Tidy model outputs

# Set seed for reproducibility
set.seed(123)

# ------------------------------------------------------------------------------
# 2. DATA IMPORT
# ------------------------------------------------------------------------------

# Import soybean parameters data
soy_data <- read_excel("data/soy_para.xlsx", sheet = 1)

# Preview data structure
cat("========== SOYBEAN PARAMETERS DATA ==========\n")
head(soy_data)
str(soy_data)
names(soy_data)

# ------------------------------------------------------------------------------
# 3. DATA PREPROCESSING
# ------------------------------------------------------------------------------

# Convert variables to appropriate types
soy_data$reps <- as.factor(soy_data$reps)
soy_data$trt <- as.factor(soy_data$trt)
soy_data$rate <- as.numeric(soy_data$rate)
soy_data$weed <- as.factor(soy_data$weed)
soy_data$location <- as.factor(soy_data$location)
soy_data$`soy.yield.bu/a` <- as.numeric(soy_data$`soy.yield.bu/a`)

# Calculate dose in g ai ha⁻¹
soy_data$dose <- soy_data$rate * 453.6 * 2.47

# Data summary
cat("\n========== DATA SUMMARY ==========\n")
cat("Total observations:", nrow(soy_data), "\n")
cat("Locations:", paste(unique(soy_data$location), collapse = ", "), "\n")
cat("Treatments:", paste(levels(soy_data$trt), collapse = ", "), "\n")

# Check available data for height and yield
cat("\n--- Data Availability ---\n")
cat("Soybean height observations:", sum(!is.na(soy_data$soy.ht.cm)), "\n")
cat("Soybean yield observations:", sum(!is.na(soy_data$`soy.yield.bu/a`)), "\n")

# ------------------------------------------------------------------------------
# 4. CHECK DATA DISTRIBUTION
# ------------------------------------------------------------------------------

# Histograms to assess data distribution
par(mfrow = c(2, 3))
hist(soy_data$ci.14dat, main = "Crop Injury 14 DAA", 
     xlab = "Injury (%)", col = "lightcoral", breaks = 20)
hist(soy_data$ci.28dat, main = "Crop Injury 28 DAA", 
     xlab = "Injury (%)", col = "lightyellow", breaks = 20)
hist(soy_data$ci.42dat, main = "Crop Injury 42 DAA", 
     xlab = "Injury (%)", col = "lightgreen", breaks = 20)
hist(soy_data$soy.ht.cm, main = "Soybean Height", 
     xlab = "Height (cm)", col = "lightblue", breaks = 20)
hist(soy_data$`soy.yield.bu/a`, main = "Soybean Yield", 
     xlab = "Yield (bu/A)", col = "plum", breaks = 20)
par(mfrow = c(1, 1))

# ------------------------------------------------------------------------------
# 5. EXPLORATORY BOXPLOTS BY LOCATION
# ------------------------------------------------------------------------------

# Crop injury by location (14 DAA)
p_ci14 <- ggplot(soy_data, aes(x = trt, y = ci.14dat)) +
  geom_boxplot() +
  facet_wrap(~ location) +
  labs(title = "Crop Injury at 14 DAA by Location",
       x = "Treatment", y = "Crop Injury (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

# Crop injury by location (28 DAA)
p_ci28 <- ggplot(soy_data, aes(x = trt, y = ci.28dat)) +
  geom_boxplot() +
  facet_wrap(~ location) +
  labs(title = "Crop Injury at 28 DAA by Location",
       x = "Treatment", y = "Crop Injury (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

# Crop injury by location (42 DAA)
p_ci42 <- ggplot(soy_data, aes(x = trt, y = ci.42dat)) +
  geom_boxplot() +
  facet_wrap(~ location) +
  labs(title = "Crop Injury at 42 DAA by Location",
       x = "Treatment", y = "Crop Injury (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

# Display diagnostic plots (uncomment to view)
# print(p_ci14)
# print(p_ci28)
# print(p_ci42)

# ------------------------------------------------------------------------------
# 6. IDENTIFY LOCATIONS WITH HEIGHT AND YIELD DATA
# ------------------------------------------------------------------------------

# Check which locations have height data
height_by_loc <- soy_data %>%
  group_by(location) %>%
  summarize(
    n_height = sum(!is.na(soy.ht.cm)),
    n_yield = sum(!is.na(`soy.yield.bu/a`)),
    .groups = "drop"
  )

cat("\n--- Data Availability by Location ---\n")
print(height_by_loc)

# Filter locations with height data
# Exclude locations with missing or insufficient height data
locations_with_height <- height_by_loc %>%
  filter(n_height > 0) %>%
  pull(location)

soy_ht_data <- soy_data %>%
  filter(location %in% locations_with_height) %>%
  filter(!is.na(soy.ht.cm))

cat("\nLocations with height data:", paste(locations_with_height, collapse = ", "), "\n")
cat("Height observations:", nrow(soy_ht_data), "\n")

# Filter locations with yield data
locations_with_yield <- height_by_loc %>%
  filter(n_yield > 0) %>%
  pull(location)

soy_yd_data <- soy_data %>%
  filter(location %in% locations_with_yield) %>%
  filter(!is.na(`soy.yield.bu/a`))

cat("\nLocations with yield data:", paste(locations_with_yield, collapse = ", "), "\n")
cat("Yield observations:", nrow(soy_yd_data), "\n")

# ------------------------------------------------------------------------------
# 7. FIGURE 8: SOYBEAN HEIGHT BY TREATMENT
# ------------------------------------------------------------------------------

cat("\n========== FIGURE 8: SOYBEAN HEIGHT ANALYSIS ==========\n")

# Linear model for soybean height
# Using -1 to get treatment means directly
model_height <- lm(soy.ht.cm ~ trt - 1, data = soy_ht_data)

cat("\n--- Model Summary ---\n")
summary(model_height)

cat("\n--- ANOVA ---\n")
print(anova(model_height))

# Estimated marginal means and Tukey HSD
res_height <- emmeans(model_height, ~ trt)
res_height_cld <- cld(res_height, Letters = letters, reversed = TRUE) %>%
  as.data.frame() %>%
  mutate(letters = trimws(.group))

cat("\n--- Tukey HSD Results ---\n")
print(res_height_cld)

# Merge Tukey letters with data for plotting
soy_ht_plot <- soy_ht_data %>%
  left_join(res_height_cld, by = "trt")

# Create Figure 8
fig8_height <- ggplot(soy_ht_plot, aes(x = factor(trt), y = soy.ht.cm)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue", alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.3, size = 1.5) +
  geom_text(data = res_height_cld, 
            aes(x = trt, y = max(soy_ht_data$soy.ht.cm, na.rm = TRUE) + 2, label = letters),
            vjust = -0.5, size = 5, fontface = "bold") +
  
  labs(
    title = "Soybean Height by Treatment",
    subtitle = "Letters indicate Tukey HSD groupings (α = 0.05)",
    x = "Treatment",
    y = "Soybean height (cm)"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.6)
  ) +
  
  ylim(0, max(soy_ht_data$soy.ht.cm, na.rm = TRUE) + 8)

print(fig8_height)

# ------------------------------------------------------------------------------
# 8. FIGURE 9: SOYBEAN YIELD BY TREATMENT
# ------------------------------------------------------------------------------

cat("\n========== FIGURE 9: SOYBEAN YIELD ANALYSIS ==========\n")

# Convert yield to kg/ha for international units
# 1 bu/A × 27.22 kg/bu × 2.47 A/ha = kg/ha
soy_yd_data$soy.yd.kg.ha <- soy_yd_data$`soy.yield.bu/a` * 27.22 * 2.47

# Note: Some locations may need to be excluded due to data quality
# Uncomment the following line if needed:
# soy_yd_data <- soy_yd_data %>% filter(!location %in% "KS'23")

# Linear Mixed Model for soybean yield
# Random effects account for location, weed species, and replicate variation
model_yield <- lmer(soy.yd.kg.ha ~ trt + (1|weed) + (1|location) + (1|reps) - 1, 
                    data = soy_yd_data)

cat("\n--- Model Summary ---\n")
summary(model_yield)

cat("\n--- ANOVA ---\n")
print(anova(model_yield))

# Estimated marginal means and Tukey HSD
res_yield <- emmeans(model_yield, ~ trt)
res_yield_cld <- cld(res_yield, Letters = letters, reversed = TRUE, alpha = 0.05) %>%
  as.data.frame() %>%
  mutate(letters = trimws(.group))

cat("\n--- Tukey HSD Results ---\n")
print(res_yield_cld)

# Merge Tukey letters with data for plotting
soy_yd_plot <- soy_yd_data %>%
  left_join(res_yield_cld, by = "trt")

# Create Figure 9
fig9_yield <- ggplot(soy_yd_plot, aes(x = factor(trt), y = soy.yd.kg.ha)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgreen", alpha = 0.7) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.4, size = 1.5) +
  geom_text(data = res_yield_cld, 
            aes(x = trt, y = max(soy_yd_data$soy.yd.kg.ha, na.rm = TRUE) * 1.05, 
                label = letters),
            vjust = -0.5, size = 5, fontface = "bold") +
  
  labs(
    title = "Soybean Yield by Treatment",
    subtitle = "Letters indicate Tukey HSD groupings (α = 0.05)",
    x = "Treatment",
    y = expression("Soybean yield (kg ha"^-1*")")
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.6)
  ) +
  
  ylim(0, max(soy_yd_data$soy.yd.kg.ha, na.rm = TRUE) * 1.15)

print(fig9_yield)

# Alternative: Yield in bu/A for US audiences
fig9_yield_bu <- ggplot(soy_yd_plot, aes(x = factor(trt), y = `soy.yield.bu/a`)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgreen", alpha = 0.7) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.4, size = 1.5) +
  geom_text(data = res_yield_cld, 
            aes(x = trt, y = max(soy_yd_data$`soy.yield.bu/a`, na.rm = TRUE) * 1.05, 
                label = letters),
            vjust = -0.5, size = 5, fontface = "bold") +
  
  labs(
    title = "Soybean Yield by Treatment",
    subtitle = "Letters indicate Tukey HSD groupings (α = 0.05)",
    x = "Treatment",
    y = "Soybean yield (bu/A)"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(linewidth = 0.8)
  ) +
  
  ylim(0, max(soy_yd_data$`soy.yield.bu/a`, na.rm = TRUE) * 1.15)

# ------------------------------------------------------------------------------
# 9. CROP INJURY GAM VISUALIZATION
# ------------------------------------------------------------------------------

cat("\n========== CROP INJURY GAM ANALYSIS ==========\n")

# Filter for metribuzin treatments
sp_mtz <- soy_data %>% 
  filter(trt %in% c(3:15))

sp_mtz$dose <- sp_mtz$rate * 453.6 * 2.47

# Create crop injury GAM plot
fig_crop_injury <- sp_mtz %>%
  ggplot(aes(x = dose)) +
  # 14 DAA injury curve (typically highest)
  geom_smooth(aes(y = ci.14dat, color = "14 DAA", fill = "14 DAA"), 
              method = "gam", formula = y ~ s(x, bs = "cs"),
              linewidth = 1, se = TRUE, alpha = 0.15) +
  # 28 DAA injury curve
  geom_smooth(aes(y = ci.28dat, color = "28 DAA", fill = "28 DAA"), 
              method = "gam", formula = y ~ s(x, bs = "cs"),
              linewidth = 1, se = TRUE, alpha = 0.15) +
  # 42 DAA injury curve (typically lowest - recovery)
  geom_smooth(aes(y = ci.42dat, color = "42 DAA", fill = "42 DAA"), 
              method = "gam", formula = y ~ s(x, bs = "cs"),
              linewidth = 1, se = TRUE, alpha = 0.15) +
  
  # 5% injury threshold line (acceptable injury level)
  geom_hline(yintercept = 5, alpha = 1, color = "black", linetype = "dotted", 
             linewidth = 0.8) +
  annotate("text", x = 280, y = 8, label = "5% acceptable injury threshold", 
           size = 3.5, hjust = 0, fontface = "italic") +
  
  # Labels and theme
  labs(
    title = "Soybean Crop Injury Response to Metribuzin Dose",
    subtitle = "Crop injury remained ≤5% even at highest dose (841 g ai/ha)",
    x = expression("Metribuzin dose (g ai ha"^-1*")"),
    y = "Crop injury (% of non-treated)",
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
  
  scale_color_manual(values = c("14 DAA" = "#E41A1C", 
                                 "28 DAA" = "#377EB8", 
                                 "42 DAA" = "#4DAF4A")) +
  scale_fill_manual(values = c("14 DAA" = "#E41A1C", 
                                "28 DAA" = "#377EB8", 
                                "42 DAA" = "#4DAF4A")) +
  
  scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 100)) +
  scale_x_continuous(
    breaks = c(210, 263, 315, 368, 420, 473, 525, 578, 630, 683, 736, 788, 841)
  )

print(fig_crop_injury)

# ------------------------------------------------------------------------------
# 10. COMPARISON: SULFENTRAZONE AND S-METOLACHLOR CROP INJURY
# ------------------------------------------------------------------------------

cat("\n========== HERBICIDE COMPARISON: CROP INJURY ==========\n")

# Compare crop injury among herbicide treatments
sp_herb <- soy_data %>% 
  filter(trt %in% c("1", "16", "17"))

cat("Observations for comparison:", nrow(sp_herb), "\n")

# 14 DAA comparison
model_ci14 <- lm(ci.14dat ~ trt, data = sp_herb)
cat("\n--- Crop Injury 14 DAA: Herbicide Comparison ---\n")
summary(model_ci14)

res_ci14 <- emmeans(model_ci14, ~ trt)
res_ci14_cld <- cld(res_ci14, Letters = letters, reversed = TRUE) %>%
  as.data.frame() %>%
  mutate(letters = trimws(.group))
print(res_ci14_cld)

# 28 DAA comparison
model_ci28 <- lm(ci.28dat ~ trt, data = sp_herb)
cat("\n--- Crop Injury 28 DAA: Herbicide Comparison ---\n")
res_ci28 <- emmeans(model_ci28, ~ trt)
res_ci28_cld <- cld(res_ci28, Letters = letters, reversed = TRUE) %>%
  as.data.frame() %>%
  mutate(letters = trimws(.group))
print(res_ci28_cld)

# 42 DAA comparison
model_ci42 <- lm(ci.42dat ~ trt, data = sp_herb)
cat("\n--- Crop Injury 42 DAA: Herbicide Comparison ---\n")
res_ci42 <- emmeans(model_ci42, ~ trt)
res_ci42_cld <- cld(res_ci42, Letters = letters, reversed = TRUE) %>%
  as.data.frame() %>%
  mutate(letters = trimws(.group))
print(res_ci42_cld)

# ------------------------------------------------------------------------------
# 11. SAVE OUTPUTS
# ------------------------------------------------------------------------------

# Create output directories
dir.create("figures", showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)

# Save Figure 8
ggsave("figures/Fig8_Soybean_Height_Boxplot.png", plot = fig8_height,
       width = 12, height = 7, dpi = 300)
ggsave("figures/Fig8_Soybean_Height_Boxplot.pdf", plot = fig8_height,
       width = 12, height = 7)

# Save Figure 9
ggsave("figures/Fig9_Soybean_Yield_Boxplot.png", plot = fig9_yield,
       width = 12, height = 7, dpi = 300)
ggsave("figures/Fig9_Soybean_Yield_Boxplot.pdf", plot = fig9_yield,
       width = 12, height = 7)

# Save Crop Injury Figure
ggsave("figures/Crop_Injury_GAM.png", plot = fig_crop_injury,
       width = 10, height = 7, dpi = 300)
ggsave("figures/Crop_Injury_GAM.pdf", plot = fig_crop_injury,
       width = 10, height = 7)

# Save Tukey results
write.csv(res_height_cld, "output/Soybean_Height_Tukey_Results.csv", row.names = FALSE)
write.csv(res_yield_cld, "output/Soybean_Yield_Tukey_Results.csv", row.names = FALSE)

cat("\n========== FILES SAVED ==========\n")
cat("Figures saved to: figures/\n")
cat("Tables saved to: output/\n")

# ------------------------------------------------------------------------------
# 12. TREATMENT KEY
# ------------------------------------------------------------------------------

cat("\n========== TREATMENT KEY ==========\n")
cat("
Treatment Descriptions:
-----------------------
trt 1  = Non-treated control
trt 2  = Weed-free control
trt 3  = Metribuzin 210 g ai/ha (0.1875 lb ai/A)
trt 4  = Metribuzin 263 g ai/ha (0.234 lb ai/A)
trt 5  = Metribuzin 315 g ai/ha (0.281 lb ai/A)
trt 6  = Metribuzin 368 g ai/ha (0.328 lb ai/A)
trt 7  = Metribuzin 420 g ai/ha (0.375 lb ai/A)
trt 8  = Metribuzin 473 g ai/ha (0.422 lb ai/A)
trt 9  = Metribuzin 525 g ai/ha (0.469 lb ai/A)
trt 10 = Metribuzin 578 g ai/ha (0.516 lb ai/A)
trt 11 = Metribuzin 630 g ai/ha (0.563 lb ai/A)
trt 12 = Metribuzin 683 g ai/ha (0.609 lb ai/A)
trt 13 = Metribuzin 736 g ai/ha (0.656 lb ai/A)
trt 14 = Metribuzin 788 g ai/ha (0.703 lb ai/A)
trt 15 = Metribuzin 841 g ai/ha (0.75 lb ai/A)
trt 16 = Sulfentrazone 420 g ai/ha
trt 17 = S-metolachlor 1790 g ai/ha
")

# ------------------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------------------

cat("\n========== SESSION INFO ==========\n")
sessionInfo()

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
