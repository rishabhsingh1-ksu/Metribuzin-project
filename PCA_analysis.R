# ==============================================================================
# Principal Component Analysis (PCA) for Metribuzin Field Trial Data
# ==============================================================================
#
# DESCRIPTION:
# This script performs Principal Component Analysis (PCA) and hierarchical
# clustering on soil and environmental data from multi-site field experiments
# evaluating metribuzin efficacy for herbicide-resistant Amaranthus weed control
# in soybean across the United States.
#
# ASSOCIATED PUBLICATION:
# Singh et al. (2025). Optimizing metribuzin rates for herbicide-resistant 
# Amaranthus weed control in soybean. Weed Technology, 39, e97.
# DOI: https://doi.org/10.1017/wet.2025.10047
#
# AUTHOR: Rishabh Singh
# CONTACT: rs81@illinois.edu
# DATE CREATED: 2022
# LAST MODIFIED: 2025
#
# DATA REQUIREMENTS:
# - Input file: PCA.xlsx (Sheet: "Final")
# - Variables: Location, soil properties, precipitation data, soil temperature
#
# OUTPUT:
# - Hierarchical clustering dendrogram
# - PCA scree plot
# - PCA biplot with location labels
# - Loading scores for principal components
#
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. SETUP AND PACKAGE LOADING
# ------------------------------------------------------------------------------

# Clear workspace
rm(list = ls())

# Load required packages
# Install packages if not already installed using:
# install.packages(c("readxl", "vegan", "ggplot2", "factoextra", "ggrepel", "dplyr"))

library(readxl)       # Read Excel files
library(vegan)        # Ecological data analysis (optional for this script)
library(ggplot2)      # Data visualization
library(factoextra)   # PCA visualization tools
library(ggrepel)      # Non-overlapping text labels
library(dplyr)        # Data manipulation

# ------------------------------------------------------------------------------
# 2. DATA IMPORT AND PREPROCESSING
# ------------------------------------------------------------------------------

# Set working directory (modify path as needed for your system)
# setwd("path/to/your/data")

# Import data from Excel file
# Sheet "Final" contains the cleaned dataset for PCA analysis
data <- read_excel("PCA.xlsx", sheet = "Final")

# Inspect data structure
str(data)
names(data)

# Extract location labels for plotting
# The "Abb" column contains abbreviated location names (e.g., "AR'23", "KS'22")
labels <- data$Abb

# ------------------------------------------------------------------------------
# 3. DATA TYPE CONVERSION
# ------------------------------------------------------------------------------

# Convert character columns to numeric for analysis
# These variables represent soil and environmental parameters collected at each site

# Soil properties
data$`Soil OM` <- as.numeric(data$`Soil OM`)       # Soil organic matter (%)
data$`Soil pH` <- as.numeric(data$`Soil pH`)       # Soil pH
data$`Clay ` <- as.numeric(data$`Clay `)           # Clay content (%)
data$Silt <- as.numeric(data$Silt)                 # Silt content (%)
data$Sand <- as.numeric(data$Sand)                 # Sand content (%)

# Precipitation and environmental variables
data$`Interval of 1st precipitation` <- as.numeric(data$`Interval of 1st precipitation`)
data$`Interval for cumulative 12.7 mm precipitation ` <- as.numeric(data$`Interval for cumulative 12.7 mm precipitation `)
data$`Cumulative Precipitation 42 DAA` <- as.numeric(data$`Cumulative Precipitation 42 DAA`)
data$`Soil temperature ` <- as.numeric(data$`Soil temperature `)
data$`Soil moisture ` <- as.numeric(data$`Soil moisture `)

# ------------------------------------------------------------------------------
# 4. PREPARE DATA MATRIX FOR PCA
# ------------------------------------------------------------------------------

# Select numeric columns for PCA (columns 3 to 13 contain the analysis variables)
# Variables included:
# - Soil OM, Soil pH, Clay, Sand, Silt
# - Interval of 1st precipitation
# - Interval for cumulative 12.7 mm precipitation
# - Cumulative precipitation 42 DAA
# - Soil temperature at application
# - Soil moisture at application

pca_data <- data[, 3:13]

# Check for missing values
cat("\nMissing values per column:\n")
print(colSums(is.na(pca_data)))

# Remove rows with missing values if necessary (optional)
# pca_data <- na.omit(pca_data)

# ------------------------------------------------------------------------------
# 5. HIERARCHICAL CLUSTERING
# ------------------------------------------------------------------------------

# Transpose data matrix for clustering by variables
loc_data_t <- t(pca_data)

# ------------------------------------------------------------------------------
# 6. PRINCIPAL COMPONENT ANALYSIS
# ------------------------------------------------------------------------------

# Perform PCA using prcomp()
# scale = TRUE: Standardizes variables to have unit variance (recommended when
#               variables are measured in different units)
# center = TRUE: Centers variables to have zero mean
pca_result <- prcomp(pca_data, scale = TRUE, center = TRUE)

# View PCA summary statistics
cat("\n========== PCA SUMMARY ==========\n")
summary(pca_result)

# Extract key PCA outputs
cat("\n--- Standard Deviations of Principal Components ---\n")
print(pca_result$sdev)

cat("\n--- Loadings (Eigenvectors) ---\n")
print(pca_result$rotation)

cat("\n--- Variable Means (Centering Values) ---\n")
print(pca_result$center)

cat("\n--- Variable Standard Deviations (Scaling Values) ---\n")
print(pca_result$scale)

cat("\n--- Principal Component Scores (First 6 rows) ---\n")
print(head(pca_result$x))

# ------------------------------------------------------------------------------
# 7. VARIANCE EXPLAINED
# ------------------------------------------------------------------------------

# Calculate variance explained by each principal component
pca_var <- pca_result$sdev^2
pca_var_percent <- round(pca_var / sum(pca_var) * 100, digits = 2)

cat("\n--- Variance Explained by Each PC (%) ---\n")
print(pca_var_percent)

cat("\n--- Cumulative Variance Explained (%) ---\n")
print(cumsum(pca_var_percent))

# ------------------------------------------------------------------------------
# 8. SCREE PLOT
# ------------------------------------------------------------------------------

# Create scree plot showing variance explained by each component
# The "elbow" in the plot helps determine the number of components to retain

fviz_eig(pca_result,
         addlabels = TRUE,
         ylim = c(0, 50),
         main = "Scree Plot - Variance Explained by Principal Components",
         xlab = "Principal Component",
         ylab = "Percentage of Variance Explained")

# Alternative base R scree plot (barplot style)
barplot(pca_var_percent,
        main = "Scree Plot",
        xlab = "Principal Component",
        ylab = "Percent Variation (%)",
        names.arg = paste0("PC", 1:length(pca_var_percent)),
        col = "steelblue")

# ------------------------------------------------------------------------------
# 9. PCA BIPLOT VISUALIZATION
# ------------------------------------------------------------------------------

# Create publication-quality biplot
# Biplot displays both sample scores (locations) and variable loadings (arrows)

# Extract PC scores for plotting
label_data <- data.frame(
  label = labels,
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2]
)

# Create biplot using factoextra
pca_plot <- fviz_pca_biplot(
  pca_result,
  geom = "point",                    # Show points only (labels added separately)
  label = "var",                     # Label variables (arrows)
  col.ind = "skyblue",               # Color for individual points
  col.var = "maroon",                # Color for variable arrows
  repel = TRUE,                      # Avoid label overlapping
  title = "PCA Biplot - Site-Year Environmental Variation"
)

# Add location labels with non-overlapping text
pca_final <- pca_plot +
  geom_text_repel(
    data = label_data,
    aes(x = PC1, y = PC2, label = label),
    size = 4,
    box.padding = 0.15,
    point.padding = 0.1,
    max.overlaps = Inf,
    color = "black"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

# Display the plot
print(pca_final)

# ------------------------------------------------------------------------------
# 10. LOADING SCORES ANALYSIS
# ------------------------------------------------------------------------------

# Identify variables contributing most to PC1
# Loading scores indicate the correlation between original variables and PCs

loading_scores_pc1 <- pca_result$rotation[, 1]
loading_magnitude <- abs(loading_scores_pc1)
loading_ranked <- sort(loading_magnitude, decreasing = TRUE)

cat("\n--- Variables Ranked by Contribution to PC1 ---\n")
print(loading_ranked)

cat("\n--- Loading Scores with Direction (+ or -) ---\n")
print(pca_result$rotation[names(loading_ranked), 1])

# Correlation between original variables and PC scores
cat("\n--- Correlation Between Variables and PC Scores ---\n")
print(round(cor(pca_data, pca_result$x), 3))

# ------------------------------------------------------------------------------
# 11. SAVE OUTPUTS (OPTIONAL)
# ------------------------------------------------------------------------------

# Save biplot as high-resolution image for publication
# ggsave("PCA_biplot.png", plot = pca_final, width = 10, height = 8, dpi = 300)
# ggsave("PCA_biplot.pdf", plot = pca_final, width = 10, height = 8)

# Save PCA results to CSV
# write.csv(pca_result$x, "PCA_scores.csv", row.names = TRUE)
# write.csv(pca_result$rotation, "PCA_loadings.csv", row.names = TRUE)

# ------------------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------------------

cat("\n========== SESSION INFO ==========\n")
sessionInfo()

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
