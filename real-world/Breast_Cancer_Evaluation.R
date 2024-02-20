# Set data path
data_file <- "data/raw_real_world_se/breast_cancer_se.rds"

# Set real-world configuration for this data set
TMT <- TRUE
NA_thr <- 0.8
DE_method <- "limma"
logFC <- TRUE
logFC_up <- 1
logFC_down <- -1
p_adj <- TRUE
alpha <- 0.05
comparisons <- NULL

# Check if DE and normalization objects already saved
if(file.exists("data/de_results_real_world/breast_de_res.rds")){
  # Read DE results and normalized object
  de_res <- readRDS("data/de_results_real_world/breast_de_res.rds")
  se_norm <- readRDS("data/preprocessed_normalized_real_world_se/breast_se_norm.rds")
} else {
  # Perform preprocessing, normalization, DE analysis
  source("real-world/Real-World_Pipeline.R")
  
}

# Quantitative and Qualitative Evaluation

# Read DE results
de_res <- readRDS("data/de_results_real_world/breast_de_res.rds")
se_norm <- readRDS("data/preprocessed_normalized_real_world_se/breast_se_norm.rds")

# DEP Comparison

plot_overview_DE_bar(de_res)
plot_overview_DE_tile(de_res)

# Intersection Analysis

plot_upset_DE(de_res)

