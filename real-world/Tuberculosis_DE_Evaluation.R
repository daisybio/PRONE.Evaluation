
# Set data path
data_file <- "data/raw_real_world_se/tuberculosis_se.rds"

# Set real-world configuration for this data set
TMT <- TRUE
NA_thr <- 0.8
DE_method <- "limma"
logFC <- FALSE
logFC_up <- 0.5
logFC_down <- 0.5
p_adj <- TRUE
alpha <- 0.01
comparisons <- c("PTB-HC", "TBL-HC", "TBL-PTB", "Rx-PTB")

# Perform preprocessing, normalization, DE analysis
source("real-world/Real-World_Pipeline.R")

# Read results
de_res <- readRDS("data/de_results_real_world/tuberculosis_de_res.rds")
se_norm <- readRDS("data/preprocessed_normalized_real_world_se/tuberculosis_se_norm.rds")

# DEP Comparison
plot_overview_DE_bar(de_res)
plot_overview_DE_tile(de_res)

de_res_1 <- apply_thresholds(de_res, logFC = TRUE, logFC_up = 1, logFC_down = -1, p_adj = TRUE, alpha = 0.05)

plot_overview_DE_bar(de_res_1)
plot_overview_DE_tile(de_res_1)


# Intersection Analysis

plot_upset_DE(de_res)


# Enrichment

#plot_intersection_enrichment(se, de_res, ain = c("Median_on_IRS", "Median_on_limBE", "RobNorm_on_IRS", "RobNorm_on_limBE", "IRS_on_GlobalMean_on_TMM"))
