# Set data path
data_name <- "tuberculosis"

# Set real-world configuration for this data set
TMT <- TRUE
NA_thr <- 0.8
DE_method <- "limma"
logFC <- TRUE
logFC_up <- 1
logFC_down <- 1
p_adj <- TRUE
alpha <- 0.05
comparisons <- NULL


# Check if DE and normalization objects already saved
if(!file.exists(paste0("data/de_results_real_world/", data_name, "_de_res.rds"))){
  # Perform preprocessing, normalization, DE analysis
  source("real-world/Real-World_Pipeline.R")
}

# Read DE results and normalized object
de_res <- readRDS(paste0("data/de_results_real_world/", data_name, "_de_res.rds"))
se_norm <- readRDS(paste0("data/preprocessed_normalized_real_world_se/", data_name, "_se_norm.rds"))

# Overview Plots

# 1. Unnormalized data

log2_pca <- plot_PCA(se_norm, ain = "log2", color_by = "Batch", label_by = "No", shape_by = "Condition") + theme
log2_boxplot <- plot_boxplots(se_norm, ain = "log2", color_by = "Batch") + theme

log2_pca / log2_boxplot

# Quantitative and Qualitative Evaluation

plot_intragroup_PMAD(se_norm) + theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "None")

# Single PCAs

assays <- names(assays(se_norm))

assays <- assays[! assays %in% c("MAD", "EigenMS", "IRS_on_MAD", "IRS_on_EigenMS", "limBE_on_MAD", "limbE_on_EigenMS")]

norm_methods_IRS <- assays[grepl("IRS_on_", assays)]
norm_methods_limBE <- assays[grepl("limBE_on_", assays)]
norm_methods_single <- assays[!assays %in% c("limBE", "IRS", "log2", "raw",norm_methods_IRS, norm_methods_limBE)]

plot_PCA(se_norm, ain = norm_methods_single, color_by = "Batch", label_by = "No", shape_by = "Condition", ncol = 4) + theme
plot_PCA(se_norm, ain = norm_methods_single, color_by = "Condition", label_by = "No", shape_by = "Batch", ncol = 4) + theme

plot_PCA(se_norm, ain = norm_methods_IRS, color_by = "Condition", label_by = "No", shape_by = "Batch", ncol = 4) + theme

plot_PCA(se_norm, ain = c("Median", "IRS_on_Median", "RobNorm", "IRS_on_RobNorm", "NormicsVSN", "IRS_on_NormicsVSN"), color_by = "Condition", label_by = "No", shape_by = "Batch", ncol = 4) + theme

plot_boxplots(se_norm, ain = norm_methods_single)
ggsave("/Users/lisiarend/Desktop/boxes.png", width = 12, height = 20)


# Complete PCA plot of all normalization methods

# find assays with _on_IRS
se_IRS <- PRONE::generate_complete_SE(se_norm, ain = norm_methods_IRS)
colData(se_IRS)$Normalization <- sapply(strsplit(colData(se_IRS)$Normalization, "_on_"), function(x) x[2])

# find assays with _on_limBE
se_limBE <- PRONE::generate_complete_SE(se_norm, ain = norm_methods_limBE)
colData(se_limBE)$Normalization <- sapply(strsplit(colData(se_limBE)$Normalization, "_on_"), function(x) x[2])

# single normalization methods
se_single <- PRONE::generate_complete_SE(se_norm, ain = norm_methods_single)

# PCA plots
pca_single <- plot_PCA(se_single, color_by = "Normalization", ain = c("all"), label_by = "No", shape_by = "Condition") + 
  theme + theme(strip.text = element_blank()) +
  ggtitle("Single Normalization Methods") + scale_color_manual(values = col_vector_norm)
pca_IRS <- plot_PCA(se_IRS, color_by = "Normalization", ain = c("all"), label_by = "No", shape_by = "Condition") + 
  theme + theme(strip.text = element_blank()) + 
  ggtitle("IRS Normalization") + scale_color_manual(values = col_vector_norm)
pca_limBE <- plot_PCA(se_limBE, color_by = "Normalization", ain = c("all"), label_by = "No", shape_by = "Condition") + 
  theme + theme(strip.text = element_blank()) +
  ggtitle("limBE Normalization") + scale_color_manual(values = col_vector_norm)

pca_single + pca_IRS + pca_limBE + patchwork::plot_layout(guides = "collect")


# DEP Comparison
de_res_IRS <- de_res[de_res$Assay %in% c(norm_methods_IRS),]
plot_overview_DE_bar(de_res_IRS, plot_type = "facet_regulation")

plot_jaccard_heatmap(de_res_IRS, plot_type = "all")

plot_upset_DE(de_res_IRS, plot_type = "stacked")
