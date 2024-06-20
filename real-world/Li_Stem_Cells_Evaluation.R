# Set data path
data_name <- "li_stem_cells"

# Set real-world configuration for this data set
TMT <- TRUE
NA_thr <- 0.8
DE_method <- "limma"
logFC <- TRUE
logFC_up <- log2(1.2)
logFC_down <- log2(0.83)
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

# Evaluation

# 1. Unnormalized data

log2_pca <- plot_PCA(se_norm, ain = "log2", color_by = "Batch", label_by = "No", shape_by = "Timepoint") + theme
log2_boxplot <- plot_boxplots(se_norm, ain = "log2", color_by = "Batch") + theme

log2_pca / log2_boxplot


# Categorize normalization methods
assays <- names(assays(se_norm))
#methods_to_remove <- c("VSN", "IRS_on_VSN", "limBE_on_VSN", "MAD", "IRS_on_MAD", "limBE_on_MAD", "EigenMS", "IRS_on_EigenMS", "limBE_on_EigenMS", "TMM", "IRS_on_TMM", "limBE_on_TMM")
#assays <- assays[!assays %in% methods_to_remove]
norm_methods_IRS <- assays[grepl("IRS_on_", assays)]
norm_methods_limBE <- assays[grepl("limBE_on_", assays)]
norm_methods_single <- assays[!assays %in% c("limBE", "IRS", "log2", norm_methods_IRS, norm_methods_limBE)]

# Quantitative and Qualitative Evaluation

#plot_intragroup_PMAD(se_norm, ain = norm_methods_single) + theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "None")

#plot_intragroup_PMAD(se_norm, ain = norm_methods_IRS) + theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "None")


# PCA of selected normalization methods

#plot_PCA(se_norm, ain = norm_methods_single, color_by = "Timepoint", shape_by = "Batch")
#plot_PCA(se_norm, ain = norm_methods_IRS, color_by = "Timepoint", shape_by = "Batch", label_by = "No")
#plot_PCA(se_norm, ain = norm_methods_limBE, color_by = "Timepoint", shape_by = "Batch", label_by = "No")


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
pca_single <- plot_PCA(se_single, color_by = "Normalization", ain = c("all"), label_by = "No", shape_by = "Batch") + 
  theme + theme(strip.text = element_blank()) + scale_color_manual(name = "Normalization Method", values = col_vector_norm)

pca_IRS <- plot_PCA(se_IRS, color_by = "Normalization", ain = c("all"), label_by = "No", shape_by = "Batch") + 
  theme + theme(strip.text = element_blank()) + scale_color_manual(name = "Normalization Method", values = col_vector_norm)

pca_IRS_timepoint <- plot_PCA(se_IRS, color_by = "Timepoint", ain = c("all"), label_by = "No", shape_by = "Batch") + 
  theme + theme(strip.text = element_blank()) 

ggarrange(pca_single, pca_IRS, labels = c("A", "B"), ncol = 2, common.legend = TRUE, legend = "right", font.label = list(size = 18))

# try other PCA plot

single_PCA_dt <- PRONE:::get_complete_pca_dt(se_norm, norm_methods_single)
coldata <- as.data.table(colData(se_norm))
single_PCA_dt <- merge(single_PCA_dt, coldata, by = "Column")
p1 <- ggplot(single_PCA_dt, aes( x= PC1, y = PC2, color = Assay, shape = Batch)) + geom_point(size = 3) + theme + scale_color_manual(name = "Normalization Method", values = col_vector_norm)
p1_timepoint <- ggplot(single_PCA_dt, aes( x= PC1, y = PC2, color = Timepoint, shape = Batch)) + geom_point(size = 3) + theme 

IRS_PCA_dt <- PRONE:::get_complete_pca_dt(se_norm, norm_methods_IRS)
IRS_PCA_dt <- merge(IRS_PCA_dt, coldata, by = "Column")
IRS_PCA_dt$Assay <- sapply(strsplit(IRS_PCA_dt$Assay, "_on_"), function(x) x[2]) 
p2 <- ggplot(IRS_PCA_dt, aes( x= PC1, y = PC2, color = Assay, shape = Batch)) + geom_point(size = 3) + theme + scale_color_manual(name = "Normalization Method", values = col_vector_norm)
p2_timepoint <- ggplot(IRS_PCA_dt, aes( x= PC1, y = PC2, color = Timepoint, shape = Batch)) + geom_point(size = 3) + theme 

p1p2 <- ggarrange(p1, p2, common.legend = TRUE, legend = "right", labels = c("A", "B"), ncol = 2)
p1p2_timepoint <- ggarrange(p1_timepoint, p2_timepoint, common.legend = TRUE, legend = "right", labels = c("A", "B"), ncol = 2)

ggarrange(p1p2, p2_timepoint, legend = "right", labels = c("", "C"), ncol = 2)

# Check densities, boxplots, normal PCA or removed methods
plot_boxplots(se_norm, ain = methods_to_remove)
plot_densities(se_norm, ain = methods_to_remove)
plot_PCA(se_norm, ain = methods_to_remove, color_by = "Batch", shape_by = "Timepoint", label_by = "No")


#pca_limBE <- plot_PCA(se_limBE, color_by = "Normalization", ain = c("all"), label_by = "No", shape_by = "Batch") + 
#  theme + theme(strip.text = element_blank()) +
#  ggtitle("limBE Normalization")  + scale_color_manual(values = col_vector_norm)

#pca_single + pca_IRS + pca_limBE + patchwork::plot_layout(guides = "collect")

# DEP Comparison

# Methods based on IRS-normalized data
de_res_IRS <- de_res[de_res$Assay %in% c(norm_methods_IRS),]
de_res_IRS$Assay <- sapply(strsplit(de_res_IRS$Assay, "_on_"), function(x) x[2])
plot_overview_DE_bar(de_res_IRS, plot_type = "facet_regulation") + ggtitle("") + theme

# remove EigenMS
de_res_IRS_no_EigenMS <- de_res_IRS[de_res_IRS$Assay != "EigenMS",]
plot_overview_DE_bar(de_res_IRS_no_EigenMS, plot_type = "facet_regulation")

# apply other thresholds
new_de_res_IRS <- PRONE::apply_thresholds(de_res_IRS)

plot_overview_DE_bar(new_de_res_IRS, plot_type = "facet_regulation") +  ggtitle("") + theme

plot_volcano_DE(new_de_res_IRS, comparisons = c("D0-D14")) + theme


new_de_res_IRS <- new_de_res_IRS[!new_de_res_IRS$Assay %in% c("TMM", "EigenMS"),]
plot_overview_DE_bar(new_de_res_IRS, plot_type = "facet_regulation") +  ggtitle("") + theme

plot_jaccard_heatmap(new_de_res_IRS, plot_type = "all") + theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# Stacked barplot of overlaps
plot_upset_DE(new_de_res_IRS, plot_type = "stacked")$upset

# Consensus
DE_candidates <- PRONE::extract_consensus_DE_candidates(new_de_res_IRS, ain = NULL, norm_thr = 0.9, per_comparison = FALSE)

