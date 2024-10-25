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
norm_methods_single <- assays[!assays %in% c("limBE", "IRS", norm_methods_IRS, norm_methods_limBE)]

# Quantitative and Qualitative Evaluation

PMAD_IRS <- get_PMAD_dt(se_norm, ain = norm_methods_IRS, condition = "Timepoint")
PMAD_IRS$BEC <- "Yes"
PMAD_IRS$Normalization <- sapply(strsplit(as.character(PMAD_IRS$Normalization), "_on_"), function(x) x[2])

PMAD_single <- get_PMAD_dt(se_norm, ain = norm_methods_single, condition = "Timepoint")
PMAD_single$BEC <- "No"

dt <- rbind(PMAD_IRS, PMAD_single)

p2_li <- ggplot(dt, aes(x = BEC, y = PMAD)) + geom_violin(aes(fill = BEC), alpha = 0.4) +
  geom_boxplot(aes(fill = BEC), alpha = 0.7, width = 0.3, outliers = FALSE) +
  geom_jitter(aes(color = Normalization, shape = get(metadata(se_norm)$condition)), size = 3) + 
  theme + scale_color_manual(values = col_vector_norm, name = "Normalization Method") + 
  scale_fill_manual(values = c("#F8766D", "#00BFC4"), name = "Batch Effect Corrected (IRS)") +
  guides(shape = guide_legend(title = metadata(se_norm)$condition)) +
  xlab("Batch Effect Corrected (IRS)")

# PCA of selected normalization methods

#plot_PCA(se_norm, ain = norm_methods_single, color_by = "Timepoint", shape_by = "Batch")
#plot_PCA(se_norm, ain = norm_methods_IRS, color_by = "Timepoint", shape_by = "Batch", label_by = "No")
#plot_PCA(se_norm, ain = norm_methods_limBE, color_by = "Timepoint", shape_by = "Batch", label_by = "No")


# Complete PCA plot of all normalization methods

colData(se_norm)$Timepoint <- factor(colData(se_norm)$Timepoint, levels = c("D0", "D3", "D7", "D14"))

# try other PCA plot

generate_individual_PCA <- function(se_norm, methods_to_exclude){
  # exclude assays
  if(!is.null(methods_to_exclude)){
    norm_methods_IRS <- norm_methods_IRS[! norm_methods_IRS %in% paste0("IRS_on_",methods_to_exclude)]
    norm_methods_single <- norm_methods_single[! norm_methods_single %in% methods_to_exclude]
  }
  
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 
                                                   "qual", ]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, 
                              rownames(qual_col_pals)))
  col_vector <- rev(col_vector)

  single_PCA_dt <- PRONE:::get_complete_pca_dt(se_norm, norm_methods_single)
  coldata <- as.data.table(colData(se_norm))
  single_PCA_dt <- merge(single_PCA_dt, coldata, by = "Column")
  p1 <- ggplot(single_PCA_dt, aes( x= PC1, y = PC2, color = Assay, shape = Batch)) + geom_point(size = 3) + theme + scale_color_manual(name = "Normalization Method", values = col_vector_norm)
  p1_timepoint <- ggplot(single_PCA_dt, aes( x= PC1, y = PC2, color = Timepoint, shape = Batch)) + 
    geom_point(size = 3) + theme + guides(color = guide_legend(order = 2), shape = guide_legend(order = 1)) +
    scale_color_manual(name = "Timepoint", 
                       values = col_vector)
  
  IRS_PCA_dt <- PRONE:::get_complete_pca_dt(se_norm, norm_methods_IRS)
  IRS_PCA_dt <- merge(IRS_PCA_dt, coldata, by = "Column")
  IRS_PCA_dt$Assay <- sapply(strsplit(IRS_PCA_dt$Assay, "_on_"), function(x) x[2]) 
  p2 <- ggplot(IRS_PCA_dt, aes( x= PC1, y = PC2, color = Assay, shape = Batch)) + geom_point(size = 3) + theme + scale_color_manual(name = "Normalization Method", values = col_vector_norm)
  p2_timepoint <- ggplot(IRS_PCA_dt, aes( x= PC1, y = PC2, color = Timepoint, shape = Batch)) + 
    geom_point(size = 3) + theme + guides(color = guide_legend(order = 2), shape = guide_legend(order = 1)) +
    scale_color_manual(name = "Timepoint", 
                                values = col_vector)
  
  # Combine plots 1 and 2 with a shared legend
  pca_idv <- ((p1 + p2) / (p1_timepoint + p2_timepoint)) + 
    plot_layout(guides = "collect") & 
    theme(legend.position = "right")
  
  return(pca_idv)
}

pca_idv <- generate_individual_PCA(se_norm, NULL)

# find assays with _on_IRS
se_IRS <- PRONE::generate_complete_SE(se_norm, ain = norm_methods_IRS)
colData(se_IRS)$Normalization <- sapply(strsplit(colData(se_IRS)$Normalization, "_on_"), function(x) x[2])
pca_all_IRS <- plot_PCA(se_IRS, color_by = "Normalization", ain = c("all"), label_by = "No", shape_by = "Batch") + 
  theme + theme(strip.text = element_blank()) + scale_color_manual(name = "Normalization Method", values = col_vector_norm)

methods_to_exclude <- c("IRS_on_VSN", "IRS_on_MAD")
se_IRS <- PRONE::generate_complete_SE(se_norm, ain = norm_methods_IRS[!norm_methods_IRS %in% methods_to_exclude])
colData(se_IRS)$Normalization <- sapply(strsplit(colData(se_IRS)$Normalization, "_on_"), function(x) x[2])
pca_not_all_IRS <- plot_PCA(se_IRS, color_by = "Normalization", ain = c("all"), label_by = "No", shape_by = "Batch") + 
  theme + theme(strip.text = element_blank()) + scale_color_manual(name = "Normalization Method", values = col_vector_norm) +
  guides(color = guide_legend(order = 2), shape = guide_legend(order = 1, title = "Batch"))

(pca_all_IRS + pca_not_all_IRS) + plot_layout(guides = "auto")




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

DE_plot_D0_D14 <- plot_overview_DE_bar(new_de_res_IRS, plot_type = "single", comparisons = c("D0-D14"))[[1]] +  ggtitle("") + theme + theme(strip.background =element_rect(fill="white")) + ylab("Normalization Method") + theme


DE_plot <- plot_overview_DE_bar(new_de_res_IRS, plot_type = "facet_regulation") +  ggtitle("") + theme + theme(strip.background =element_rect(fill="white")) + ylab("Normalization Method")

plot_volcano_DE(new_de_res_IRS, comparisons = c("D0-D14")) + theme


new_de_res_IRS <- new_de_res_IRS[!new_de_res_IRS$Assay %in% c("TMM", "EigenMS"),]
plot_overview_DE_bar(new_de_res_IRS, plot_type = "facet_regulation") +  ggtitle("") + theme

jaccard_plot <- plot_jaccard_heatmap(new_de_res_IRS, plot_type = "all") + theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


jaccard_plot_D0_D14 <- plot_jaccard_heatmap(new_de_res_IRS, plot_type = "all", comparisons = c("D0-D14")) + 
  theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + ggplot2::scale_fill_gradient(low = "white", 
                                                                                                    high = "#0072B2")


((p1 + p2) / (p1_timepoint + p2_timepoint)) + 
  plot_layout(axis_titles = "collect", guides = "collect") + 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))
  
  
  (DE_plot_D0_D14 / jaccard_plot_D0_D14) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))
  
pca_idv / DE_plot_D0_D14 / jaccard_plot_D0_D14 + plot_layout(heights = c(2,1,1))


# p <- ggplot2::ggplot(melted_m, ggplot2::aes(x = Var1, 
#                                             y = Var2, fill = value)) + ggplot2::geom_tile(color = "black") + 
#   ggplot2::geom_text(ggplot2::aes(label = round(value, 
#                                                 digits = 2), color = value > 0.75)) + ggplot2::scale_fill_gradient(low = "white", 
#                                                                                                                    high = "#0072B2") + ggplot2::labs(x = "Normalization Method", 
#                                                                                                                                                      y = "Normalization Method", fill = "Jaccard Similarity") + theme +
#   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
#                                                      vjust = 0.5)) + scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none")


t <- ggarrange(panel_PCA, DE_plot, jaccard_plot, ncol = 1, common.legend = FALSE, legend = "right", labels = c("","D", "E"), heights = c(0.7,0.5, 0.5))
ggsave("figures/test_li_et_al.png",t, width = 14, height = 16, dpi = 300)

# Stacked barplot of overlaps
plot_upset_DE(new_de_res_IRS, plot_type = "stacked")$upset

# Consensus
DE_candidates <- PRONE::extract_consensus_DE_candidates(new_de_res_IRS, ain = NULL, norm_thr = 0.9, per_comparison = FALSE)

