# DE Analysis Paper Evaluation of Real-World Data Sets

# Read data


## Li Stem Cells
data_name <- "li_stem_cells"
li_de_res <- readRDS(paste0("data/de_results_real_world/", data_name, "_de_res.rds"))
li_se_norm <- readRDS(paste0("data/preprocessed_normalized_real_world_se/", data_name, "_se_norm.rds"))
li_condition <- "Timepoint"
li_batch <- "Batch"

## Mouse P450 
data_name <- "mouse_liver_cytochrome_P450"
p450_de_res <- readRDS(paste0("data/de_results_real_world/", data_name, "_de_res.rds"))
p450_se_norm <- readRDS(paste0("data/preprocessed_normalized_real_world_se/", data_name, "_se_norm.rds"))

## CPTAC
data_name <- "CPTAC_OV"
cptac_de_res <- readRDS(paste0("data/de_results_real_world/", data_name, "_de_res.rds"))
cptac_se_norm <- readRDS(paste0("data/preprocessed_normalized_real_world_se/", data_name, "_se_norm.rds"))
cptac_condition <- "Pathological_Status"
cptac_batch <- "Pool"

# Define Norm Methods
assays <- names(assays(li_se_norm))
norm_methods_IRS <- assays[grepl("IRS_on_", assays)]
norm_methods_limBE <- assays[grepl("limBE_on_", assays)]
norm_methods_single <- assays[!assays %in% c("limBE", "IRS", norm_methods_IRS, norm_methods_limBE)]


# Intragroup Variation Plots
## -> see limBE_IRS.R file

# DE Analysis

li_de_res_new <- PRONE::apply_thresholds(li_de_res) # logFC =1, p.adj < 0.05
li_de_res_new <- li_de_res_new[li_de_res_new$Assay %in% c(norm_methods_limBE),]
li_de_res_new$Assay <- sapply(strsplit(li_de_res_new$Assay, "_on_"), function(x) x[2]) 

cptac_de_res <- cptac_de_res[cptac_de_res$Assay %in% c(norm_methods_limBE),]
cptac_de_res$Assay <- sapply(strsplit(cptac_de_res$Assay, "_on_"), function(x) x[2])

li_overview_de_plot <- plot_overview_DE_bar(li_de_res_new, plot_type = "single", comparisons = c("D0-D14"))[[1]] +  ggtitle("") + theme + theme(strip.background =element_rect(fill="white"), legend.position = "none") + ylab("Normalization Method") + ggtitle("DE Results of Cell Culture Dataset dR1")
li_jaccard_heatmap <- plot_jaccard_heatmap(li_de_res_new, plot_type = "all", comparisons = c("D0-D14")) + theme + theme(axis.text.x = element_text(angle = 90, vjust =0.5), legend.position = "none")

cptac_overview_de_plot <- plot_overview_DE_bar(cptac_de_res, plot_type = "single")[[1]]  +  ggtitle("") + theme + theme(strip.background =element_rect(fill="white"), legend.position = c(0.85,0.7), legend.box.background = element_rect(colour = "black")) + ylab("Normalization Method") + ggtitle("DE Results of Clinical Cancer Dataset dR2")
cptac_jaccard_heatmap <- plot_jaccard_heatmap(cptac_de_res, plot_type = "all") + theme + theme(axis.text.x = element_text(angle = 90, vjust =0.5), legend.position = c(0.85,0.8), legend.box.background = element_rect(colour = "black")) 

cptac_jaccard_heatmap <- cptac_jaccard_heatmap + theme(axis.text.x = element_text(hjust = 1)) + scale_x_discrete(limits = rev) + ggtitle("Intersection Analysis of DE Results of Clincal Cancer Dataset dR2")
li_jaccard_heatmap <- li_jaccard_heatmap + theme(axis.text.x = element_text(hjust = 1)) + scale_x_discrete(limits = rev) + ggtitle("Intersection Analysis of DE Results of Cell Culture Dataset dR1")


li_jaccard_heatmap$layers[[2]] <- NULL
li_jaccard_heatmap <- li_jaccard_heatmap + geom_text(aes(label = round(value, digits = 1), color = value > 0.5))
cptac_jaccard_heatmap$layers[[2]] <- NULL
cptac_jaccard_heatmap <- cptac_jaccard_heatmap + geom_text(aes(label = round(value, digits = 1), color = value > 0.5))

(li_overview_de_plot + cptac_overview_de_plot + plot_layout(axis_titles = "collect", guides = "keep", axes = "collect"))  / 
  (li_jaccard_heatmap + cptac_jaccard_heatmap + plot_layout(axis_titles = "collect", guides = "keep", axes = "collect")) + 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 23))

(li_overview_de_plot + cptac_overview_de_plot + plot_layout(guides = "keep", axes = "collect", axis_titles = "keep"))  / 
  (li_jaccard_heatmap + cptac_jaccard_heatmap + plot_layout(guides = "keep", axes = "collect", axis_titles = "keep")) + 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 23))

ggsave("figures/de_results_real_world_paper.png", width = 14, height = 12, dpi = 300)

# Poster

#li_overview_de_plot <- li_overview_de_plot + scale_fill_manual(name = "Regulation", values = c("#D55E00", "#0072B2"), labels = c("Up", "Down")) + ggtitle(NULL)
#cptac_overview_de_plot <- cptac_overview_de_plot + scale_fill_manual(name = "Regulation", values = c("#D55E00", "#0072B2"), labels = c("Up", "Down")) + ggtitle(NULL)

#li_jaccard_heatmap$layers[[2]] <- NULL
#li_jaccard_heatmap <- li_jaccard_heatmap + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + labs(fill = "Jaccard") + scale_x_discrete(limits = rev) + labs(x = "Normalization Method *")
#cptac_jaccard_heatmap$layers[[2]] <- NULL
#cptac_jaccard_heatmap <- cptac_jaccard_heatmap + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + labs(fill = "Jaccard") + scale_x_discrete(limits = rev) + labs( x = "Normalization Method *")

#(li_overview_de_plot + cptac_overview_de_plot + plot_layout(axis_titles = "collect", guides = "collect", axes = "collect"))  / 
#  (li_jaccard_heatmap + cptac_jaccard_heatmap + plot_layout(axis_titles = "collect", guides = "collect", axes = "collect")) +
#  plot_annotation(tag_levels = list(c("E", "F", "G", "H"))) &
#  theme(plot.tag = element_text(face = "bold", size = 22, margin = unit(c(0,0.2,-0.5,0), "cm")))


(li_overview_de_plot + cptac_overview_de_plot + li_jaccard_heatmap+ cptac_jaccard_heatmap + plot_layout(widths = c(1.5,1.5,2,2), ncol = 4, axis_titles = "collect", guides = "collect", axes = "collect")) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 22, margin = margin(-10, 0, -25, 0)))

ggsave("figures/de_results_poster.png", width = 16, height = 11, dpi = 300)


# Other comparisons Li et. al dataset

plot_overview_DE_bar(li_de_res_new, plot_type = "facet_regulation") + ggtitle("") + theme + theme(strip.background =element_rect(fill="white")) + ylab("Normalization Method") 
ggsave("figures/de_results_li_other_comps.png", width = 15, height = 8, dpi = 300)

li_jaccard_heatmap_other <- plot_jaccard_heatmap(li_de_res_new, plot_type = "facet_comp") + facet_wrap(~Comparison, ncol = 2) + theme + theme(axis.text.x = element_text(angle = 90, vjust =0.5), strip.background =element_rect(fill="white")) + ylab("Normalization Method") + theme(axis.text.x = element_text(hjust = 1)) + scale_x_discrete(limits = rev)
li_jaccard_heatmap_other$layers[[2]] <- NULL
li_jaccard_heatmap_other <- li_jaccard_heatmap_other + geom_text(aes(label = round(value, digits = 1), color = value > 0.5))


ggsave("figures/de_results_jaccard_li_other_comps.png", width = 15, height = 18, dpi = 300)

# Volcano plot Li et. al dataset

plot_volcano_DE(li_de_res_new[!is.na(li_de_res_new$Change),], comparisons = c("D0-D14"), facet_norm = TRUE)[[1]] + theme + theme(strip.background =element_rect(fill="white")) + ylab("Normalization Method") + ggtitle("")
ggsave("figures/de_volcano_jaccard_li.png", width = 18, height = 15, dpi = 300)


# Mouse Dataset


## PMAD

p450_pmad <- plot_intragroup_PMAD(p450_se_norm) + scale_fill_manual(name = "Normalization Method", values = col_vector_norm) + theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + theme(legend.position = "none")

p450_pca <- plot_PCA(p450_se_norm, label_by = "No") + theme

plot_densities(p450_se_norm)

single_PCA_dt <- PRONE:::get_complete_pca_dt(p450_se_norm)
coldata <- as.data.table(colData(p450_se_norm))
single_PCA_dt <- merge(single_PCA_dt, coldata, by = "Column")
p450_pca_methods <- ggplot(single_PCA_dt, aes( x= PC1, y = PC2, color = Assay, shape = Condition)) + 
  geom_point(size = 2) + theme  +
  scale_color_manual(name = "Normalization Method", values = col_vector_norm) +
  guides(shape = guide_legend(order = 0), color = guide_legend(order = 1))

p450_pca_condition <- ggplot(single_PCA_dt, aes( x= PC1, y = PC2, color = Condition)) + 
  geom_point(size = 2) + theme 

p450_pmad / (p450_pca_methods + p450_pca_condition) + plot_layout(axis_titles = "collect", guides = "collect", axes = "collect") + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold"))
ggsave("figures/pca_pmad_mouse.png", width = 12, height = 8, dpi = 300)




## DE Results

p450_de_res <- PRONE::apply_thresholds(p450_de_res)

p450_de_res$Assay <- factor(p450_de_res$Assay, levels = sort(unique(p450_de_res$Assay)))

p450_overview_de_plot <- plot_overview_DE_bar(p450_de_res, plot_type = "single")[[1]] +  ggtitle("") + theme + theme(strip.background =element_rect(fill="white"), legend.position = c(0.85,0.8)) + ylab("Normalization Method") 
p450_jaccard_heatmap <- plot_jaccard_heatmap(p450_de_res, plot_type = "all") + theme + theme(axis.text.x = element_text(angle = 90, vjust =0.5), legend.position = c(0.85,0.8)) 

p450_jaccard_heatmap <- p450_jaccard_heatmap + theme(axis.text.x = element_text(hjust = 1)) + scale_x_discrete(limits = rev)


p450_jaccard_heatmap$layers[[2]] <- NULL
p450_jaccard_heatmap <- p450_jaccard_heatmap + geom_text(aes(label = round(value, digits = 1), color = value > 0.5))

p450_overview_de_plot + p450_jaccard_heatmap + plot_layout(axis_titles = "collect", guides = "keep", axes = "collect") + 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 20)) & theme(legend.position = "top")

ggsave("figures/de_results_p450.png", width = 16, height = 10, dpi = 300)

