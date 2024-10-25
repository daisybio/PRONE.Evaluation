# Intragroup Variation Paper Evaluation of Real-World Data Sets

get_complete_PCA_new <- function (se, ain = NULL){
  complete_pca_dt <- NULL
  if (is.null(ain)) {
    assays <- names(SummarizedExperiment::assays(se))
  }
  else {
    assays <- c(ain)
  }
  for (ain in assays) {
    if (ain != "raw") {
      dt <- SummarizedExperiment::assays(se)[[ain]]
      dt <- dt[!apply(is.na(dt), 1, any), ]
      rm <- which(apply(t(dt), 2, stats::var) == 0)
      if (!identical(rm, integer(0))) {
        dt <- dt[-rm, ]
      }
      pca <- stats::prcomp(t(dt), scale. = TRUE)
      pc1 <- round(summary(pca)$importance["Proportion of Variance", 
                                           "PC1"] * 100, 2)
      pc2 <- round(summary(pca)$importance["Proportion of Variance", 
                                           "PC2"] * 100, 2)
      pc3 <- round(summary(pca)$importance["Proportion of Variance", 
                                           "PC3"] * 100, 2)
      pca_dt <- data.table::data.table(PC1 = pca$x[, "PC1"], 
                                       PC2 = pca$x[, "PC2"], 
                                       PC3 = pca$x[, "PC3"], Column = row.names(pca$x))
      pca_dt$Assay <- ain
      pca_dt$PC1_Perc <- pc1
      pca_dt$PC2_Perc <- pc2
      pca_dt$PC3_Perc <- pc3
      if (is.null(complete_pca_dt)) {
        complete_pca_dt <- pca_dt
      }
      else {
        complete_pca_dt <- rbind(complete_pca_dt, pca_dt)
      }
    }
  }
  return(complete_pca_dt)
}


# Silhouette coefficient
# function that calculates the silhouette coefficient based on a known cluster (i.e. batch or treatment)
# calculates silhouette width average
calc.sil = function(
    x, # the PC variates
    y1, y2 = NULL, # factor of interest, e.g. known batch info or known treatment info
    name.y1, name.y2 = NULL # character of the factor of interest
){
  library(cluster)
  # calculate the distance, here euclidean is appropriate for PCA, NOT for t-SNE
  dist.res = daisy(x, metric = 'euclidean')
  # for factor 1
  sil.batch.res1 = silhouette(x = as.numeric(y1), dist = dist.res)
  # if factor 2 is provided
  if(!is.null(y2)) sil.batch.res2 = silhouette(x = as.numeric(y2), dist = dist.res)
  
  # extract average width silhouette per level
  res1 = c(summary(sil.batch.res1)['clus.avg.widths']$clus.avg.widths)
  names(res1) = levels(y1)
  
  
  if(!is.null(y2)){
    res2 = c(summary(sil.batch.res2)['clus.avg.widths']$clus.avg.widths)
    names(res2) = levels(y2)
  }
  
  # output data for plotting
  if(!is.null(y2)){
    silh.coeff = c(res1, res2)
    Cluster = c(levels(y1), levels (y2))
    Type = c(rep(name.y1, nlevels(y1)), rep(name.y2, nlevels(y2)))
    data.plot = data.frame(silh.coeff, Cluster, Type)
    
  }else{
    silh.coeff = c(res1)
    Cluster = c(levels(y1))
    Type = rep(name.y1, nlevels(y1))
    data.plot = data.frame(silh.coeff, Cluster, Type)
  }
  
  return(invisible(data.plot))
}


calculate_silhouette_on_PCA <- function(se, condition, batch, shapes){
  
  pc_dt <- get_complete_PCA_new(se, NULL)
  
  # replace values by numbers
  # create mapping of column names to numbers
  condition_num <- as.factor(as.numeric(as.factor(se[[condition]])))
  names(condition_num) <- se[["Column"]]
  
  batch_num <- as.factor(as.numeric(as.factor(se[[batch]])))
  names(batch_num) <- se[["Column"]]
  
  # mapping
  mapping_cond <- unique(data.table(Name = se[[condition]], Value = condition_num))
  mapping_batch <- unique(data.table(Name = se[[batch]], Value = batch_num))
  
  
  all_silh <- NULL
  for(ain in unique(pc_dt$Assay)){
    print(ain)
    pc_ain <- pc_dt[pc_dt$Assay == ain, c("PC1", "PC2","PC3", "Column")]
    # order pc_ain by Column of colData
    pc_ain <- pc_ain[match(rownames(colData(se)), pc_ain$Column),]
    pc_ain$Column <- NULL
    pc_ain <- as.matrix(pc_ain)
    rownames(pc_ain) <- se[["Column"]]
    silh <- calc.sil(as.matrix(pc_ain), y1 = condition_num, y2 = batch_num, name.y1 = condition, name.y2 = batch)
    silh$Assay <- ain
    if(is.null(all_silh)){
      all_silh <- silh
    } else {
      all_silh <- rbind(all_silh, silh)
    }
  }
  
  # remap cluster and batch numbers to values
  silh_condition <- all_silh[all_silh$Type == condition,]
  silh_condition$Cluster <- mapping_cond[match(silh_condition$Cluster, mapping_cond$Value), "Name"][[1]]
  
  silh_batch <- all_silh[all_silh$Type == batch,]
  silh_batch$Cluster <- mapping_batch[match(silh_batch$Cluster, mapping_batch$Value), "Name"][[1]]
  
  all_silh <- rbind(silh_condition, silh_batch)
  
  all_silh$BEC <- "No"
  all_silh$BEC[grepl("IRS_on_", all_silh$Assay)] <- "IRS"
  all_silh$Assay[grepl("IRS_on_", all_silh$Assay)] <- sapply(strsplit(as.character(all_silh$Assay[grepl("IRS_on_", all_silh$Assay)]), "_on_"), function(x) x[2])
  
  all_silh$BEC[grepl("limBE_on_", all_silh$Assay)] <- "limBE"
  all_silh$Assay[grepl("limBE_on_", all_silh$Assay)] <- sapply(strsplit(as.character(all_silh$Assay[grepl("limBE_on_", all_silh$Assay)]), "_on_"), function(x) x[2])
  
  # change log2
  all_silh$BEC[all_silh$Assay == "IRS"] <- "IRS"
  all_silh$Assay[all_silh$Assay == "IRS"] <- "log2"
  
  all_silh$BEC[all_silh$Assay == "limBE"] <- "limBE"
  all_silh$Assay[all_silh$Assay == "limBE"] <- "log2"
  
  all_silh$BEC <- factor(all_silh$BEC, levels = c("No", "IRS", "limBE"))
  all_silh$Type <- factor(all_silh$Type, levels = c(batch, condition))
  
  return(all_silh)
}

# Helper functions

get_PMAD_dt <- function(se, ain, condition){
  ain <- PRONE:::check_input_assays(se, ain)
  if (is.null(ain)) {
    return(NULL)
  }
  assays <- SummarizedExperiment::assays(se)
  exclude <- names(assays)[!names(assays) %in% ain]
  for (e in exclude) {
    assays[e] <- NULL
  }
  assays["raw"] <- NULL
  for (a in names(assays)) {
    assays[[a]] <- as.matrix(assays[[a]])
  }
  condition <- PRONE:::get_condition_value(se, condition)
  avgmadmem <- NormalyzerDE:::calculateAvgMadMem(assays, data.table::as.data.table(SummarizedExperiment::colData(se))[[condition]])
  # Add rownames as column named by condition
  conds <- rownames(avgmadmem)
  avgmadmem <- data.table::as.data.table(avgmadmem)
  avgmadmem[[condition]] <- conds
  pmad_long <- melt(avgmadmem, id.vars = condition, variable.name = "Normalization", value.name = "PMAD")
  return(pmad_long)
}



## Li Stem Cells
data_name <- "li_stem_cells"
li_de_res <- readRDS(paste0("data/de_results_real_world/", data_name, "_de_res.rds"))
li_se_norm <- readRDS(paste0("data/preprocessed_normalized_real_world_se/", data_name, "_se_norm.rds"))
li_condition <- "Timepoint"
li_batch <- "Batch"
li_shapes <- c(0, 1, 2, 3, 15, 16, 17)

## CPTAC
data_name <- "CPTAC_OV"
cptac_de_res <- readRDS(paste0("data/de_results_real_world/", data_name, "_de_res.rds"))
cptac_se_norm <- readRDS(paste0("data/preprocessed_normalized_real_world_se/", data_name, "_se_norm.rds"))
cptac_condition <- "Pathological_Status"
cptac_batch <- "Pool"
cptac_shapes <- c(15, 16, seq(0,12))

# Define Norm Methods
assays <- names(assays(li_se_norm))
norm_methods_IRS <- assays[grepl("IRS_on_", assays)]
norm_methods_limBE <- assays[grepl("limBE_on_", assays)]
norm_methods_single <- assays[!assays %in% c("IRS", "limBE", norm_methods_IRS, norm_methods_limBE)]


## Silhouette 

# Li Silhouette
tmp_li <- calculate_silhouette_on_PCA(li_se_norm, li_condition, li_batch, li_shapes)
tmp_li$Cluster <- factor(tmp_li$Cluster, levels = c("D0", "D3", "D7", "D14", "Set1", "Set2", "Set3"))

comparisons <- list(c("IRS", "limBE"), c("limBE", "No"), c("IRS", "No"))

# Summary
plt_li_condition <- ggplot(tmp_li[tmp_li$Type == "Timepoint",], aes(x = BEC, y = silh.coeff)) +
  geom_violin(aes(fill = BEC), alpha = 0.4, width = 1) +
  geom_boxplot(aes(fill = BEC), alpha = 0.7, width = 0.1, outliers = FALSE) +
  stat_compare_means(data = tmp_li[tmp_li$Type == "Timepoint",], comparisons = comparisons, method = "wilcox.test", paired = TRUE) +
  theme +
  scale_fill_manual(values = c("#F8766D", "#00BFC4", "#619CFF"), name = "Batch Effect Correction") +
  guides(shape = guide_legend(title =  "Timepoint", order = 1, position = "top", nrow = 1), color = FALSE, fill = FALSE) +
  xlab("Batch Effect Correction") +
  scale_shape_manual(values = c(16, 17, 18, 19)) + labs(y = "Silhouette Coefficient") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

plt_li_batch <- ggplot(tmp_li[tmp_li$Type == "Batch",], aes(x = BEC, y = silh.coeff)) +
  geom_violin(aes(fill = BEC), alpha = 0.4, width = 1) +
  geom_boxplot(aes(fill = BEC), alpha = 0.7, width = 0.1, outliers = FALSE) +
  stat_compare_means(data = tmp_li[tmp_li$Type == "Batch",], comparisons = comparisons, method = "wilcox.test", paired = TRUE) +
  theme +
  scale_fill_manual(values = c("#F8766D", "#00BFC4", "#619CFF"), name = "Batch Effect Correction") +
  guides(shape = guide_legend(title =  "Batch", order = 1, position = "top", nrow = 1), color = FALSE, fill = FALSE) +
  xlab("") +
  scale_shape_manual(values = c(1,2,3)) + labs(y = "Silhouette Coefficient") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# CPTAC Silhouette
library(plyr)
tmp_cptac <- calculate_silhouette_on_PCA(cptac_se_norm, cptac_condition, cptac_batch, cptac_shapes)

nbatch <- length(unique(tmp_cptac[tmp_cptac$Type == "Pool",]$Cluster))
levels(tmp_cptac$Type) <- c(levels(tmp_cptac$Type), "Batch", "Pathological Status")
tmp_cptac$Type[tmp_cptac$Type == "Pool"] <- "Batch"
tmp_cptac$Type[tmp_cptac$Type == "Pathological_Status"] <- "Pathological Status"

tmp_cptac_condition <- tmp_cptac[tmp_cptac$Type == "Pathological Status",]
tmp_cptac_batch <- tmp_cptac[tmp_cptac$Type == "Batch",]

from <- c("01CPTAC","02CPTAC","03CPTAC","04CPTAC","05CPTAC","06CPTAC","07CPTAC","08CPTAC","09CPTAC","10CPTAC","11CPTAC","12CPTAC","13CPTAC")
to <- c("Set1", "Set2", "Set3", "Set4", "Set5", "Set6", "Set7", "Set8", "Set9", "Set10", "Set11", "Set12", "Set13")
tmp_cptac_batch$Cluster <- plyr::mapvalues(tmp_cptac_batch$Cluster, from, to)

tmp_cptac_batch$Cluster <- factor(tmp_cptac_batch$Cluster, levels = to)

comparisons <- list(c("IRS", "limBE"), c("limBE", "No"), c("IRS", "No"))

plt_cptac_condition <- ggplot(tmp_cptac_condition, aes(x = BEC, y = silh.coeff)) + 
  geom_violin(aes(fill = BEC), alpha = 0.4, width = 1) +
  geom_boxplot(aes(fill = BEC), alpha = 0.7, width = 0.1, outliers = FALSE) +
  theme_bw() + theme +
  stat_compare_means(data = tmp_cptac_condition, comparisons = comparisons, method = "wilcox.test", paired = TRUE) +
  theme + 
  scale_fill_manual(values = c("#F8766D", "#00BFC4", "#619CFF"), name = "Batch Effect Correction") +
  guides(shape = guide_legend(title =  "Pathological Status", order = 1, position = "top", nrow = 1), color = FALSE, fill = FALSE) +
  xlab("Batch Effect Correction") + scale_shape_manual(values = c(15, 16)) + labs(y = "Silhouette Coefficient")

plt_cptac_batch <- ggplot(tmp_cptac_batch, aes(x = BEC, y = silh.coeff)) + 
  geom_violin(aes(fill = BEC), alpha = 0.4, width = 1) +
  geom_boxplot(aes(fill = BEC), alpha = 0.7, width = 0.1, outliers = FALSE) +
  theme_bw() + theme +
  stat_compare_means(data = tmp_cptac_batch, comparisons = comparisons, method = "wilcox.test", paired = TRUE) +
  theme + 
  scale_fill_manual(values = c("#F8766D", "#00BFC4", "#619CFF"), name = "Batch Effect Correction") +
  guides(shape = guide_legend(title =  "Batch", order = 1, position = "top", nrow = 2), fill = FALSE) +
  xlab("Batch Effect Correction") + scale_shape_manual(values = seq(0,12)) + labs(y = "Silhouette Coefficient")

(plt_li_batch + plt_li_condition + plt_cptac_batch + plt_cptac_condition) + 
  plot_layout(guides = "collect", axis = "collect", axes = "collect", ncol = 2) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 22, margin = unit(c(0,0.2,-0.5,0), "cm")))


ggsave("figures/Silhouette_Coefficient_IRS_limBE_all.png", width = 12, height = 8)


# Supplemental Figures (per normalization method)

# Li
tmp_li$Cluster <- factor(tmp_li$Cluster, levels = c("Set1", "Set2", "Set3", "D0", "D3", "D7", "D14"))
ggplot(tmp_li, aes(x = Type, y = silh.coeff, fill = BEC), aes(x = Type, y = silh.coeff)) + 
  geom_boxplot(aes(fill = BEC), alpha = 0.7, outlier.shape = NA, position = position_dodge(width = 1)) + 
  geom_point(aes(group = BEC, shape = Cluster), position = position_jitterdodge(jitter.width = 0.8, dodge.width = 1)) +
  facet_wrap(~Assay, ncol = 4) +
  theme + scale_fill_manual(name = "Batch Effect Correction", values = c("#F8766D", "#00BFC4", "#619CFF")) + 
  labs(y = "Silhouette Coefficient", x = "Sample Group") + scale_shape_manual(values = c(1,2,3, 16, 17, 18, 19), name = "Sample Group \n(Batch or Timepoint)")

ggsave("figures/Silhouette_Coefficient_IRS_limBE_li.png", width = 12, height = 12)


# CPTAC
ggplot(rbind(tmp_cptac_batch, tmp_cptac_condition), aes(x = Type, y = silh.coeff)) + 
  geom_boxplot(aes(fill = BEC), alpha = 0.7, outlier.shape = NA, position = position_dodge(width = 1)) + 
  geom_point(aes(group = BEC, shape = Cluster), position = position_jitterdodge(jitter.width = 0.8, dodge.width = 1)) +
  facet_wrap(~Assay, ncol = 4) +
  theme + scale_fill_manual(name = "Batch Effect Correction", values = c("#F8766D", "#00BFC4", "#619CFF")) + 
  labs(y = "Silhouette Coefficient", x = "Sample Group") + scale_shape_manual(values = c(seq(0,12), 16, 17), name = "Sample Group \n(Batch or Pathological Status)")

ggsave("figures/Silhouette_Coefficient_IRS_limBE_cptac.png", width = 12, height = 12)


# Intragroup Variation (Group) -> PMAD

PMAD_IRS <- get_PMAD_dt(li_se_norm, ain = c("IRS", norm_methods_IRS), condition = li_condition)
PMAD_IRS$BEC <- "IRS"
PMAD_IRS$Normalization <- factor(PMAD_IRS$Normalization, levels = c(levels(PMAD_IRS$Normalization), "IRS_on_log2"))
PMAD_IRS$Normalization[PMAD_IRS$Normalization == "IRS"] <- "IRS_on_log2"
PMAD_IRS$Normalization <- sapply(strsplit(as.character(PMAD_IRS$Normalization), "_on_"), function(x) x[2])

PMAD_single <- get_PMAD_dt(li_se_norm, ain = norm_methods_single, condition = li_condition)
PMAD_single$BEC <- "No"

PMAD_limBE <- get_PMAD_dt(li_se_norm, ain = c("limBE", norm_methods_limBE), condition = li_condition)
PMAD_limBE$BEC <- "limBE"
PMAD_limBE$Normalization <- factor(PMAD_limBE$Normalization, levels = c(levels(PMAD_limBE$Normalization), "limBE_on_log2"))
PMAD_limBE$Normalization[PMAD_limBE$Normalization == "IRS"] <- "limBE_on_log2"
PMAD_limBE$Normalization <- sapply(strsplit(as.character(PMAD_limBE$Normalization), "_on_"), function(x) x[2])

li_pmad_dt <- rbind(PMAD_IRS, PMAD_single, PMAD_limBE)
li_pmad_dt$Normalization <- factor(li_pmad_dt$Normalization, levels = sort(as.character(unique(li_pmad_dt$Normalization))))

PMAD_IRS <- get_PMAD_dt(cptac_se_norm, ain = c(norm_methods_IRS, "IRS"), condition = cptac_condition)
PMAD_IRS$BEC <- "IRS"
PMAD_IRS$Normalization <- factor(PMAD_IRS$Normalization, levels = c(levels(PMAD_IRS$Normalization), "IRS_on_log2"))
PMAD_IRS$Normalization[PMAD_IRS$Normalization == "IRS"] <- "IRS_on_log2"
PMAD_IRS$Normalization <- sapply(strsplit(as.character(PMAD_IRS$Normalization), "_on_"), function(x) x[2])

PMAD_single <- get_PMAD_dt(cptac_se_norm, ain = norm_methods_single, condition = cptac_condition)
PMAD_single$BEC <- "No"

PMAD_limBE <- get_PMAD_dt(cptac_se_norm, ain = c("limBE", norm_methods_limBE), condition = cptac_condition)
PMAD_limBE$BEC <- "limBE"
PMAD_limBE$Normalization <- factor(PMAD_limBE$Normalization, levels = c(levels(PMAD_limBE$Normalization), "limBE_on_log2"))
PMAD_limBE$Normalization[PMAD_limBE$Normalization == "IRS"] <- "limBE_on_log2"
PMAD_limBE$Normalization <- sapply(strsplit(as.character(PMAD_limBE$Normalization), "_on_"), function(x) x[2])

cptac_pmad_dt <- rbind(PMAD_IRS, PMAD_single, PMAD_limBE)
cptac_pmad_dt$Normalization <- factor(cptac_pmad_dt$Normalization, levels = sort(as.character(unique(li_pmad_dt$Normalization))))


comparisons <- list(c("IRS", "limBE"), c("IRS", "No"), c("limBE", "No"))



# Li Stem Cells PMAD Plot
li_pmad_dt$Timepoint <- factor(li_pmad_dt$Timepoint, levels = c("D0", "D3", "D7", "D14"))
li_pmad_dt$BEC <- factor(li_pmad_dt$BEC, levels = c("No", "IRS", "limBE"))

li_pmad_plot <- ggplot(li_pmad_dt, aes(x = BEC, y = PMAD)) + geom_violin(aes(fill = BEC), alpha = 0.4, width = 1) +
  geom_boxplot(aes(fill = BEC), alpha = 0.7, width = 0.1, outliers = FALSE) +
  geom_jitter(aes(color = Normalization, shape = get(li_condition)), size = 2) + 
  stat_compare_means(data = li_pmad_dt, comparisons = comparisons, method = "wilcox.test", paired = TRUE, label.x = 1.5, label.y = c(0.49, 0.52, 0.54)) +
  theme + scale_color_manual(values = col_vector_norm, name = "Normalization Method") + 
  scale_fill_manual(values = c("#F8766D", "#00BFC4", "#619CFF"), name = "Batch Effect Correction") +
  scale_shape_manual(values = c(16, 17, 18, 19)) + 
  guides(shape = guide_legend(title = li_condition, order = 1, position = "top"), 
         fill = "none", 
         color = "none") +
  xlab("Batch Effect Correction")

# CPTAC PMAD Plot
cptac_pmad_dt$BEC <- factor(cptac_pmad_dt$BEC, levels = c("No", "IRS", "limBE"))
cptac_pmad_plot <- ggplot(cptac_pmad_dt, aes(x = BEC, y = PMAD)) + geom_violin(aes(fill = BEC), alpha = 0.4, width = 1) +
  geom_boxplot(aes(fill = BEC), alpha = 0.7, width = 0.1, outliers = FALSE) +
  geom_jitter(aes(color = Normalization, shape = get(cptac_condition)), size = 2) + 
  stat_compare_means(data = cptac_pmad_dt, comparisons = comparisons, method = "wilcox.test", paired = TRUE, label.x = 1.5, label.y = c(0.95, 0.98, 1.01)) +
  theme + scale_color_manual(values = col_vector_norm, name = "Normalization Method") + 
  scale_fill_manual(values = c("#F8766D", "#00BFC4", "#619CFF"), name = "Batch Effect Correction") +
  scale_shape_manual(values = c(16, 17)) +
  guides(shape = guide_legend(title = "Pathological Status", order = 1, position = "top"), 
         fill = "none", 
         color = guide_legend(order = 3, ncol = 1)) +
  xlab("Batch Effect Correction")


li_pmad_plot + cptac_pmad_plot + plot_layout(axis_titles = "collect") + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 20)) 

ggsave("figures/PMAD_real_world_plot.png", width = 13, height = 7, dpi = 300)





# Hierarchical Clustering

# CPTAC Heatmap of Log2

se <- cptac_se_norm

clustering <- sapply(names(assays(se)), function(ain){
  dt <- as.data.table(assays(se)[[ain]])
  clusters <- dendsort::dendsort(stats::hclust(stats::dist(t(dt))))
  return(clusters$labels[clusters$order])
})

clustering <- t(clustering)
clustering <- as.data.frame(clustering)
clustering$Normalization <- row.names(clustering)
clustering <- clustering[clustering$Normalization != "raw",]

cd <- as.data.table(colData(se))

clustering_long <- melt(clustering, value.name = "Column", variable.name = "Position", id.vars = c("Normalization"))
clustering_long <- merge(clustering_long, cd, by = "Column")

clustering_long_IRS <- clustering_long[clustering_long$Normalization %in% c("IRS", norm_methods_IRS),]
clustering_long_IRS$BEC <- "IRS"
clustering_long_IRS$Normalization[clustering_long_IRS$Normalization == "IRS"] <- "IRS_on_log2"
clustering_long_IRS$Normalization <- sapply(strsplit(as.character(clustering_long_IRS$Normalization), "_on_"), function(x) x[2])

clustering_long_single <- clustering_long[clustering_long$Normalization %in% norm_methods_single,]
clustering_long_single$BEC <- "No"

clustering_long_limBE <- clustering_long[clustering_long$Normalization %in% c("limBE", norm_methods_limBE),]
clustering_long_limBE$BEC <- "limBE"
clustering_long_limBE$Normalization[clustering_long_limBE$Normalization == "limBE"] <- "limBE_on_log2"
clustering_long_limBE$Normalization <- sapply(strsplit(as.character(clustering_long_limBE$Normalization), "_on_"), function(x) x[2])

clustering_new <- rbind(clustering_long_IRS, clustering_long_single, clustering_long_limBE)

clustering_new$BEC <- factor(clustering_new$BEC, levels = c("No", "IRS", "limBE"))
clustering_new$Pool <- plyr::mapvalues(clustering_new$Pool, from, to)
clustering_new$Pool <- factor(clustering_new$Pool, levels = to)

qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual",]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[17:33]

heatmap_batch <- ggplot(clustering_new, aes(x = Position, y = Normalization, fill = Pool)) + geom_tile(color = "white") +
  scale_fill_manual(values = col_vector, name = "Batch") +
  facet_wrap(~BEC) + 
  theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + labs(x = "Samples", y = "Normalization Method")

heatmap_condition <- ggplot(clustering_new, aes(x = Position, y = Normalization, fill = Pathological_Status)) + geom_tile(color = "white") +
  scale_fill_brewer(palette = "Paired", name = "Pathological Status") +
  facet_wrap(~BEC) + 
  theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + labs(x = "Samples", y = "Normalization Method")

(heatmap_batch / heatmap_condition) + plot_layout(axis_titles = "collect") + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 20)) 

ggsave("figures/hierarchical_clustering_CPTAC.png", width = 11, height = 9, dpi = 300)





# PCA


# PCA Plots

## Li Stem Cells

qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == "qual", ]
col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- rev(col_vector)

### Single 
single_PCA_dt <- PRONE:::get_complete_pca_dt(li_se_norm, norm_methods_single)
coldata <- as.data.table(colData(li_se_norm))
single_PCA_dt <- merge(single_PCA_dt, coldata, by = "Column")

p1 <- ggplot(single_PCA_dt, aes( x= PC1, y = PC2, color = Assay, shape = Batch)) + 
  geom_point(size = 2) + theme  +
  scale_color_manual(name = "Normalization Method", values = col_vector_norm) +
  guides(shape = "none")

p1_timepoint <- ggplot(single_PCA_dt, aes( x= PC1, y = PC2, color = Timepoint, shape = Batch)) + 
  geom_point(size = 2) + theme + guides(color = "none") +
  scale_color_manual(name = "Timepoint", 
                     values = col_vector)
### IRS
IRS_PCA_dt <- PRONE:::get_complete_pca_dt(li_se_norm, c("IRS", norm_methods_IRS))
IRS_PCA_dt <- merge(IRS_PCA_dt, coldata, by = "Column")
IRS_PCA_dt$Assay[IRS_PCA_dt$Assay == "IRS"] <- "IRS_on_log2"
IRS_PCA_dt$Assay <- sapply(strsplit(IRS_PCA_dt$Assay, "_on_"), function(x) x[2]) 

p2 <- ggplot(IRS_PCA_dt, aes( x= PC1, y = PC2, color = Assay, shape = Batch)) + geom_point(size = 2) + theme + 
  scale_color_manual(name = "Normalization Method", values = col_vector_norm) + guides(color = "none")

p2_timepoint <- ggplot(IRS_PCA_dt, aes( x= PC1, y = PC2, color = Timepoint, shape = Batch)) + 
  geom_point(size = 2) + theme  +
  scale_color_manual(name = "Timepoint", 
                     values = col_vector)

### LimBE
limBE_PCA_dt <- PRONE:::get_complete_pca_dt(li_se_norm, c("limBE", norm_methods_limBE))
limBE_PCA_dt <- merge(limBE_PCA_dt, coldata, by = "Column")
limBE_PCA_dt$Assay[limBE_PCA_dt$Assay == "limBE"] <- "limBE_on_log2"
limBE_PCA_dt$Assay <- sapply(strsplit(limBE_PCA_dt$Assay, "_on_"), function(x) x[2]) 

p3 <- ggplot(limBE_PCA_dt, aes( x= PC1, y = PC2, color = Assay, shape = Batch)) + geom_point(size = 2) + theme + 
  scale_color_manual(name = "Normalization Method", values = col_vector_norm) + guides(color = "none")

p3_timepoint <- ggplot(limBE_PCA_dt, aes( x= PC1, y = PC2, color = Timepoint, shape = Batch)) + 
  geom_point(size = 2) + theme  +
  scale_color_manual(name = "Timepoint", 
                     values = col_vector)



((p1 + p1_timepoint) / (p2 + p2_timepoint) / (p3 + p3_timepoint) + plot_layout(guides = "collect")) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15)) 

ggsave("figures/li_pca_plots.png", width = 12, height = 12, dpi = 300)


# CPTAC PCA Plots

single_PCA_dt <- PRONE:::get_complete_pca_dt(cptac_se_norm, norm_methods_single)
coldata <- as.data.table(colData(cptac_se_norm))
single_PCA_dt <- merge(single_PCA_dt, coldata, by = "Column")

p1 <- ggplot(single_PCA_dt, aes( x= PC1, y = PC2, color = Assay)) + 
  geom_point(size = 2) + theme  +
  scale_color_manual(name = "Normalization Method", values = col_vector_norm)

p1_condition <- ggplot(single_PCA_dt, aes( x= PC1, y = PC2, color = Pathological_Status)) + 
  geom_point(size = 2) + theme +
  scale_color_manual(name = "Pathological Status", 
                     values = col_vector)

p1_batch <- ggplot(single_PCA_dt, aes( x= PC1, y = PC2, color = Pool)) + 
  geom_point(size = 2) + theme + scale_color_manual(name = "Batch", values = col_vector)

IRS_PCA_dt <- PRONE:::get_complete_pca_dt(cptac_se_norm, norm_methods_IRS)
IRS_PCA_dt <- merge(IRS_PCA_dt, coldata, by = "Column")
IRS_PCA_dt$Assay <- sapply(strsplit(IRS_PCA_dt$Assay, "_on_"), function(x) x[2]) 

p2 <- ggplot(IRS_PCA_dt, aes( x= PC1, y = PC2, color = Assay)) + geom_point(size = 2) + theme + 
  scale_color_manual(name = "Normalization Method", values = col_vector_norm) + guides(color = "none")

p2_condition <- ggplot(IRS_PCA_dt, aes( x= PC1, y = PC2, color = Pathological_Status)) + 
  geom_point(size = 2) + theme  +
  scale_color_manual(name = "Pathological Status", 
                     values = col_vector)

p2_batch <- ggplot(IRS_PCA_dt, aes( x= PC1, y = PC2, color = Pool)) + 
  geom_point(size = 2) + theme + scale_color_manual(name = "Batch", values = col_vector)

((p1 + p2 + plot_layout(guides = "collect")) / (p1_condition + p2_condition + plot_layout(guides = "collect")) / (p1_batch + p2_batch + plot_layout(guides = "collect"))) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold")) 

ggsave("figures/cptac_pca_plots.png", width = 15, height = 18, dpi = 300)

pca_single_color_pool <- plot_PCA(cptac_se_norm, ain = norm_methods_single, color_by = "Pool", label_by = "No", shape_by = "Pathological_Status") + theme
pca_IRS_color_pool <- plot_PCA(cptac_se_norm, ain = norm_methods_IRS, color_by = "Pool", label_by = "No", shape_by = "Pathological_Status") + theme






# # Check variance of the data before and after batch effect correction (within-batch, within-group, between-batch, between-group)
# library(variancePartition)
# 
# calculate_group_batch_variance <- function(se, ain, condition, batch) {
#   form <- as.formula(paste0("~ (1|", condition, ") + (1|", batch, ")"))
#   m <- as.matrix(assays(se)[[ain]])
#   m[is.na(m)] <- NA
#   #m <- na.omit(m)
#   md <- as.data.frame(colData(se))
#   md <- md[, c("Column", condition, batch), ]
#   # make factors of the categorical variables
#   md[, condition] <- factor(md[, condition])
#   md[, batch] <- factor(md[, batch])
#   fit <- fitExtractVarPartModel(m, form, data = md)
#   variances <- data.frame(n1 = fit[[1]], n2 = fit[[2]])
#   colnames(variances) <- colnames(fit)[1:2]
#   return(fit)
# }
# 
# # Cell Culture (Li) Dataset
# full_variances <- NULL
# all_assays <- names(assays(li_se_norm))
# all_assays <- all_assays[all_assays != "raw"]
# for (ain in all_assays) {
#   print(ain)
#   variances <- calculate_group_batch_variance(li_se_norm, ain, li_condition, li_batch)
#   variances$Assay <- ain
#   if (is.null(full_variances)) {
#     full_variances <- variances
#   } else {
#     full_variances <- rbind(full_variances, variances)
#   }
# }
# 
# saveRDS(full_variances, "data/variances_li_stem_cells.rds")
# 
# # Cancer (Hu) Dataset
# full_variances <- NULL
# all_assays <- names(assays(cptac_se_norm))
# all_assays <- all_assays[all_assays != "raw"]
# for (ain in all_assays) {
#   print(ain)
#   variances <- calculate_group_batch_variance(cptac_se_norm, ain, cptac_condition, cptac_batch)
#   variances$Assay <- ain
#   if (is.null(full_variances)) {
#     full_variances <- variances
#   } else {
#     full_variances <- rbind(full_variances, variances)
#   }
# }
# 
# saveRDS(full_variances, "data/variances_cptac_ovarian.rds")
# full_variances <- readRDS("data/variances_cptac_ovarian.rds")
# 
# condition <- cptac_condition
# batch <- cptac_batch
# 
# full_variances$BEC <- "No"
# full_variances$BEC[grepl("IRS_on_", full_variances$Assay)] <- "IRS"
# full_variances$Assay[grepl("IRS_on_", full_variances$Assay)] <- sapply(strsplit(as.character(full_variances$Assay[grepl("IRS_on_", full_variances$Assay)]), "_on_"), function(x)
#   x[2])
# 
# full_variances$BEC[grepl("limBE_on_", full_variances$Assay)] <- "limBE"
# full_variances$Assay[grepl("limBE_on_", full_variances$Assay)] <- sapply(strsplit(as.character(full_variances$Assay[grepl("limBE_on_", full_variances$Assay)]), "_on_"), function(x)
#   x[2])
# 
# full_variances$Residuals <- NULL
# full_variances$BEC <- factor(full_variances$BEC, levels = c("No", "IRS", "limBE"))
# full_variances <- melt(
#   full_variances,
#   id.vars = c("Assay", "BEC"),
#   measure.vars = c(batch, condition),
#   value.name = "Variance",
#   variable.name = "Variable"
# )
# ggplot(full_variances[!full_variances$Assay %in% c("log2", "IRS", "limBE"), ], aes(x = Variable, y = Variance, fill = BEC)) + geom_boxplot() + facet_wrap( ~
#                                                                                                                                                              Assay) +
#   labs(x = "Variable", y = "Proportion Variance (%)") + theme_bw()
# 
# 
# ggplot(full_variances[!full_variances$Assay %in% c("log2", "IRS", "limBE"), ], aes(x = BEC, y = Variance, fill = Variable)) + geom_boxplot() + facet_wrap( ~
#                                                                                                                                                              Assay) +
#   labs(x = "Batch Effect Correction", y = "Proportion Variance (%)") + theme_bw()
# 
# # remove cases with variable = Pool AND BEC = NO or IRS
# test <- full_variances[full_variances$Variable != "Pool" |
#                          full_variances$BEC == "limBE", ]
# ggplot(test[!test$Assay %in% c("log2", "IRS", "limBE"), ], aes(x = BEC, y = Variance, fill = Variable)) + geom_boxplot(outliers = FALSE) + facet_wrap( ~
#                                                                                                                                                          Assay) +
#   labs(x = "Batch Effect Correction", y = "Proportion Variance (%)") + theme_bw()
# 
# 
# # fitExtractVarPartModel (categorical variables modeled as random effects -> (1|c), standard linear model with + c)
# # - exprObj: matrix with intenstiies (proteins x samples)
# # - formula: specified variables for the linear mixed model (only specify covariates since rows of exprObj are automatically used as a response (a + b + (1|c)))
# # - data: data.frame with columns corresponding to formula


# Mouse P450 Dataset
data_name <- "mouse_liver_cytochrome_P450"
p450_de_res <- readRDS(paste0("data/de_results_real_world/", data_name, "_de_res.rds"))
p450_se_norm <- readRDS(paste0("data/preprocessed_normalized_real_world_se/", data_name, "_se_norm.rds"))
se <- p450_se_norm
condition <- "Condition"

pc_dt <- get_complete_PCA_new(se, NULL)

# replace values by numbers
# create mapping of column names to numbers
condition_num <- as.factor(as.numeric(as.factor(se[[condition]])))
names(condition_num) <- se[["Column"]]
  
# mapping
mapping_cond <- unique(data.table(Name = se[[condition]], Value = condition_num))

  
all_silh <- NULL
for(ain in unique(pc_dt$Assay)){
  print(ain)
  pc_ain <- pc_dt[pc_dt$Assay == ain, c("PC1", "PC2","PC3", "Column")]
  # order pc_ain by Column of colData
  pc_ain <- pc_ain[match(rownames(colData(se)), pc_ain$Column),]
  pc_ain$Column <- NULL
  pc_ain <- as.matrix(pc_ain)
  rownames(pc_ain) <- se[["Column"]]
  silh <- calc.sil(as.matrix(pc_ain), y1 = condition_num, y2 = NULL, name.y1 = condition, name.y2 = batch)
  silh$Assay <- ain
  if(is.null(all_silh)){
    all_silh <- silh
  } else {
    all_silh <- rbind(all_silh, silh)
  }
}

# remap cluster and batch numbers to values
all_silh$Cluster <- mapping_cond[match(all_silh$Cluster, mapping_cond$Value), "Name"][[1]]

# plot
ggplot(all_silh, aes(x = Assay, y = silh.coeff, color = Assay)) + 
  geom_point(aes(shape = Cluster), size = 3) +
  theme + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "None") +
  labs(x = "Normalization Method", y = "Silhouette Coefficient") +
  scale_color_manual(values = col_vector_norm)
  