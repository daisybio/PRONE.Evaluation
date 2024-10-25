# Set data path
data_name <- "CPTAC_OV"

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
if(!file.exists(paste0("data/de_results_real_world/", data_name, "_de_res.rds"))){
  # Perform preprocessing, normalization, DE analysis
  source("real-world/Real-World_Pipeline.R")
}

# Read DE results and normalized object
de_res <- readRDS(paste0("data/de_results_real_world/", data_name, "_de_res.rds"))
se_norm <- readRDS(paste0("data/preprocessed_normalized_real_world_se/", data_name, "_se_norm.rds"))

# Evaluation

# 1. Unnormalized data

log2_pca <- plot_PCA(se_norm, ain = "log2", color_by = "Pool", label_by = "No", shape_by = "Pathological_Status") + theme
log2_boxplot <- plot_boxplots(se_norm, ain = "log2", color_by = "Pathological_Status") + theme

log2_pca + log2_boxplot


# Define normalization methods

assays <- names(assays(se_norm))
norm_methods_IRS <- assays[grepl("IRS_on_", assays)]
norm_methods_limBE <- assays[grepl("limBE_on_", assays)]
norm_methods_single <- assays[!assays %in% c("limBE", "IRS", norm_methods_IRS, norm_methods_limBE)]


# Quantitative and Qualitative Evaluation

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

PMAD_IRS <- get_PMAD_dt(se_norm, ain = norm_methods_IRS, condition = "Pathological_Status")
PMAD_IRS$BEC <- "Yes"
PMAD_IRS$Normalization <- sapply(strsplit(as.character(PMAD_IRS$Normalization), "_on_"), function(x) x[2])

PMAD_single <- get_PMAD_dt(se_norm, ain = norm_methods_single, condition = "Pathological_Status")
PMAD_single$BEC <- "No"

dt <- rbind(PMAD_IRS, PMAD_single)

p1 <- ggplot(dt, aes(x = BEC, y = PMAD)) + geom_boxplot(aes(fill = BEC), alpha = 0.4, outliers = FALSE) + geom_jitter(aes(color = Normalization, shape = get(metadata(se_norm)$condition)), size = 3) + 
  theme + scale_color_manual(values = col_vector_norm, name = "Normalization Method") + 
  scale_fill_manual(values = c("#F8766D", "#00BFC4"), name = "Batch Effect Corrected (IRS)") +
  guides(shape = guide_legend(title = metadata(se_norm)$condition))

p2_li <- ggplot(dt, aes(x = BEC, y = PMAD)) + geom_violin(aes(fill = BEC), alpha = 0.4) +
  geom_boxplot(aes(fill = BEC), alpha = 0.7, width = 0.3, outliers = FALSE) +
  geom_jitter(aes(color = Normalization, shape = get(metadata(se_norm)$condition)), size = 3) + 
  theme + scale_color_manual(values = col_vector_norm, name = "Normalization Method") + 
  scale_fill_manual(values = c("#F8766D", "#00BFC4"), name = "Batch Effect Corrected (IRS)") +
  guides(shape = guide_legend(title = metadata(se_norm)$condition)) +
  xlab("Batch Effect Corrected (IRS)")

p3 <- ggplot(dt, aes(x = BEC, y = PMAD, fill = get(metadata(se_norm)$condition))) + geom_split_violin(alpha = 0.4) +
  theme + scale_color_manual(values = col_vector_norm, name = "Normalization Method") + 
  scale_fill_manual(values = c("#F8766D", "#00BFC4"), name = "Batch Effect Corrected (IRS)")

p4 <- ggplot(dt, aes(x = BEC, y = PMAD, fill = Pathological_Status)) + geom_violin(alpha = 0.4) + 
  geom_jitter(aes(color = Normalization, group = Pathological_Status), size = 3, position=position_jitterdodge(jitter.width = 2)) + 
  theme + scale_color_manual(values = col_vector_norm, name = "Normalization Method") + 
  scale_fill_manual(values = c("#F8766D", "#00BFC4"), name = "Pathological Status") +
  guides(shape = guide_legend(title = metadata(se_norm)$condition)) + xlab("Batch Effect Corrected With IRS")


(p1 + p2) / (p3 + p4)

plot_PCA(se_norm, ain = c("RobNorm", "Median", "NormicsVSN"), shape_by = "Pathological_Status", color_by = "Pool", label_by = "No")
plot_PCA(se_norm, ain = c("RobNorm", "Median", "NormicsVSN"), color_by = "Pathological_Status", label_by = "No")

plot_PCA(se_norm, ain = c("IRS_on_RobNorm", "IRS_on_Median", "IRS_on_NormicsVSN"), shape_by = "Pathological_Status", color_by = "Pool", label_by = "No")
plot_PCA(se_norm, ain = c("IRS_on_RobNorm", "IRS_on_Median", "IRS_on_NormicsVSN"), color_by = "Pathological_Status", label_by = "No")

ain <- c("IRS_on_RobNorm", "IRS_on_Median", "IRS_on_NormicsVSN", "IRS_on_MAD", "IRS_on_EigenMS")
plot_PCA(se_norm, ain = ain, color_by = "Pathological_Status", label_by = "No")
plot_PCA(se_norm, ain = ain, color_by = "Anatomic_Site_Tumor", label_by = "No")
plot_PCA(se_norm, ain = ain, color_by = "Origin_Site_Disease", label_by = "No")

plot_PCA(se_norm, ain = norm_methods_single, color_by = "Pool", label_by = "No", shape_by = "Pathological_Status", ellipse = FALSE, ncol = 3) + theme
plot_PCA(se_norm, ain = norm_methods_single, color_by = "Pathological_Status", label_by = "No", shape_by = "No", ellipse = FALSE, ncol = 3) + theme
ggsave("figures/CPTAC_OV/PCA_single_methods.png", width = 13, height = 15)
plot_PCA(se_norm, ain = norm_methods_IRS, color_by = "Pathological_Status", label_by = "No", shape_by = "No", ellipse = FALSE, ncol = 3) + theme
ggsave("figures/CPTAC_OV/PCA_IRS_methods.png", width = 13, height = 15)


# PCA Overall
methods_to_exclude <- NULL

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
p1 <- ggplot(single_PCA_dt, aes( x= PC1, y = PC2, color = Assay, shape = get("Pathological_Status"))) + geom_point(size = 3) + theme + 
  scale_color_manual(name = "Normalization Method", values = col_vector_norm) + guides(shape = guide_legend(title = "Pathological Status"))
p1_timepoint <- ggplot(single_PCA_dt, aes( x= PC1, y = PC2, color = get("Pathological_Status"))) + 
  geom_point(size = 3) + theme + guides(color = guide_legend(order = 2)) +
  scale_color_manual(name = "Pathological Status", 
                     values = col_vector)

IRS_PCA_dt <- PRONE:::get_complete_pca_dt(se_norm, norm_methods_IRS)
IRS_PCA_dt <- merge(IRS_PCA_dt, coldata, by = "Column")
IRS_PCA_dt$Assay <- sapply(strsplit(IRS_PCA_dt$Assay, "_on_"), function(x) x[2]) 
p2 <- ggplot(IRS_PCA_dt, aes( x= PC1, y = PC2, color = Assay, shape = get("Pathological_Status"))) + geom_point(size = 3) + theme + 
  scale_color_manual(name = "Normalization Method", values = col_vector_norm) + guides(shape = guide_legend(title = "Pathological Status"))
p2_timepoint <- ggplot(IRS_PCA_dt, aes( x= PC1, y = PC2, color = get("Pathological_Status"))) + 
  geom_point(size = 3) + theme + guides(color = guide_legend(order = 2)) +
  scale_color_manual(name = "Pathological Status", 
                     values = col_vector)


# Combine plots 1 and 2 with a shared legend
pca_idv <- ((p1 + p2) / (p1_timepoint + p2_timepoint)) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")






# PCA of all methods
assays <- names(assays(se_norm))
#methods_to_remove <- c("VSN", "IRS_on_VSN", "limBE_on_VSN", "MAD", "IRS_on_MAD", "limBE_on_MAD", "EigenMS", "IRS_on_EigenMS", "limBE_on_EigenMS", "TMM", "IRS_on_TMM", "limBE_on_TMM")
#assays <- assays[!assays %in% methods_to_remove]
norm_methods_IRS <- assays[grepl("IRS_on_", assays)]
norm_methods_limBE <- assays[grepl("limBE_on_", assays)]
norm_methods_single <- assays[!assays %in% c("limBE", "IRS", "log2", norm_methods_IRS, norm_methods_limBE)]


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



plot_densities(se_norm, ain = ain, color_by = "Pathological_Status", ncol = 3)

ain <- c("IRS_on_RobNorm", "IRS_on_Median", "IRS_on_NormicsVSN")
plot_boxplots(se_norm, ain = ain, color_by = "Pathological_Status")

# Complete PCA plot of all normalization methods
assays <- names(assays(se_norm))

#assays <- assays[!assays %in% c("EigenMS", "MAD", "IRS_on_EigenMS", "limBE_on_EigenMS", "IRS_on_MAD", "limBE_on_MAD")]#, "VSN", "IRS_on_VSN", "limBE_on_VSN", "TMM", "IRS_on_TMM", "limBE_on_TMM")]

# find assays with _on_IRS
norm_methods_IRS <- assays[grepl("IRS_on_", assays)]
se_IRS <- PRONE::generate_complete_SE(se_norm, ain = norm_methods_IRS)
colData(se_IRS)$Normalization <- sapply(strsplit(colData(se_IRS)$Normalization, "_on_"), function(x) x[2])

# find assays with _on_limBE
se_limBE <- PRONE::generate_complete_SE(se_norm, ain = norm_methods_limBE)
colData(se_limBE)$Normalization <- sapply(strsplit(colData(se_limBE)$Normalization, "_on_"), function(x) x[2])

# single normalization methods
norm_methods_single <- assays[!assays %in% c("limBE", "IRS", "log2", norm_methods_IRS, norm_methods_limBE)]
se_single <- PRONE::generate_complete_SE(se_norm, ain = norm_methods_single)

# PCA plots
pca_single <- plot_PCA(se_single, color_by = "Normalization", ain = c("all"), label_by = "No", shape_by = "Pathological_Status") + 
  theme + theme(strip.text = element_blank()) +
  ggtitle("Single Normalization Methods")# + scale_color_manual(values = col_vector_norm)

pca_IRS <- plot_PCA(se_IRS, color_by = "Normalization", ain = c("all"), label_by = "No", shape_by = "Pathological_Status") + 
  theme + theme(strip.text = element_blank()) + 
  ggtitle("IRS Normalization") # + scale_color_manual(values = col_vector_norm)

pca_limBE <- plot_PCA(se_limBE, color_by = "Normalization", ain = c("all"), label_by = "No", shape_by = "Pool") + 
  theme + theme(strip.text = element_blank()) +
  ggtitle("limBE Normalization") # + scale_color_manual(values = col_vector_norm)

pca_single + pca_IRS + pca_limBE + patchwork::plot_layout(guides = "collect")

# DEP Comparison

# Methods based on IRS-normalized data
de_res_IRS <- de_res[de_res$Assay %in% c(norm_methods_IRS),]
de_res_IRS$Assay <- sapply(strsplit(de_res_IRS$Assay, "_on_"), function(x) x[2]) 

plot_overview_DE_bar(de_res_IRS, plot_type = "single")[[1]]  +  ggtitle("") + theme + theme(strip.background =element_rect(fill="white")) + ylab("Normalization Method") + theme



new_de_res_IRS <- PRONE::apply_thresholds(de_res_IRS, logFC_up = 2, logFC_down = -2)

overview_de_bar <- plot_overview_DE_bar(new_de_res_IRS, plot_type = "single")[[1]]  +  ggtitle("") + theme + theme(strip.background =element_rect(fill="white")) + ylab("Normalization Method") + theme

jaccard_heatmap <- plot_jaccard_heatmap(new_de_res_IRS, plot_type = "all")

(overview_de_bar / jaccard_heatmap) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

plot_volcano_DE(new_de_res_IRS)[[1]] + theme + theme(strip.background =element_rect(fill="white"))

# Stacked barplot of overlaps
plot_upset_DE(new_de_res_IRS, plot_type = "single")

upset_res <- plot_upset_DE(new_de_res_IRS, plot_type = "stacked", min_degree = 15)
upset_res_table <- upset_res$table
unique_res <- upset_res_table[upset_res_table$`Number of Intersected Assays` == 1,]

# Statistical enrichment analysis


dt <- new_de_res_IRS[!new_de_res_IRS$Assay %in% c("IRS_on_TMM", "IRS_on_VSN"),]
ain <- unique(dt$Assay)
id_column <- "Gene.Names"
queries <- lapply(ain, function(method) {
  tmp <- dt[dt$Assay == method, ]
  query <- tmp[[id_column]]
  query <- unique(query[query != ""])
  query
})
names(queries) <- ain
gprofiler_res <- gprofiler2::gost(query = queries, organism = "mmusculus", sources = c("KEGG"), correction_method = "fdr")
gprofiler_res <- gprofiler_res$result

ggplot(gprofiler_res, aes(x = query, y = term_name, fill = p_value)) + geom_tile() +
  theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + labs(x = "Normalization Method", y = "KEGG Patheay", fill = "P-Value")


# Consensus
dt <- new_de_res_IRS[!new_de_res_IRS$Assay %in% c("IRS_on_TMM", "IRS_on_VSN"),]
DE_candidates <- extract_consensus_DE_candidates(dt, ain = NULL, norm_thr = 0.9, per_comparison = TRUE)

# check distribution of p.adjust values in the different normalization methods

ggplot(new_de_res_IRS, aes(x = adj.P.Val)) +
  geom_histogram(alpha = 0.5) +
  facet_wrap(~Assay, scales = "free") +
  labs(x = "p.adjust", fill = "Normalization Method") +
  theme

ggplot(new_de_res_IRS, aes(x = logFC)) +
  geom_histogram(alpha = 0.5) +
  facet_wrap(~Assay, scales = "free") +
  labs(x = "logFC", fill = "Normalization Method") +
  theme

ggplot(new_de_res_IRS, aes(x = P.Value)) +
  geom_histogram(alpha = 0.5) +
  facet_wrap(~Assay, scales = "free") +
  labs(x = "P.Value", fill = "Normalization Method") +
  theme


# VSN Normalization

# VSN 
se_VSN <- subset_SE_by_norm(se_norm, ain = c("raw", "log2"))#, "VSN", "RobNorm", "Median"))#, "limBE_on_VSN", "IRS_on_VSN", "IRS_on_RobNorm", "IRS_on_Median", "IRS_on_Median", "limBE_on_Median"))

se_VSN <- vsnNorm(se_VSN, ain = "raw", aout = "VSN_0.9", VSN_quantile = 0.9)
se_VSN <- vsnNorm(se_VSN, ain = "raw", aout = "VSN_0.8", VSN_quantile = 0.8)
se_VSN <- vsnNorm(se_VSN, ain = "raw", aout = "VSN_0.5", VSN_quantile = 0.5)
se_VSN <- vsnNorm(se_VSN, ain = "raw", aout = "VSN_1", VSN_quantile = 1)


se_VSN_big <- generate_complete_SE(se_VSN)
plot_PCA(se_VSN_big, color_by = "Normalization", label_by = "No", shape_by = "Batch")

se_VSN <- normalize_se_custom(se, methods = c("VSN", "IRS_on_VSN", "limBE_on_VSN"))



normalize_se_custom <- function(se, methods, combination_pattern = "_on_", gamma.0 = 0.5, reduce_correlation_by = 1, NormicsVSN_quantile = 0.8, top_x = 50, VSN_quantile = 0.9){
  browser()
  # extract combination of methods
  if(!is.null(combination_pattern)){
    comb_methods <- methods[stringr::str_detect(methods, combination_pattern)] #  combined methods
    sing_methods <- methods[!methods %in% comb_methods] # single methods
  } else {
    sing_methods <- methods
  }
  
  # single normalization
  if(length(sing_methods) > 0){
    se <- normalize_se_single(se, sing_methods, gamma.0 = gamma.0, reduce_correlation_by = reduce_correlation_by, NormicsVSN_quantile = NormicsVSN_quantile, top_x = top_x, VSN_quantile = VSN_quantile)
  }
  # combined normalization
  if(!is.null(combination_pattern)){
    if(length(comb_methods) > 0){
      for(m in comb_methods){
        methods_split <- strsplit(m, combination_pattern)[[1]]
        
        # Initialize method and ain with the last two elements
        length_split <- length(methods_split)
        ain <- methods_split[length_split]
        method <- methods_split[length_split - 1]
        se <- normalize_se_combination(se, c(method), c(ain), combination_pattern, gamma.0 = gamma.0, reduce_correlation_by = reduce_correlation_by, NormicsVSN_quantile = NormicsVSN_quantile, top_x = top_x, VSN_quantile = VSN_quantile)
        
        # If there are more than two methods in the combination, process the rest
        if (length(methods_split) > 2) {
          for (i in (length_split - 2):1) {
            ain <- paste(method, ain, sep = combination_pattern)
            method <- methods_split[i]
            
            # Apply the subsequent combination
            se <- normalize_se_combination(se, c(method), c(ain), combination_pattern, gamma.0 = gamma.0, reduce_correlation_by = reduce_correlation_by, NormicsVSN_quantile = NormicsVSN_quantile, top_x = top_x, VSN_quantile = VSN_quantile)
          }
        }
      }
    }
  }
  return(se)
}



# Load raw se

source("real-world/data_preparation/li_stem_cells_prep.R")

se <- medianNorm(se, ain = "raw", aout = "Median")

se <- irsNorm(se, ain = "Median", aout = "IRS_on_Median")

se <- limmaNorm(se, ain = "Median", aout = "limBE_on_Median")

se <- limmaNorm(se, ain = "IRS_on_Median", aout = "limBE_on_IRS_on_Median")


plot_PCA(se, color_by = "Timepoint", shape_by = "Batch") + theme
