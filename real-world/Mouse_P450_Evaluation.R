# Set data path
data_name <- "mouse_liver_cytochrome_P450"

# Set real-world configuration for this data set
TMT <- FALSE
NA_thr <- 0.8
DE_method <- "limma"
logFC <- TRUE
logFC_up <- 0
logFC_down <- 0
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

se_norm <- vsnNorm(se_norm, ain = "raw", aout = "VSN_0.5", VSN_quantile = 0.5)
se_norm <- normicsNorm(se_norm, ain = "raw", aout = "NormicsVSN_0.5", method = "NormicsVSN", NormicsVSN_quantile = 0.5)

# Evaluation

col_vector_norm <- c(col_vector_norm, "log2" = "yellow")

# 1. Unnormalized data

log2_pca <- plot_PCA(se_norm, ain = "log2", color_by = "Condition", label_by = "No") + theme
log2_boxplot <- plot_boxplots(se_norm, ain = "log2", color_by = "Condition") + theme

log2_pca / log2_boxplot

#plot_heatmap(se_norm, ain = "log2", color_by = c("Batch", "Timepoint"))


# Quantitative and Qualitative Evaluation

plot_intragroup_PMAD(se_norm) + theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "None") +
  ggtitle("Intragroup PMAD") + scale_fill_manual(values = col_vector_norm)

# PCA, Boxplots, densities of single normalization methods

plot_PCA(se_norm, label_by = "No")
plot_boxplots(se_norm, ncol = 4)
plot_densities(se_norm, ncol = 4)

# Complete PCA plot of all normalization methods

# All normalization methods
se_big_norms <- generate_complete_SE(se_norm, ain =  NULL)
plot_PCA(se_big_norms, color_by = "Normalization", ain = c("all"), label_by = "No", shape_by = "Condition") + 
  theme + theme(strip.text = element_blank()) +
  ggtitle("Single Normalization Methods") + scale_color_manual(values = col_vector_norm)


single_PCA_dt <- PRONE:::get_complete_pca_dt(se_norm, ain = NULL)
coldata <- as.data.table(colData(se_norm))
single_PCA_dt <- merge(single_PCA_dt, coldata, by = "Column")
p1 <- ggplot(single_PCA_dt, aes( x= PC1, y = PC2, color = Assay, shape = Condition)) + geom_point(size = 3) + theme + scale_color_manual(name = "Normalization Method", values = col_vector_norm)



assays <- names(assays(se_norm))
assays <- assays[! assays %in% c("MAD")]

se_big_norms <- generate_complete_SE(se_norm, ain =  assays)
plot_PCA(se_big_norms, color_by = "Normalization", ain = c("all"), label_by = "No", shape_by = "Condition") + 
  theme + theme(strip.text = element_blank()) +
  ggtitle("Single Normalization Methods") + scale_color_manual(values = col_vector_norm)

# DEP Comparison

plot_overview_DE_bar(de_res)

plot_volcano_DE(de_res)

de_res <- de_res[!de_res$Assay %in% c("EigenMS", "MAD", "TMM", "log2")]


# Intersection Analysis

plot_jaccard_heatmap(de_res)

plot_upset_DE(de_res, plot_type = "single", min_degree = 15)



# Biomarker Research

# convert proteins to genes
signif_res <- de_res[de_res$Change != "No Change",]
signif_res$Gene.Names <- NULL
gprofiler_res <- gprofiler2::gconvert(query = unique(signif_res$Protein.IDs), organism = "mmusculus", target = "ENTREZGENE")
gprofiler_res <- gprofiler_res[, c("input", "target")]
colnames(gprofiler_res) <- c("Protein.IDs", "Gene.Names")

signif_res <- merge(signif_res, gprofiler_res, by = "Protein.IDs", all.x = TRUE)

# Biomarker mRNA Extraction (Original Publication)
mRNA_DEGs <- fread("data/raw_real_world_files/mouse_liver_cytochrome_P450/supp_mRNA_DEGs.txt")
mRNA_DEGs <- unique(mRNA_DEGs)
# Run gprofiler to get a better coverage
gprofiler_res_DEGs <- gprofiler2::gconvert(query = mRNA_DEGs$Gene.Names, organism = "mmusculus", target = "ENTREZGENE")
gprofiler_res_DEGs <- gprofiler_res_DEGs[, c("target")]

tmp <- signif_res[signif_res$Gene.Names %in% gprofiler_res_DEGs,]
ggplot(tmp, aes(x = Assay, fill = Assay)) + geom_bar(stat = "count") + theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "None") +
  ggtitle("Biomarker mRNA Extraction (Original Publication)") + scale_fill_manual(values = col_vector_norm)

# Statistical Enrichment
dt <- signif_res[!signif_res$Assay %in% c("EigenMS", "TMM", "log2", "MAD"),]
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
ggsave("figures/P450_enrichment.png", width = 10, height = 14)

# Only include those that are not found by all methods

gprof_res_wide <- dcast(gprofiler_res[, c("query", "term_name", "p_value")], term_name ~ query, value.var = "p_value")
gprof_res_wide_not_complete <- gprof_res_wide[apply(gprof_res_wide, 1, function(r) any(is.na(r))),]

# hierarchical clustering of gprof_res_wide
terms <- gprof_res_wide_not_complete$term_name
gprof_res_wide_not_complete$term_name <- NULL
gprof_res_wide_not_complete[!is.na(gprof_res_wide_not_complete)] <- 1
gprof_res_wide_not_complete[is.na(gprof_res_wide_not_complete)] <- 0
gprof_res_not_complete_matrix <- as.matrix(gprof_res_wide_not_complete)
rownames(gprof_res_not_complete_matrix) <- terms
hclust_obj <- hclust(dist(gprof_res_not_complete_matrix))

gprof_shared <- gprofiler_res[gprofiler_res$term_name %in% terms,]
gprof_shared$term_name <- factor(gprof_shared$term_name, levels = hclust_obj$labels[hclust_obj$order])

ggplot(gprof_shared, aes(x = query, y = term_name, fill = p_value)) + geom_tile() +
  theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + labs(x = "Normalization Method", y = "KEGG Patheay", fill = "P-Value")

# check number of shared, unique, partially shared terms

gprof_shared_dt <- data.table("term_name" = gprof_res_wide$term_name, "intersections" = apply(gprof_res_wide[,-1], 1, function(r) sum(!is.na(r))))
gprof_shared_dt$class <- "partially shared"
gprof_shared_dt$class[gprof_shared_dt$intersections == 1] <- "unique"
gprof_shared_dt$class[gprof_shared_dt$intersections == ncol(gprof_res_wide[, -1])] <- "shared"

ggplot(gprof_shared_dt, aes(x = class)) + geom_bar(stat = "count")

# Consensus Method

extract_consensus_DE_candidates <- function(de_res, ain = NULL, norm_thr = 0.8, per_comparison = FALSE){
  # extract only significant changes
  dt <- de_res[de_res$Change != "No Change",]
  if(!is.null(ain)){
    dt <- dt[dt$Assay %in% ain,]
  } else {
    ain <- unique(dt$Assay)
  }
  if(per_comparison){
    dt <- dt[, c("Protein.IDs", "Assay", "Comparison")]
    count_dt <- dt %>% dplyr::group_by(Protein.IDs, Comparison) %>% dplyr::summarise("Intersections" = n())
    count_dt$Percentage <- count_dt$Intersections / length(ain)
    consensus_de <- count_dt[count_dt$Percentage >= norm_thr,]
    consensus_list <- lapply(unique(consensus_de$Comparison), function(x) consensus_de[consensus_de$Comparison == x,]$Protein.IDs)
    names(consensus_list) <- unique(consensus_de$Comparison)
    return(consensus_list)
  } else {
    dt <- dt[, c("Protein.IDs", "Assay")]
    count_dt <- dt %>% dplyr::group_by(Protein.IDs) %>% dplyr::summarise("Intersections" = n())
    count_dt$Percentage <- count_dt$Intersections / length(ain)
    consensus_de <- count_dt[count_dt$Percentage >= norm_thr,]
    return(consensus_de$Protein.IDs)
  }
}

consensus_de <- extract_consensus_DE_candidates(de_res, ain = NULL, norm_thr = 0.8)


# Check standard deviation of all proteins vs. proteins only identified as potential non-DE using Normics

# Extract log2 data
dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[["raw"]])
coldata <- data.table::as.data.table(SummarizedExperiment::colData(se))
rowdata <- data.table::as.data.table(SummarizedExperiment::rowData(se))


# Extract only proteins identified as potential non-DE using Normics
TMT_ratio <- FALSE
top_x <- 50
cols <- c("Protein_ID", "Rank_sum_RS", "Mean_correlation_MC", 
  "Coefficient_of_variation_CV", "Variance_correlation_VC")
dt_matrix <- as.matrix(dt)
means <- apply(dt_matrix, 1, mean, na.rm = TRUE)
variances <- apply(dt_matrix, 1, stats::var, na.rm = TRUE)
longlist <- data.frame(matrix(ncol = length(cols), nrow = 0))
colnames(longlist) <- cols
cor_matrix <- stats::cor(t(dt_matrix), method = "spearman", 
  use = "pairwise.complete.obs")
cor_df <- as.data.frame(cor_matrix)
for (i in 1:nrow(cor_df)) {
  temp <- as.numeric(cor_df[i, -i])
  CV <- ifelse(!TMT_ratio, sqrt(variances[i])/means[i], 
    sqrt(variances[i]))
  longlist <- rbind(longlist, data.frame(Protein_ID = rowdata$Protein.IDs[i], 
    Rank_sum_RS = NA, Mean_correlation_MC = mean(temp, 
      na.rm = TRUE), Coefficient_of_variation_CV = CV, 
    Variance_correlation_VC = stats::var(temp, na.rm = TRUE)))
}
longlist <- longlist[order(longlist[[cols[4]]]), ]
longlist$Rank_sum_RS <- seq_len(nrow(longlist)) - 1
longlist <- longlist[order(-longlist[[cols[3]]]), ]
longlist$Rank_sum_RS <- longlist$Rank_sum_RS + seq_len(nrow(longlist)) - 1
longlist <- longlist[order(longlist[[cols[2]]]), ]
shortlist <- longlist[1:top_x, ]
dt_ids <- cbind(dt, rowdata[, "Protein.IDs"])
colnames(dt_ids)[length(dt) + 1] <- "Protein.IDs"
dt_shortlist <- dt_ids[dt_ids$Protein.IDs %in% shortlist$Protein_ID, 
]
dt_shortlist$Protein.IDs <- NULL


sd_Normics_proteins <- data.frame("Std" = apply(dt_shortlist, 1, function(r) sd(r, na.rm = TRUE)), "Type" = "Normics")
# Per protein -> standard deviation
sd_all_proteins <- data.frame("Std" = apply(dt, 1, function(r) sd(r, na.rm = TRUE)), "Type" = "All")

sd_dt <- rbind(sd_Normics_proteins, sd_all_proteins)

ggplot(sd_dt, aes(x = Std, color = Type)) + geom_density()

