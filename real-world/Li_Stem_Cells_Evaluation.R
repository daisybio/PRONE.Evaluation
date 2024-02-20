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

# Quantitative and Qualitative Evaluation

plot_intragroup_PMAD(se_norm) + theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "None")


# Complete PCA plot of all normalization methods

generate_complete_SE <- function(se, ain = NULL){
  # create one big data table
  big_dt <- NULL
  big_cd <- NULL
  # check if assays all in se
  ain <- PRONE.R::check_input_assays(se, ain)
  for(a in ain){
    if(!a %in% c("raw")){
      dt <- as.data.table(assays(se)[[a]])
      colnames(dt) <- paste0(a, "_", colnames(dt))
      cd <- as.data.table(colData(se))
      cd$Column <- paste0(a,"_",cd$Column)
      cd$Normalization <- a
      if(is.null(big_dt)){
        big_dt <- dt
        big_cd <- cd
      } else {
        big_dt <- cbind(big_dt, dt)
        big_cd <- rbind(big_cd, cd)
      }
    }
  }
  se_big <- SummarizedExperiment::SummarizedExperiment(assays = list(all = big_dt), colData = big_cd, rowData = as.data.table(rowData(se)))
  metadata(se_big) <- metadata(se)
  return(se_big)
}

se_big_norms <- generate_complete_SE(se_norm, ain = c("EigenMS", "GlobalMean", "GlobalMedian", "LoessCyc", "RobNorm", "LoessF", "MAD", "Mean", "Median", "NormicsMedian", "NormicsVSN", "Quantile", "RLR", "RlrMA", "RlrMACyc", "TMM", "VSN"))
plot_PCA(se_big_norms, color_by = "Normalization", ain = c("all"), label_by = "No", shape_by = "Batch") + theme

se_big_batch <- generate_complete_SE(se_norm, ain = c("IRS", "limBE"))
plot_PCA(se_big_batch, color_by = "Normalization", ain = c("all"), label_by = "No", shape_by = "Batch") + theme

se_big_batch_norms <- generate_complete_SE(se_norm, ain = c("Median_on_IRS", "Median_on_limBE", "RobNorm_on_limBE", "RobNorm_on_IRS", "IRS_on_GlobalMean_on_TMM"))
plot_PCA(se_big_batch_norms, color_by = "Normalization", ain = c("all"), label_by = "No", shape_by = "Batch") + theme


# DEP Comparison

plot_overview_DE_bar(de_res)
plot_overview_DE_tile(de_res)

# Intersection Analysis

plot_upset_DE(de_res)

