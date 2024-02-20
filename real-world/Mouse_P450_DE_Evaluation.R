# Set data path
data_name <- "mouse_liver_cytochrome_P450"

# Set real-world configuration for this data set
TMT <- FALSE
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

# All normalization methods
se_big_norms <- generate_complete_SE(se_norm, ain =  NULL)
plot_PCA(se_big_norms, color_by = "Normalization", ain = c("all"), label_by = "No", shape_by = "Condition") + theme

# Remove MAD
norm_methods <- names(assays(se_norm))
norm_methods <- norm_methods[norm_methods != "MAD"]
se_big_norms <- generate_complete_SE(se_norm, ain =  norm_methods)
plot_PCA(se_big_norms, color_by = "Normalization", ain = c("all"), label_by = "No", shape_by = "Condition") + theme


# DEP Comparison

plot_overview_DE_bar(de_res)

# Intersection Analysis

plot_upset_DE(de_res)


