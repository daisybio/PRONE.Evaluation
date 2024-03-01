####### -------- Pipeline for Evaluation of a Real-World Data Set -------- #######

# Load all packages
message("Loading of packages...")
source("config.R")

# Create directories
dir.create("data/raw_real_world_se/")
dir.create("data/preprocessed_normalized_real_world_se/")
dir.create("data/de_results_real_world/")

# Load data set in SummarizedExperiment objects if this has not been done before
data_file <- paste0("data/raw_real_world_se/", data_name, "_se.rds")
if(!file.exists(data_file)){
  message("Load data sets into SE...")
  sapply(list.files("real-world/data_preparation/", pattern = ".R", full.names = TRUE), source)
}

# Load data into SE
se <- readRDS(data_file)

dim(se)

# Pre-filtering

se <- filter_complete_NA_proteins(se)
se <- filter_proteins_by_value(se, column_name = "Reverse", values = c("+"))
se <- filter_proteins_by_value(se, column_name = "Potential.contaminant", values = c("+"))
se <- filter_proteins_by_value(se, column_name = "Only.identified.by.site", values = c("+"))

# Overview Plots

#plot_condition_overview(se)

#plot_nr_prot_samples(se)

#plot_tot_int_samples(se)

# NA Quality Control

#get_NA_overview(se)

#plot_NA_frequency(se)

#plot_NA_density(se)

#plot_NA_heatmap(se)

se <- filter_NA_proteins_by_threshold(se, thr = NA_thr)

#plot_NA_heatmap(se)

# Outlier Detection

poma_results <- detect_outliers_POMA(se, ain = "raw")

#poma_results$polygon_plot

#poma_results$distance_boxplot

#poma_results$outliers

se <- remove_POMA_outliers(se, poma_results$outliers)

#saveRDS(se, "data/li_et_al_preprocessed_se.rds")

# Normalization

norm_methods <- get_normalization_methods()
norm_methods <- norm_methods[norm_methods != "HarmonizR"] # this does not work yet!

if(!TMT){
  norm_methods <- norm_methods[!norm_methods %in% c("limBE", "IRS")]
}


se <- normalize_se(se, methods = norm_methods, reduce_correlation_by = 1, gamma.0 = 0.1)

if(TMT){
  se <- normalize_se_combination(se, methods = norm_methods[!norm_methods %in% c("IRS", "limBE")], ains = c("IRS", "limBE"), "_on_", reduce_correlation_by = 1, gamma.0 = 0.1)
}

# Evaluation

# Visual quality control
#plot_boxplots(se)


#plot_PCA(se, color_by = "Batch", shape_by = "Condition", ncol = 4, label_by = "No")

#plot_densities(se, color_by = "Column")

# Intragroup variation
#plot_intragroup_PMAD(se)
#plot_intragroup_correlation(se)
#plot_intragroup_PEV(se)

# Differential expression analysis
if(TMT){
  se <- remove_reference_samples(se)
}

out_file <- paste0("data/preprocessed_normalized_real_world_se/", data_name, "_se_norm.rds")
saveRDS(se, file = out_file)

if(is.null(comparisons)){
  comparisons <- specify_comparisons(se)
}

de_res <- run_DE(se, 
                 comparisons = comparisons,
                 ain = NULL,
                 condition = NULL,
                 DE_method = DE_method,
                 covariate = NULL,
                 logFC = logFC,
                 logFC_up = logFC_up,
                 logFC_down = logFC_down,
                 p_adj = p_adj,
                 alpha = alpha
                 )

out_file <- paste0("data/de_results_real_world/", data_name, "_de_res.rds")
saveRDS(de_res, file = out_file)

