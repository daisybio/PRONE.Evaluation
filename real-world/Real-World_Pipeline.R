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

se <- filter_out_complete_NA_proteins(se)
se <- filter_out_proteins_by_value(se, column_name = "Reverse", values = c("+"))
se <- filter_out_proteins_by_value(se, column_name = "Potential.contaminant", values = c("+"))
se <- filter_out_proteins_by_value(se, column_name = "Only.identified.by.site", values = c("+"))

# Overview Plots

#plot_condition_overview(se)

#plot_nr_prot_samples(se)

#plot_tot_int_samples(se)

# NA Quality Control

#get_NA_overview(se)

#plot_NA_frequency(se)

#plot_NA_density(se)


#plot_NA_heatmap(se)

#plot_heatmap(se, ain = "log2", color_by = "Batch")

se <- filter_out_NA_proteins_by_threshold(se, thr = NA_thr)

#plot_NA_heatmap(se)

# Outlier Detection

poma_results <- detect_outliers_POMA(se, ain = "log2")

#poma_results$polygon_plot

#poma_results$distance_boxplot

#poma_results$outliers

#se <- remove_POMA_outliers(se, poma_results$outliers)

saveRDS(se, "data/li_et_al_preprocessed_se.rds")

# Normalization

norm_methods <- get_normalization_methods()

if(!TMT){
  norm_methods <- norm_methods[!norm_methods %in% c("limBE", "IRS")]
}


se <- normalize_se(se, methods = norm_methods, reduce_correlation_by = 1, gamma.0 = 0.1)

se <- vsnNorm(se, ain = "raw", aout = "VSN_0.5", VSN_quantile = 0.5)
se <- normicsNorm(se, ain = "raw", aout = "NormicsVSN_0.5", NormicsVSN_quantile = 0.5)


if(TMT){
  se <- normalize_se_combination(se, methods = c("IRS", "limBE"), ains = norm_methods[!norm_methods %in% c("IRS", "limBE")], combination_pattern = "_on_", reduce_correlation_by = 1, gamma.0 = 0.1)
  se <- irsNorm(se, ain = "VSN_0.5", aout = "IRS_on_VSN_0.5")
  se <- limmaNorm(se, ain = "VSN_0.5", aout = "limBE_on_VSN_0.5")
  se <- irsNorm(se, ain = "NormicsVSN_0.5", aout = "IRS_on_NormicsVSN_0.5")
  se <- limmaNorm(se, ain = "NormicsVSN_0.5", aout = "limBE_on_NormicsVSN_0.5")
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






method <- "NormicsVSN"
on_raw = TRUE
reduce_correlation_by = 1
NormicsVSN_quantile = 0.8
TMT_ratio = FALSE
top_x = 50

stopifnot(method %in% c("NormicsVSN", "NormicsMedian"))
stopifnot(is.numeric(top_x))
dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
stopifnot(top_x <= nrow(dt))
coldata <- data.table::as.data.table(SummarizedExperiment::colData(se))
rowdata <- data.table::as.data.table(SummarizedExperiment::rowData(se))
if (on_raw & ain != "raw") {
  if (ain == "log2") {
    warning("Log2 data specified as ain but on_raw is set to TRUE. On_raw=TRUE forces data to be transformed to raw-scale.")
  }
  dt <- 2^dt
}
if (ain == "raw" & on_raw == FALSE) {
  warning("Raw data specified as ain but on_raw is set to FALSE. On_raw=FALSE forces data to be transformed to log2-scale.")
  dt <- log2(dt)
}

dt_reduced <- dt[seq(1, nrow(dt), by = reduce_correlation_by), ]
rowdata_reduced <- rowdata[seq(1, nrow(rowdata), by = reduce_correlation_by),]
cols <- c("Protein_ID", "Rank_sum_RS", "Mean_correlation_MC", 
          "Coefficient_of_variation_CV", "Variance_correlation_VC")
dt_matrix <- as.matrix(dt_reduced)
means <- apply(dt_matrix, 1, mean, na.rm = TRUE)
variances <- apply(dt_matrix, 1, stats::var, na.rm = TRUE)
longlist <- data.frame(matrix(ncol = length(cols), nrow = 0))
colnames(longlist) <- cols
cor_matrix <- stats::cor(t(dt_matrix), method = "spearman", 
                         use = "pairwise.complete.obs")
cor_df <- as.data.frame(cor_matrix)
for (i in seq_len(nrow(cor_df))) {
  temp <- as.numeric(cor_df[i, -i])
  CV <- ifelse(!TMT_ratio, sqrt(variances[i])/means[i], 
               sqrt(variances[i]))
  longlist <- rbind(longlist, data.frame(Protein_ID = rowdata_reduced$Protein.IDs[i], 
                                         Rank_sum_RS = NA, Mean_correlation_MC = mean(temp, 
                                                                                      na.rm = TRUE), Coefficient_of_variation_CV = CV, 
                                         Variance_correlation_VC = stats::var(temp, na.rm = TRUE)))
}
longlist <- longlist[order(longlist[[cols[4]]]), ]
longlist$Rank_sum_RS <- seq_len(nrow(longlist)) - 1
longlist <- longlist[order(-longlist[[cols[3]]]), ]
longlist$Rank_sum_RS <- longlist$Rank_sum_RS + seq_len(nrow(longlist)) - 
  1
longlist <- longlist[order(longlist[[cols[2]]]), ]
shortlist <- longlist[seq_len(top_x), ]
if (method == "NormicsVSN") {
  dt_ids <- cbind(dt, rowdata[, "Protein.IDs"])
  colnames(dt_ids)[length(dt) + 1] <- "Protein.IDs"
  dt_shortlist <- dt_ids[dt_ids$Protein.IDs %in% shortlist$Protein_ID, 
  ]
  dt_shortlist$Protein.IDs <- NULL
  dt_to_norm <- Biobase::ExpressionSet(assayData = as.matrix(dt_shortlist))
  fit <- vsn::vsn2(dt_to_norm, lts.quantile = NormicsVSN_quantile)
  dt_to_norm <- Biobase::ExpressionSet(assayData = as.matrix(dt))
  norm_data <- vsn::predict(fit, dt_to_norm, log2scale = TRUE)
  norm_dt <- express_to_DT(expr_data = norm_data, column_names = colnames(dt), 
                           row_names = rownames(dt))
} else {
  dt_ids <- cbind(dt, rowdata[["Protein.IDs"]])
  colnames(dt_ids)[length(dt) + 1] <- "Protein.IDs"
  dt_shortlist <- dt_ids[dt_ids$Protein.IDs %in% shortlist$Protein_ID, 
  ]
  dt_shortlist$Protein.IDs <- NULL
  ratios <- as.data.frame(lapply(dt_shortlist, stats::median, 
                                 na.rm = TRUE), col.names = colnames(dt_shortlist))
  dt <- as.data.frame(dt)
  norm_dt <- dt
  for (j in seq_along(dt)) {
    ri <- ratios[1, j]
    norm_dt[, j] <- dt[, j]/ri * mean(as.numeric(ratios[1, 
    ]), na.rm = TRUE)
  }
  if (on_raw) {
    norm_dt <- log2(norm_dt)
  }
  norm_dt <- tib_to_DF(norm_dt, colnames(norm_dt), rownames(norm_dt))
}
SummarizedExperiment::assay(se, aout, FALSE) <- norm_dt


complete_dt <- NULL
for(ain in names(assays(se_norm))){
  dt <- as.data.table(assays(se_norm)[[ain]])
  dt_long <- melt(dt, value.name = "Intensity", variable.name = "Column")
  cd <- as.data.table(colData(se_norm))
  merged_dt <- merge(dt_long, cd, by = "Column")
  merged_dt$assay <- ain
  if(is.null(complete_dt)){
    complete_dt <- merged_dt
  } else {
    complete_dt <- rbind(complete_dt, merged_dt)
  }
}



ggplot(merged_dt, aes(x = Pathological_Status, y = Intensity, fill = Assay)) + geom_boxplot()
