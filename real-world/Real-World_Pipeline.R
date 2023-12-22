####### -------- Pipeline for Evaluation of a Real-World Data Set -------- #######

# Load all packages
message("Loading of packages...")
source("config.R")

# Set real-world configuration

data_file <- "data/raw_real_world_se/tuberculosis_se.rds"

NA_thr <- 0.8

DE_method <- "limma"
logFC <- FALSE
logFC_up <- 1
logFC_down <- -1
p_adj <- TRUE
alpha <- 0.05

# Load data into SE

se <- readRDS(data_file)

dim(se)

# Pre-filtering

se <- filter_complete_NA_proteins(se)
se <- filter_proteins_by_value(se, column_name = "Reverse", values = c("+"))
se <- filter_proteins_by_value(se, column_name = "Potential.contaminant", values = c("+"))
se <- filter_proteins_by_value(se, column_name = "Only.identified.by.site", values = c("+"))

# Overview Plots

plot_condition_overview(se)

plot_nr_prot_samples(se)

plot_tot_int_samples(se)

# NA Quality Control

get_NA_overview(se)

plot_NA_frequency(se)

plot_NA_density(se)

plot_NA_heatmap(se)

se <- filter_NA_proteins_by_threshold(se, thr = NA_thr)

plot_NA_heatmap(se)

# Outlier Detection

poma_results <- detect_outliers_POMA(se, ain = "raw")

poma_results$polygon_plot

poma_results$distance_boxplot

poma_results$outliers

se <- remove_POMA_outliers(se, poma_results$outliers)

# Normalization

norm_methods <- get_normalization_methods()

se <- normalize_se(se, methods = norm_methods)

# Add combinations of methods

# TODO

# Evaluation

# Visual quality control
plot_boxplots(se)
plot_PCA(se, shape_by = "Batch")
plot_densities(se, color_by = "Column")


# Intragroup variation
plot_intragroup_PMAD(se)
plot_intragroup_correlation(se)
plot_intragroup_PEV(se)

# Differential expression analysis
se <- remove_reference_samples(se)

comparisons <- specify_comparisons(se)
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

