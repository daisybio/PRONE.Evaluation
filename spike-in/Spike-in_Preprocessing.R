####### -------- Preprocessing and Normalization Script -------- #######
nr_initial <- nrow(se)
nr_samples <- ncol(se)
nr_conditions <- length(unique(as.data.table(colData(se))[[metadata(se)$condition]]))

# Remove proteins with missing values in all samples
se <- filter_out_complete_NA_proteins(se)

# Remove reverse hits, contamininants, and only identified by site
se <- filter_out_proteins_by_value(se, "Reverse", "+")
se <- filter_out_proteins_by_value(se, "Potential.contaminant", "+")
se <- filter_out_proteins_by_value(se, "Only.identified.by.site", "+")

nas <- get_NA_overview(se, ain = "log2")

# Overview Plots
if(plot){
  print(nas)
  
  plot_condition_overview(se, condition = NULL)
  
  plot_identified_spiked_proteins(se, color_by = "Condition")
  
  plot_histogram_spiked(se, condition = "Condition")
  
  plot_profiles_spiked(se)
}

nr_prefilter <- nrow(se)
na_percentage <- paste0(round(nas$NA.Percentage, 2), " %")

# Filter Missing Values

if(plot){
  plot_NA_heatmap(se, color_by = NULL, label_by = NULL, cluster_samples = TRUE, cluster_proteins = TRUE)
  
  plot_NA_density(se)
  
  plot_NA_frequency(se)
}

se <- filter_out_NA_proteins_by_threshold(se, thr = NA_thr)

nr_final <- nrow(se)

if(plot){
  plot_NA_heatmap(se, color_by = NULL, label_by = NULL, cluster_samples = TRUE, cluster_proteins = TRUE)
}

# Quality Control of Samples

if(plot){
  plot_nr_prot_samples(se, color_by = NULL, label_by = NULL)
  
  plot_tot_int_samples(se, color_by = NULL, label_by = NULL)
}

if(performPOMA){
  
  poma_res <- detect_outliers_POMA(se, ain = "log2")
  
  poma_res$polygon_plot
  
  poma_res$distance_boxplot
  
  poma_res$outliers
}

if(removePOMAoutliers){
  se <- remove_outliers_POMA(se, ain = "log2")
}


# Normalization
se_norm <- normalize_se(se, norm_methods, combination_pattern = NULL, gamma.0 = 0.1)

if(plot){
  plot_boxplots(se_norm, ain = NULL, color_by = NULL, label_by = NULL, ncol = 3, facet_norm = TRUE)
  plot_densities(se_norm, ain = NULL, color_by = NULL, facet_norm = TRUE)
  plot_PCA(se_norm, ain = NULL, color_by = NULL, label_by = "No", shape_by = NULL, facet_norm = TRUE, facet_by = NULL)
  plot_intragroup_PMAD(se_norm, ain = NULL, condition = NULL, diff = TRUE)
}


