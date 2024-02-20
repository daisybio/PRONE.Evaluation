# Script for the evaluation of NORMICS

se_norm <- readRDS("data/preprocessed_normalized_spike_in_se/Ecoli_human_MaxLFQ_se_norm.rds")

# Add NORMICS normalization methods
se_norm <- normicsNorm(se_norm, ain = "log2", aout = "NormicsVSN_top50", method = "NormicsVSN", reduce_correlation_by = 1, NormicsVSN_quantile = 0.9, top_x = 50)
se_norm <- normicsNorm(se_norm, ain = "log2", aout = "NormicsMedian_top50", method = "NormicsMedian", reduce_correlation_by = 1, NormicsVSN_quantile = 0.9, top_x = 50)

se_norm <- normicsNorm(se_norm, ain = "log2", aout = "NormicsVSN_top200", method = "NormicsVSN", reduce_correlation_by = 1, NormicsVSN_quantile = 0.9, top_x = 200)
se_norm <- normicsNorm(se_norm, ain = "log2", aout = "NormicsMedian_top200", method = "NormicsMedian", reduce_correlation_by = 1, NormicsVSN_quantile = 0.9, top_x = 200)


# PMAD
plot_intragroup_PMAD(se_norm, ain = c("Median", "NormicsVSN_top50", "NormicsMedian_top50","NormicsVSN_top200", "NormicsMedian_top200", "VSN", "log2", "RobNorm"))


# Run DE 
de_res <- run_DE(se_norm,
                 comparisons = specify_comparisons(se_norm))

# Visualize DE Results
stats <- get_spiked_stats_DE(se_norm, de_res)

plot_TP_FP_spiked_bar(stats, ain = NULL, comparisons = NULL)
plot_stats_spiked_heatmap(stats, ain = NULL, metrics = c("Precision", "F1Score"))


