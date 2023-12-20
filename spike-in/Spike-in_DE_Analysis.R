
# Specify comparisons
comparisons <- specify_comparisons(se_norm, condition = "Condition", sep = NULL, control = NULL)

# Run DE
de_results <- run_DE(se = se_norm, 
                    comparisons = comparisons,
                    ain = NULL, 
                    condition = NULL, 
                    DE_method = DE_method, 
                    covariate = NULL, 
                    logFC = logFC, 
                    logFC_up = logFC_up, 
                    logFC_down = logFC_down, 
                    p_adj = p_adj,
                    alpha = alpha, 
                    B = 100, 
                    K = 500)
                     

stats <- get_spiked_stats_DE(se_norm, de_results)

if(plot){
  plot_TP_FP_spiked_bar(stats, ain = NULL, comparisons = NULL)
  
  plot_stats_spiked_heatmap(stats, ain = NULL, metrics = c("Precision", "F1Score"))
  
  plot_TP_FP_spiked_box(stats)

}

# AUC
tmp <- plot_ROC_AUC_spiked(se_norm, de_results)

if(plot){
  tmp$ROC
  tmp$AUC_box + theme + ggtitle("") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
}

auc <- tmp$AUC_dt

# TODO: add fold changes spiked plot, pvalues spiked plot, logFC_thresholds spiked plot

