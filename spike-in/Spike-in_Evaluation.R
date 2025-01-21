####### -------- Script for Evaluation of Performance of Normalization Approaches -------- #######

dataset_naming <- c("dS1", "dS2", "dS3", "dS4", "dS5", "dS6")
names(dataset_naming) <- c("CPTAC6_UPS1_Valikangas", "yeast_UPS1_Ramus", "Ecoli_human_Ionstar", "Ecoli_human_MaxLFQ", "Ecoli_human_DEqMS", "yeast_human_OConnell")

####### -------- Read normalized SE -------- #######
se_norm_list <- list()
for(se_rds in list.files(file.path("data/preprocessed_normalized_spike_in_se/"), pattern = "_se_norm.rds")){
  # read RDS
  se_norm <- readRDS(file.path("data/preprocessed_normalized_spike_in_se/", se_rds))
  data_name <- strsplit(se_rds, "_se_norm.rds")[1][[1]]
  se_norm_list[[data_name]] <- se_norm
}

####### -------- Proportion of spike-in proteins -------- #######
sapply(se_norm_list, function(se){
  spike_col <- metadata(se)$spike_column
  rd <- as.data.table(rowData(se))
  return(table(rd[[spike_col]]))
})

####### -------- Calculate PMAD and Plot -------- #######

PMAD_intra_overall <- NULL
PMAD_diff_overall <- NULL

for(dataset in names(se_norm_list)){
  se <- se_norm_list[[dataset]]
  assays <- SummarizedExperiment::assays(se)
  assays["raw"] <- NULL
  condition <- S4Vectors::metadata(se)$condition
  # PMAD intragroup
  avgmadmem <- calculateAvgMadMem(assays, data.table::as.data.table(SummarizedExperiment::colData(se))[,get(condition)])
  PMAD_intra <- data.table::as.data.table(avgmadmem)
  # set column named by condition variable to rownames of avgmadmem
  PMAD_intra$Condition <- rownames(avgmadmem)
  PMAD_intra <- data.table::melt(PMAD_intra, measure.vars = colnames(PMAD_intra)[colnames(PMAD_intra) != "Condition"], variable.name = "Normalization", value.name = "PMAD")
  PMAD_intra$Normalization <- factor(PMAD_intra$Normalization, levels = sort(as.character(unique(PMAD_intra$Normalization))))
  PMAD_intra$Dataset <- dataset
  if(is.null(PMAD_intra_overall)){
    PMAD_intra_overall <- PMAD_intra
  } else {
    PMAD_intra_overall <- rbind(PMAD_intra_overall, PMAD_intra)
  }
}

# take median over all datasets
PMAD_per_norm <- PMAD_intra_overall %>% dplyr::group_by(Normalization) %>% dplyr::summarize(Median = median(PMAD, na.rm = TRUE)) %>% as.data.table()
# sort norm methods by median PMAD
PMAD_per_norm <- PMAD_per_norm[order(PMAD_per_norm$Median),]

# Calculate median and median absolute deviation over sample groups per normalization method and dataset
PMAD_table <- PMAD_intra_overall %>% dplyr::group_by(Normalization, Dataset, Condition) %>% dplyr::summarize(Median = median(PMAD, na.rm = TRUE), MAD = mad(PMAD, na.rm = TRUE))
PMAD_table$Med_MAD <- paste0(round(PMAD_table$Median, 3), " (", round(PMAD_table$MAD, 3), ")")
PMAD_table$Median <- NULL
PMAD_table$MAD <- NULL
for(dataset in names(dataset_naming)){
  PMAD_table[PMAD_table$Dataset == dataset,]$Dataset <- dataset_naming[dataset]
}
PMAD_table <- dcast(PMAD_table, Normalization~Dataset, value.var = "Med_MAD")

PMAD_table <- PMAD_table[order(match(PMAD_table$Normalization, PMAD_per_norm$Normalization)),]
rownames(PMAD_table) <- PMAD_table$Normalization
PMAD_table$Normalization <- NULL
write.csv(PMAD_table, file = "tables/PMAD_median.csv", col.names = TRUE, row.names = TRUE)


# plot
PMAD_intra_overall$Normalization <- factor(PMAD_intra_overall$Normalization, levels = rev(PMAD_per_norm$Normalization))
PMAD_intra_overall$Pairs <- paste0(PMAD_intra_overall$Dataset, "_", PMAD_intra_overall$Condition)
PMAD_intra_overall[order(PMAD_intra_overall$Pairs, PMAD_intra_overall$Normalization), ]
PMAD_plot <- ggplot(PMAD_intra_overall, aes(y=PMAD, x = Normalization, fill = Normalization)) + theme + 
  geom_boxplot() +
  labs(y = "PMAD", x="Normalization Method") + guides(fill="none") + scale_fill_manual(values = col_vector_norm) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, color = "black") +
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "MAD", paired = TRUE, vjust = 0.7, hjust = 1, label.y = 0.48) +
  #stat_compare_means(label.y = 0.3) +
  coord_flip()
  
####### -------- Calculate Pearson Correlation and Plot -------- #######

cor_intra_overall <- NULL
cor_diff_overall <- NULL
for(dataset in names(se_norm_list)){
  se <- se_norm_list[[dataset]]
  assays <- SummarizedExperiment::assays(se)
  assays["raw"] <- NULL
  condition <- S4Vectors::metadata(se)$condition
  condition_vector <- data.table::as.data.table(SummarizedExperiment::colData(se))[,get(condition)]
  # corr intragroup
  sampleGroupsWithReplicates <- names(table(condition_vector)[table(condition_vector) > 1])
  repCor <- calculateSummarizedCorrelationVector(
    assays,
    condition_vector,
    sampleGroupsWithReplicates,
    "pearson"
  )
  cor_intra <- data.table::as.data.table(repCor)
  cor_intra <- data.table::melt(cor_intra, measure.vars = colnames(cor_intra), variable.name = "Normalization", value.name = "Correlation")
  cor_intra$Normalization <- factor(cor_intra$Normalization, levels = sort(as.character(unique(cor_intra$Normalization))))
  cor_intra$Dataset <- dataset
  if(is.null(cor_intra_overall)){
    cor_intra_overall <- cor_intra
  } else {
    cor_intra_overall <- rbind(cor_intra_overall, cor_intra)
  }
}

# median over all datasets
cor_per_norm <- cor_intra_overall %>% dplyr::group_by(Normalization) %>% dplyr::summarise(Median = median(Correlation, na.rm = TRUE))
cor_per_norm <- cor_per_norm[order(cor_per_norm$Median),]

# median and median absolute deviation for each normalization method and each dataset
cor_intra <- cor_intra_overall %>% dplyr::group_by(Dataset, Normalization) %>% dplyr::summarise(Median = median(Correlation, na.rm = TRUE), MAD = mad(Correlation, na.rm = TRUE))
cor_intra$Correlation <- paste0(round(cor_intra$Median, 3), " (", round(cor_intra$MAD, 3), ")")

for(dataset in names(dataset_naming)){
  cor_intra[cor_intra$Dataset == dataset,]$Dataset <- dataset_naming[dataset]
}

final_cor_intra <- cor_intra[, c("Dataset", "Normalization", "Correlation")]
final_cor_intra <- data.table::dcast(final_cor_intra, Normalization ~ Dataset, value.var = "Correlation")
final_cor_intra <- final_cor_intra[order(match(final_cor_intra$Normalization, cor_per_norm$Normalization), decreasing = TRUE),]
rownames(final_cor_intra) <- final_cor_intra$Normalization
final_cor_intra$Normalization <- NULL
write.csv(final_cor_intra, file = "tables/pearson_correlation_median.csv", col.names = TRUE, row.names = TRUE)

cor_intra_overall$Normalization <- factor(cor_intra_overall$Normalization, levels = rev(PMAD_per_norm$Normalization))

corr_plot <- ggplot(cor_intra_overall, aes(y=Correlation, x = Normalization, fill = Normalization)) + theme + 
  geom_boxplot() +
  labs(y = "Pearson Correlation", x="Normalization Method") + guides(fill="none") + scale_fill_manual(values = col_vector_norm) +
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2, color = "black") +
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "EigenMS", paired = TRUE, vjust = 0.7, hjust = 1, label.y = 1.03) +
  #stat_compare_means(label.y = 0.85) +
  coord_flip()


# Combine PMAD and Correlation plots
PMAD_plot + corr_plot + plot_layout(guides = "collect", axis_titles = "collect", axes = "collect") + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 20))

ggsave("figures/intragroup_variation_overall.png", width = 12, height = 6)
ggsave("figures/paper_figures/Figure3.jpeg", width = 12, height = 6, dpi = 600)

####### -------- Read DE Result Statistics -------- #######

stats_list <- list()
for(stats_file in list.files(file.path("data/de_results_spike_in/"), pattern = "de_stats.rds")){
  # read RDS
  stats <- readRDS(file.path("data/de_results_spike_in/", stats_file))
  data_name <- strsplit(stats_file, "_de_stats.rds")[1][[1]]
  stats$Dataset <- data_name
  stats_list[[data_name]] <- stats
}

stats <- rbindlist(stats_list)
stats$FDR <- stats$FP / (stats$FP + stats$TP)

stats$F1Score <- (2 * stats$TP) / (2 * stats$TP + stats$FP + stats$FN)


####### -------- Read AUCs -------- #######

auc_list <- list()
for(auc_file in list.files(file.path("data/de_results_spike_in/"), pattern = "de_auc.rds")){
  # read RDS
  auc <- readRDS(file.path("data/de_results_spike_in/", auc_file))
  data_name <- strsplit(auc_file, "_de_auc.rds")[1][[1]]
  auc$Dataset <- data_name
  auc_list[[data_name]] <- auc
}

auc <- rbindlist(auc_list)
auc <- auc[, c("Assay", "Dataset", "Comparison", "AUC")]


####### -------- Plot DE Plots -------- #######

# sort methods according to median F1 score
medians <- stats[,c("Assay", "F1Score")] %>% group_by(Assay) %>% dplyr::summarise(Median = median(F1Score, na.rm=TRUE))
medians <- medians[order(medians$Median),]

stats$Assay <- factor(stats$Assay, levels = medians$Assay)
auc$Assay <- factor(auc$Assay, levels = medians$Assay)

for(dataset in names(dataset_naming)){
  stats[stats$Dataset == dataset,]$Dataset <- dataset_naming[dataset]
}



# F1 Score Violin 
dt <- stats[, c("Assay", "F1Score", "Dataset", "Comparison"), with = FALSE]
dt$Assay <- factor(dt$Assay, levels = medians$Assay)

f1_violins <- ggplot(dt, aes (x = Assay, y = F1Score, fill = Assay)) + 
  geom_violin(width = 1, alpha = 0.5) +
  geom_boxplot(width = 0.3) +
  #geom_jitter(aes(shape = Dataset), width = 0.2, alpha = 0.5) +
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2, color = "black", show.legend = FALSE) +
  stat_compare_means(label = "p.format", method = "wilcox.test", ref.group = "RobNorm", paired = TRUE, vjust = 0.7, hjust = 0.5) +
  coord_flip() + 
  scale_fill_manual(name = "Normalization Method", values = col_vector_norm[unique(stats$Assay)]) + 
  labs( x="Normalization Methods", y = "F1 Score") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1)) +
  guides(fill = "none") +
  theme

ggsave("figures/f1_scores_violins.png", width = 12, height = 10)

# F1 Score Boxplot
f1score_final <- ggplot(dt, aes (x = Assay, y = F1Score, fill = Assay)) + 
  geom_boxplot() +
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2, color = "black", show.legend = FALSE) +
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "RobNorm", paired = TRUE, vjust = 0.7, hjust = 0.5) +
  #stat_compare_means(label.y = 0.1, label.x = 18.1) +
  coord_flip() + 
  scale_fill_manual(name = "Normalization Method", values = col_vector_norm[unique(stats$Assay)]) + 
  labs( x="", y = "F1 Score") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1)) +
  theme + guides(fill = "none") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# FPR Boxplot
fpr_box <- ggplot(stats, aes (y = Assay, x = FPR, fill = Assay)) +
  geom_boxplot() + 
  scale_fill_manual(name = "Normalization Method", values = col_vector_norm[unique(stats$Assay)]) + 
  theme + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs( y="", x = "False Positive Rate") + guides(fill = "none") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1))

# AUC Boxplot
auc_box <- ggplot(auc, aes(y = Assay, x = AUC, fill = Assay)) + 
  geom_boxplot() + 
  scale_fill_manual(name = "Normalization Method", values = col_vector_norm[unique(stats$Assay)]) + 
  theme + guides(fill = "none") + 
  labs( y="Normalization Method", x = "AUC") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1))

# F1 Score Heatmap (By Dataset and Normalization Method)
dt <- stats[, c("Assay", "F1Score", "Dataset"), with = FALSE]
median_MAD <- dt %>% dplyr::group_by(Dataset, Assay) %>% dplyr::summarize(Median = median(F1Score, na.rm = TRUE), MAD = mad(F1Score, na.rm = TRUE))
median_MAD$Assay <- factor(median_MAD$Assay, levels = medians$Assay)
median_MAD$Text <- paste0(round(median_MAD$Median, 2), " (", round(median_MAD$MAD, 2), ")")

for(dataset in names(dataset_naming)){
  median_MAD[median_MAD$Dataset == dataset,]$Dataset <- dataset_naming[dataset]
}

heatmap_med <- ggplot(median_MAD, aes(x = Dataset, y = Assay, fill = Median)) + geom_tile(colour="white") + labs(fill = "Median F1Score") + 
  theme + theme(legend.position = "right") + ggplot2::scale_fill_gradient(low = "white", 
                                                                           high = "steelblue")  +
  geom_text(aes(label = Text), size = 5, color = "black") + labs(y = "Normalization Method", x = "Dataset") +
  guides(fill = guide_colorbar(position = "top", barwidth = 10)) +
  theme(legend.title = element_text(margin = margin(-15, 0, 0, 0, "pt")))

# Combine all plots
(auc_box + fpr_box + f1score_final) / heatmap_med + plot_layout(heights = c(1,1)) & plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 22, margin = unit(c(0,0,-0.5,0), "cm")) )
ggsave("figures/spike_in_DE_results.png", width = 14, height = 12)
ggsave("figures/paper_figures/Figure4.jpeg", width = 14, height = 12, dpi = 600)


####### -------- Read DE Results -------- #######

de_list <- list()
for(de_file in list.files(file.path("data/de_results_spike_in/"), pattern = "de_results.rds")){
  # read RDS
  de_res <- readRDS(file.path("data/de_results_spike_in/", de_file))
  data_name <- strsplit(de_file, "_de_results.rds")[1][[1]]
  de_res$Dataset <- data_name
  # add Truth values
  de_res <- calculate_expected_logFC(se_norm_list[[data_name]], de_res)
  de_list[[data_name]] <- de_res
}

de_res <- rbindlist(de_list)
de_res$Difference <- de_res$logFC - de_res$Truth

col_vector_norm <- c(col_vector_norm, "log2" = "yellow")

####### -------- Area under the Precision Recall Curve ------- #######

library(PRROC)

aupr_values <- data.frame("Dataset" = c(), "Assay" = c(), "Comparison" = c(), "AUPR" = c(), "Nr_Spike" = c(), "Nr_Background" = c(), "Nr_All" = c())
for(assay in unique(de_res$Assay)){
  dt <- de_res[de_res$Assay == assay,]
  for(dataset in unique(dt$Dataset)){
    tmp_1 <- dt[dt$Dataset == dataset,]
    for(comp in unique(tmp_1$Comparison)){
      tmp_2 <- tmp_1[tmp_1$Comparison == comp,]
      # Remove NA values
      tmp_2 <- tmp_2[!is.na(tmp_2$adj.P.Val),]
      nr_spike <- nrow(tmp_2[tmp_2$Spiked == "Spike-In",])
      nr_bg <- nrow(tmp_2[tmp_2$Spiked == "BG",])
      pr <- pr.curve(scores.class0 = 1-tmp_2[tmp_2$Spiked == "Spike-In",]$adj.P.Val, 
                     scores.class1 = 1-tmp_2[tmp_2$Spiked == "BG",]$adj.P.Val,
                     curve = TRUE, sorted = FALSE, rand.compute = TRUE)
      r <- data.frame("Dataset" = c(dataset), "Assay" = c(assay), "Comparison" = c(comp),
                      "AUPR" = c(pr$auc.integral), "Nr_Spike" = nr_spike, 
                      "Nr_Background" = nr_bg, "Nr_All" = nr_spike + nr_bg)
      aupr_values <- rbind(aupr_values, r)
    }
    
  }
}

aupr_values$Random <- aupr_values$Nr_Spike / aupr_values$Nr_All
aupr_values$Ratio <- aupr_values$AUPR / aupr_values$Random

median_aupr <- aupr_values %>% group_by(Assay) %>% summarize(Median = median(AUPR))
median_aupr <- median_aupr[order(median_aupr$Median, decreasing = FALSE),]

aupr_values$Assay <- factor(aupr_values$Assay, levels = median_aupr$Assay)

auprc <- ggplot(aupr_values, aes (x = Assay, y = AUPR, fill = Assay)) + 
  geom_boxplot() +
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2, color = "black", show.legend = FALSE) +
  coord_flip() + 
  scale_fill_manual(name = "Normalization Method", values = col_vector_norm[unique(stats$Assay)]) + 
  labs( x="", y = "AUPR Ratio (AUPR / Random)") +  theme +
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "RobNorm", paired = TRUE, vjust = 0.7, hjust = 0.5) +
  stat_compare_means(label.y = 0.1, label.x = 18.1)

mean_std_aupr <- aupr_values %>% group_by(Dataset, Assay) %>% summarize(Mean = mean(AUPR), SD = sd(AUPR))
mean_std_aupr$Text <- paste0(round(mean_std_aupr$Mean, 2), " (", round(mean_std_aupr$SD, 2), ")")
mean_std_aupr$Text[is.na(mean_std_aupr$SD)] <- paste0(round(mean_std_aupr$Mean[is.na(mean_std_aupr$SD)], 2))
mean_std_aupr$Facet <- "Normalization Method"

random_aupr <- unique(aupr_values[, c("Dataset", "Random")])
colnames(random_aupr) <- c("Dataset", "Text")
random_aupr$Assay <- "Random"
random_aupr$Mean <- random_aupr$Text
random_aupr$SD <- 0
random_aupr$Text <- as.character(round(random_aupr$Text, 2))
random_aupr$Facet <- "Random"

mean_std_aupr <- rbind(mean_std_aupr, random_aupr)

for(dataset in names(dataset_naming)){
  mean_std_aupr[mean_std_aupr$Dataset == dataset,]$Dataset <- dataset_naming[dataset]
}

overall_mean <- aupr_values %>% group_by(Assay) %>% summarize(Mean = mean(AUPR), SD = sd(AUPR))
overall_mean <- overall_mean[order(overall_mean$Mean, decreasing = FALSE),]

mean_std_aupr$Assay <- factor(mean_std_aupr$Assay, levels = c(as.character(overall_mean$Assay), "Random"))

ggplot(mean_std_aupr, aes (x = Dataset, y = Assay, fill = Mean, label = Text)) + geom_tile(color = "white") + labs(fill = "Mean AUPR") + facet_grid(Facet~., scales = "free", space = "free") +
  theme + theme(legend.position = "top", strip.text = element_blank(), strip.background = element_blank()) + scale_fill_distiller(palette = "YlOrRd", direction = 1, limits = c(0,1)) + 
  geom_text() + guides(fill = guide_colorbar(barwidth = 10, barheight = 1.5)) + labs(y = "Normalization Method", x = "Data Set")

####### -------- Plot LogFC Observed vs. Expected -------- #######

bg <- ggplot(de_res[de_res$Spiked == "BG",], aes(x = logFC, color = Assay)) + geom_density() + geom_vline(xintercept = 0, linetype="dotted") + 
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) + scale_color_manual(name = "Normalization Method", values = col_vector_norm) +
  theme + labs(x = "LogFC", y = "Density") 

spike <- ggplot(de_res[de_res$Spiked != "BG",], aes( x = Assay, y = Difference, fill = Assay)) + geom_violin(width = 1.0) + 
  geom_boxplot(width = 0.2, outlier.shape = NA) +  stat_boxplot(geom="errorbar", width = 0.2)+ 
  theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + geom_hline(yintercept = 0, linetype="dotted") +
  scale_fill_manual(name = "Normalization Method", values = col_vector_norm) + 
  labs(x = "Normalization Method", y="Observed - Theoretical LogFCs") + ggtitle("Spike-In Proteins")


# helper function
replace_outliers <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  x[(x <= Q1 - 1.5 * IQR) | (x >= Q3 + 1.5 * IQR)] <- NA
  x
}

bg_de_res <- de_res[de_res$Spiked == "BG",]
bg_de_res <- bg_de_res %>% group_by(Assay) %>%
  mutate(logFC_with_NA = replace_outliers(logFC)) 


bg <- ggplot(bg_de_res, aes(x = Assay, y = logFC_with_NA, fill = Assay)) +  geom_violin(width = 1.0, trim = TRUE) + 
  geom_boxplot(aes(y = logFC), width = 0.2, outliers = FALSE) +  stat_boxplot(geom="errorbar", width = 0.2)+ 
  theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + geom_hline(yintercept = 0, linetype="dotted") +
  scale_fill_manual(name = "Normalization Method", values = col_vector_norm) + 
  labs(x = "Normalization Method", y="LogFCs") + ggtitle("Background Proteins")

tmp <- bg_de_res %>% group_by(Assay) %>% summarize(Median = median(logFC, na.rm = TRUE), SD = sd(logFC, na.rm = TRUE)) %>% arrange(Median)
tmp[order(tmp$Median),]

spike / bg + plot_layout(guides = "collect", axis_titles = "collect", axes = "collect") + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

#ggarrange(spike, bg, ncol = 1, labels = c("A", "B"), common.legend = TRUE, legend = "right")
ggsave("figures/expected_observed_logFCs_overall.png", width = 10, height = 10)


