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
  PMAD_intra <- data.table::melt(PMAD_intra, measure.vars = colnames(PMAD_intra), variable.name = "Normalization", value.name = "PMAD")
  PMAD_intra$Normalization <- factor(PMAD_intra$Normalization, levels = sort(as.character(unique(PMAD_intra$Normalization))))
  PMAD_intra$Dataset <- dataset
  avgmadmempdiff<- calculatePercentageAvgDiffInMat(avgmadmem)
  PMAD_diff <- data.table(Normalization = names(assays), PMAD = avgmadmempdiff)
  sort_methods <- sort(as.character(unique(PMAD_diff$Normalization)))
  PMAD_diff$Normalization <- factor(PMAD_diff$Normalization, levels=sort_methods)
  PMAD_diff$Dataset <- dataset
  if(is.null(PMAD_intra_overall)){
    PMAD_intra_overall <- PMAD_intra
    PMAD_diff_overall <- PMAD_diff
  } else {
    PMAD_intra_overall <- rbind(PMAD_intra_overall, PMAD_intra)
    PMAD_diff_overall <- rbind(PMAD_diff_overall, PMAD_diff)
  }
}

PMAD_intra <- PMAD_intra_overall %>% dplyr::group_by(Dataset, Normalization) %>% dplyr::summarise(Average = mean(PMAD, na.rm = TRUE), SD = sd(PMAD, na.rm = TRUE))
PMAD_intra$PMAD <- paste0(round(PMAD_intra$Average, 2), " (±", round(PMAD_intra$SD, 2), ")")
PMAD_intra <- PMAD_intra[order(PMAD_intra$Average),]

# table with percent reduction per dataset and per normalization method
PMAD_diff_overall <- PMAD_diff_overall[PMAD_diff_overall$Normalization != "log2",]
perc_reduction <- PMAD_diff_overall
perc_reduction$PMAD <- round(perc_reduction$PMAD, 2)
perc_reduction <- dcast(perc_reduction, Normalization~Dataset, value.var = "PMAD")
perc_reduction <- perc_reduction[, c("Normalization", names(dataset_naming)), with = FALSE]
colnames(perc_reduction) <- c("Normalization", unname(dataset_naming))
write.csv(perc_reduction, file = "tables/percent_reduction_PMAD_individual.csv", col.names = TRUE, row.names = TRUE)


final_PMAD_intra <- PMAD_intra[, c("Dataset", "Normalization", "PMAD")]
final_PMAD_intra <- data.table::dcast(final_PMAD_intra, Normalization ~ Dataset, value.var = "PMAD")
final_PMAD_intra <- final_PMAD_intra[final_PMAD_intra$Normalization != "log2",]

PMAD_diff_ranks <- PMAD_diff_overall %>% dplyr::group_by(Dataset) %>% dplyr::mutate(rank = dplyr::dense_rank(PMAD))
final_PMAD_diff <- PMAD_diff_ranks %>% dplyr::group_by(Normalization) %>% dplyr::summarise(Average = mean(rank, na.rm = TRUE), SD = sd(rank, na.rm = TRUE))
final_PMAD_diff$`Mean of Ranks` <- paste0(round(final_PMAD_diff$Average, 2), " (±", round(final_PMAD_diff$SD, 2), ")")
final_PMAD_diff <- final_PMAD_diff[order(-final_PMAD_diff$Average),]
final_PMAD_diff$Normalization <- factor(final_PMAD_diff$Normalization, levels = final_PMAD_diff$Normalization)

PMAD_table <- merge(final_PMAD_intra, final_PMAD_diff, by = "Normalization")
PMAD_table <- PMAD_table[order(PMAD_table$Average),]
PMAD_table$Average <- NULL
PMAD_table$SD <- NULL

pPMAD <- ggplot(final_PMAD_diff, aes(x = Average, y = Normalization)) + geom_point(size = 5) + geom_errorbar(aes(xmin = Average - SD, xmax = Average + SD), width = .4, position = position_dodge(.9)) +
  theme + labs(x = "Mean Rank of Percent Reduction in PMAD", y="Normalization Method") 
ggsave("figures/PMAD_overall_ranks.png", width = 12, height = 8)


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

cor_intra <- cor_intra_overall %>% dplyr::group_by(Dataset, Normalization) %>% dplyr::summarise(Average = mean(Correlation, na.rm = TRUE), SD = sd(Correlation, na.rm = TRUE))
cor_intra$Correlation <- paste0(round(cor_intra$Average, 2), " (±", round(cor_intra$SD, 2), ")")
cor_intra <- cor_intra[order(cor_intra$Average),]

final_cor_intra <- cor_intra[, c("Dataset", "Normalization", "Correlation")]
final_cor_intra <- data.table::dcast(final_cor_intra, Normalization ~ Dataset, value.var = "Correlation")
final_cor_intra <- final_cor_intra[final_cor_intra$Normalization != "log2",]

cor_overall <- cor_intra_overall %>% dplyr::group_by(Normalization) %>% dplyr::summarise(Average = mean(Correlation, na.rm = TRUE), SD = sd(Correlation, na.rm = TRUE), Median = median(Correlation, na.rm = TRUE))
cor_overall <- cor_overall[order(cor_overall$Median),]
methods <- cor_overall$Normalization
ggplot(cor_overall, aes(x = Average, y = Normalization)) + geom_point(size = 5) + geom_errorbar(aes(xmin = Average - SD, xmax = Average + SD), width = .2, position = position_dodge(.9)) +
  theme + labs(x = "Mean Pearson Correlation Coefficient", y = "Normalization Method") 

cor_intra_overall <- cor_intra_overall[cor_intra_overall$Normalization != "log2",]
cor_intra_overall$Normalization <- factor(cor_intra_overall$Normalization, levels = final_PMAD_diff$Normalization)
pcor <- ggplot(cor_intra_overall, aes( x= Correlation, y=Normalization, fill = Normalization)) + 
  geom_boxplot() + theme + labs(x = "Pearson Correlation Coefficient", y="Normalization Method") + guides(fill="none") + scale_fill_manual(values = col_vector_norm)
ggsave("figures/pearson_correlation_overall.png", width = 12, height = 8)

ggpubr::ggarrange(pPMAD, pcor, ncol = 2, labels = c("A", "B"))
ggsave("figures/intragroup_variation_overall.png", width = 12, height = 4)


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
auc <- auc[, c("Assay", "Dataset", "AUC")]


####### -------- Plot DE Plots -------- #######

# sort methods according to median F1 score
medians <- stats[,c("Assay", "F1Score")] %>% group_by(Assay) %>% dplyr::summarise(Median = median(F1Score, na.rm=TRUE))
medians <- medians[order(medians$Median),]

stats$Assay <- factor(stats$Assay, levels = medians$Assay)
auc$Assay <- factor(auc$Assay, levels = medians$Assay)

f1_box <- ggplot(stats[stats$Assay != "log2"], aes (y = Assay, x = F1Score, fill = Assay)) + 
  geom_boxplot() + 
  scale_fill_manual(name = "Normalization", values = col_vector_norm[unique(stats$Assay)]) + 
  theme + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs( y="", x = "F1 Score") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1))

fpr_box <- ggplot(stats[stats$Assay != "log2"], aes (y = Assay, x = FPR, fill = Assay)) +
  geom_boxplot() + 
  scale_fill_manual(name = "Normalization", values = col_vector_norm[unique(stats$Assay)]) + 
  theme + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs( y="", x = "False Positive Rate") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1))

auc_box <- ggplot(auc[auc$Assay != "log2",], aes(y = Assay, x = AUC, fill = Assay)) + 
  geom_boxplot() + 
  scale_fill_manual(name = "Normalization", values = col_vector_norm[unique(stats$Assay)]) + 
  theme + 
  labs( y="Normalization Method", x = "AUC") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1))

t <- ggarrange(auc_box, fpr_box, f1_box, ncol = 3, common.legend = TRUE, legend = "none", labels = c("A", "B", "C"), font.label = list(size = 18), vjust = 0.01, hjust = 0.01, widths = c(0.65, 0.5, 0.5))
ggarrange(NULL,t, heights = c(0.1, 1), ncol = 1)
ggsave("figures/DE_panel_spike_in.png", width = 12, height = 6)


# F1 Score Heatmap (By Dataset and Normalization Method)
dt <- stats[, c("Assay", "F1Score", "Dataset"), with = FALSE]
averages <- dt %>% dplyr::group_by(Dataset, Assay) %>% dplyr::summarize(Average = mean(F1Score, na.rm = TRUE), SD = sd(F1Score, na.rm = TRUE))
averages$Assay <- factor(averages$Assay, levels = medians$Assay)

for(dataset in names(dataset_naming)){
  averages[averages$Dataset == dataset,]$Dataset <- dataset_naming[dataset]
}

ggplot(averages[averages$Assay != "log2",], aes(x = Dataset, y = Assay, fill = Average)) + geom_tile(colour="white") + labs(fill = "Mean F1Score") + 
  theme + theme(legend.position = "top") + scale_fill_distiller(palette = "YlOrRd", direction = 1, limits = c(0,1)) + 
  geom_text(aes(label = round(Average, 2))) + guides(fill = guide_colorbar(barwidth = 10, barheight = 1.5)) + labs(y = "Normalization Method", x = "Data Set")
ggsave("figures/DE_F1score_heatmap_spike_in.png", width = 8, height = 6)


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

####### -------- Plot LogFC Observed vs. Expected -------- #######

bg <- ggplot(de_res[de_res$Spiked == "BG",], aes(x = logFC, color = Assay)) + geom_density() + geom_vline(xintercept = 0, linetype="dotted") + 
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) + scale_color_manual(name = "Normalization Method", values = col_vector_norm) +
  theme + labs(x = "LogFC", y = "Density") 

spike <- ggplot(de_res[de_res$Spiked != "BG",], aes( x = Assay, y = Difference, fill = Assay)) + geom_violin(width = 1.0) + 
  geom_boxplot(width = 0.2, outlier.shape = NA) +  stat_boxplot(geom="errorbar", width = 0.2)+ 
  theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + geom_hline(yintercept = 0, linetype="dotted") +
  scale_fill_manual(name = "Normalization Method", values = col_vector_norm) + 
  labs(x = "Normalization Method", y="Observed - Theoretical LogFCs") 

ggarrange(spike, bg, ncol = 1, labels = c("A", "B"), common.legend = TRUE, legend = "right")
ggsave("figures/expected_observed_logFCs_overall.png", width = 8, height = 10)


