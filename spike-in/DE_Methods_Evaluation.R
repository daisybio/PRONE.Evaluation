# Script for evaluating the performance of the differential expression methods (limma vs. ROTS)

library(data.table)
library(PRONE)
library(ggplot2)
library(dplyr)

# Load ROTS results
de_list_ROTS <- list()
stats_list_ROTS <- list()
for(de_file in list.files(file.path("data/de_results_spike_in_ROTS/"), pattern = "de_results_ROTS.rds")){
  # read RDS
  de_res <- readRDS(file.path("data/de_results_spike_in_ROTS/", de_file))
  data_name <- strsplit(de_file, "_de_results_ROTS.rds")[1][[1]]
  de_res$Dataset <- data_name
  se_norm <- readRDS(file.path("data/preprocessed_normalized_spike_in_se/", paste0(data_name, "_se_norm.rds")))
  # calculate stats
  stats <- PRONE::get_spiked_stats_DE(se_norm, de_res)
  stats$Dataset <- data_name
  de_list_ROTS[[data_name]] <- de_res
  stats_list_ROTS[[data_name]] <- stats
}

de_res_ROTS <- rbindlist(de_list_ROTS)
stats_ROTS <- rbindlist(stats_list_ROTS)
stats_ROTS$Method <- "ROTS"
stats_ROTS$F1Score <- (2 * stats_ROTS$TP) / (2 * stats_ROTS$TP + stats_ROTS$FP + stats_ROTS$FN)


# Load limma results

stats_list_limma <- list()
for(stats_file in list.files(file.path("data/de_results_spike_in/"), pattern = "de_stats.rds")){
  # read RDS
  stats <- readRDS(file.path("data/de_results_spike_in/", stats_file))
  data_name <- strsplit(stats_file, "_de_stats.rds")[1][[1]]
  stats$Dataset <- data_name
  stats_list_limma[[data_name]] <- stats
}

stats_limma <- rbindlist(stats_list_limma)
stats_limma$Method <- "limma"
stats_limma$F1Score <- (2 * stats_limma$TP) / (2 * stats_limma$TP + stats_limma$FP + stats_limma$FN)


# Combine ROTS and limma stats
stats <- rbind(stats_ROTS, stats_limma)

# Plot the results for a single data set
dt <- data.table::as.data.table(stats[, c("Method", "Assay", "TP", "FP", "Dataset")])
dt <- melt(dt, id.vars = c("Method", "Assay", "Dataset"), value.name = "Count", variable.name = "Type")

# sort methods according to median F1 score
medians <- stats_ROTS[,c("Assay", "Method", "F1Score")] %>% group_by(Assay, Method) %>% dplyr::summarise(Median = median(F1Score, na.rm=TRUE))
medians <- medians[order(medians$Median),]

stats_ROTS$Assay <- factor(stats_ROTS$Assay, levels = medians$Assay)

f1score_ROTS <- ggplot(stats_ROTS, aes (x = Assay, y = F1Score, fill = Assay)) + 
  geom_boxplot() +
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2, color = "black", show.legend = FALSE) +
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "LoessF", paired = TRUE, vjust = 0.7, hjust = 0.5) +
  #stat_compare_means(label.y = 0.1, label.x = 18.1) +
  coord_flip() + 
  scale_fill_manual(name = "Normalization Method", values = col_vector_norm[unique(stats$Assay)]) + 
  labs( x="", y = "F1 Score of ROTS") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1)) +
  theme +
  guides(fill = "none")

# tps of limma / tps of ROTS for each assay and each data set
dt_rots <- stats_ROTS[, c("Assay","Comparison", "Dataset", "TP", "FP")]
colnames(dt_rots) <- c("Assay","Comparison", "Dataset", "TP_ROTS", "FP_ROTS")

dt_limma <- stats_limma[, c("Assay", "Comparison","Dataset" , "TP", "FP")]
colnames(dt_limma) <- c("Assay", "Comparison", "Dataset", "TP_limma", "FP_limma")

dt <- merge(dt_rots, dt_limma, by = c("Assay", "Comparison", "Dataset"))

dt$DivTP <-  dt$TP_limma / dt$TP_ROTS
dt$DivFP <- dt$FP_limma / dt$FP_ROTS

# Division
dt <- dt[, c("Assay", "Dataset", "DivTP", "DivFP")]
colnames(dt) <- c("Assay", "Dataset", "TP", "FP")
dt <- melt(dt, value.name = "Value", variable.name = "Class", measure.vars = c("TP", "FP"))

dt$Assay <- factor(dt$Assay, levels = medians$Assay)
tps_fps_div <- ggplot(dt, aes(x = Value, y = Assay, fill = Class)) + geom_boxplot() + 
  theme + scale_x_log10() + guides(fill = guide_legend(position = "top")) +
  scale_fill_manual(values = c("#D55E00", "#0072B2")) + labs(x = "limma DEP Count / ROTS DEP Count \n (log10)", y = "Normalization Method")

# Paired test of ROTS and limma F1 Scores

f1_ROTS <- stats_ROTS[, c("Assay", "Comparison", "Dataset", "F1Score")]
colnames(f1_ROTS) <- c("Assay", "Comparison", "Dataset", "ROTS")

f1_limma <- stats_limma[, c("Assay", "Comparison", "Dataset", "F1Score")]
colnames(f1_limma) <- c("Assay", "Comparison", "Dataset", "limma")

# Merge f1 scores based on assay, comparison and dataset
f1_both <- merge(f1_ROTS, f1_limma, by = c("Assay", "Comparison", "Dataset"))

# Wilcox paired test
f1_wilcox <- sapply(unique(f1_both$Assay), function(x){
  return(wilcox.test(f1_both[f1_both$Assay == x,]$ROTS, f1_both[f1_both$Assay == x,]$limma, paired = TRUE)$p.value)
})

f1_wilcox <- as.data.frame(f1_wilcox)
f1_wilcox$Assay <- rownames(f1_wilcox)

f1_wilcox$Significance <- ifelse(f1_wilcox$f1_wilcox <= 0.05, "*", "")
f1_wilcox[f1_wilcox$f1_wilcox <= 0.01,]$Significance <- "**"
f1_wilcox[f1_wilcox$f1_wilcox <= 0.001,]$Significance <- "***"
#f1_wilcox[f1_wilcox$f1_wilcox <= 0.0001,]$Significance <- "****"

f1_wilcox$Assay <- factor(f1_wilcox$Assay, levels = medians$Assay)

f1_wilcox_plot <- ggplot(f1_wilcox, aes(x = f1_wilcox, y = Assay, fill = Assay)) + 
  geom_bar(stat = "identity", color = "black") + 
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") + 
  theme + labs(x = "P-Value of Paired Wilcox-Test - \n F1 Score ROTS vs. Limma", y = "Normalization Method") + 
  scale_fill_manual(values = col_vector_norm) + 
  geom_text(aes(label = Significance), size = 6, nudge_x = 0.006, nudge_y = -0.15) +
  guides(fill = "none")


(tps_fps_div + f1score_ROTS + f1_wilcox_plot) + plot_layout(axes = "collect") + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18))

ggsave("figures/DE_ROTS_spike_in.png", width = 12, height = 6)


 