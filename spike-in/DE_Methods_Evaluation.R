# Script for evaluating the performance of the differential expression methods (limma vs. ROTS)

library(data.table)
library(PRONE.R)
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
  stats <- PRONE.R::get_spiked_stats_DE(se_norm, de_res)
  stats$Dataset <- data_name
  de_list_ROTS[[data_name]] <- de_res
  stats_list_ROTS[[data_name]] <- stats
}

de_res_ROTS <- rbindlist(de_list_ROTS)
stats_ROTS <- rbindlist(stats_list_ROTS)
stats_ROTS$Method <- "ROTS"


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

# Combine ROTS and limma stats
stats <- rbind(stats_ROTS, stats_limma)

# Plot the results for a single data set
dt <- data.table::as.data.table(stats[, c("Method", "Assay", "TP", "FP", "Dataset")])
dt <- melt(dt, id.vars = c("Method", "Assay", "Dataset"), value.name = "Count", variable.name = "Type")

tps <- ggplot(dt[dt$Type == "TP",], aes(x = Assay, y = Count, fill = Method)) +
  geom_boxplot() + facet_wrap(~Dataset, scales = "free_y", ncol = 1) + theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ggtitle("True Positives")

fps <- ggplot(dt[dt$Type == "FP",], aes(x = Assay, y = Count, fill = Method)) +
  geom_boxplot() + facet_wrap(~Dataset, scales = "free_y", ncol = 1) + theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ggtitle("False Positives")

library(patchwork)
tps + fps + plot_layout(guides = "collect")
ggsave("figures/DE_Methods_Evaluation.png", width = 10, height = 12, dpi = 300)

# Compare mean F1 score for limma and ROTS across different data sets

f1_df <- stats[, c("Method", "Dataset", "Assay", "F1Score")]
f1_df <- f1_df %>% group_by(Method, Dataset, Assay) %>% summarize(mean_f1 = mean(F1Score, na.rm = TRUE))

ggplot(f1_df, aes(x = Assay, y = mean_f1, fill = Method)) + 
  geom_bar(stat="identity", position = position_dodge()) + 
  facet_wrap(~Dataset, scales = "free_y", ncol = 2) +
  theme +
  coord_flip()
ggsave("figures/DE_Methods_Evaluation_Barplot.png", width = 10, height = 12, dpi = 300)


# Compare ranking of methods

# sort methods according to median F1 score
medians <- stats_ROTS[,c("Assay", "Method", "F1Score")] %>% group_by(Assay, Method) %>% dplyr::summarise(Median = median(F1Score, na.rm=TRUE))
medians <- medians[order(medians$Median),]

stats_ROTS$Assay <- factor(stats_ROTS$Assay, levels = medians$Assay)

f1_box <- ggplot(stats_ROTS[stats_ROTS$Assay != "log2"], aes (y = Assay, x = F1Score, fill = Assay)) + 
  geom_boxplot() + 
  scale_fill_manual(name = "Normalization", values = col_vector_norm[unique(stats_ROTS$Assay)]) + 
  theme + theme(legend.position = "none") +
  labs( y="Normalization Method", x = "F1 Score") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1))

# Ecoli MaxLFQ
dt <- data.table::as.data.table(stats[, c("Method", "Assay", "TP", "FP", "Dataset")])
dt <- dt[dt$Dataset == "Ecoli_human_MaxLFQ",]
dt$Dataset <- NULL
dt <- melt(dt, id.vars = c("Method", "Assay"), value.name = "Count", variable.name = "Type")
dt$Assay <- factor(dt$Assay, levels = medians$Assay)

tps <- ggplot(dt[dt$Type == "TP",], aes(x = Assay, y = Count, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) + theme + labs(y = "True Positives", x = "Normalization Method", fill = "DE Method") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) + coord_flip()

fps <- ggplot(dt[dt$Type == "FP",], aes(x = Assay, y = Count, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) + theme + labs(y = "False Positives", x = "Normalization Method", fill = "DE Method") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) + coord_flip()

(tps + fps) / f1_box + plot_layout(guides = "collect")


# tps of limma / tps of ROTS for each assay and each data set
dt_rots <- stats_ROTS[, c("Assay","Comparison", "Dataset", "TP", "FP")]
colnames(dt_rots) <- c("Assay","Comparison", "Dataset", "TP_ROTS", "FP_ROTS")

dt_limma <- stats_limma[, c("Assay", "Comparison","Dataset" , "TP", "FP")]
colnames(dt_limma) <- c("Assay", "Comparison", "Dataset", "TP_limma", "FP_limma")

dt <- merge(dt_rots, dt_limma, by = c("Assay", "Comparison", "Dataset"))

dt <- dt[dt$Assay != "log2",]

dt$DivTP <-  dt$TP_limma / dt$TP_ROTS
dt$DivFP <- dt$FP_limma / dt$FP_ROTS

# Division
dt <- dt[, c("Assay", "Dataset", "DivTP", "DivFP")]
colnames(dt) <- c("Assay", "Dataset", "TP", "FP")
dt <- melt(dt, value.name = "Value", variable.name = "Class", measure.vars = c("TP", "FP"))

dt$Assay <- factor(dt$Assay, levels = medians$Assay)
tps_fps_div <- ggplot(dt, aes(x = Value, y = Assay, fill = Class)) + geom_boxplot() + 
  theme + scale_x_log10() +
  scale_fill_manual(values = c("#D55E00", "#0072B2")) + labs(x = "limma DEP Count / ROTS DEP Count (log10)", y = "Normalization Method")

ggarrange(tps_fps_div, f1_box, labels = c("A", "B"))
ggsave("figures/DE_ROTS_spike_in.png", width = 12, height = 6)
