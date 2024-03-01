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
medians <- stats[,c("Assay", "F1Score", "Method")] %>% group_by(Assay, Method) %>% dplyr::summarise(Median = median(F1Score, na.rm=TRUE))

medians_limma <- medians[medians$Method == "limma",]
medians_ROTS <- medians[medians$Method == "ROTS",]

medians_limma <- medians_limma[order(medians_limma$Median),]
medians_ROTS <- medians_ROTS[order(medians_ROTS$Median),]
