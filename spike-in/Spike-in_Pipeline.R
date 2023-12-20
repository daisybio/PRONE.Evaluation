####### -------- Pipeline for Evaluation of Spike-in Data Sets -------- #######

# Set spike-in configuration
plot <- FALSE
NA_thr <- 0.8
performPOMA <- FALSE
norm_methods <- PRONE.R::get_normalization_methods()
norm_methods <- norm_methods[!norm_methods %in% c("IRS", "HarmonizR", "limBE")]

DE_method <- "limma"
logFC <- FALSE
logFC_up <- 1
logFC_down <- -1
p_adj <- TRUE
alpha <- 0.05

# Load all packages
message("Loading of packages...")
source("config.R")

# Create directories
dir.create("data/raw_spike_in_se/")
dir.create("data/preprocessed_normalized_spike_in_se/")
dir.create("data/de_results_spike_in/")

# Load all data sets in SummarizedExperiment objects
message("Load data sets into SE...")
sapply(list.files("spike-in/data_preparation/", pattern = ".R", full.names = TRUE), source)


# Preprocess and Normalize the individual data sets
message("Preprocess and normalize data sets...")
for(se_rds in list.files(file.path("data/raw_spike_in_se/"), pattern = "_se.rds")){
  message(paste0("...", se_rds))
  # read RDS
  se <- readRDS(file.path("data/raw_spike_in_se/", se_rds))
  # preprocess & normalize SummarizedExperiment
  source("spike-in/Spike-in_Preprocessing.R")
  # save RDS
  se_norm_rds <- paste0(strsplit(se_rds, ".rds")[1][[1]], "_norm.rds")
  saveRDS(se_norm, file.path("data/preprocessed_normalized_spike_in_se/", se_norm_rds))
}

# DE Analysis
message("Run DE analysis...")

for(se_rds in list.files(file.path("data/preprocessed_normalized_spike_in_se/"), pattern = "_se_norm.rds")){
  message(paste0("...", se_rds))
  # read RDS
  se_norm <- readRDS(file.path("data/preprocessed_normalized_spike_in_se/", se_rds))
  # preprocess & normalize SummarizedExperiment
  source("spike-in/Spike-in_DE_Analysis.R")
  # save RDS
  de_results_file <- paste0(strsplit(se_rds, "_se_norm.rds")[1][[1]], "_de_results.rds")
  saveRDS(de_results, file.path("data/de_results_spike_in/", de_results_file))
  stats_file <- paste0(strsplit(se_rds, "_se_norm.rds")[1][[1]], "_de_stats.rds")
  saveRDS(stats, file.path("data/de_results_spike_in/", stats_file))
  auc_file <- paste0(strsplit(se_rds, "_se_norm.rds")[1][[1]], "_de_auc.rds")
  saveRDS(auc, file.path("data/de_results_spike_in/", auc_file))
}
