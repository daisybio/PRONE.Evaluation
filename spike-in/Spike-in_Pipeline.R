####### -------- Pipeline for Evaluation of Spike-in Data Sets -------- #######

# Set spike-in configuration
plot <- FALSE
NA_thr <- 0.8
performPOMA <- FALSE

DE_method <- "limma"
logFC <- FALSE
logFC_up <- 1
logFC_down <- -1
p_adj <- TRUE
alpha <- 0.05

# Load all packages
message("Loading of packages...")
source("config.R")

norm_methods <- norm_methods[!norm_methods %in% c("IRS", "limBE")]

# Create directories
dir.create("data/raw_spike_in_se/")
dir.create("data/preprocessed_normalized_spike_in_se/")
dir.create("data/de_results_spike_in/")

# Load all data sets in SummarizedExperiment objects
message("Load data sets into SE...")
sapply(list.files("spike-in/data_preparation/", pattern = ".R", full.names = TRUE), source)
                                      
# Preprocess and normalize the individual data sets
message("Preprocess and normalize data sets...")
overview_list <- list()
for(se_rds in list.files(file.path("data/raw_spike_in_se/"), pattern = "_se.rds")){
  message(paste0("...", se_rds))
  # read RDS
  se <- readRDS(file.path("data/raw_spike_in_se/", se_rds))
  # preprocess & normalize SummarizedExperiment
  source("spike-in/Spike-in_Preprocessing.R")
  # save overview table
  # Create overview table
  overview_table <- data.table("Dataset" = c(se_rds), "Samples" = c(nr_samples),"Conditions" = c(nr_conditions), "Initial" = c(nr_initial), "Prefilter" = c(nr_prefilter), "MV rate" = c(na_percentage), "Final" = c(nr_final))
  overview_list[[se_rds]] <- overview_table
  # save RDS
  se_norm_rds <- paste0(strsplit(se_rds, ".rds")[1][[1]], "_norm.rds")
  saveRDS(se_norm, file.path("data/preprocessed_normalized_spike_in_se/", se_norm_rds))
}

# Construct overview table
overview_table <- rbindlist(overview_list)
overview_table$Dataset <- sapply(strsplit(overview_table$Dataset,"_se.rds"), "[", 1) 
dataset_naming <- c("dS1", "dS2", "dS3", "dS4", "dS5", "dS6", "dS7")
names(dataset_naming) <- c("CPTAC6_UPS1_Valikangas", "yeast_UPS1_Ramus", "yeast_UPS1_Pursiheimo_Valikangas", "Ecoli_human_Ionstar", "Ecoli_human_MaxLFQ", "Ecoli_human_DEqMS", "yeast_human_OConnell")
overview_table <- overview_table[match(names(dataset_naming), overview_table$Dataset),]
overview_table$Dataset <- as.vector(dataset_naming)
write.csv(overview_table, file = "tables/overview_spike_in_table.csv", col.names = TRUE, row.names = FALSE)

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