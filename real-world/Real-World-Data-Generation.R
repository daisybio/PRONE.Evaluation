# Load all packages
message("Loading of packages...")
source("config.R")

# Create directories
dir.create("data/raw_real_world_se/")

# Load all data sets in SummarizedExperiment objects
message("Load data sets into SE...")
sapply(list.files("real-world/data_preparation/", pattern = ".R", full.names = TRUE), source)
