# Install PRONE.R from github and build vignettes
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("lisiarend/PRONE.R", build_vignettes = TRUE, dependencies = TRUE, force = TRUE)
# Load and attach PRONE.R 
library("PRONE.R")

# Set parameters
DE_method <- "ROTS"
logFC <- FALSE
logFC_up <- 1
logFC_down <- -1
p_adj <- TRUE
alpha <- 0.05

for(se_rds in list.files(file.path("in/"), pattern = "_se.rds")){
  # Read RDS
  se_norm <- readRDS(file.path("in/", se_rds))

  # Specify comparisons
  comparisons <- PRONE::specify_comparisons(se_norm, condition = "Condition", sep = NULL, control = NULL)

  # Run DE
  de_results <- PRONE::run_DE(se = se_norm, 
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
  
  # Save DE as RDS
  de_results_file <- paste0(strsplit(se_rds, "_se_norm.rds")[1][[1]], "_de_results.rds")
  saveRDS(de_results, file.path("out/", de_results_file))
}

