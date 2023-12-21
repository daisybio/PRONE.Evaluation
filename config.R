## ----- Load required libraries -----
required_packages <- c("data.table", "readxl", "SummarizedExperiment", "stringr", "ggplot2", "RColorBrewer", "ggpubr")
for(package in required_packages){
  if(!require(package,character.only = TRUE, quietly = TRUE)) install.packages(package, dependencies = TRUE, quietly = TRUE)
  library(package, character.only = TRUE, quietly = TRUE)
}

# Install PRONE.R from github and build vignettes
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("lisiarend/PRONE.R", build_vignettes = TRUE, dependencies = TRUE, force = TRUE)
# Load and attach PRONE.R 
library("PRONE.R")

# Source external functions
source("spike-in/functions.R")


# Theme for ggplot2
old_theme <- theme_get()
new_theme <- old_theme + theme(axis.text.x = element_text(size = 12), 
                               axis.text.y = element_text(size = 12), 
                               strip.text = element_text(size = 12), 
                               axis.title.x = element_text(size = 12), 
                               axis.title.y = element_text(size = 12),
                               legend.text = element_text(size = 12),
                               legend.title = element_text(size = 12))
theme_set(new_theme)

theme <- theme_bw() + theme(axis.text.x = element_text(size = 12), 
                            axis.text.y = element_text(size = 12), 
                            strip.text = element_text(size = 12), 
                            axis.title.x = element_text(size = 12), 
                            axis.title.y = element_text(size = 12),
                            legend.text = element_text(size = 12),
                            legend.title = element_text(size = 12)
)

# Color vector
norm_methods <- get_normalization_methods()
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual",]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector_norm <- col_vector[10:(length(norm_methods)+9)]
names(col_vector_norm) <- norm_methods
