## ----- Load required libraries -----
required_packages <- c("data.table", "readxl", "stringr", "ggplot2", "RColorBrewer", "ggpubr", "BiocManager", "patchwork", "dplyr")
for(package in required_packages){
  if(!require(package,character.only = TRUE, quietly = TRUE)) install.packages(package, dependencies = TRUE, quietly = TRUE)
  library(package, character.only = TRUE, quietly = TRUE)
}

if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) BiocManager::install("SummarizedExperiment")


# Install PRONE.R from github and build vignettes
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("lisiarend/PRONE", build_vignettes = FALSE, dependencies = TRUE)

# Load and attach PRONE.R 
library(PRONE)
library(SummarizedExperiment)

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
#norm_methods <- get_normalization_methods()
#norm_methods <- norm_methods[!norm_methods %in% c("IRS", "limBE")]
#qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual",]
#col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#col_vector_norm <- col_vector[10:(length(norm_methods)+9)]
#names(col_vector_norm) <- norm_methods

col_vector_norm <- c("#E6AB02","#A6761D","#E7298A","#FDDAEC","#B15928", "#E41A1C", "#A6CEE3", "#1F78B4", "#B2DF8A", "#4DAF4A","#355E3B", "#D9D9D9",
                    "grey40","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A", "#eb7f78","#ece4f1")
names(col_vector_norm) <- c("GlobalMean", "GlobalMedian", "Median", "Mean", "Quantile", "VSN", "LoessF", "LoessCyc", "RLR", "RlrMA", "RlrMACyc", "EigenMS", "MAD", "RobNorm", "TMM", "NormicsVSN", "NormicsMedian", "VSN_0.5", "NormicsVSN_0.5")
col_vector_norm <- c("log2" = "yellow", col_vector_norm)

