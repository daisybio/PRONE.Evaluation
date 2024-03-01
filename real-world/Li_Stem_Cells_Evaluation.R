# Set data path
data_name <- "li_stem_cells"

# Set real-world configuration for this data set
TMT <- TRUE
NA_thr <- 0.8
DE_method <- "limma"
logFC <- TRUE
logFC_up <- log2(1.2)
logFC_down <- log2(0.83)
p_adj <- TRUE
alpha <- 0.05
comparisons <- NULL


# Check if DE and normalization objects already saved
if(!file.exists(paste0("data/de_results_real_world/", data_name, "_se_norm.rds"))){
# Perform preprocessing, normalization, DE analysis
  source("real-world/Real-World_Pipeline.R")
}

# Read DE results and normalized object
de_res <- readRDS(paste0("data/de_results_real_world/", data_name, "_de_res.rds"))
se_norm <- readRDS(paste0("data/preprocessed_normalized_real_world_se/", data_name, "_se_norm.rds"))

# Evaluation

# 1. Unnormalized data

log2_pca <- plot_PCA(se_norm, ain = "log2", color_by = "Batch", label_by = "No", shape_by = "Timepoint") + theme
log2_boxplot <- plot_boxplots(se_norm, ain = "log2", color_by = "Batch") + theme

log2_pca / log2_boxplot

# Quantitative and Qualitative Evaluation

plot_intragroup_PMAD(se_norm) + theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "None")


# Complete PCA plot of all normalization methods

generate_complete_SE <- function(se, ain = NULL){
  # create one big data table
  big_dt <- NULL
  big_cd <- NULL
  # check if assays all in se
  ain <- PRONE.R::check_input_assays(se, ain)
  for(a in ain){
    if(!a %in% c("raw")){
      dt <- as.data.table(assays(se)[[a]])
      colnames(dt) <- paste0(a, "_", colnames(dt))
      cd <- as.data.table(colData(se))
      cd$Column <- paste0(a,"_",cd$Column)
      cd$Normalization <- a
      if(is.null(big_dt)){
        big_dt <- dt
        big_cd <- cd
      } else {
        big_dt <- cbind(big_dt, dt)
        big_cd <- rbind(big_cd, cd)
      }
    }
  }
  se_big <- SummarizedExperiment::SummarizedExperiment(assays = list(all = big_dt), colData = big_cd, rowData = as.data.table(rowData(se)))
  metadata(se_big) <- metadata(se)
  return(se_big)
}



assays <- names(assays(se_norm))

assays <- assays[!assays %in% c("EigenMS", "MAD", "EigenMS_on_IRS", "EigenMS_on_limBE", "MAD_on_IRS", "MAD_on_limBE")]

# find assays with _on_IRS
norm_methods_IRS <- assays[grepl("_on_IRS", assays)]
se_IRS <- generate_complete_SE(se_norm, ain = norm_methods_IRS)
colData(se_IRS)$Normalization <- sapply(strsplit(colData(se_IRS)$Normalization, "_on_"), function(x) x[1])

# find assays with _on_limBE
norm_methods_limBE <- assays[grepl("_on_limBE", assays)]
se_limBE <- generate_complete_SE(se_norm, ain = norm_methods_limBE)
colData(se_limBE)$Normalization <- sapply(strsplit(colData(se_limBE)$Normalization, "_on_"), function(x) x[1])

# single normalization methods
norm_methods_single <- assays[!assays %in% c("limBE", "IRS", "log2", norm_methods_IRS, norm_methods_limBE)]
se_single <- generate_complete_SE(se_norm, ain = norm_methods_single)

# PCA plots
pca_single <- plot_PCA(se_single, color_by = "Normalization", ain = c("all"), label_by = "No", shape_by = "Batch") + 
  theme + theme(strip.text = element_blank()) +
  ggtitle("Single Normalization Methods")
pca_IRS <- plot_PCA(se_IRS, color_by = "Normalization", ain = c("all"), label_by = "No", shape_by = "Batch") + 
  theme + theme(strip.text = element_blank()) +
  ggtitle("IRS Normalization")
pca_limBE <- plot_PCA(se_limBE, color_by = "Normalization", ain = c("all"), label_by = "No", shape_by = "Batch") + 
  theme + theme(strip.text = element_blank()) +
  ggtitle("limBE Normalization")

pca_single + pca_IRS + pca_limBE + patchwork::plot_layout(guides = "collect")

# DEP Comparison

# Methods based on IRS-normalized data
de_res_IRS <- de_res[de_res$Assay %in% norm_methods_IRS,]
plot_overview_DE_bar(de_res_IRS, plot_type = "facet_regulation")
plot_overview_DE_bar(de_res_IRS[de_res_IRS$Assay != "EigenMS_on_IRS",], plot_type = "facet_regulation")

plot_overview_DE_tile(de_res_IRS)
plot_overview_DE_tile(de_res_IRS[de_res_IRS$Assay != "EigenMS_on_IRS",])

# Intersection Analysis

plot_upset_DE(de_res_IRS[de_res_IRS$Assay != "EigenMS_on_IRS",], plot_type = "stacked", min_degree = 12)


new_de_res_IRS <- PRONE.R::apply_thresholds(de_res_IRS)
plot_overview_DE_bar(new_de_res_IRS, plot_type = "facet_regulation")

plot_overview_DE_bar(new_de_res_IRS[!new_de_res_IRS$Assay %in% c("EigenMS_on_IRS", "TMM_on_IRS")], plot_type = "facet_regulation")

# New Heatmap (x -> normalization methods, y -> normalization methods, shared DEGs)
# get matrix with rows and columns being the normalization methods and the values the number of shared protein IDs
dt <- new_de_res_IRS
#dt <- dt[dt$Comparison == "D0-D14",]
dt <- dt[dt$Change != "No Change",]
dt <- dcast(dt, Assay ~ Protein.IDs, fun.aggregate = length)
assays <- dt$Assay
intersections <- as.data.frame(lapply(dt[, -1], function(x) as.integer(x > 0)))
matrix_dt <- as.matrix(intersections)
shared_ids_matrix <- matrix_dt %*% t(matrix_dt)
rownames(shared_ids_matrix) <- assays
colnames(shared_ids_matrix) <- assays
shared_ids_matrix[upper.tri(shared_ids_matrix)] <- NA
melted_dt <- melt(shared_ids_matrix, na.rm = TRUE)
ggplot(melted_dt, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "#f2f569", high = "#4daabd") +
  labs(x = "Normalization Method", y = "Normalization Method", fill = "Number of Shared DEGs") +
  theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# Jaccard index
dt <- t(matrix_dt)
dt <- as.data.frame(dt)
colnames(dt) <- assays

jaccard <- function(x, y){
  M.11 <- sum(x == 1 & y == 1)
  M.10 <- sum(x == 1 & y == 0)
  M.01 <- sum(x == 0 & y == 1)
  return(M.11 / (M.11 + M.10 + M.01))
}

m <- matrix(data = NA, nrow = ncol(dt), ncol = ncol(dt))
for(i in 1:ncol(dt)){
  for(j in 1:ncol(dt)){
    col1 <- colnames(dt)[i]
    col2 <- colnames(dt)[j]
    if(col1 == col2){
      m[i, j] <- 1
    } else if(i > j){
      m[i, j] <- jaccard(dt[, col1], dt[, col2])
    }
  }
}
colnames(m) <- colnames(dt)
rownames(m) <- colnames(dt)
melted_m <- melt(m, na.rm = TRUE)

ggplot(melted_m, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, digits = 2))) +
  scale_fill_gradient(low = "white", high = "#4daabd") +
  labs(x = "Normalization Method", y = "Normalization Method", fill = "Jaccard Similarity") +
  theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

plot_volcano_DE(new_de_res_IRS, comparisons = c("D0-D14"))

# check distribution of p.adjust values in the different normalization methods

ggplot(new_de_res_IRS, aes(x = adj.P.Val)) +
  geom_histogram(alpha = 0.5) +
  facet_wrap(~Assay, scales = "free") +
  labs(x = "p.adjust", fill = "Normalization Method") +
  theme

ggplot(new_de_res_IRS, aes(x = logFC)) +
  geom_histogram(alpha = 0.5) +
  facet_wrap(~Assay, scales = "free") +
  labs(x = "logFC", fill = "Normalization Method") +
  theme

ggplot(new_de_res_IRS, aes(x = P.Value)) +
  geom_histogram(alpha = 0.5) +
  facet_wrap(~Assay, scales = "free") +
  labs(x = "P.Value", fill = "Normalization Method") +
  theme

# Focus on one comparison

new_de_res_IRS_D0_D14 <- new_de_res_IRS[new_de_res_IRS$Comparison == "D0-D14",]

low_regulated_protein_de <- new_de_res_IRS_D0_D14[new_de_res_IRS_D0_D14$Protein.IDs == "A0A087WYH7;O15230",]

# Check data distribution of VSN, NormicsVSN, NormicsMedian

plot_densities(se_norm, ain = c("VSN", "Median", "RobNorm", "NormicsMedian", "log2", "NormicsVSN", "LoessF"), color_by = "Timepoint") + theme
plot_boxplots(se_norm, ain = c("VSN", "Median", "RobNorm", "NormicsMedian", "log2", "NormicsVSN", "LoessF"), color_by = "Timepoint") + theme


# Check performing log2 before NormicsMedian
se_norm <- vsnNorm(se, aout = "VSN_new")
plot_densities(se_norm, ain = c("VSN", "Median", "VSN_new", "LoessF"), color_by = "Timepoint") + theme


# try other vsn normalization

se <- readRDS("data/li_et_al_preprocessed_se.rds")

# Run IRS on raw data (but result = log2)
se <- irsNorm(se, ain = "raw", aout = "IRS")

# Run VSN on IRS (IRS -> raw data -> glog)
se <- vsnNorm(se, ain = "IRS", aout = "VSN_on_IRS")
se <- medianNorm(se, ain = "IRS", aout = "Median_on_IRS")
se <- robNorm(se, ain = "IRS", aout = "RobNorm_on_IRS")

se <- normicsNorm(se, ain = "IRS", aout = "NormicsVSN_on_IRS", method = "NormicsVSN", NormicsVSN_quantile = 0.9, reduce_correlation_by = 3)
se <- normicsNorm(se, ain = "IRS", aout = "NormicsVSN_on_IRS_0.8", method = "NormicsVSN", NormicsVSN_quantile = 0.8, reduce_correlation_by = 3)
se <- normicsNorm(se, ain = "IRS", aout = "NormicsMedian_on_IRS", method = "NormicsMedian", NormicsVSN_quantile = 0.8, reduce_correlation_by = 3)


SummarizedExperiment::assays(se)$IRS <- NULL
se <- remove_reference_samples(se)

comparisons <- specify_comparisons(se)
de_res <- run_DE(se, comparisons = comparisons)

plot_overview_DE_bar(de_res, plot_type = "facet_regulation")
