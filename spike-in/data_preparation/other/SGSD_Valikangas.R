###### ------- Spiked SGSD Data Preparation Script (Tommy Valikangas) ------- ######

## ----- Load Data and Convert To SummarizedExperiment -----
data <- load("data/raw_spike_in_files/Spike_in_datas_Valikangas.RData")
data <- sgsd.data

md <- data.table(Column = colnames(data))
md$Replicate <- sapply(strsplit(as.character(md$Column),"_"), "[", 4)
md$Condition <- c(rep("A", 3), rep("B", 3), rep("C", 3), rep("D", 3), rep("E", 3), rep("F", 3), rep("G",3), rep("H", 3))
md$Label <- paste0(md$Condition, "_", md$Replicate)

data$Protein.IDs <- rownames(data)
data$Gene.names <- ""
data$Spiked <- rep("BG", nrow(data))
data$Spiked[data$Protein.IDs %in% sgsds.spike] <- "Spiked"

se <- load_spike_data(data, 
                      md, 
                      spike_column = "Spiked",
                      spike_value = "Spiked",
                      spike_concentration = NULL,
                      protein_column = "Protein.IDs", 
                      gene_column = "Gene.names",
                      ref_samples = NULL,
                      batch_column = NULL,
                      condition_column = "Condition",
                      label_column = "Label")
colData(se)$Condition <- factor(colData(se)$Condition, levels = unique(colData(se)$Condition))

saveRDS(se, "data/raw_spike_in_se/yeast_UPS1_Pursiheimo_Valikangas_se.rds")


