###### ------- Spiked CPTAC6 Data Preparation Script (Tommy Valikangas) ------- ######

## ----- Load Data and Convert To SummarizedExperiment -----
data <- load("data/raw_spike_in_files/Spike_in_datas_Valikangas.RData")
data <- cptac.data

md <- data.table(Column = colnames(data))
md$Replicate <- sapply(strsplit(as.character(md$Column),"_"), "[", 1)
md$Condition <- sapply(strsplit(as.character(md$Column),"_"), "[", 4)
md$Condition <- substr(md$Condition, start = 2, stop = 2)
md$Label <- paste0(md$Condition, "_", md$Replicate)

concentrations <- data.table(Condition = c("A","B", "C", "D"), Concentration = c(0.25, 0.74, 2.2,6.7))
md <- merge(md, concentrations, by="Condition")

data$Gene.names <- ""
data$Protein.IDs <- rownames(data)
data$Spiked <- rep("BG", nrow(data))
data$Spiked[data$Protein.IDs %in% cptac.spike] <- "UPS1"

se <- load_spike_data(data, 
                      md, 
                      spike_column = "Spiked",
                      spike_value = "UPS1",
                      spike_concentration = "Concentration",
                      protein_column = "Protein.IDs", 
                      gene_column = "Gene.names",
                      ref_samples = NULL,
                      batch_column = NULL,
                      condition_column = "Condition",
                      label_column = "Label")
colData(se)$Condition <- factor(colData(se)$Condition, levels = unique(colData(se)$Condition))

saveRDS(se, "data/raw_spike_in_se/CPTAC6_UPS1_Valikangas_se.rds")
