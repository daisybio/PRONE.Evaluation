###### ------- Spiked Yeast Lysate UPS1 (Pursiheimo) Data Preparation Script (Tommy Valikangas) ------- ######


## ----- Load Data and Convert To SummarizedExperiment -----
data <- load("data/raw_spike_in_files/Spike_in_datas_Valikangas.RData")
data <- ups1.data

md <- data.table(Column = colnames(data))
md[md$Column == "110616_yeast_ups_10fmol"]$Column <- "110616_yeast_ups_10fmol_r1"
colnames(data) <- md$Column
md$Replicate <- sapply(strsplit(as.character(md$Column),"_"), "[", 5)
md$Condition <- sapply(strsplit(as.character(md$Column),"_"), "[", 4)
md$Label <- paste0(md$Condition, "_", md$Replicate)
md$Concentration <- as.numeric(sapply(strsplit(as.character(md$Condition),"fmol"), "[", 1))

concentrations <- data.table(Tmp = c("A","B", "C", "D", "E"), Concentration = c(2, 4, 10,25,50))
md <- merge(md, concentrations, by="Concentration")
md$Condition <- NULL
colnames(md)[colnames(md) == "Tmp"] <- "Condition"

data$Spiked <- rep("BG", nrow(data))
data$Spiked[grepl("ups", rownames(data))] <- "UPS1"
data$Gene.names <- ""
data$Protein.IDs <- rownames(data)

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

saveRDS(se, "data/raw_spike_in_se/yeast_UPS1_Pursiheimo_Valikangas_se.rds")
