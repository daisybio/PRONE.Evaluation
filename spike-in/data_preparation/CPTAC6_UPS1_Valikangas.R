###### ------- Spiked CPTAC6 Data Preparation Script (Tommy Valikangas) ------- ######

data_file <- "data/raw_spike_in_files/CPTAC6_UPS1/CPTAC6_UPS1_Valikangas_protein_intensities.csv"
md_file <- "data/raw_spike_in_files/CPTAC6_UPS1/CPTAC6_UPS1_Valikangas_metadata.csv"

if(file.exists(data_file)){
  message("File already exists. Skipping data preparation.")
  data <- fread(data_file, check.names = FALSE, integer64 = "numeric")
  md <- fread(md_file, check.names = FALSE)
} else {
  message("File does not exist. Proceeding with data preparation.")
  ## ----- Load Data and Convert To SummarizedExperiment -----
  data <- load("data/raw_spike_in_files/CPTAC6_UPS1/original_data/Spike_in_datas_Valikangas.RData")
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
  
  # write md and data file for zenodo
  write.csv(md, md_file, row.names = FALSE)
  write.csv(data, data_file, row.names = FALSE)
}

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
