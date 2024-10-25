###### ------- E.coli Added to Human Background (MaxLFQ) Data Preparation Script ------- ######

data_file <- "data/raw_spike_in_files/Ecoli_human_DEqMS/Ecoli_human_DEqMS_protein_intensities.csv"
md_file <- "data/raw_spike_in_files/Ecoli_human_DEqMS/Ecoli_human_DEqMS_metadata.csv"

if(file.exists(data_file)){
  message("File already exists. Skipping data preparation.")
  data <- fread(data_file, check.names = FALSE, integer64 = "numeric")
  md <- fread(md_file, check.names = FALSE)
} else {
  message("File does not exist. Proceeding with data preparation.")
  ## ----- Load Data and Convert To SummarizedExperiment -----
  data <- fread("data/raw_spike_in_files/Ecoli_human_DEqMS/original_data/input_pwilmart.txt", integer64 = "numeric")
  
  columns <- colnames(data)[3:length(colnames(data))]
  md <- data.table(Column = columns)
  md$Condition <- sapply(strsplit(as.character(md$Column),"_"), "[", 3)
  md$Replicate <- sapply(strsplit(as.character(md$Column),"_"), "[", 1)
  md$Replicate <- ifelse(md$Replicate == "A", 1, ifelse(md$Replicate == "B", 2, ifelse(md$Replicate == "C", 3, 4)))
  md[md$Condition == "7pt5",]$Condition <- "7.5"
  md$Label <- paste0("Rep_", md$Replicate, "_Spiked_", md$Condition)
  md$Concentration <- as.numeric(md$Condition)
  
  concentrations <- data.table(Tmp = c("A","B","C"), Concentration = c(7.5, 15, 45))
  md <- merge(md, concentrations, by="Concentration")
  md$Condition <- NULL
  colnames(md)[colnames(md) == "Tmp"] <- "Condition"
  
  data$`Gene names` <- NA
  data$HorE <- ifelse(data$HorE == "E.coli", "ECOLI", "HUMAN")
  colnames(data)[colnames(data) == "HorE"] <- "Spiked"
  
  # write md and data file for zenodo
  write.csv(data, data_file, row.names = FALSE)
  write.csv(md, md_file, row.names = FALSE)
}

se <- load_spike_data(data, 
                      md, 
                      spike_column = "Spiked",
                      spike_value = "ECOLI",
                      spike_concentration = "Concentration",
                      protein_column = "Accession", 
                      gene_column = "Gene names",
                      ref_samples = NULL,
                      batch_column = NULL,
                      condition_column = "Condition",
                      label_column = "Label")
colData(se)$Condition <- factor(colData(se)$Condition, levels = unique(colData(se)$Condition))

saveRDS(se, "data/raw_spike_in_se/Ecoli_human_DEqMS_se.rds")
