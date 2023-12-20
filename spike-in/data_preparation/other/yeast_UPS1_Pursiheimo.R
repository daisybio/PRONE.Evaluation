###### ------- Spiked Yeast Lysate UPS1 (Pursiheimo) Data Preparation Script ------- ######

## ----- Load Data and Convert To SummarizedExperiment -----
data <- fread("data/raw_spike_in_files/yeast_UPS1_Pursiheimo/YEAST-Data-NonNormalized.csv", check.names = FALSE, integer64 = "numeric", dec=",")


# Exclude proteins that were only quantified using 1 peptide
data <- data[data$`Peptides used for quantitation` >= 2,]


md <- data.table(Column = colnames(data)[4:length(colnames(data))])
md[7,]$Column <- "110616_yeast_ups1_10fmol_r1"
colnames(data)[4:length(colnames(data))] <- md$Column
md$Replicate <- sapply(strsplit(as.character(md$Column),"_"), "[", 5)
md$Condition <- sapply(strsplit(as.character(md$Column),"_"), "[", 4)
md$Label <- paste0(md$Condition, "_", md$Replicate)
md$Concentration <- as.numeric(sapply(strsplit(as.character(md$Condition),"fmol"), "[", 1))

concentrations <- data.table(Tmp = c("A","B", "C", "D", "E"), Concentration = c(2, 4, 10,25,50))
md <- merge(md, concentrations, by="Concentration")
md$Condition <- NULL
colnames(md)[colnames(md) == "Tmp"] <- "Condition"

data$Spiked <- rep("BG", nrow(data))
data$Spiked[grepl("ups", data$Accession)] <- "UPS1"
data$Gene.names <- ""

se <- load_spike_data(data, 
                      md, 
                      spike_column = "Spiked",
                      spike_value = "UPS1",
                      spike_concentration = "Concentration",
                      protein_column = "Accession", 
                      gene_column = "Gene.names",
                      ref_samples = NULL,
                      batch_column = NULL,
                      condition_column = "Condition",
                      label_column = "Label")
colData(se)$Condition <- factor(colData(se)$Condition, levels = unique(colData(se)$Condition))

saveRDS(se, "data/raw_spike_in_se/yeast_UPS1_Pursiheimo_se.rds")
