###### ------- E.coli Added to Human Background (Ionstar) Data Preparation Script ------- ######

## ----- Load Data and Convert To SummarizedExperiment -----
data <- fread("data/raw_spike_in_files/Ecoli_human_Ionstar/proteinGroups.txt", check.names = TRUE, integer64 = "numeric", dec=",")
data <- as.data.table(data)

# Check if some protein groups are mixed
mixed <- grepl("Homo sapiens.*Escherichia|Escherichia.*Homo sapiens", data$Fasta.headers)
data <- data[!mixed,]

data$Spiked <- rep("HUMAN", nrow(data))
data$Spiked[grepl("Escherichia coli", data$Fasta.headers)] <- "ECOLI"

columns <- colnames(data)[grepl("Intensity", colnames(data))]
columns <- columns[2:length(columns)]
md <- data.table(Column = columns)
md$Label <- sapply(strsplit(as.character(md$Column),"[.]"), "[", 2)
md$Replicate <- sapply(strsplit(as.character(md$Label),""), "[", 2)
md$Condition <- toupper(sapply(strsplit(as.character(md$Label),""), "[", 1))

concentrations <- data.table(Condition = c("A","B", "C", "D", "E"), Concentration = c(3, 4.5, 6, 7.5, 9))
md <- merge(md, concentrations, on="Condition")
md$Label <- paste0("Rep_", md$Replicate, "_Spiked_", md$Condition, "(", md$Concentration, "%)")

se <- load_spike_data(data, 
                      md, 
                      spike_column = "Spiked",
                      spike_value = "ECOLI",
                      spike_concentration = "Concentration",
                      protein_column = "Protein.IDs", 
                      gene_column = "Gene.names",
                      ref_samples = NULL,
                      batch_column = NULL,
                      condition_column = "Condition",
                      label_column = "Label")
colData(se)$Condition <- factor(colData(se)$Condition, levels = unique(colData(se)$Condition))

saveRDS(se, "data/raw_spike_in_se/Ecoli_human_Ionstar_se.rds")
