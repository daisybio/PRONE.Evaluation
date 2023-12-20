###### ------- E.coli Added to Human Background (MaxLFQ) Data Preparation Script ------- ######

## ----- Load Data and Convert To SummarizedExperiment -----
data <- fread("data/raw_spike_in_files/Ecoli_human_MaxLFQ/proteomebenchmark/proteinGroups.txt", check.names = TRUE, integer64 = "numeric", dec=",")


# Check if some protein groups are mixed
mixed <- grepl("Homo sapiens.*Escherichia|Escherichia.*Homo sapiens", data$Fasta.headers)
data <- data[!mixed,]

data$Spiked <- rep("HUMAN", nrow(data))
data$Spiked[grepl("ECOLI", data$Fasta.headers)] <- "ECOLI"

columns <- colnames(data)[grepl("LFQ.intensity", colnames(data))]
md <- data.table(Column = columns)
md$Label <- sapply(strsplit(as.character(md$Column),"[.]"), "[", 3)
md$Replicate <- sapply(strsplit(as.character(md$Label),""), "[", 2)
md$Condition <- sapply(strsplit(as.character(md$Label),""), "[", 1)
md$Condtion <- factor(md$Condition, level = c("L", "H"))
md$Concentration <- c(30,30,30,10,10,10)

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

saveRDS(se, "data/raw_spike_in_se/Ecoli_human_MaxLFQ_se.rds")

