###### ------- Spiked Yeast Lysate UPS1 (Ramus) Data Preparation Script ------- ######

## ----- Load Data and Convert To SummarizedExperiment -----
data <- fread("data/raw_spike_in_files/yeast_UPS1_Ramus/proteinGroups.txt", check.names = FALSE, integer64 = "numeric")

# Check if some protein groups are mixed
mixed <- grepl("Homo sapiens.*Saccharomyces|Saccharomyces.*Homo sapiens", data$`Fasta headers`)
data <- data[!mixed,]

colnames(data)[colnames(data) == "Potential contaminant"] <- "Potential.contaminant"
colnames(data)[colnames(data) == "Only identified by site"] <- "Only.identified.by.site"

data$Spiked <- rep("BG", nrow(data))
data$Spiked[grepl("ups", data$`Protein IDs`)] <- "UPS1"

md <- read_excel("data/raw_spike_in_files/yeast_UPS1_Ramus/metadata.xlsx", col_names = FALSE)
md <- as.data.table(md)
colnames(md) <- c("Intensity", "Column", "Label", "Condition", "Batch")
md$Intensity <- NULL
md$Batch <- NULL
md$Concentration <- as.numeric(md$Condition)
md$Condition <- NULL
md$Concentration <- md$Concentration / 1000
md$Replicate <- sapply(strsplit(as.character(md$Label),"_"), "[", 2)
conditions <- data.table("Concentration" = sort(unique(md$Concentration)), "Condition" = c("A", "B", "C", "D", "E", "F", "G", "H", "I"))
md <- merge(md, conditions, by="Concentration")
md$Label <- paste0(md$Condition , "_", md$Replicate)

comparisons <- c("I-D", "I-F", "H-G")



se <- load_spike_data(data, 
                      md, 
                      spike_column = "Spiked",
                      spike_value = "UPS1",
                      spike_concentration = "Concentration",
                      protein_column = "Protein IDs", 
                      gene_column = "Gene names",
                      ref_samples = NULL,
                      batch_column = NULL,
                      condition_column = "Condition",
                      label_column = "Label")
colData(se)$Condition <- factor(colData(se)$Condition, levels = unique(colData(se)$Condition))

se <- se[,!se$Condition %in% c("A", "B", "C", "D")]
se$Condition <- droplevels(se$Condition)

saveRDS(se, "data/raw_spike_in_se/yeast_UPS1_Ramus_se.rds")

