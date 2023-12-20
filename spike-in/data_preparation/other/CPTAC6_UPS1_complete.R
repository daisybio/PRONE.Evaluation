###### ------- CPTAC Study 6 (UPS1 Spiked) Data Preparation Script ------- ######


## ----- Load Data and Convert To SummarizedExperiment -----
data <- fread("data/raw_spike_in_files/CPTAC6_UPS1/complete/proteinGroups.txt", check.names = FALSE, integer64 = "numeric")

colnames(data)[colnames(data) == "Potential contaminant"] <- "Potential.contaminant"
colnames(data)[colnames(data) == "Only identified by site"] <- "Only.identified.by.site"

# Check if some protein groups are mixed
mixed <- grepl("ups.*YEAST|YEAST.*ups", data$`Protein IDs`)
data <- data[!mixed,]

data$Spiked <- rep("BG", nrow(data))
data$Spiked[grepl("ups", data$`Protein IDs`)] <- "UPS1"
data$`Gene names` <- NA


column <- colnames(data)[grepl("Intensity", colnames(data))]
column <- column[2:length(column)]
md <- data.table(Column = column)
md$Label <- sapply(strsplit(as.character(md$Column)," "), "[", 2)
md$Replicate <- sapply(strsplit(as.character(md$Label),"_"), "[", 2)
md$Condition <- sapply(strsplit(as.character(md$Label),"_"), "[", 1)
md$Condition <- sapply(strsplit(as.character(md$Condition),""), "[", 2)
md$Lab <- ifelse(md$Replicate %in% c(1,2,3), "Lab1", ifelse(md$Replicate %in% c(4,5,6), "Lab2", "Lab3"))

concentrations <- data.table(Condition = c("A","B", "C", "D", "E"), Concentration = c(0.25, 0.74, 2.2,6.7,20))
md <- merge(md, concentrations, by="Condition")

# Select One Lab to Analyze
exclude <- md[md$Lab %in% c("Lab2", "Lab3"),]$Column
md <- md[md$Lab %in% c("Lab1"),]
data <- data[, !colnames(data) %in% exclude, with=FALSE]

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

# remove condition E since Valikangas did not included E in their analysis
se <- se[,se$Condition != "E"]
se$Condition <- droplevels(se$Condition)

saveRDS(se, "data/raw_spike_in_se/CPTAC6_UPS1_complete_se.rds")

