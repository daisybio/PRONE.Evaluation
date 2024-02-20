###### ------- Breast Cancer (proteiNorm) Data Preparation Script ------- ######

## ----- Load Data and Convert To SummarizedExperiment -----

data <- fread("data/raw_real_world_files/breast_cancer_proteiNorm/proteinGroups.txt", check.names = TRUE)

filtered_data <- fread("data/raw_real_world_files/breast_cancer_proteiNorm/proteinGroups_filtered.txt")
sample_names <- colnames(filtered_data)[-1]

md <- read_excel("data/raw_real_world_files/breast_cancer_proteiNorm/meta.xlsx")
colnames(md) <- c("Sample", "Label", "Condition", "Batch")
md$Column <- sample_names

ref_samples <- md[md$Condition == "pool", ]$Column

se <- load_data(data, 
                md, 
                protein_column = "Protein.IDs", 
                gene_column = "Gene.names",
                ref_samples = ref_samples,
                batch_column = "Batch",
                condition_column = "Condition",
                label_column = "Label")

saveRDS(se, "data/raw_real_world_se/breast_cancer_se.rds")
