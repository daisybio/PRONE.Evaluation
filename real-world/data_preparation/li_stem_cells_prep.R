###### ------- Li Stem Cells Data Preparation Script ------- ######

## ----- Load Data and Convert To SummarizedExperiment -----

data <- fread("data/raw_real_world_files/li_stem_cells/proteinGroups_CoSy.Bio_filtered_remapped.txt")

# create a metadata data.table
md <- data.frame(Column = colnames(data)[5:19])
md$Batch <- sapply(strsplit(as.character(md$Column),"_"), "[", 1)
tmp <- sapply(strsplit(as.character(md$Column),"_"), "[", 2)
md$Timepoint <- sapply(strsplit(as.character(tmp),"[-]"), "[", 1)
md$Replicate <- sapply(strsplit(as.character(tmp),"[-]"), "[", 2)
md$Replicate[is.na(md$Replicate)] <- "reference"
ref_samples <- md[md$Replicate == "reference",]$Column

se <- load_data(data, 
                md, 
                protein_column = "Protein IDs", 
                gene_column = "Gene names",
                ref_samples = ref_samples,
                batch_column = "Batch",
                condition_column = "Timepoint",
                label_column = "Column")

saveRDS(se, "data/raw_real_world_se/li_stem_cells_se.rds")
