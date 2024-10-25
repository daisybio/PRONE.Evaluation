###### ------- Li Stem Cells Data Preparation Script ------- ######


data_file <- "data/raw_real_world_files/li_stem_cells/li_stem_cells_protein_intensities.csv"
md_file <- "data/raw_real_world_files/li_stem_cells/li_stem_cells_metadata.csv"

if(file.exists(data_file)){
  message("File already exists. Skipping data preparation.")
  data <- fread(data_file, check.names = FALSE, integer64 = "numeric")
  md <- fread(md_file, check.names = FALSE)
} else {
  message("File does not exist. Proceeding with data preparation.")
  ## ----- Load Data and Convert To SummarizedExperiment -----
  
  data <- fread("data/raw_real_world_files/li_stem_cells/original_data/proteinGroups_CoSy.Bio_filtered_remapped.txt")
  
  # create a metadata data.table
  md <- data.frame(Column = colnames(data)[5:19])
  md$Batch <- sapply(strsplit(as.character(md$Column),"_"), "[", 1)
  tmp <- sapply(strsplit(as.character(md$Column),"_"), "[", 2)
  md$Timepoint <- sapply(strsplit(as.character(tmp),"[-]"), "[", 1)
  md$Replicate <- sapply(strsplit(as.character(tmp),"[-]"), "[", 2)
  md$Replicate[is.na(md$Replicate)] <- "reference"
  ref_samples <- md[md$Replicate == "reference",]$Column
  
  # write md and data file for zenodo
  write.csv(data, data_file, row.names = FALSE)
  write.csv(md, md_file, row.names = FALSE)
}

se <- load_data(data, 
                md, 
                protein_column = "Protein IDs", 
                gene_column = "Gene names",
                ref_samples = ref_samples,
                batch_column = "Batch",
                condition_column = "Timepoint",
                label_column = "Column")

saveRDS(se, "data/raw_real_world_se/li_stem_cells_se.rds")
