###### ------- CPTAC OV Data Preparation Script ------- ######


data_file <- "data/raw_real_world_files/CPTAC_OV/CPTAC_OV_protein_intensities.csv"
md_file <- "data/raw_real_world_files/CPTAC_OV/CPTAC_OV_metadata.csv"

if(file.exists(data_file)){
  message("File already exists. Skipping data preparation.")
  data <- fread(data_file, check.names = FALSE, integer64 = "numeric")
  md <- fread(md_file, check.names = FALSE)
} else {
  message("File does not exist. Proceeding with data preparation.")
  ## ----- Load Data and Convert To SummarizedExperiment -----
  data <- fread("data/raw_real_world_files/CPTAC_OV/original_data/TMT_CPTAC_MaxQuant/JHU_output/proteinGroups.txt")

  # Original publication -> design
  design <- read_excel("data/raw_real_world_files/CPTAC_OV/original_data/original_publication/NIHMS1668244-supplement-Table_S1.xlsx", sheet = "Experimental Design")
  design <- as.data.table(design)
  design <- design[1:13, 2:11]
  
  # Make a long vector out of the tibble
  sample_vector <- unlist(lapply(seq_len(nrow(design)), function(i) as.list(design[i,])))
  
  # Get Column names
  column_vector <- colnames(data)[grep("Reporter intensity corrected", colnames(data))]
  
  # Create meta-data
  md <- data.table("Column" = column_vector, "Sample_ID" = sample_vector)
  
  # Original publication -> conditions
  clinical <- read_excel("data/raw_real_world_files/CPTAC_OV/original_data/original_publication/NIHMS1668244-supplement-Table_S1.xlsx", sheet = "Clinical Information")
  clinical <- as.data.table(clinical)
  
  # QC (9), pooled samples (13), + SPL050 not included in clinical information
  
  # Remove samples SPL039 and SPL047
  md <- md[!md$Sample_ID %in% c("SPL 039", "SPL 047")] # no metadata for these samples
  
  # Merge md with clinical data
  md <- merge(md, clinical, by = "Sample_ID", all.x = TRUE, all.y = FALSE)
  
  md$Pool <- sapply(strsplit(md$Column," "), "[", 5)
  
  ref_samples <- md[md$Sample_ID == "pooled sample",]$Column
  
  # exclude qc samples from meta data
  md <- md[!md$Sample_ID %in% c("QC1", "QC2", "QC3", "QC4", "QC5", "QC6", "QC7", "QC8", "QC9")]
  
  # add ref to condition
  md[md$Sample_ID == "pooled sample"]$Pathological_Status <- "pooled sample"
  
  # add numbers to pooled samples
  md[md$Sample_ID == "pooled sample"]$Sample_ID <- paste(md[md$Sample_ID == "pooled sample"]$Sample_ID, md[md$Sample_ID == "pooled sample"]$Pool)
  
  # replace NAs with "pooled samples" in all columns for pooled sample except Sample_ID, Column, Pool, Pathological_Status
  md[is.na(md)] <- "NA"
  
  # write md and data file for zenodo
  write.csv(data, data_file, row.names = FALSE)
  write.csv(md, md_file, row.names = FALSE)
}

# Load SummarizedExperiment
se <- load_data(data, 
                md, 
                protein_column = "Protein IDs", 
                gene_column = "Gene names",
                ref_samples = ref_samples,
                batch_column = "Pool",
                condition_column = "Pathological_Status",
                label_column = "Sample_ID")

cd <- as.data.frame(colData(se))
cd$Pathological_Status[cd$Pathological_Status == "Non-Tumor"] <- "NonTumor"
colData(se) <- DataFrame(cd)
saveRDS(se, "data/raw_real_world_se/CPTAC_OV_se.rds")
