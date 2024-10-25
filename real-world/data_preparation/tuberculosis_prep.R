###### ------- Tuberculosis Data Preparation Script ------- ######

data_file <- "data/raw_real_world_files/tuberculosis/tuberculosis_protein_intensities.csv"
md_file <- "data/raw_real_world_files/tuberculosis/tuberculosis_metadata.csv"

if(file.exists(data_file)){
  message("File already exists. Skipping data preparation.")
  data <- fread(data_file, check.names = FALSE, integer64 = "numeric")
  md <- fread(md_file, check.names = FALSE)
} else {
  message("File does not exist. Proceeding with data preparation.")
  ## ----- Load Data and Convert To SummarizedExperiment -----
  data <- fread("data/raw_real_world_files/tuberculosis/original_data/proteinGroups.txt", integer64 = "numeric")
  md <- fread("data/raw_real_world_files/tuberculosis/original_data/quant.columns.txt")
  colnames(md) <- c("Column", "Label", "Condition", "Batch")
  md$Condition <- plyr::mapvalues(md$Condition, from=c("Common.reference","Healthy.control", "Pulmonary.tuberculosis", "Tuberculous.lymphadenitis", "Treated.pulmonary.tuberculosis"), to = c("ref", "HC", "PTB", "TBL", "Rx"))
  
  ref_samples <- md[md$Condition == "ref",]$Column
  
  md <- md[, c("Column", "Condition", "Batch", "Label"), with = FALSE]
  
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
                condition_column = "Condition",
                label_column = "Label")
se$Condition <- factor(se$Condition, levels = c("ref", "HC", "PTB", "Rx", "TBL"))

se <- remove_samples_manually(se, column = "Label", values = c("1.HC_Pool1", "1.HC_Pool2"))

saveRDS(se, "data/raw_real_world_se/tuberculosis_se.rds")
