###### ------- Mouse Liver Cytochrome P450 ------- ######

data_file <- "data/raw_real_world_files/mouse_liver_cytochrome_P450/mouse_liver_cytochrome_P450_protein_intensities.csv"
md_file <- "data/raw_real_world_files/mouse_liver_cytochrome_P450/mouse_liver_cytochrome_P450_metadata.csv"

if(file.exists(data_file)){
  message("File already exists. Skipping data preparation.")
  data <- fread(data_file, check.names = FALSE, integer64 = "numeric", check.names = FALSE)
  md <- fread(md_file, check.names = FALSE)
} else {
  message("File does not exist. Proceeding with data preparation.")
  ## ----- Load Data and Convert To SummarizedExperiment -----
  data <- fread("data/raw_real_world_files/mouse_liver_cytochrome_P450/original_data/MOUSE-Data-NonNormalized.csv")
  
  # create a metadata data.table
  md <- data.frame(Column = colnames(data)[3:ncol(data)])
  md$Animal <- sapply(strsplit(as.character(md$Column),"_"), "[", 1)
  md$Condition <- sapply(strsplit(as.character(md$Column),"_"), "[", 2)
  
  accessions <- data$Accession
  used_for_quant <- data$`Peptides used for quantitation`
  
  data <- sapply(md$Column, function(x){
    column <- data[, get(x)]
    return(as.numeric(gsub(",", ".",column)))
  })
  data <- as.data.table(data)
  data$`Gene names` <- NA
  data$Accession <- accessions
  data$`Peptides used for quantitation` <- used_for_quant
  
  
  # write md and data file for zenodo
  write.csv(data, data_file, row.names = FALSE)
  write.csv(md, md_file, row.names = FALSE)
}


se <- load_data(data, 
                md, 
                protein_column = "Accession", 
                gene_column = "Gene names",
                ref_samples = NULL,
                batch_column = NULL,
                condition_column = "Condition",
                label_column = "Column")

saveRDS(se, "data/raw_real_world_se/mouse_liver_cytochrome_P450_se.rds")
