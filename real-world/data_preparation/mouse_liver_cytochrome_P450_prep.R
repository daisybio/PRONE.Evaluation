###### ------- Mouse Liver Cytochrome P450 ------- ######

## ----- Load Data and Convert To SummarizedExperiment -----
data <- fread("data/raw_real_world_files/mouse_liver_cytochrome_P450/MOUSE-Data-NonNormalized.csv")

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

se <- load_data(data, 
                md, 
                protein_column = "Accession", 
                gene_column = "Gene names",
                ref_samples = NULL,
                batch_column = NULL,
                condition_column = "Condition",
                label_column = "Column")

saveRDS(se, "data/raw_real_world_se/mouse_liver_cytochrome_P450_se.rds")
