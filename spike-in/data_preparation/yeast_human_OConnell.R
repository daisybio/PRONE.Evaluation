###### ------- Yeast Added to Human Background (O'Connell/EmpiRe) Data Preparation Script ------- ######


data_file <- "data/raw_spike_in_files/yeast_human_Oconnell/yeast_human_Oconnell_protein_intensities.csv"
md_file <- "data/raw_spike_in_files/yeast_human_Oconnell/yeast_human_Oconnell_metadata.csv"

if(file.exists(data_file)){
  message("File already exists. Skipping data preparation.")
  data <- fread(data_file, check.names = FALSE, integer64 = "numeric")
  md <- fread(md_file, check.names = FALSE)
} else {
  ## ----- Load Data and Convert To SummarizedExperiment -----
  data <- fread("data/raw_spike_in_files/yeast_human_Oconnell/original_data/proteinGroups.txt", integer64 = "numeric")
  
  colnames(data)[colnames(data) == "Potential contaminant"] <- "Potential.contaminant"
  colnames(data)[colnames(data) == "Only identified by site"] <- "Only.identified.by.site"
  
  md <- fread("data/raw_spike_in_files/yeast_human_Oconnell/original_data/labels.txt", header = FALSE)
  colnames(md) <- c("Column", "Condition")
  md$Concentration <- c(rep(10,3), rep(5, 4), rep(3.3, 4))
  md$Condition <- factor(c(rep("C",3), rep("B", 4), rep("A", 4)), levels = c("A", "B", "C"))
  md <- md[order(md$Condition),]
    
  # get human and yeast ids
  human <- fread("data/raw_spike_in_files/yeast_human_Oconnell/original_data/human_ids.txt", header = FALSE)
  human <- human$V1
  
  yeast <- fread("data/raw_spike_in_files/yeast_human_Oconnell/original_data/yeast_ids.txt", header = FALSE)
  yeast <- yeast$V1
  
  get_organism <- function(ids){
    ids <- str_split_1(ids, ";")
    org <- sapply(ids, function(id){
      if(id %in% human){
        return("HUMAN")
      } else if (id %in% yeast){
        return("YEAST")
      } else return("NA")
    })
    if(all(org == "HUMAN")){
      return("HUMAN")
    } else if (all(org == "YEAST")){
      return("YEAST")
    } else return(NA)
  }
  
  organisms <- sapply(data$`Protein IDs`, get_organism)
  orgs <- data.table(`Protein IDs` = names(organisms), "Spiked" = unlist(organisms))
  
  data <- merge(data, orgs, by = "Protein IDs")
  
  data <- data[!is.na(data$Spiked),]
  
  # write md and data file for zenodo
  write.csv(data, data_file, row.names = FALSE)
  write.csv(md, md_file, row.names = FALSE)
}

se <- load_spike_data(data, 
                      md, 
                      spike_column = "Spiked",
                      spike_value = "YEAST",
                      spike_concentration = "Concentration",
                      protein_column = "Protein IDs", 
                      gene_column = "Gene names",
                      ref_samples = NULL,
                      batch_column = NULL,
                      condition_column = "Condition",
                      label_column = NULL)
colData(se)$Condition <- factor(colData(se)$Condition, levels = unique(colData(se)$Condition))

saveRDS(se, "data/raw_spike_in_se/yeast_human_OConnell_se.rds")





