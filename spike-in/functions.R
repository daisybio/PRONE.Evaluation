#' Return list containing vector positions of values in string
#'
#' @param targetVector vector of sample groups
#' @return indexList List where key is condition level and values are indices
#'   for the condition
#' @keywords internal
getIndexList <- function(targetVector) {
  
  indexList <- list()
  uniqVals <- unique(targetVector)
  for (val in uniqVals) {
    indexList[[toString(val)]] <- which(targetVector == val)
  }
  indexList
}


#' Calculate average MAD (Median Absolute Deviation) for each feature in
#' each condition and then calculates the average for each replicate group
#'
#' @param methodList List containing normalized matrices.
#' @param sampleReplicateGroups Condition header.
#' @return condAvgMadMat Matrix with average MAD for each biological condition.
#' @keywords internal
calculateAvgMadMem <- function(methodList, sampleReplicateGroups) {
  
  groupIndexList <- getIndexList(sampleReplicateGroups)
  
  calculateAvgFeatureMadForGroup <- function(groupIndices, methodData) {
    groupData <- methodData[,..groupIndices] # .. added
    featureMAD <- matrixStats::rowMads(as.matrix(groupData), na.rm=TRUE) # as.matrix added
    featureMAD
  }
  
  calculateGroupMadForMethod <- function(methodData, groups, indexList) {
    
    # Extracts groups of replicates and calculate MAD for each feature
    featureMADMat <- vapply(
      indexList,
      calculateAvgFeatureMadForGroup,
      rep(0, nrow(methodData)),
      methodData=methodData
    )
    
    methodRepGroupMADMean <- colMeans(featureMADMat, na.rm=TRUE)
    methodRepGroupMADMean
  }
  
  condAvgMadMat <- vapply(
    methodList,
    calculateGroupMadForMethod,
    rep(0, length(unique(sampleReplicateGroups))),
    groups=sampleReplicateGroups,
    indexList=groupIndexList
  )
  
  condAvgMadMat
}

#' General function for calculating percentage difference of average column
#' means in matrix
#'
#' @param targetMat Matrix for which column means should be compared
#' @return percDiffVector Vector with percentage difference, where first element
#'   always will be 100
#' @keywords internal
calculatePercentageAvgDiffInMat <- function(targetMat) {
  
  calculatePercDiff <- function (sampleIndex, mat) {
    mean(mat[, sampleIndex]) * 100 / mean(mat[,"log2"])
  }
  
  percDiffVector <- vapply(
    seq_len(ncol(targetMat)),
    calculatePercDiff,
    0,
    mat=targetMat)
  
  percDiffVector
}

#' Calculates correlation values between replicates for each condition matrix.
#' Finally returns a matrix containing the results for all dataset
#'
#' @param methodlist List containing normalized matrices for each normalization
#'   method
#' @param allReplicateGroups Vector with condition groups matching the columns
#'   found in the normalization methods
#' @param sampleGroupsWithReplicates Unique vector with condition groups
#'   present in two or more samples
#' @param corrType Type of correlation (Pearson or Spearman)
#' @return avgCorSum Matrix with column corresponding to normalization
#' approaches and rows corresponding to replicate group
#' @keywords internal
calculateSummarizedCorrelationVector <- function(
    methodlist, allReplicateGroups, sampleGroupsWithReplicates, corrType) {
  
  validCorrTypes <- c("pearson", "spearman")
  if (!corrType %in% validCorrTypes) {
    stop("Unknown correlation type: ",
         corrType,
         " valid are: ",
         paste(validCorrTypes, collapse=", "))
  }
  
  corr_combination_count <- function(allReplicateGroups) {
    replicate_counts <- table(allReplicateGroups)
    sum(vapply(
      replicate_counts,
      function(count) { (count * (count-1)) / 2 },
      0))
  }
  
  avgCorSum <- vapply(
    methodlist,
    calculateCorrSum,
    rep(0, corr_combination_count(allReplicateGroups)),
    allReplicateGroups=allReplicateGroups,
    sampleGroupsWithReplicates=sampleGroupsWithReplicates,
    corrType=corrType
  )
  
  avgCorSum
}

#' Calculates internal correlations for each condition having at least two
#' samples and returns a vector with correlation values corresponding to each
#' condition
#'
#' @param methodData Expression data matrix
#' @param allReplicateGroups Full condition header corresponding to data tables
#'   columns
#' @param sampleGroupsWithReplicates Unique conditions where number of
#'   replicates exceeds one
#' @param corrType Type of correlation (Pearson or Spearman)
#' @return corSums
#' @keywords internal
calculateCorrSum <- function(methodData, allReplicateGroups,
                             sampleGroupsWithReplicates, corrType) {
  methodData <- as.matrix(methodData) # line added
  corSums <- vector()
  for (groupNbr in seq_along(sampleGroupsWithReplicates)) {
    
    specificReplicateVals <- as.matrix(
      methodData[, which(allReplicateGroups == sampleGroupsWithReplicates[groupNbr])])
    class(specificReplicateVals) <- "numeric"
    corVals <- stats::cor(
      specificReplicateVals ,
      use="pairwise.complete.obs",
      method=corrType)
    
    for (index in seq_len(ncol(specificReplicateVals) - 1)) {
      corSums <- c(
        corSums,
        corVals[index, -(seq_len(index)), drop="FALSE"]
      )
    }
  }
  
  corSums
}

calculate_expected_logFC <- function(se, de_results){
  coldata <- data.table::as.data.table(SummarizedExperiment::colData(se))
  condition <- S4Vectors::metadata(se)$condition
  # retrieve concentrations to calculate ground truth logFC per comparison
  concentration <- S4Vectors::metadata(se)$spike_concentration
  concentrations <- unique(coldata[, c(condition, concentration), with=FALSE])
  comps <- data.table::data.table(Comparison = unique(de_results$Comparison))
  comps$SampleA <- sapply(strsplit(as.character(comps$Comparison),"-"), "[", 1)
  comps$SampleB <- sapply(strsplit(as.character(comps$Comparison),"-"), "[", 2)
  comps <- merge(comps, concentrations, by.x = "SampleA", by.y = condition)
  comps <- comps %>% dplyr::rename("ConcentrationA" = concentration)
  comps <- merge(comps, concentrations, by.x = "SampleB", by.y = condition)
  comps <- comps %>% dplyr::rename("ConcentrationB" = concentration)
  comps$Truth <- log2(comps$ConcentrationA / comps$ConcentrationB)
  de_results <- merge(de_results, comps[, c("Comparison", "Truth")], by="Comparison")
  spike <- S4Vectors::metadata(se)$spike_column
  spike_value <- S4Vectors::metadata(se)$spike_value
  de_results$Spiked <- ifelse(de_results[[spike]] == spike_value, "Spike-In", "BG")
  return(de_results)
}
