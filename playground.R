run_DE_test <- function (se, comparisons, ain = NULL, condition = NULL, DE_method = "limma", 
                         covariate = NULL, logFC = TRUE, logFC_up = 1, logFC_down = -1, 
                         p_adj = TRUE, alpha = 0.05, B = 100, K = 500) 
{
  params <- check_DE_parameters(se, ain = ain, condition = condition, 
                                comparisons = comparisons, DE_method = DE_method, covariate = covariate, 
                                logFC = logFC, logFC_up = logFC_up, logFC_down = logFC_down, 
                                p_adj = p_adj, alpha = alpha, B = B, K = K)
  ain <- params[["ain"]]
  condition <- params[["condition"]]
  if ("raw" %in% ain) {
    message("DE Analysis will not be performed on raw data.")
    ain <- ain[ain != "raw"]
  }
  de_res <- NULL
  for (method in c(ain)) {
    if(method == "LoessCyc"){
      browser()
    }
    de_chunk <- run_DE_single(se, method = method, condition = condition, 
                              comparisons = comparisons, DE_method = DE_method, 
                              covariate = covariate, logFC = logFC, logFC_up = logFC_up, 
                              logFC_down = logFC_down, p_adj = p_adj, alpha = alpha, 
                              B = B, K = K)
    print(colnames(de_chunk))
    if (is.null(de_res)) {
      de_res <- de_chunk
    }
    else {
      de_res <- rbind(de_res, de_chunk)
    }
    message(paste0(method, ": DE analysis completed."))
  }
  return(de_res)
}
