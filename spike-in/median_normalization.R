
dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[["log2"]])

# Lisi
sample_med <- apply(dt, 2, stats::median, na.rm = TRUE)
mean_med <- mean(sample_med, na.rm = TRUE) # = target_median
norm_dt <- t(t(dt)/sample_med)
norm_dt <- norm_dt * mean_med
norm_dt <- data.table::as.data.table(norm_dt)

# Vivi

dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[["raw"]])

target_median <- mean(apply(dt, 2, median, na.rm = TRUE)) 
norm_factors <- target_median / apply(dt, 2, median, na.rm = TRUE)
norm_dt_vivi <- sweep(dt, 2, norm_factors, FUN = "*")

norm_dt_vivi <- log2(norm_dt_vivi)

# Compare

norm_dt_NormalyzerDE <- 