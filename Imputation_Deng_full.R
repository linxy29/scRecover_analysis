library(ImputeSingle)
setwd("/data/miaozhun/analysis/imputation/Deng")
counts_original <- readRDS(file = "/data/miaozhun/analysis/imputation/Deng/readCounts_GSE45719.rds")
counts_original <- as.matrix(counts_original)

# Run ImputeSingle
for(Kcluster in c(2, 5, 10)){
  for(depth in c(2, 5, 10, 20, 50)){
    outputDir <- paste0("./Deng_ImputeSingle_full_K_", Kcluster, "_d_", depth, "/")
    ImputeSingle(counts = counts_original, Kcluster = Kcluster, outputDir = outputDir, depth = depth, SAVER = TRUE, MAGIC = TRUE, parallel = TRUE)
    gc()
  }
}




