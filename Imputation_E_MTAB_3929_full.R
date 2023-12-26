library(ImputeSingle)
setwd("/data/miaozhun/analysis/imputation/E_MTAB_3929")
load("//data/miaozhun/analysis/imputation/E_MTAB_3929/E_MTAB_3929.Rdata")
counts_original <- E_MTAB_3929[rowSums(E_MTAB_3929 != 0) > 0,]
counts_original <- as.matrix(counts_original)

# Run ImputeSingle
for(Kcluster in c(2, 5, 10)){
  for(depth in c(2, 5, 10, 20, 50)){
    outputDir <- paste0("./E_MTAB_3929_ImputeSingle_full_K_", Kcluster, "_d_", depth, "/")
    ImputeSingle(counts = counts_original, Kcluster = Kcluster, outputDir = outputDir, depth = depth, SAVER = TRUE, MAGIC = TRUE, parallel = TRUE)
    gc()
  }
}




