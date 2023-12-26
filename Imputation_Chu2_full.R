library(ImputeSingle)
setwd("/data/miaozhun/analysis/imputation/Chu2")
counts_original <- readRDS(file = "/data/miaozhun/analysis/imputation/Chu2/GSE75748_sc_time_course_ec.rds")
counts_original <- as.matrix(counts_original)

# Run ImputeSingle
for(Kcluster in c(2, 5, 10)){
  for(depth in c(2, 5, 10, 20, 50)){
    outputDir <- paste0("./Chu2_ImputeSingle_full_K_", Kcluster, "_d_", depth, "/")
    ImputeSingle(counts = counts_original, Kcluster = Kcluster, outputDir = outputDir, depth = depth, SAVER = TRUE, MAGIC = TRUE, parallel = TRUE)
    gc()
  }
}




