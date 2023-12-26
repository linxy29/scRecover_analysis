library(ImputeSingle)
setwd("/data/miaozhun/analysis/imputation/10x_neuron_1k_v3")
counts_original <- readRDS(file = "/data/miaozhun/analysis/imputation/10x_neuron_1k_v3/10x_neuron_1k_v3.rds")
counts_original <- as.matrix(counts_original)

# Run ImputeSingle
for(Kcluster in c(2, 5, 10)){
  for(depth in c(2, 5, 10, 20, 50)){
    if(Kcluster == 2 & (depth == 2 | depth == 5)){
      print(paste0("Jump Kcluster = ", Kcluster, ", depth = ", depth))
      next;
    }
    outputDir <- paste0("./10x_neuron_1k_v3_ImputeSingle_full_K_", Kcluster, "_d_", depth, "/")
    ImputeSingle(counts = counts_original, Kcluster = Kcluster, outputDir = outputDir, depth = depth, SAVER = TRUE, MAGIC = TRUE, parallel = TRUE)
    gc()
  }
}




