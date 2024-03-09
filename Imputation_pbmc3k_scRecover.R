library(scRecover)
setwd("~/Documents/Data/scRecover/pbmc3k")
pbmc <- readRDS(file = "./pbmc3k_seuratLable.rds")
counts_original <- pbmc@assays$RNA@counts
counts_original <- as.matrix(counts_original)

# Run ImputeSingle
for(Kcluster in c(2, 5, 10)){
  for(depth in c(10, 20, 50)){
    outputDir <- paste0("./pbmc3k_K_", Kcluster, "_d_", depth, "/")
    scRecover(counts = counts_original, Kcluster = Kcluster, outputDir = outputDir, depth = depth, SAVER = TRUE, MAGIC = TRUE, parallel = TRUE)
    gc()
  }
}
