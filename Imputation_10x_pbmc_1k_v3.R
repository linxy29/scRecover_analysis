# Prepare data
library(Matrix)
counts <- readMM(file = "F:/data/10x_pbmc_1k_v3/filtered_feature_bc_matrix/matrix.mtx")
geneName <- read.table(file = "F:/data/10x_pbmc_1k_v3/filtered_feature_bc_matrix/features.tsv", header = FALSE, stringsAsFactors = FALSE)
geneName <- geneName[,1]
cellName <- read.table(file = "F:/data/10x_pbmc_1k_v3/filtered_feature_bc_matrix/barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
cellName <- cellName[,1]
row.names(counts) <- geneName
colnames(counts) <- cellName
counts <- counts[rowSums(counts != 0) > 0,]
counts <- as.matrix(counts)
saveRDS(counts, file = "F:/data/10x_pbmc_1k_v3/10x_pbmc_1k_v3.rds")





# Plot imputation results
# PLot Accuracy
load("results_10x_pbmc_1k_v3_ImputeSingle.Rdata")
pdf("F:/Plots/ImputeSingle/10x_pbmc_1k_v3/Accuracy.pdf", width = 14, height = 9)
par(mfrow=c(3,4))
for(Kcluster in c(2, 5, 10)){
  for(percentage in c(0.2, 0.4, 0.6, 0.8)){
    accuracy_whether_impute <- NULL
    accuracy_scImpute_original <- NULL
    accuracy_scImpute_filter <- NULL
    accuracy_SAVER_original <- NULL
    accuracy_SAVER_filter <- NULL
    accuracy_MAGIC_original <- NULL
    accuracy_MAGIC_filter <- NULL
    for(depth in c(2, 5, 10, 20, 50)){
      resultsName <- paste0("`10x_pbmc_1k_v3_ImputeSingle_p_", percentage, "_K_", Kcluster, "_d_", depth, "`")
      results <- eval(parse(text = resultsName))
      accuracy_whether_impute <- c(accuracy_whether_impute, results["Accuracy = (TP + TN) / (TP + FP + TN + FN)", "whether_impute"]*100)
      accuracy_scImpute_original <- c(accuracy_scImpute_original, results["Accuracy = (TP + TN) / (TP + FP + TN + FN)", "scImpute_original"]*100)
      accuracy_scImpute_filter <- c(accuracy_scImpute_filter, results["Accuracy = (TP + TN) / (TP + FP + TN + FN)", "scImpute_filter"]*100)
      accuracy_SAVER_original <- c(accuracy_SAVER_original, results["Accuracy = (TP + TN) / (TP + FP + TN + FN)", "SAVER_original"]*100)
      accuracy_SAVER_filter <- c(accuracy_SAVER_filter, results["Accuracy = (TP + TN) / (TP + FP + TN + FN)", "SAVER_filter"]*100)
      accuracy_MAGIC_original <- c(accuracy_MAGIC_original, results["Accuracy = (TP + TN) / (TP + FP + TN + FN)", "MAGIC_original"]*100)
      accuracy_MAGIC_filter <- c(accuracy_MAGIC_filter, results["Accuracy = (TP + TN) / (TP + FP + TN + FN)", "MAGIC_filter"]*100)
    }
    
    depth <- c(2, 5, 10, 20, 50)
    cex <- 1.8
    lwd <- 2
    allNum <- c(accuracy_SAVER_original, accuracy_scImpute_original, accuracy_SAVER_filter, accuracy_scImpute_filter)
    plot(accuracy_SAVER_original~depth, col="darkolivegreen3", pch=15, lwd=lwd, cex = cex, ylim=c(min(allNum) - 5, max(allNum) + 5), type="b", bty="l", xlab="Depth", ylab="Accuracy / %", xaxt="n")
    axis(1, at=depth, labels=depth)
    lines(accuracy_scImpute_original~depth, col="darkorange", pch=17, lwd=lwd, type="b", cex = cex)
    lines(accuracy_SAVER_filter~depth, col="darkolivegreen3", pch=22, lwd=lwd, type="b", cex = cex)
    lines(accuracy_scImpute_filter~depth, col="darkorange", pch=24, lwd=lwd, type="b", cex = cex)
    title(paste0("percentage = ", percentage, ", Kcluster = ", Kcluster))
  }
}
dev.off()

pdf("F:/Plots/ImputeSingle/10x_pbmc_1k_v3/Accuracy_legend.pdf", width = 6, height = 3)
plot.new()
legend("center", 
       legend = c("scImpute original", "scImpute - ImputeSingle", "SAVER / MAGIC original", "SAVER / MAGIC - ImputeSingle"), 
       col = c("darkorange", "darkorange", "darkolivegreen3", "darkolivegreen3"), 
       pch = c(17, 24, 15 ,22), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.3,
       lwd = 2,
       pt.bg = 'white',
       text.col = "black",
       horiz = F)
dev.off()



# PLot predicted dropout number
load("results_10x_pbmc_1k_v3_ImputeSingle.Rdata")
pdf("F:/Plots/ImputeSingle/10x_pbmc_1k_v3/Dropout.pdf", width = 14, height = 9)
par(mfrow=c(3,4))
for(Kcluster in c(2, 5, 10)){
  for(percentage in c(0.2, 0.4, 0.6, 0.8)){
    dropout_downsampling <- NULL
    dropout_whether_impute <- NULL
    dropout_scImpute_original <- NULL
    dropout_scImpute_filter <- NULL
    dropout_SAVER_original <- NULL
    dropout_SAVER_filter <- NULL
    dropout_MAGIC_original <- NULL
    dropout_MAGIC_filter <- NULL
    for(depth in c(2, 5, 10, 20, 50)){
      resultsName <- paste0("`10x_pbmc_1k_v3_ImputeSingle_p_", percentage, "_K_", Kcluster, "_d_", depth, "`")
      results <- eval(parse(text = resultsName))
      dropout_downsampling <- c(dropout_downsampling, results["Dropout number per cell", "downsampling"])
      dropout_whether_impute <- c(dropout_whether_impute, results["Imputed zero to non-zero number", "whether_impute"])
      dropout_scImpute_original <- c(dropout_scImpute_original, results["Imputed zero to non-zero number", "scImpute_original"])
      dropout_scImpute_filter <- c(dropout_scImpute_filter, results["Imputed zero to non-zero number", "scImpute_filter"])
      dropout_SAVER_original <- c(dropout_SAVER_original, results["Imputed zero to non-zero number", "SAVER_original"])
      dropout_SAVER_filter <- c(dropout_SAVER_filter, results["Imputed zero to non-zero number", "SAVER_filter"])
      dropout_MAGIC_original <- c(dropout_MAGIC_original, results["Imputed zero to non-zero number", "MAGIC_original"])
      dropout_MAGIC_filter <- c(dropout_MAGIC_filter, results["Imputed zero to non-zero number", "MAGIC_filter"])
    }
    
    depth <- c(2, 5, 10, 20, 50)
    cex <- 1.8
    lwd <- 2
    allNum <- c(dropout_downsampling, dropout_SAVER_original, dropout_scImpute_original, dropout_SAVER_filter, dropout_scImpute_filter)
    plot(dropout_downsampling~depth, col="royalblue", pch=16, lwd=lwd, cex = cex, ylim=c(min(allNum) - 10, max(allNum) + 10), type="b", bty="l", xlab="Depth", ylab="Predicted dropout number", xaxt="n")
    axis(1, at=depth, labels=depth)
    lines(dropout_SAVER_original~depth, col="darkolivegreen3", pch=15, lwd=lwd, type="b", cex = cex)
    lines(dropout_scImpute_original~depth, col="darkorange", pch=17, lwd=lwd, type="b", cex = cex)
    lines(dropout_SAVER_filter~depth, col="darkolivegreen3", pch=22, lwd=lwd, type="b", cex = cex)
    lines(dropout_scImpute_filter~depth, col="darkorange", pch=24, lwd=lwd, type="b", cex = cex)
    title(paste0("percentage = ", percentage, ", Kcluster = ", Kcluster))
  }
}
dev.off()

pdf("F:/Plots/ImputeSingle/10x_pbmc_1k_v3/Dropout_legend.pdf", width = 6, height = 4)
plot.new()
legend("center", 
       legend = c("True dropout number", "scImpute original", "scImpute - ImputeSingle", "SAVER / MAGIC original", "SAVER / MAGIC - ImputeSingle"), 
       col = c("royalblue", "darkorange", "darkorange", "darkolivegreen3", "darkolivegreen3"), 
       pch = c(16, 17, 24, 15 ,22), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.3,
       lwd = 2,
       pt.bg = 'white',
       text.col = "black",
       horiz = F)
dev.off()





