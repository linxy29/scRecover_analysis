# Plot Imputation Correlation
library(RColorBrewer)
corImputation <- function(path, plotFolder, counts_original){
  print(Sys.time())
  print(paste("Processing files of", path))
  dir.create(plotFolder, showWarnings = FALSE)
  
  print(paste(Sys.time(), "Reading files..."))
  counts_down <- read.csv(file = paste0(path, "raw_data.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  whether_impute_inz <- readRDS(file = paste0(path, "tempFile/whether_impute_inz.rds"))
  scImpute_count <- read.csv(file = paste0(path, "tempFile/scimpute_count.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  SAVER_count <- read.csv(file = paste0(path, "tempFile/SAVER_count.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  MAGIC_count <- read.csv(file = paste0(path, "tempFile/MAGIC_count.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  scImpute_filter <- read.csv(file = paste0(path, "scImpute_filter.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  SAVER_filter <- read.csv(file = paste0(path, "SAVER_filter.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  MAGIC_filter <- read.csv(file = paste0(path, "MAGIC_filter.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  
  index <- which(rowSums(counts_down != 0) > 0)
  counts_original <- counts_original[index,]
  counts_down <- counts_down[index,]
  whether_impute_inz <- whether_impute_inz[index,]
  scImpute_count <- scImpute_count[index,]
  SAVER_count <- SAVER_count[index,]
  MAGIC_count <- MAGIC_count[index,]
  scImpute_filter <- scImpute_filter[index,]
  SAVER_filter <- SAVER_filter[index,]
  MAGIC_filter <- MAGIC_filter[index,]
  
  # Calculate gene correlation and cell correlation
  print(paste(Sys.time(), "Calculating correlation..."))
  cor.gene_counts_down <- mapply(cor, as.data.frame(t(counts_original)), as.data.frame(t(counts_down)))
  cor.cell_counts_down <- mapply(cor, as.data.frame(counts_original), as.data.frame(counts_down))
  cor.gene_scImpute_count <- mapply(cor, as.data.frame(t(counts_original)), as.data.frame(t(scImpute_count)))
  cor.cell_scImpute_count <- mapply(cor, as.data.frame(counts_original), as.data.frame(scImpute_count))
  cor.gene_scImpute_filter <- mapply(cor, as.data.frame(t(counts_original)), as.data.frame(t(scImpute_filter)))
  cor.cell_scImpute_filter <- mapply(cor, as.data.frame(counts_original), as.data.frame(scImpute_filter))
  cor.gene_SAVER_count <- mapply(cor, as.data.frame(t(counts_original)), as.data.frame(t(SAVER_count)))
  cor.cell_SAVER_count <- mapply(cor, as.data.frame(counts_original), as.data.frame(SAVER_count))
  cor.gene_SAVER_filter <- mapply(cor, as.data.frame(t(counts_original)), as.data.frame(t(SAVER_filter)))
  cor.cell_SAVER_filter <- mapply(cor, as.data.frame(counts_original), as.data.frame(SAVER_filter))
  cor.gene_MAGIC_count <- mapply(cor, as.data.frame(t(counts_original)), as.data.frame(t(MAGIC_count)))
  cor.cell_MAGIC_count <- mapply(cor, as.data.frame(counts_original), as.data.frame(MAGIC_count))
  cor.gene_MAGIC_filter <- mapply(cor, as.data.frame(t(counts_original)), as.data.frame(t(MAGIC_filter)))
  cor.cell_MAGIC_filter <- mapply(cor, as.data.frame(counts_original), as.data.frame(MAGIC_filter))
  
  cor.gene <- list(Observed = cor.gene_counts_down, scImpute = cor.gene_scImpute_count, scImpute_filter = cor.gene_scImpute_filter, SAVER = cor.gene_SAVER_count, SAVER_filter = cor.gene_SAVER_filter, MAGIC = cor.gene_MAGIC_count, MAGIC_filter = cor.gene_MAGIC_filter)
  names(cor.gene) <- c("Observed", "scImpute", "scImpute+", "SAVER", "SAVER+", "MAGIC", "MAGIC+")
  cor.cell <- list(Observed = cor.cell_counts_down, scImpute = cor.cell_scImpute_count, scImpute_filter = cor.cell_scImpute_filter, SAVER = cor.cell_SAVER_count, SAVER_filter = cor.cell_SAVER_filter, MAGIC = cor.cell_MAGIC_count, MAGIC_filter = cor.cell_MAGIC_filter)
  names(cor.cell) <- c("Observed", "scImpute", "scImpute+", "SAVER", "SAVER+", "MAGIC", "MAGIC+")
  saveRDS(cor.gene, file = paste0(plotFolder, "cor.gene.rds"))
  saveRDS(cor.cell, file = paste0(plotFolder, "cor.cell.rds"))
  
  png(paste0(plotFolder, "Gene and cell correlation.png"), width = 600, height = 350)
  par(mfrow=c(1,2))
  boxplot(cor.gene, col=c("white", brewer.pal(6, "Paired")), ylab = "Gene correlation with reference", xaxt = "n")
  text(x = seq_along(cor.gene) + 0.5, y = par("usr")[3], labels = names(cor.gene), srt = 45, adj = 1.3, xpd = TRUE)
  boxplot(cor.cell, col=c("white", brewer.pal(6, "Paired")), ylab = "Cell correlation with reference", xaxt = "n")
  text(x = seq_along(cor.cell) + 0.5, y = par("usr")[3], labels = names(cor.cell), srt = 45, adj = 1.3, xpd = TRUE)
  dev.off()
  
  
  # Calculate gene-to-gene CMD and cell-to-cell CMD
  normalizeData <- function(x, y = x) {
    sf <- colSums(y)/mean(colSums(y))
    return(sweep(x, 2, sf, "/"))
  }
  
  calc_cmd <- function(R1, R2) {
    traceR1R2 <- sum(mapply(crossprod, as.data.frame(R1), as.data.frame(R2)))
    R1.norm <- norm(R1, type = "F")
    R2.norm <- norm(R2, type = "F")
    return(1-traceR1R2/(R1.norm*R2.norm))
  }
  
  print(paste(Sys.time(), "Normalization..."))
  norm.counts_original <- normalizeData(as.matrix(counts_original))
  norm.counts_down <- normalizeData(as.matrix(counts_down))
  norm.scImpute_count <- normalizeData(as.matrix(scImpute_count))
  norm.scImpute_filter <- normalizeData(as.matrix(scImpute_filter))
  norm.SAVER_count <- normalizeData(as.matrix(SAVER_count))
  norm.SAVER_filter <- normalizeData(as.matrix(SAVER_filter))
  norm.MAGIC_count <- normalizeData(as.matrix(MAGIC_count))
  norm.MAGIC_filter <- normalizeData(as.matrix(MAGIC_filter))
  
  print(paste(Sys.time(), "Calculating cor.g2g..."))
  cor.g2g_counts_original <- cor(t(norm.counts_original))
  cor.g2g_counts_down <- cor(t(norm.counts_down))
  cor.g2g_scImpute_count <- cor(t(norm.scImpute_count))
  cor.g2g_scImpute_filter <- cor(t(norm.scImpute_filter))
  cor.g2g_SAVER_count <- cor(t(norm.SAVER_count))
  cor.g2g_SAVER_filter <- cor(t(norm.SAVER_filter))
  cor.g2g_MAGIC_count <- cor(t(norm.MAGIC_count))
  cor.g2g_MAGIC_filter <- cor(t(norm.MAGIC_filter))
  
  print(paste(Sys.time(), "Calculating cmd.g2g..."))
  cmd.g2g_counts_down <- calc_cmd(cor.g2g_counts_original, cor.g2g_counts_down)
  cmd.g2g_scImpute_count <- calc_cmd(cor.g2g_counts_original, cor.g2g_scImpute_count)
  cmd.g2g_scImpute_filter <- calc_cmd(cor.g2g_counts_original, cor.g2g_scImpute_filter)
  cmd.g2g_SAVER_count <- calc_cmd(cor.g2g_counts_original, cor.g2g_SAVER_count)
  cmd.g2g_SAVER_filter <- calc_cmd(cor.g2g_counts_original, cor.g2g_SAVER_filter)
  cmd.g2g_MAGIC_count <- calc_cmd(cor.g2g_counts_original, cor.g2g_MAGIC_count)
  cmd.g2g_MAGIC_filter <- calc_cmd(cor.g2g_counts_original, cor.g2g_MAGIC_filter)
  
  print(paste(Sys.time(), "Calculating cor.c2c..."))
  cor.c2c_counts_original <- cor(norm.counts_original)
  cor.c2c_counts_down <- cor(norm.counts_down)
  cor.c2c_scImpute_count <- cor(norm.scImpute_count)
  cor.c2c_scImpute_filter <- cor(norm.scImpute_filter)
  cor.c2c_SAVER_count <- cor(norm.SAVER_count)
  cor.c2c_SAVER_filter <- cor(norm.SAVER_filter)
  cor.c2c_MAGIC_count <- cor(norm.MAGIC_count)
  cor.c2c_MAGIC_filter <- cor(norm.MAGIC_filter)
  
  print(paste(Sys.time(), "Calculating cmd.c2c..."))
  cmd.c2c_counts_down <- calc_cmd(cor.c2c_counts_original, cor.c2c_counts_down)
  cmd.c2c_scImpute_count <- calc_cmd(cor.c2c_counts_original, cor.c2c_scImpute_count)
  cmd.c2c_scImpute_filter <- calc_cmd(cor.c2c_counts_original, cor.c2c_scImpute_filter)
  cmd.c2c_SAVER_count <- calc_cmd(cor.c2c_counts_original, cor.c2c_SAVER_count)
  cmd.c2c_SAVER_filter <- calc_cmd(cor.c2c_counts_original, cor.c2c_SAVER_filter)
  cmd.c2c_MAGIC_count <- calc_cmd(cor.c2c_counts_original, cor.c2c_MAGIC_count)
  cmd.c2c_MAGIC_filter <- calc_cmd(cor.c2c_counts_original, cor.c2c_MAGIC_filter)
  
  cmd.g2g <- c(cmd.g2g_counts_down, cmd.g2g_scImpute_count, cmd.g2g_scImpute_filter, cmd.g2g_SAVER_count, cmd.g2g_SAVER_filter, cmd.g2g_MAGIC_count, cmd.g2g_MAGIC_filter)
  names(cmd.g2g) <- c("Observed", "scImpute", "scImpute+", "SAVER", "SAVER+", "MAGIC", "MAGIC+")
  cmd.c2c <- c(cmd.c2c_counts_down, cmd.c2c_scImpute_count, cmd.c2c_scImpute_filter, cmd.c2c_SAVER_count, cmd.c2c_SAVER_filter, cmd.c2c_MAGIC_count, cmd.c2c_MAGIC_filter)
  names(cmd.c2c) <- c("Observed", "scImpute", "scImpute+", "SAVER", "SAVER+", "MAGIC", "MAGIC+")
  saveRDS(cmd.g2g, file = paste0(plotFolder, "cmd.g2g.rds"))
  saveRDS(cmd.c2c, file = paste0(plotFolder, "cmd.c2c.rds"))
  
  png(paste0(plotFolder, "Gene-to-gene and cell-to-cell CMD.png"), width = 600, height = 350)
  par(mfrow=c(1,2))
  x <- barplot(cmd.g2g, col = c("white", brewer.pal(6, "Paired")), ylab = "Gene-to-gene CMD", xaxt = "n", ylim = c(0, 1.1*max(cmd.g2g)))
  text(x = x + 0.5, y = par("usr")[3], labels = names(cmd.g2g), srt = 45, adj = 1.3, xpd = TRUE)
  x <- barplot(cmd.c2c, col = c("white", brewer.pal(6, "Paired")), ylab = "Cell-to-cell CMD", xaxt = "n")
  text(x = x + 0.5, y = par("usr")[3], labels = names(cmd.c2c), srt = 45, adj = 1.3, xpd = TRUE)
  dev.off()
  
}



library(Matrix)
setwd("/data/miaozhun/analysis/imputation/10x_pbmc_1k_v3")
dir.create("Correlation_Plot", showWarnings = FALSE)
counts_original <- readRDS(file = "/data/miaozhun/analysis/imputation/10x_pbmc_1k_v3/10x_pbmc_1k_v3.rds")
counts_original <- as.matrix(counts_original)

for(percentage in c(0.2, 0.4, 0.6, 0.8)){
  for(Kcluster in c(2, 5, 10)){
    for(depth in c(2, 5, 10, 20, 50)){
      resultsName <- paste0("10x_pbmc_1k_v3_ImputeSingle_p_", percentage, "_K_", Kcluster, "_d_", depth)
      path <- paste0("./", resultsName, "/")
      plotFolder <- paste0("./Correlation_Plot/", resultsName, "/")
      corImputation(path, plotFolder, counts_original)
    }
  }
}




