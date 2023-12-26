# Counts sampling function
countsSampling <- function(counts, fraction){
  n <- floor(sum(counts) * fraction)
  readsGet <- sort(sample(1:sum(counts), n))
  cumCounts <-  c(0, cumsum(counts))
  counts_New <- hist(readsGet, breaks = cumCounts, plot=FALSE)$count
  return(counts_New)
}

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

set.seed(999)
load("F:/data/E-MTAB-3929/E_MTAB_3929.Rdata")
E_MTAB_3929 <- E_MTAB_3929[rowSums(E_MTAB_3929 != 0) > 0,]
counts_original <- E_MTAB_3929[,sample(1:ncol(E_MTAB_3929), 1000)]
counts_original <- counts_original[rowSums(counts_original != 0) > 0,]
counts_original <- counts_original[sample(1:nrow(counts_original), 2000),]

# set.seed(999)
# counts_original <- readRDS(file = "F:/data/10x_heart_1k_v3/10x_heart_1k_v3.rds")
# counts_original <- counts_original[,sample(1:ncol(counts_original), 1001)]
# counts_original <- counts_original[rowSums(counts_original != 0) > 0,]
# counts_original <- counts_original[sample(1:nrow(counts_original), 2000),]
# counts_original <- counts_original[,sample(1:ncol(counts_original), 1000)]
# counts_original <- as.matrix(counts_original)

percentage <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3)

cor.dat <- vector("list", 2)
names(cor.dat) <- c("gene", "cell")
for (i in 1:2) {
  cor.dat[[i]] <- vector("list", length(percentage))
  names(cor.dat[[i]]) <- percentage
}

cmd.g2g <- rep(NA, length(percentage))
cmd.c2c <- rep(NA, length(percentage))
names(cmd.g2g) <- percentage
names(cmd.c2c) <- percentage

counts_full <- counts_original
for(i in 1:length(percentage)){
  cat("\r",paste0("p = ", percentage[i]))
  counts_down <- counts_full
  for(j in 1:ncol(counts_full))
    counts_down[,j] <- countsSampling(counts_full[,j], percentage[i])
  
  counts_original <- counts_full[rowSums(counts_down != 0) > 0, colSums(counts_down != 0) > 0]
  counts_down <- counts_down[rowSums(counts_down != 0) > 0, colSums(counts_down != 0) > 0]
  
  cor.dat[["gene"]][[i]] <- mapply(cor, as.data.frame(t(counts_original)), as.data.frame(t(counts_down)))
  cor.dat[["cell"]][[i]] <- mapply(cor, as.data.frame(counts_original), as.data.frame(counts_down))

  norm.counts_original <- normalizeData(as.matrix(counts_original))
  norm.counts_down <- normalizeData(as.matrix(counts_down))

  cor.g2g_counts_original <- cor(t(norm.counts_original))
  cor.g2g_counts_down <- cor(t(norm.counts_down))
  cmd.g2g[i] <- calc_cmd(cor.g2g_counts_original, cor.g2g_counts_down)

  cor.c2c_counts_original <- cor(norm.counts_original)
  cor.c2c_counts_down <- cor(norm.counts_down)
  cmd.c2c[i] <- calc_cmd(cor.c2c_counts_original, cor.c2c_counts_down)
  
}

plotFolder <- "F:/Plots/ImputeSingle/Correlation Plot/test/"

png(paste0(plotFolder, "Gene and cell correlation.png"), width = 600, height = 350)
par(mfrow=c(1,2))
boxplot(cor.dat[["gene"]], xlab = "Downsampling efficiency", ylab = "Gene correlation with reference", xaxt = "n")
text(x = seq_along(cor.dat[["gene"]]) + 0.5, y = par("usr")[3], labels = names(cor.dat[["gene"]]), srt = 45, adj = 1.3, xpd = TRUE)
boxplot(cor.dat[["cell"]], xlab = "Downsampling efficiency", ylab = "Cell correlation with reference", xaxt = "n")
text(x = seq_along(cor.dat[["gene"]]) + 0.5, y = par("usr")[3], labels = names(cor.dat[["gene"]]), srt = 45, adj = 1.3, xpd = TRUE)
dev.off()

png(paste0(plotFolder, "Gene-to-gene and cell-to-cell CMD.png"), width = 600, height = 350)
par(mfrow=c(1,2))
x <- barplot(cmd.g2g, xlab = "Downsampling efficiency", ylab = "Gene-to-gene CMD", xaxt = "n", ylim = c(0, 1.1*max(cmd.g2g)))
text(x = x + 0.5, y = par("usr")[3], labels = names(cmd.g2g), srt = 45, adj = 1.3, xpd = TRUE)
x <- barplot(cmd.c2c, xlab = "Downsampling efficiency", ylab = "Cell-to-cell CMD", xaxt = "n")
text(x = x + 0.5, y = par("usr")[3], labels = names(cmd.c2c), srt = 45, adj = 1.3, xpd = TRUE)
dev.off()











