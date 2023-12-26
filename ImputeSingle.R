# Estimate dropout gene number in a cell
estDropoutNum <- function(sample){
  library(preseqR)
  
  sample <- as.numeric(sample)
  sample <- ceiling(sample)
  histCounts <- cbind(as.numeric(names(table(sample))), as.numeric(unname(table(sample))))
  histCounts <- histCounts[-1,]
  colnames(histCounts) <- c("j", "n_j")
  preseqR <- ztnb.mincount(histCounts)
  geneNumPredict <- round(preseqR(1e30))
  dropoutNum <- geneNumPredict - sum(sample != 0)
}

# Functions for barcode sample complete read counts matrix
complete_matrix <- read.csv("outxls1Ex.csv", row.names = 1, stringsAsFactors = FALSE)
UMI_counts <- get_UMI_counts(complete_matrix)
hist_raw_counts_list <- get_hist_raw_counts(complete_matrix)
hist_RUG_counts_list <- get_hist_RUG_counts(complete_matrix)
save(UMI_counts, hist_raw_counts_list, hist_RUG_counts_list, file = "barcode_scRNAseq.Rdata")

get_UMI_counts <- function(complete_matrix){
  complete_matrix <- as.matrix(complete_matrix)
  system.time(UMI_counts <- structure(vapply(complete_matrix, UMI_string_to_UMI_counts, as.integer(1)), dim=dim(complete_matrix)))
  dimnames(UMI_counts) <- list(row.names(complete_matrix), colnames(complete_matrix))
  return(UMI_counts)
}

UMI_string_to_UMI_counts <- function(UMI_string){
  if(UMI_string == "0")
    return(as.integer(0))
  UMI_elements <- strsplit(UMI_string, split = ";")[[1]]
  return(length(UMI_elements))
}

get_hist_raw_counts <- function(complete_matrix){
  complete_matrix <- as.matrix(complete_matrix)
  raw_counts <- structure(vapply(complete_matrix, UMI_string_to_raw_counts, as.integer(1)), dim=dim(complete_matrix))
  dimnames(raw_counts) <- list(row.names(complete_matrix), colnames(complete_matrix))
  
  hist_raw_counts_list <- list()
  for(i in 1:ncol(raw_counts)){
    raw_counts_one <- raw_counts[,i]
    hist_raw_counts <- cbind(as.numeric(names(table(raw_counts_one))), as.numeric(unname(table(raw_counts_one))))
    hist_raw_counts <- hist_raw_counts[-1,]
    colnames(hist_raw_counts) <- c("j", "n_j")
    hist_raw_counts_list[[i]] <- hist_raw_counts
  }
  names(hist_raw_counts_list) <- colnames(complete_matrix)
  return(hist_raw_counts_list)
}

UMI_string_to_raw_counts <- function(UMI_string){
  if(UMI_string == "0")
    return(as.integer(0))
  UMI_elements <- strsplit(UMI_string, split = ";")[[1]]
  RUG_counts <- as.integer(sapply(strsplit(UMI_elements, split = ":"), "[[", 2))
  raw_counts <- sum(RUG_counts)
  return(raw_counts)
}

get_hist_RUG_counts <- function(complete_matrix){
  hist_RUG_counts_list <- list()
  for(i in 1:ncol(complete_matrix)){
    oneCell <- complete_matrix[,i]
    RUG_counts <- unlist(sapply(oneCell, UMI_string_to_RUG_counts))
    hist_RUG_counts <- cbind(as.numeric(names(table(RUG_counts))), as.numeric(unname(table(RUG_counts))))
    hist_RUG_counts <- hist_RUG_counts[-1,]
    colnames(hist_RUG_counts) <- c("j", "n_j")
    hist_RUG_counts_list[[i]] <- hist_RUG_counts
  }
  names(hist_RUG_counts_list) <- colnames(complete_matrix)
  return(hist_RUG_counts_list)
}

UMI_string_to_RUG_counts <- function(UMI_string){
  if(UMI_string == "0")
    return(as.integer(0))
  UMI_elements <- strsplit(UMI_string, split = ";")[[1]]
  RUG_counts <- as.integer(sapply(strsplit(UMI_elements, split = ":"), "[[", 2))
  return(RUG_counts)
}

# Normalization function
normalization <- function(counts){
  # Filter genes with all zero counts
  counts <- as.matrix(counts)
  counts <- ceiling(counts)
  sampleNum <- ncol(counts)
  if(any(rowSums(counts) == 0))
    message("Removing ", sum(rowSums(counts) == 0), " rows of genes with all zero counts")
  counts_NAZ <- counts[rowSums(counts) != 0,]
  geneNum_NAZ <- nrow(counts_NAZ)
  
  # Normalization
  GEOmean <- rep(NA,geneNum_NAZ)
  for (i in 1:geneNum_NAZ)
  {
    gene_NZ <- counts_NAZ[i,counts_NAZ[i,] > 0]
    GEOmean[i] <- exp(sum(log(gene_NZ), na.rm=TRUE) / length(gene_NZ))
  }
  S <- rep(NA, sampleNum)
  counts_norm <- counts_NAZ
  for (j in 1:sampleNum)
  {
    sample_j <- counts_NAZ[,j]/GEOmean
    S[j] <- median(sample_j[which(sample_j != 0)])
    counts_norm[,j] <- counts_NAZ[,j]/S[j]
  }
  counts_norm <- ceiling(counts_norm)
  return(counts_norm)
}



# MLE of parameters of ZINB
mleZINB <- function(counts_1){
  suppressPackageStartupMessages({
    library(stats)
    library(MASS)
    library(pscl)
    library(bbmle)
    library(gamlss)
  })
  counts_1 <- as.numeric(counts_1)
  if(sum(counts_1 == 0) > 0){
    if(sum(counts_1 == 0) == length(counts_1)){
      theta_1 <- 1
      mu_1 <- 0
      size_1 <- 1
      prob_1 <- size_1/(size_1 + mu_1)
    }else{
      options(show.error.messages = FALSE)
      zinb_try <- try(gamlssML(counts_1, family="ZINBI"), silent=TRUE)
      options(show.error.messages = TRUE)
      if('try-error' %in% class(zinb_try)){
        zinb_try_twice <- try(zeroinfl(formula = counts_1 ~ 1 | 1, dist = "negbin"), silent=TRUE)
        if('try-error' %in% class(zinb_try_twice)){
          # print("MLE of ZINB failed!");
          parameters <- c(NA, NA, NA, NA)
          names(parameters) <- c("theta", "mu", "size", "prob")
          return(parameters)
        }else{
          zinb_1 <- zinb_try_twice
          theta_1 <- plogis(zinb_1$coefficients$zero);names(theta_1) <- NULL
          mu_1 <- exp(zinb_1$coefficients$count);names(mu_1) <- NULL
          size_1 <- zinb_1$theta;names(size_1) <- NULL
          prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
        }
      }else{
        zinb_1 <- zinb_try
        theta_1 <- zinb_1$nu;names(theta_1) <- NULL
        mu_1 <- zinb_1$mu;names(mu_1) <- NULL
        size_1 <- 1/zinb_1$sigma;names(size_1) <- NULL
        prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
      }
    }
  }else{
    op <- options(warn=2)
    nb_try <- try(glm.nb(formula = counts_1 ~ 1), silent=TRUE)
    options(op)
    if('try-error' %in% class(nb_try)){
      nb_try_twice <- try(fitdistr(counts_1, "Negative Binomial"), silent=TRUE)
      if('try-error' %in% class(nb_try_twice)){
        nb_try_again <- try(mle2(counts_1~dnbinom(mu=exp(logmu),size=1/invk), data=data.frame(counts_1), start=list(logmu=0,invk=1), method="L-BFGS-B", lower=c(logmu=-Inf,invk=1e-8)), silent=TRUE)
        if('try-error' %in% class(nb_try_again)){
          nb_try_fourth <- try(glm.nb(formula = counts_1 ~ 1), silent=TRUE)
          if('try-error' %in% class(nb_try_fourth)){
            # print("MLE of NB failed!");
            parameters <- c(NA, NA, NA, NA)
            names(parameters) <- c("theta", "mu", "size", "prob")
            return(parameters)
          }else{
            nb_1 <- nb_try_fourth
            theta_1 <- 0
            mu_1 <- exp(nb_1$coefficients);names(mu_1) <- NULL
            size_1 <- nb_1$theta;names(size_1) <- NULL
            prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
          }
        }else{
          nb_1 <- nb_try_again
          theta_1 <- 0
          mu_1 <- exp(nb_1@coef["logmu"]);names(mu_1) <- NULL
          size_1 <- 1/nb_1@coef["invk"];names(size_1) <- NULL
          prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
        }
      }else{
        nb_1 <- nb_try_twice
        theta_1 <- 0
        mu_1 <- nb_1$estimate["mu"];names(mu_1) <- NULL
        size_1 <- nb_1$estimate["size"];names(size_1) <- NULL
        prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
      }
    }else{
      nb_1 <- nb_try
      theta_1 <- 0
      mu_1 <- exp(nb_1$coefficients);names(mu_1) <- NULL
      size_1 <- nb_1$theta;names(size_1) <- NULL
      prob_1 <- size_1/(size_1 + mu_1);names(prob_1) <- NULL
    }
  }
  parameters <- c(theta_1, mu_1, size_1, prob_1)
  names(parameters) <- c("theta", "mu", "size", "prob")
  return(parameters)
}



# Main
library(stats)
library(MASS)
library(pscl)
library(bbmle)
library(gamlss)
library(preseqR)
library(scImpute)
library(SAVER)
library(Rmagic)
library(BiocParallel)

count_path <- "D:/Research/R_workspace/rawdata_scImpute/Fig4/raw_data.csv"
counts_raw <- read.csv(count_path, header = TRUE, row.names = 1)
# cell_type <- read.csv("D:/Research/R_workspace/rawdata_scImpute/Fig2/cell_type.csv", header = TRUE, row.names = 1)
tempFile_path <-"D:/Research/R_workspace/temp/"
dir.create(tempFile_path, showWarnings = FALSE)

scimpute(# full path to raw count matrix
  count_path = count_path, 
  infile = "csv",           # format of input file
  outfile = "csv",          # format of output file
  out_dir = tempFile_path,  # full path to output directory
  labeled = FALSE,          # cell type labels not available
  drop_thre = 0.5,          # threshold set on dropout probability
  Kcluster = 3,             # 2 cell subpopulations
  labels = NULL,            # Each cell type should have at least two cells for imputation
  ncores = 1)               # number of cores used in parallel computation

library(doParallel)
cl <- makeCluster(4, outfile = "")
registerDoParallel(cl)
counts_SAVER <- saver(counts_raw)
stopCluster(cl)
counts_SAVER <- counts_SAVER$estimate
counts_MAGIC <- run_magic(counts_raw, t = 6)
counts_scImpute <- read.csv(file = paste0(tempFile_path, "scimpute_count.csv"), header = TRUE, row.names = 1)
write.csv(counts_SAVER, file = paste0(tempFile_path, "SAVER_count.csv"))
write.csv(counts_MAGIC, file = paste0(tempFile_path, "MAGIC_count.csv"))

clust <- readRDS(paste0(tempFile_path, "clust.rds"))
names(clust) <- colnames(counts_raw)
nclust <- sum(!is.na(unique(clust)))

counts_raw_used <- counts_raw[,!is.na(clust)]
counts_norm <- normalization(counts_raw_used)
whether_impute <- counts_raw
whether_impute[,] <- FALSE

parallel <- TRUE
BPPARAM <- SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)
register(BPPARAM)

# # Estimate dropout gene number in each cell
dropoutNum <- NULL
if(!parallel){
  for(i in 1:ncol(counts_raw_used)){
    cat("\r",paste0("ImputeSingle is estimating dropout gene number in ", i, " of ", ncol(counts_raw_used), " non-outlier cells"))
    dropoutNum <- c(dropoutNum, estDropoutNum(counts_raw_used[,i]))
  }
  names(dropoutNum) <- colnames(counts_raw_used)
}else{
  message("ImputeSingle is estimating dropout gene number in ", ncol(counts_raw_used), " non-outlier cells")
  dropoutNum <- do.call(c, bplapply(as.data.frame(counts_raw_used), FUN = estDropoutNum, BPPARAM = BPPARAM))
}

# ZINB MLE of each cluster gene by gene
ZINB_parameters_list <- list()
P_dropout_cc_list <- list()
P_dropout_mat <- NULL
for(cc in 1:nclust){
  message("Processing ", cc, " of ", nclust, " cell clusters")
  cells <- names(clust[clust %in% cc])
  counts_norm_cc <- counts_norm[, cells, drop = FALSE]
  
  # ZINB MLE of each cluster gene by gene
  ZINB_parameters <- NULL
  if(!parallel){
    for(i in 1:nrow(counts_norm_cc)){
      cat("\r",paste0("ImputeSingle is analyzing ", i, " of ", nrow(counts_norm), " genes in cluster ", cc))
      ZINB_parameters <- rbind(ZINB_parameters, mleZINB(counts_norm_cc[i,]))
    }
    row.names(ZINB_parameters) <- row.names(counts_norm_cc)
    message("\r")
  }else{
    message("ImputeSingle is analyzing ", nrow(counts_norm), " genes in cluster ", cc)
    ZINB_parameters <- do.call(rbind, bplapply(as.data.frame(t(counts_norm_cc)), FUN = mleZINB, BPPARAM = BPPARAM))
  }
  ZINB_parameters_list[[cc]] <- ZINB_parameters
  
  # Estimate dropout probability of each gene
  P_dropout <- NULL
  for(i in 1:nrow(counts_norm_cc)){
    if(any(is.na(ZINB_parameters[i,])))
      P_dropout <- c(P_dropout, 0)
    else
      P_dropout <- c(P_dropout, (1 - ZINB_parameters[i, "theta"]) * dnbinom(x=0, size=ZINB_parameters[i, "size"], prob=ZINB_parameters[i, "prob"]) / (ZINB_parameters[i, "theta"] + (1 - ZINB_parameters[i, "theta"]) * dnbinom(x=0, size=ZINB_parameters[i, "size"], prob=ZINB_parameters[i, "prob"])))
  }
  names(P_dropout) <- row.names(counts_norm_cc)
  P_dropout_mat <- cbind(P_dropout_mat, P_dropout)
  
  # Get dropout probability for each gene in each cell
  P_dropout_cc <- counts_norm_cc
  P_dropout_cc[,1:ncol(counts_norm_cc)] <- P_dropout
  P_dropout_cc[counts_norm_cc != 0] <- 0
  P_dropout_cc_list[[cc]] <- P_dropout_cc
  
  # Determine the position need to be imputed
  P_dropout_rank <- apply(-P_dropout_cc, 2, rank)
  whether_impute_cc <- sweep(P_dropout_rank, MARGIN = 2, dropoutNum[cells], FUN = "<=")
  whether_impute[row.names(whether_impute_cc), cells] <- whether_impute_cc
}
colnames(P_dropout_mat) <- paste0("CellCluster_", 1:nclust)

whether_impute_iz <- whether_impute
whether_impute_inz <- whether_impute
whether_impute_inz[counts_raw != 0] <- TRUE

counts_scImpute_iz <- counts_raw + counts_scImpute * whether_impute_iz
counts_scImpute_inz <- counts_scImpute * whether_impute_inz

counts_SAVER_iz <- counts_raw + counts_SAVER * whether_impute_iz
counts_SAVER_inz <- counts_SAVER * whether_impute_inz

counts_MAGIC_iz <- counts_raw + counts_MAGIC * whether_impute_iz
counts_MAGIC_inz <- counts_MAGIC * whether_impute_inz

write.csv(counts_scImpute_iz, file = paste0(tempFile_path, "counts_scImpute_iz.csv"))
write.csv(counts_scImpute_inz, file = paste0(tempFile_path, "counts_scImpute_inz.csv"))
write.csv(counts_SAVER_iz, file = paste0(tempFile_path, "counts_SAVER_iz.csv"))
write.csv(counts_SAVER_inz, file = paste0(tempFile_path, "counts_SAVER_inz.csv"))
write.csv(counts_MAGIC_iz, file = paste0(tempFile_path, "counts_MAGIC_iz.csv"))
write.csv(counts_MAGIC_inz, file = paste0(tempFile_path, "counts_MAGIC_inz.csv"))























