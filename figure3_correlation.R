library(scRecover)
library(ggplot2)

####### Calculate correlation
calculateCor <- function(path, counts_original){
  print(Sys.time())
  print(paste("Processing files of", path))
  
  counts_down <- read.csv(file = paste0(path, "/raw_data.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  whether_impute_inz <- readRDS(file = paste0(path, "/tempFile/whether_impute_inz.rds"))
  scImpute_count <- read.csv(file = paste0(path, "/tempFile/scimpute_count.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  SAVER_count <- read.csv(file = paste0(path, "/tempFile/SAVER_count.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  MAGIC_count <- read.csv(file = paste0(path, "/tempFile/MAGIC_count.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  scImpute_filter <- read.csv(file = paste0(path, "/scRecover+scImpute.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  SAVER_filter <- read.csv(file = paste0(path, "/scRecover+SAVER.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  MAGIC_filter <- read.csv(file = paste0(path, "/scRecover+MAGIC.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  
  cor_counts_down <- sapply(1:nrow(counts_original), function(i) {
    cor(unlist(counts_original[i, ]), unlist(counts_down[i, ]), use = "complete.obs")
  })
  cor_scImpute_count <- sapply(1:nrow(counts_original), function(i) {
    cor(unlist(counts_original[i, ]), unlist(scImpute_count[i, ]), use = "complete.obs")
  })
  cor_scImpute_filter <- sapply(1:nrow(counts_original), function(i) {
    cor(unlist(counts_original[i, ]), unlist(scImpute_filter[i, ]), use = "complete.obs")
  })
  cor_SAVER_count <- sapply(1:nrow(counts_original), function(i) {
    cor(unlist(counts_original[i, ]), unlist(SAVER_count[i, ]), use = "complete.obs")
  })
  cor_SAVER_filter <- sapply(1:nrow(counts_original), function(i) {
    cor(unlist(counts_original[i, ]), unlist(SAVER_filter[i, ]), use = "complete.obs")
  })
  #cor_MAGIC_count <- analysisImputation(counts_original, counts_down, MAGIC_count)
  #cor_MAGIC_filter <- analysisImputation(counts_original, counts_down, MAGIC_filter)
  results <- data.frame(cor_counts_down = cor_counts_down, cor_scImpute_count = cor_scImpute_count, cor_scImpute_filter = cor_scImpute_filter, cor_SAVER_count = cor_SAVER_count, cor_SAVER_filter = cor_SAVER_filter)
  rownames(results) = rownames(SAVER_count)
  return(results)
}

setwd("~/Documents/Data/scRecover/E-MTAB-3929")
load("~/Documents/Data/scRecover/E-MTAB-3929/E_MTAB_3929.Rdata")
counts_original <- E_MTAB_3929[rowSums(E_MTAB_3929 != 0) > 0,]
## save the conuts_original to csv file
#write.csv(counts_original, file = "./counts_original.csv")
to_be_saved_obj <- NULL

#for(percentage in c(0.01, 0.05, 0.1, 0.25, 0.3, 0.35)){
for(percentage in c(0.25, 0.3, 0.35)){
  for(Kcluster in c(2, 5, 10)){
    for(depth in c(2, 5, 10, 20, 50)){
      resultsName <- paste0("E_MTAB_3929_ImputeSingle_p_", percentage, "_K_", Kcluster, "_d_", depth)
      path <- paste0("./", resultsName)
      results <- calculateCor(path, counts_original)
      assign(resultsName, results)
      to_be_saved_obj[[resultsName]] <- results
    }
  }
}

#save(to_be_saved_obj, file = "~/Documents/Data/scRecover/E-MTAB-3929/results_E_MTAB_3929_correlation_v2.Rdata")

####### Plot correlation

load("~/Documents/Data/scRecover/E-MTAB-3929/results_E_MTAB_3929_correlation_v2.Rdata")

process_data = function(df, name) {
  data <- df %>%
    rownames_to_column("Gene") %>% 
    mutate(
      Name = name,
      p = as.numeric(str_extract(Name, "(?<=_p_)[0-9.]+")),
      K = as.numeric(str_extract(Name, "(?<=_K_)[0-9]+")),
      d = as.numeric(str_extract(Name, "(?<=_d_)[0-9]+"))) %>%
    select(-Name) %>%
    pivot_longer(
      cols = starts_with("cor"),
      names_to = "Method",
      values_to = "Correlation"
    ) %>%
    mutate(Method = str_remove(Method, "cor_")) %>% 
    drop_na(Correlation)
  return(data)
}

combined_df = mapply(process_data, to_be_saved_obj, names(to_be_saved_obj), SIMPLIFY = FALSE) %>% do.call(rbind, .)

# First, create a summarized dataset with mean values
processed_df = combined_df %>%
  filter(Method != "counts_down") %>%
  group_by(p, d, Method, Gene) %>%
  summarise(Correlation = max(Correlation, na.rm = TRUE)) %>%
  ungroup()

saveRDS(processed_df, "correlation.rds")

summary_values <- processed_df %>% 
  # Calculate mean Correlation for each combination of d and Method
  group_by(p, d, Method) %>%
  summarise(MeanCorrelation = mean(Correlation, na.rm = TRUE),
            MinCorrelation = min(Correlation, na.rm = TRUE),
            MaxCorrelation = max(Correlation, na.rm = TRUE)) %>%
  ungroup()

write_csv(summary_values, file = "summary_correlation.csv")

# Now, plot the data with boxplots and add mean Correlation values as text
processed_df %>% 
  #filter(p %in% p_values_of_interest) %>% 
  mutate(p = str_c("p = ", p)) %>% 
  filter(d >= 10) %>% 
  mutate(d = as.factor(d),
         Method = as.factor(Method)) %>% 
  mutate(Method = str_replace(Method, "_filter", " + scRecover")) %>% 
  ggplot(aes(x = d, y = Correlation, color = Method, group = interaction(d, Method))) +
  geom_boxplot() +
  facet_wrap(~p, scales = "free", ncol = 2) +
  labs(title = paste("Change of correlation for different methods and d values"),
       x = "Prediction Depth (d)",
       y = "Correlation") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figure3.pdf")
