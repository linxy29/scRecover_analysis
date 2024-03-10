setwd("~/Documents/Data/scRecover/E-MTAB-3929")
load("~/Documents/Data/scRecover/E-MTAB-3929/E_MTAB_3929.Rdata")
counts_original <- E_MTAB_3929[rowSums(E_MTAB_3929 != 0) > 0,]

addTrueZero <- function(path, counts_original){
  print(Sys.time())
  print(paste("Processing files of", path))
  
  scImpute_count <- read.csv(file = paste0(path, "/tempFile/scimpute_count.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  SAVER_count <- read.csv(file = paste0(path, "/tempFile/SAVER_count.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  ## whether_impute is the matrix got from ZINB indicating which item should be imputed. A matrix contains true and false.
  ## whether_impute_iz: for non zero value in the original matrix, turn it to false (True are those imputed items); 
  ## whether_impute_inz: for non zero value in the original matrix, turn it to true (True are all non zero values, either it has values in the original matrix or is imputed);
  true_zero = counts_original != 0
  scImute_truezero = scImpute_count * true_zero
  SAVER_truezero = SAVER_count * true_zero
  tempFileDir <- paste0(path, "/tempFile/")
  write.csv(scImute_truezero, file = paste0(tempFileDir, "scImute_truezero.csv"))
  write.csv(SAVER_truezero, file = paste0(tempFileDir, "SAVER_truezero.csv"))
}

summaryImputationSub <- function(path, counts_original, results){
  print(Sys.time())
  print(paste("Processing files of", path))
  
  counts_down <- read.csv(file = paste0(path, "/raw_data.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  scImpute_truezero <- read.csv(file = paste0(path, "/tempFile/scImute_truezero.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  SAVER_truezero <- read.csv(file = paste0(path, "/tempFile/SAVER_truezero.csv"), stringsAsFactors = FALSE, header = T, row.names = 1)
  
  analysisImputation <- function(counts_original, counts_down, counts_imputation){
    feature <- c("Expressed gene number", "Non-zero number per cell", "Zero number per cell")
    results <- data.frame(original = rep(NA, times = length(feature)), downsampling = rep(NA, times = length(feature)), imputation = rep(NA, times = length(feature)), row.names = feature)
    results["Expressed gene number", ] <- c(sum(rowSums(counts_original != 0) > 0), sum(rowSums(counts_down != 0) > 0), sum(rowSums(counts_imputation != 0) > 0))
    results["Non-zero number per cell", ] <- c(mean(colSums(counts_original != 0)), mean(colSums(counts_down != 0)), mean(colSums(counts_imputation != 0)))
    results["Zero number per cell", ] <- c(mean(colSums(counts_original == 0)), mean(colSums(counts_down == 0)), mean(colSums(counts_imputation == 0)))
    results["Dropout number per cell",] <- c(0, mean(colSums(counts_down == 0 & counts_original != 0)), mean(colSums(counts_imputation == 0 & counts_original != 0)))
    results["Imputed zero to non-zero number",] <- c(NA, NA, mean(colSums(counts_imputation != 0 & counts_down == 0)))
    results["True Positive (Recovered dropout)",] <- c(NA, NA, mean(colSums(counts_imputation != 0 & counts_down == 0 & counts_original != 0)))
    results["False Positive (Wrong imputed true zero)",] <- c(NA, NA, mean(colSums(counts_imputation != 0 & counts_down == 0 & counts_original == 0)))
    results["False Negative (Remained dropout)",] <- c(NA, NA, mean(colSums(counts_imputation == 0 & counts_original != 0)))
    results["True Negative (Remained true zero)",] <- c(NA, NA, mean(colSums(counts_imputation == 0 & counts_original == 0)))
    results["Precision = TP / (TP + FP)",] <- c(NA, NA, results["True Positive (Recovered dropout)",3]/(results["True Positive (Recovered dropout)",3] + results["False Positive (Wrong imputed true zero)",3]))
    results["Recall = TP / (TP + FN)",] <- c(NA, NA, results["True Positive (Recovered dropout)",3]/(results["True Positive (Recovered dropout)",3] + results["False Negative (Remained dropout)",3]))
    results["Accuracy = (TP + TN) / (TP + FP + TN + FN)",] <- c(NA, NA, (results["True Positive (Recovered dropout)",3] + results["True Negative (Remained true zero)",3])/(results["True Positive (Recovered dropout)",3] + results["False Positive (Wrong imputed true zero)",3] + results["False Negative (Remained dropout)",3] + results["True Negative (Remained true zero)",3]))
    results <- round(results, digits = 3)
    return(results)
  }
  results_scImpute_truezero <- analysisImputation(counts_original, counts_down, scImpute_truezero)
  results_SAVER_truezero <- analysisImputation(counts_original, counts_down, SAVER_truezero)
  results$scImpute_truezero <- results_scImpute_truezero[,3]
  results$SAVER_truezero <- results_SAVER_truezero[,3]
  return(results)
}

load("~/Documents/Data/scRecover/E-MTAB-3929/results_E_MTAB_3929_ImputeSingle_v2.Rdata")
#for(percentage in c(0.01, 0.05, 0.1, 0.2, 0.25, 0.3, 0.35)){
for(percentage in c(0.25, 0.3, 0.35)){         ## add new results
  for(Kcluster in c(2, 5, 10)){
    for(depth in c(2, 5, 10, 20, 50)){
      resultsName <- paste0("E_MTAB_3929_ImputeSingle_p_", percentage, "_K_", Kcluster, "_d_", depth)
      path <- paste0("./", resultsName)
      #addTrueZero(path, counts_original)
      to_be_saved_obj[[resultsName]] <- summaryImputationSub(path, counts_original, to_be_saved_obj[[resultsName]])
      #assign(resultsName, results)
      #to_be_saved_obj[[resultsName]] <- results
    }
  }
}

save(to_be_saved_obj, file = "~/Documents/Data/scRecover/E-MTAB-3929/results_E_MTAB_3929_ImputeSingle_v3.Rdata")

############### add new statistics
# Function to calculate and add metrics to the data frame
add_metrics_to_df <- function(df) {
  TP <- df["True Positive (Recovered dropout)", ]
  FP <- df["False Positive (Wrong imputed true zero)", ]
  FN <- df["False Negative (Remained dropout)", ]
  TN <- df["True Negative (Remained true zero)", ]
  
  # Calculate the statistics
  df["Precision",] <- TP / (TP + FP)
  #df["Recall",] <- TP / (TP + FN)
  df["Accuracy",] <- (TP + TN) / (TP + FP + FN + TN)
  #df["F1_Score",] <- 2 * (df["Precision",] * df["Recall",]) / (df["Precision",] + df["Recall",])
  df["Specificity",] <- TN / (TN + FP)
  df["False_Positive_Rate",] <- FP / (FP + TN)
  #df["NPV",] <- TN / (TN + FN)
  df["FDR",] <- FP / (TP + FP)
  
  # Calculate MCC
  #mcc_numerator <- TP * TN - FP * FN
  #mcc_denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  #df["MCC",] <- mcc_numerator / mcc_denominator
  
  return(df)
}

# Iterate over each data frame in the list and add the metrics
for (name in names(to_be_saved_obj)) {
  to_be_saved_obj[[name]] <- add_metrics_to_df(to_be_saved_obj[[name]])
}

print("New statsitics like Precision, Accuracy, Specificity, False_Positive_Rate, FDR are added.")

############### Gather results together
items_to_remove <- grep("p_0.1|p_0.2", names(to_be_saved_obj), value = TRUE)

# Remove these items from the list
filtered_list <- to_be_saved_obj[!names(to_be_saved_obj) %in% items_to_remove]

# Define the metrics names to extract
metrics <- c("Precision", "Accuracy", "Specificity", "False_Positive_Rate", "FDR")

# Create a list to hold each metric's data frames
metrics_data <- setNames(vector("list", length(metrics)), metrics)

# Iterate through each metric and extract the corresponding data across all data frames
for (metric in metrics) {
  # Extract the metric data and the name of the data frame and store in a list
  metrics_data[[metric]] <- mapply(function(df, name) {
    if (metric %in% rownames(df)) {
      return(data.frame(Name = name, Value = df[metric, ], stringsAsFactors = FALSE))
    } else {
      return(data.frame(Name = name, Value = NA, stringsAsFactors = FALSE)) # In case the metric is not found
    }
  }, filtered_list, names(filtered_list), SIMPLIFY = FALSE)
  
  # Combine the list of data frames for the current metric into one data frame
  metrics_data[[metric]] <- do.call(rbind, metrics_data[[metric]])
}

# At this point, 'metrics_data' is a list where each element is a data frame for a specific metric
# Each data frame contains the metric values and the names of the original data frames

# plot 
library(ggplot2)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

# Assuming each data frame in metrics_data list has the same structure as your example for Precision

# Define a function to process each data frame
process <- function(df, metric_name) {
  # 1. Extract p, K, d from the Name column
  data <- df %>%
    mutate(
      p = as.numeric(str_extract(Name, "(?<=_p_)[0-9.]+")),
      K = as.numeric(str_extract(Name, "(?<=_K_)[0-9]+")),
      d = as.numeric(str_extract(Name, "(?<=_d_)[0-9]+"))) %>%
    select(-Name, -Value.original, -Value.downsampling) %>%
    pivot_longer(
      cols = starts_with("Value"),
      names_to = "Method",
      values_to = "Value"
    ) %>%
    drop_na(Value)
  
  # 2. Group by p and d, and select the best K for each metric
  # Determine the function to use based on the metric type
  best_value_function <- case_when(
    metric_name %in% c("Precision", "Accuracy", "Specificity") ~ "max",
    metric_name %in% c("False_Positive_Rate", "FDR") ~ "min",
    TRUE ~ NA_character_  # This should not happen but provides a fallback
  )
  
  # Add a column to indicate the best value row
  data <- data %>%
    group_by(p, d, Method) %>%
    mutate(
      Best_Value = case_when(
        best_value_function == "max" & Value == max(Value, na.rm = TRUE) ~ TRUE,
        best_value_function == "min" & Value == min(Value, na.rm = TRUE) ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    filter(Best_Value) %>%
    filter(d >= 10) %>% 
    mutate(Method = str_remove(Method, "Value.")) %>%
    filter(Method %in% c("scImpute_truezero", "SAVER_truezero")) %>%
    mutate(Method = str_replace(Method, "_filter", " + scRecover")) %>% 
    #filter(p != 0.01) %>% 
    ungroup()
  
  return(data)
}

# Initialize an empty list to store processed data frames
processed_dfs <- list()

# Iterate over each metric, process the data, and store the result with metric_name
for(metric_name in names(metrics_data)) {
  print(paste("Processing", metric_name))
  df <- metrics_data[[metric_name]]
  processed_dfs[[metric_name]] <- process(df, metric_name)
}

metrics_of_interest <- c("Accuracy", "Precision", "Specificity")

# Combine data frames into one, preserving metric names
combined_df <- bind_rows(processed_dfs[metrics_of_interest], .id = "Metric")
write_csv(combined_df, file = "truezero.csv")

combined_df %>% 
  mutate(p = str_c("p = ", p)) %>% 
  ggplot(aes(x = Metric, y = Value, color = Method)) +
  geom_boxplot() +
  #facet_grid(Metric ~ .) +  # Changed to facet_grid for better control
  labs(title = "Metrics Comparison Across Different Scenarios",
       x = "Prediction Depth (d)",
       y = "Metric Value") +
  theme_classic() +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_text(angle = 0))  # Adjust text angle if needed

ggsave("truezero.pdf")


combined_df <- bind_rows(processed_dfs[metrics_of_interest], .id = "Metric")