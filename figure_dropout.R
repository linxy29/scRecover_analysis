# Install and load necessary libraries
library(tidyverse)

setwd("~/Documents/Data/scRecover/E-MTAB-3929")

# load result file
load("~/Documents/Data/scRecover/E-MTAB-3929/results_E_MTAB_3929_ImputeSingle_v2.Rdata")

process_data = function(df, name) {
  data <- df %>%
    mutate(Name = name,
      p = as.numeric(str_extract(Name, "(?<=_p_)[0-9.]+")),
      K = as.numeric(str_extract(Name, "(?<=_K_)[0-9]+")),
      d = as.numeric(str_extract(Name, "(?<=_d_)[0-9]+"))) %>%
    mutate(truth = df[4,2]) %>% 
    select(-original, -whether_impute, -MAGIC_original, -MAGIC_filter, -Name, -downsampling)
  return(data[5,])
}

combined_df = mapply(process_data, to_be_saved_obj, names(to_be_saved_obj), SIMPLIFY = FALSE) %>% do.call(rbind, .)



pvals_of_interest = c(0.05, 0.1, 0.15, 0.25, 0.3, 0.35)

processed_df = combined_df %>% 
  filter(p %in% pvals_of_interest) %>% 
  select(truth, everything()) %>% 
  pivot_longer(truth:SAVER_filter, names_to = "Method", values_to = "Value") %>% 
  group_by(p, d, Method) %>%
  summarise(Value = max(Value, na.rm = TRUE)) %>%
  ungroup() %>% 
  mutate(p = str_c("p = ", p),
         Method = str_replace(Method, "_filter", " + scRecover"))

write_csv(processed_df, file = "predicted_dropout_num.csv")

ggplot(processed_df, aes(x = d, y = Value, color = Method)) +
  geom_line() +
  geom_point() +
  facet_wrap(~p, scales = "free") +
  labs(x = "Prediction Depth (d)",
       y = "Predicted dropout number") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("dropout.pdf")





