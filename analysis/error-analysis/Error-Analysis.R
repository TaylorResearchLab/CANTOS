library(readxl)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(viridis)

sheet_names <- excel_sheets("categorized_tumor_errors_v7.xlsx")

# Read all sheets into a named list of data frames
all_sheets <- lapply(sheet_names, read_excel, path = "categorized_tumor_errors_v7.xlsx")
names(all_sheets) <- sheet_names

for (sheet in sheet_names) {
  assign(sheet, read_excel(path = "categorized_tumor_errors_v7.xlsx", sheet = sheet), envir = .GlobalEnv)
}

for(sheet in sheet_names){
  df <- get(sheet)
  df$method_name <- sheet
  colnames(df)[c(5,6)]<-c("method_prediction","method_eval")
  assign(sheet, df, envir = .GlobalEnv)
}

combined_df <- do.call(rbind, lapply(sheet_names, get))

combined_df<- combined_df%>%filter(ground_truth!="NF")
combined_df<-combined_df%>%filter(method_eval==0)

ind<-which(combined_df$method_name=="LTE3")
combined_df$method_name[ind]<-"LTE3+Euclidean Distance"

ind<-which(combined_df$method_name=="ada2")
combined_df$method_name[ind]<-"ADA-002+Euclidean distance"

ind<-which(combined_df$method_name=="LTE3_AP")
combined_df$method_name[ind]<-"LTE-3+AP"

ind<-which(combined_df$method_name=="ada2_ap")
combined_df$method_name[ind]<-"ADA-002+AP"

ind<-which(combined_df$method_name=="LTE3_kmeans")
combined_df$method_name[ind]<-"LTE-3+K-means"

ind<-which(combined_df$method_name=="ada2_kmeans")
combined_df$method_name[ind]<-"ADA-002+K-means"

ind<-which(combined_df$method_name=="cosine_ap")
combined_df$method_name[ind]<-"Cosine+AP"

ind<-which(combined_df$method_name=="jw_ap")
combined_df$method_name[ind]<-"Jaro Winkler+AP"

ind<-which(combined_df$method_name=="lv_ap")
combined_df$method_name[ind]<-"Levenshtein+AP"

ind<-which(combined_df$method_name=="cosine")
combined_df$method_name[ind]<-"Cosine"

ind<-which(combined_df$method_name=="jw")
combined_df$method_name[ind]<-"Jaro Winkler"

ind<-which(combined_df$method_name=="lv")
combined_df$method_name[ind]<-"Levenshtein"

ind<-which(combined_df$method_name=="llama")
combined_df$method_name[ind]<-"Llama3.0+Euclidean Distance"

ind<-which(combined_df$method_name=="biobert")
combined_df$method_name[ind]<-"BioBERT+Euclidean Distance"

ind<-which(combined_df$method_name=="medllama")
combined_df$method_name[ind]<-"MedLlama_13B+Euclidean Distance"

ind<-which(combined_df$method_name=="pubmedbert")
combined_df$method_name[ind]<-"PubMedBERT+Euclidean Distance"

ind<-which(combined_df$method_name=="modernbert")
combined_df$method_name[ind]<-"ModernBERT+Euclidean Distance"

ind<-which(combined_df$method_name=="medllama_7b")
combined_df$method_name[ind]<-"MedLlama2+Euclidean Distance"

ind<-which(combined_df$method_name=="llama_32_3b")
combined_df$method_name[ind]<-"Llama3.2_3B+Euclidean Distance"

ind<-which(combined_df$method_name=="phi4")
combined_df$method_name[ind]<-"Phi-4+Euclidean Distance"

ind<-which(combined_df$method_name=="llama_33_70b")
combined_df$method_name[ind]<-"Llama3.3_70B+Euclidean Distance"

ind<-which(combined_df$method_name=="MiniLM_L6_v2")
combined_df$method_name[ind]<-"all-MiniLM-L6-v2+Euclidean Distance"

ind<-which(combined_df$method_name=="mpnet_base")
combined_df$method_name[ind]<-"all-mpnet-base-v2+Euclidean Distance"

ind<-which(combined_df$method_name=="roberta")
combined_df$method_name[ind]<-"all-roberta-large-v1+Euclidean Distance"

ind<-which(combined_df$method_name=="MiniLM_L12_v2")
combined_df$method_name[ind]<-"all-MiniLM-L12-v2+Euclidean Distance"

ind<-which(combined_df$method_name=="e5_large")
combined_df$method_name[ind]<-"e5-large+Euclidean Distance"

ind<-which(combined_df$method_name=="gtr_t5_large")
combined_df$method_name[ind]<-"gtr-t5-large+Euclidean Distance"

ind<-which(combined_df$method_name=="labse")
combined_df$method_name[ind]<-"LaBSE+Euclidean Distance"

ind<-which(combined_df$method_name=="sciBERT")
combined_df$method_name[ind]<-"SciBERT+Euclidean Distance"

ind<-which(combined_df$method_name=="sapBERT")
combined_df$method_name[ind]<-"SapBERT+Euclidean Distance"

ind<-which(combined_df$method_name=="cohere")
combined_df$method_name[ind]<-"embed-english-v2.0+Euclidean Distance"

ind<-which(combined_df$method_name=="deepseek")
combined_df$method_name[ind]<-"DeepSeek_8B+Euclidean Distance"

ind<-which(combined_df$method_name=="BioGPT")
combined_df$method_name[ind]<-"BioGPT+Euclidean Distance"

ind<-which(combined_df$method_name=="clincalBERT")
combined_df$method_name[ind]<-"ClinicalBERT+Euclidean Distance"

ind<-which(combined_df$method_name=="e5large_v2")
combined_df$method_name[ind]<-"e5-large-v2+Euclidean Distance"

ind<-which(combined_df$method_name=="nomic")
combined_df$method_name[ind]<-"nomic+Euclidean Distance"




# Compute normalized error % per method
error_freq <- combined_df %>%
  group_by(method_name, error_category) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(method_name) %>%
  mutate(percent = 100 * count / sum(count)) %>%
  ungroup()
#sum(error_summary$percent[error_summary$method_name=="BioGPT"])



# Optionally: reorder methods by total error volume (for cleaner heatmap)
method_order <- error_freq %>%
  group_by(method_name) %>%
  summarise(total = sum(percent)) %>%
  arrange(desc(total)) %>%
  pull(method_name)

# Plot heatmap
ggplot(error_freq, aes(x = error_category, y = factor(method_name, levels = method_order), fill = percent)) +
  geom_tile(color = "grey90") +
  scale_fill_viridis(name = "% of Errors", option = "plasma", direction = -1, begin = 0.2, end = 0.9) +
  labs(title = "Heatmap of Error Type Frequency by Method, WHO-5th",
       x = "Error Category",
       y = "Method Name") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        legend.position = "right")




#plot of stacked bar plot
ggplot(error_freq, aes(x = factor(method_name, levels = method_order),
                       y = percent,
                       fill = error_category)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_brewer(palette = "Set2")+
  labs(title = "Error Category Distribution per Method, WHO-5th",
       x = "Method",
       y = "Percentage of Errors",
       fill = "Error Category") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major.x = element_blank())





