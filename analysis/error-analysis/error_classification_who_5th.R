library(readxl)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(openxlsx)

WHO_5th_Error_Categorization_Output_Corrected <- read_excel("~/Desktop/MTP_Paper/CT_Embedding_Storage_9thJan2025/Error-Ananlysis/WHO_5th_Error_Categorization_Output_Corrected.xlsx")
WHO_5th_Error_Categorization_Output_Corrected<-WHO_5th_Error_Categorization_Output_Corrected%>%filter(method_eval==0)
WHO_5th_Error_Categorization_Output_Corrected<-WHO_5th_Error_Categorization_Output_Corrected%>%filter(ground_truth!="NF")

combined_df<-WHO_5th_Error_Categorization_Output_Corrected%>%separate_rows(error_category,sep=";\\s*")



ind<-which(combined_df$method_name=="euclidean_dist_v3")
combined_df$method_name[ind]<-"LTE3+Euclidean distance"

ind<-which(combined_df$method_name=="euclidean_dist_ada2")
combined_df$method_name[ind]<-"ADA-002+Euclidean distance"

ind<-which(combined_df$method_name=="LTE-3+AP")
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
combined_df$method_name[ind]<-"Llama3.0+Euclidean distance"

ind<-which(combined_df$method_name=="biobert")
combined_df$method_name[ind]<-"BioBERT+Euclidean distance"

ind<-which(combined_df$method_name=="medllama")
combined_df$method_name[ind]<-"MedLlama_13B+Euclidean distance"

ind<-which(combined_df$method_name=="pubmedbert")
combined_df$method_name[ind]<-"PubMedBERT+Euclidean distance"

ind<-which(combined_df$method_name=="modernbert")
combined_df$method_name[ind]<-"ModernBERT+Euclidean distance"

ind<-which(combined_df$method_name=="medllama_7b")
combined_df$method_name[ind]<-"MedLlama2+Euclidean Distance"

ind<-which(combined_df$method_name=="llama_32_3b")
combined_df$method_name[ind]<-"Llama3.2_3B+Euclidean distance"

ind<-which(combined_df$method_name=="phi4")
combined_df$method_name[ind]<-"Phi-4+Euclidean distance"

ind<-which(combined_df$method_name=="llama_33_70b")
combined_df$method_name[ind]<-"Llama3.3_70B+Euclidean distance"

ind<-which(combined_df$method_name=="MiniLM_L6_v2")
combined_df$method_name[ind]<-"all-MiniLM-L6-v2+Euclidean distance"

ind<-which(combined_df$method_name=="mpnet_base")
combined_df$method_name[ind]<-"all-mpnet-base-v2+Euclidean distance"

ind<-which(combined_df$method_name=="roberta")
combined_df$method_name[ind]<-"all-roberta-large-v1+Euclidean distance"

ind<-which(combined_df$method_name=="MiniLM_L12_v2")
combined_df$method_name[ind]<-"all-MiniLM-L12-v2+Euclidean distance"

ind<-which(combined_df$method_name=="e5_large")
combined_df$method_name[ind]<-"e5-large+Euclidean distance"

ind<-which(combined_df$method_name=="gtr_t5_large")
combined_df$method_name[ind]<-"gtr-t5-large+Euclidean distance"

ind<-which(combined_df$method_name=="labse")
combined_df$method_name[ind]<-"LaBSE+Euclidean distance"

ind<-which(combined_df$method_name=="sciBERT")
combined_df$method_name[ind]<-"SciBERT+Euclidean distance"

ind<-which(combined_df$method_name=="sapBERT")
combined_df$method_name[ind]<-"SapBERT+Euclidean distance"

ind<-which(combined_df$method_name=="cohere")
combined_df$method_name[ind]<-"embed-english-v2.0+Euclidean distance"

ind<-which(combined_df$method_name=="deepseek")
combined_df$method_name[ind]<-"DeepSeek_8B+Euclidean distance"

ind<-which(combined_df$method_name=="BioGPT")
combined_df$method_name[ind]<-"BioGPT+Euclidean distance"

ind<-which(combined_df$method_name=="clincalBERT")
combined_df$method_name[ind]<-"ClinicalBERT+Euclidean distance"

ind<-which(combined_df$method_name=="e5large_v2")
combined_df$method_name[ind]<-"e5-large-v2+Euclidean distance"

ind<-which(combined_df$method_name=="nomic")
combined_df$method_name[ind]<-"nomic+Euclidean distance"

# Compute normalized error % per method
error_freq <- combined_df %>%
  group_by(method_name, error_category) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(method_name) %>%
  mutate(percent = 100 * count / sum(count)) %>%
  ungroup()

# Optionally: reorder methods by total error volume (for cleaner heatmap)
method_order <- error_freq %>%
  group_by(method_name) %>%
  summarise(total = sum(percent)) %>%
  arrange(desc(total)) %>%
  pull(method_name)


my_colors <- c(
  "#1f77b4",  # blue
  "#ff7f0e",  # orange
  "#2ca02c",  # green
  "#7f7f7f",  # gray
  "#9467bd",  # purple
  "#8c564b",  # brown
  "#e377c2",  # pink
  "#d62728",  # red
  "#bcbd22",  # lime
  "#aec7e8",   # light blue (contrast with first blue)
  "#17becf"  # cyan
  
)

ggplot(error_freq, aes(x = factor(method_name, levels = method_order),
                       y = percent,
                       fill = error_category)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = my_colors)+
  labs(title = "Error Category Distribution per Method, WHO-5th",
       x = "Method",
       y = "Percentage of Errors",
       fill = "Error Category") +
  geom_bar(stat = "identity", color = "white")+
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
        panel.grid.major.x = element_blank())+
  geom_text(aes(label = ifelse(percent > 3, paste0(round(percent), "%"), "")),
            position = position_stack(vjust = 0.5),
            color = "black",
            size = 3)+ coord_flip()

