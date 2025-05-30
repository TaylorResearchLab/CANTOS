library(readxl)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(openxlsx)

setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
error_analysis_dir <- file.path(analysis_dir,"error-analysis")

combined_df_5th <- read_excel(paste(error_analysis_dir,"/WHO_5th_Error_Categorization_Output.xlsx",sep=""))
combined_df_5th<-combined_df_5th%>%filter(method_eval==0)
combined_df_5th<-combined_df_5th%>%filter(ground_truth!="NF")

combined_df_5th<-combined_df_5th%>%separate_rows(error_category,sep=";\\s*")



ind<-which(combined_df_5th$method_name=="euclidean_dist_v3")
combined_df_5th$method_name[ind]<-"LTE-3+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_dist_ada2")
combined_df_5th$method_name[ind]<-"ADA-002+Euclidean distance"

ind<-which(combined_df_5th$method_name=="af_v3")
combined_df_5th$method_name[ind]<-"LTE-3+AP"

ind<-which(combined_df_5th$method_name=="af_ada2")
combined_df_5th$method_name[ind]<-"ADA-002+AP"

ind<-which(combined_df_5th$method_name=="kmeans_v3")
combined_df_5th$method_name[ind]<-"LTE-3+K-means"

ind<-which(combined_df_5th$method_name=="kmeans_ada2")
combined_df_5th$method_name[ind]<-"ADA-002+K-means"

ind<-which(combined_df_5th$method_name=="af_cosine")
combined_df_5th$method_name[ind]<-"Cosine+AP"

ind<-which(combined_df_5th$method_name=="af_jw")
combined_df_5th$method_name[ind]<-"Jaro Winkler+AP"

ind<-which(combined_df_5th$method_name=="af_lv")
combined_df_5th$method_name[ind]<-"Levenshtein+AP"

ind<-which(combined_df_5th$method_name=="cosine_match")
combined_df_5th$method_name[ind]<-"Cosine"

ind<-which(combined_df_5th$method_name=="jw_match")
combined_df_5th$method_name[ind]<-"Jaro Winkler"

ind<-which(combined_df_5th$method_name=="lv_match")
combined_df_5th$method_name[ind]<-"Levenshtein"

ind<-which(combined_df_5th$method_name=="euclidean_dist_llama")
combined_df_5th$method_name[ind]<-"Llama3.0+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_dist_biobert")
combined_df_5th$method_name[ind]<-"BioBERT+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_dist_medllama")
combined_df_5th$method_name[ind]<-"MedLlama_13B+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_dist_pubmedbert")
combined_df_5th$method_name[ind]<-"PubMedBERT+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_dist_modernbert")
combined_df_5th$method_name[ind]<-"ModernBERT+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_dist_medllama_7b")
combined_df_5th$method_name[ind]<-"MedLlama2+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_dist_llama_32_3b")
combined_df_5th$method_name[ind]<-"Llama3.2_3B+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_dist_phi4")
combined_df_5th$method_name[ind]<-"Phi-4+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_dist_llama_33_70b")
combined_df_5th$method_name[ind]<-"Llama3.3_70B+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_dist_MiniLM_L6_v2")
combined_df_5th$method_name[ind]<-"all-MiniLM-L6-v2+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_mpnet_base")
combined_df_5th$method_name[ind]<-"all-mpnet-base-v2+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_dist_roberta")
combined_df_5th$method_name[ind]<-"all-roberta-large-v1+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_dist_MiniLM_L12_v2")
combined_df_5th$method_name[ind]<-"all-MiniLM-L12-v2+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_dist_e5_large")
combined_df_5th$method_name[ind]<-"e5-large+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_dist_gtr_t5_large")
combined_df_5th$method_name[ind]<-"gtr-t5-large+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_dist_labse")
combined_df_5th$method_name[ind]<-"LaBSE+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_dist_sciBERT")
combined_df_5th$method_name[ind]<-"sciBERT+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_dist_sapBERT")
combined_df_5th$method_name[ind]<-"sapBERT+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_dist_cohere")
combined_df_5th$method_name[ind]<-"embed-english-v2.0+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_dist_deepseek")
combined_df_5th$method_name[ind]<-"DeepSeek_8B+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_dist_BioGPT")
combined_df_5th$method_name[ind]<-"BioGPT+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_dist_clincalBERT")
combined_df_5th$method_name[ind]<-"ClinicalBERT+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_dist_e5large_v2")
combined_df_5th$method_name[ind]<-"e5-large-v2+Euclidean distance"

ind<-which(combined_df_5th$method_name=="euclidean_dist_nomic")
combined_df_5th$method_name[ind]<-"nomic+Euclidean distance"

# Compute normalized error % per method
error_freq <- combined_df_5th %>%
  group_by(method_name, error_category) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(method_name) %>%
  mutate(percent = 100 * count / sum(count)) %>%
  ungroup()



method_ranking_5th <- c("all-MiniLM-L12-v2+Euclidean distance", "e5-large+Euclidean distance", "all-mpnet-base-v2+Euclidean distance", "LTE-3+Euclidean distance", "all-MiniLM-L6-v2+Euclidean distance",
                        "LTE-3+AP", "LTE-3+K-means", "ADA-002+Euclidean distance", "ADA-002+AP", "ADA-002+K-means", 
                        "nomic+Euclidean distance", "all-roberta-large-v1+Euclidean distance", "e5-large-v2+Euclidean distance",
                        "gtr-t5-large+Euclidean distance", "embed-english-v2.0+Euclidean distance", "sapBERT+Euclidean distance",
                        "BioGPT+Euclidean distance", "ClinicalBERT+Euclidean distance", "LaBSE+Euclidean distance",
                        "MedLlama_13B+Euclidean distance", "Levenshtein", "Levenshtein+AP", "Llama3.3_70B+Euclidean distance", 
                        "Jaro Winkler", "Jaro Winkler+AP", "Cosine", "Cosine+AP", "DeepSeek_8B+Euclidean distance", "PubMedBERT+Euclidean distance",
                        "Llama3.2_3B+Euclidean distance", "BioBERT+Euclidean distance", "Llama3.0+Euclidean distance",
                        "sciBERT+Euclidean distance", "ModernBERT+Euclidean distance", "Phi-4+Euclidean distance", 
                        "MedLlama2+Euclidean distance")

method_ranking_5th<-rev(method_ranking_5th)



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

plt<-ggplot(error_freq, aes(x = factor(method_name, levels = method_ranking_5th),
                       y = percent,
                       fill = error_category)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = my_colors)+
  labs(title = "Error category distribution per method, WHO-5th",
       x = "Method",
       y = "Percentage of Errors",
       fill = "Error Category") +
  geom_bar(stat = "identity", color = "white")+
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))+
  geom_text(aes(label = ifelse(percent > 3, paste0(round(percent,1), "%"), "")),
            position = position_stack(vjust = 0.5),
            color = "black",
            size = 3)+ coord_flip()

plt


# count bar plots

unique_categories <- unique(error_freq$error_category)
pdf(paste(error_analysis_dir,"/SD8_error_counts_by_category_5th.pdf",sep=""), width = 10, height = 6)

for (category in unique_categories) {
  # Filter for the current error category
  subset_df <- error_freq %>% filter(error_category == category)
  
  # Create the plot
  p <- ggplot(subset_df, aes(x = method_name, y = count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = paste("Count of", category, "by Method,WHO-5th"),
         x = "Method Name",
         y = "Error Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  # Print the plot to the current PDF page
  print(p)
}
# Close the PDF device
dev.off()