library(readxl)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)

setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
error_analysis_dir <- file.path(analysis_dir,"error-analysis")

sheet_names_who_5th <- excel_sheets(paste(error_analysis_dir,"/error_category_who_5th.xlsx",sep=""))

all_sheets <- lapply(sheet_names_who_5th, read_excel, path = paste(error_analysis_dir,"/error_category_who_5th.xlsx",sep=""))
names(all_sheets) <- sheet_names_who_5th

for (sheet in sheet_names_who_5th) {
  assign(sheet, read_excel(path = paste(error_analysis_dir,"/error_category_who_5th.xlsx",sep=""), sheet = sheet), envir = .GlobalEnv)
}

for(sheet in sheet_names_who_5th){
  df <- get(sheet)
  df$method_name <- sheet
  colnames(df)[c(5,6)]<-c("method_prediction","method_eval")
  assign(sheet, df, envir = .GlobalEnv)
}

combined_df <- do.call(rbind, lapply(sheet_names_who_5th, get))



combined_df<- combined_df%>%filter(ground_truth!="NF")
combined_df<-combined_df%>%filter(method_eval==0)

ind<-which(combined_df$method_name=="LTE3")
combined_df$method_name[ind]<-"LTE3+Euclidean distance"

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

combined_df_expanded<-combined_df%>%separate_rows(error_category,sep=";\\s*")


# Compute normalized error % per method
error_freq <- combined_df_expanded %>%
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



unique_categories <- unique(error_freq$error_category)

for (category in unique_categories) {
  subset_df <- error_freq %>% filter(error_category == category)
  
  # Create the plot
  p <- ggplot(subset_df, aes(x = method_name, y = count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = paste("Method wise error counts ", category, ", WHO-5th "),
         x = "Method Name",
         y = "Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  print(p)
}

