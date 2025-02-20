#This script is used evaluate the accuracy of each clustering+standardization methods. 

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(tidyverse)
  library(stringi)
  library(qdapRegex)
  library(ghql)
  library(readxl)
  library(dbscan)
  library(isotree)
  library(infotheo)
  library(combinat)
  library(ggplot2)
  library(reshape2)
})

# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
intermediate_dir <- file.path(analysis_dir,"intermediate")
intermediate_dir_5th <- file.path(analysis_dir,"intermediate_5th")
result_dir <- file.path(analysis_dir,"results")
result_dir_5th <- file.path(analysis_dir,"results_5th")

source(paste(util_dir,"/compute_mi_pairwise.R",sep = ""))

tumor_all_gt<-read.csv(paste(result_dir,"/tumor_manually_validated_all.csv",sep = ""))
tumor_5thed_gt<-read.csv(paste(result_dir_5th,"/tumor_manually_validated_5th.csv",sep = ""))


tumor_all_gt<-tumor_all_gt[,c(-1)]
tumor_5thed_gt<-tumor_5thed_gt[,c(-1)]

tumor_all_gt<-tumor_all_gt%>%filter(ground_truth %in% c("G","MG","G-Manual"))
tumor_5thed_gt<-tumor_5thed_gt%>%filter(ground_truth %in% c("G","MG","G-Manual"))

data_5th<- tumor_5thed_gt[,seq(6,76,2)]
data_all<- tumor_all_gt[,seq(6,76,2)]

col_select <- c("valid_euclidean_dist_v3","valid_euclidean_dist_ada2","valid_af_ad2","valid_af_v3","valid_kmeans_v3",
                "valid_kmeans_ada2","valid_euclidean_dist_MiniLM_L6_v2","valid_euclidean_mpnet_base",
                "valid_euclidean_dist_e5_large","valid_euclidean_dist_MiniLM_L12_v2","valid_euclidean_dist_nomic")


data_5th<-data_5th%>%dplyr::select(col_select)
data_all<-data_all%>%dplyr::select(col_select)


# Compute pairwise Mutual Information
MI_mat_5th <-compute_mi_pairwise(data_5th)
MI_mat_all <-compute_mi_pairwise(data_all)

method_names<-c("LTE3+Euclidean_dist","ADA002+Euclidean_dist","ADA002+AP",                      
                "LTE3+ADA002","LTE3+kmeans","ADA002+kmeans",                 
                "all_MiniLM_L6_v2+euclidean_dist",  "all_mpnet_base+euclidean_dist","e5_large+euclidean_dist",     
                "all_MiniLM_L12_v2+euclidean_dist", "nomic+euclidean_dist")

colnames(MI_mat_5th)<-method_names
colnames(MI_mat_all)<-method_names
rownames(MI_mat_5th)<-method_names
rownames(MI_mat_all)<-method_names



which(MI_mat_5th==min(MI_mat_5th),arr.ind = TRUE)
which(MI_mat_all==min(MI_mat_all),arr.ind = TRUE)  
method_combinations <- combn(method_names, 3, simplify = FALSE)

compute_avg_mi <- function(methods, mi_matrix) {
  sub_matrix <- mi_matrix[methods, methods]  # Extract sub-matrix for selected methods
  avg_mi <- mean(sub_matrix[lower.tri(sub_matrix)])  # Average off-diagonal mutual information
  return(avg_mi)
}

mi_scores_5th <- sapply(method_combinations, function(combo) compute_avg_mi(combo, MI_mat_5th))
mi_scores_all <- sapply(method_combinations, function(combo) compute_avg_mi(combo, MI_mat_all))

collapsed_colnames <- sapply(method_combinations, function(x) paste(x, collapse = "+"))


mi_scores_5th<-as.data.frame(cbind(collapsed_colnames,mi_scores_5th))
mi_scores_all<-as.data.frame(cbind(collapsed_colnames,mi_scores_all))

colnames(mi_scores_5th)<- c("methods","mi_score")
colnames(mi_scores_all)<- c("methods","mi_score")

mi_scores_5th<-mi_scores_5th[order(mi_scores_5th$mi_score),]
mi_scores_all<-mi_scores_all[order(mi_scores_all$mi_score),]
rownames(mi_scores_5th)<-NULL
rownames(mi_scores_all)<-NULL

print(paste("Best diverse combination (lowest mutual information):", mi_scores_5th$methods[1]))
print(paste("Best diverse combination (lowest mutual information):", mi_scores_all$methods[1]))



################
combined_pred_5th<-tumor_5thed_gt%>%dplyr::select(nct_id,Tumor_Names,ground_truth,ground_truth_val,af_ada2,
                                                                     valid_af_ad2,kmeans_v3,valid_kmeans_v3,
                                                  euclidean_dist_MiniLM_L12_v2,valid_euclidean_dist_MiniLM_L12_v2)



combined_pred_5th$combined_predictions<-NA
combined_pred_5th$valid_combined_predictions<-NA


for(iter in 1:nrow(combined_pred_5th)){
  ground_truths<- combined_pred_5th$ground_truth_val[iter]
  ground_truths<-unique(unlist(strsplit(ground_truths,";")))
  
  kmeans_v3_pred<- combined_pred_5th$kmeans_v3[iter]
  af_ada2_pred<-combined_pred_5th$af_ada2[iter]
  minilm_pred <-combined_pred_5th$euclidean_dist_MiniLM_L12_v2[iter]
  
  table_pred<-as.data.frame(table(c(kmeans_v3_pred,af_ada2_pred,minilm_pred)))
  most_frequent <- as.character(table_pred$Var1[table_pred$Freq== max(table_pred$Freq)])
  
  if(length(most_frequent)==1){
    combined_pred_5th$combined_predictions[iter]=most_frequent
  }else{
    combined_pred_5th$combined_predictions[iter]=combined_pred_5th$euclidean_dist_MiniLM_L12_v2[iter]
  }
  
  if(combined_pred_5th$combined_predictions[iter] %in% ground_truths){
    combined_pred_5th$valid_combined_predictions[iter]=1
  }else{
    combined_pred_5th$valid_combined_predictions[iter]=0
  }
  
}

print(sum(combined_pred_5th$valid_combined_predictions)/nrow(combined_pred_5th))

################################

combined_pred_all<-tumor_all_gt%>%dplyr::select(nct_id,Tumor_Names,ground_truth,ground_truth_val,kmeans_v3,
                                                valid_kmeans_v3,kmeans_ada2, valid_kmeans_ada2,euclidean_dist_MiniLM_L12_v2,
                                                valid_euclidean_dist_MiniLM_L12_v2)

combined_pred_all$combined_predictions<-NA
combined_pred_all$valid_combined_predictions<-NA
for(iter in 1:nrow(combined_pred_all)){
  ground_truths<- combined_pred_all$ground_truth_val[iter]
  ground_truths<-unique(unlist(strsplit(ground_truths,";")))
  
  kmeans_v3_pred<- combined_pred_all$kmeans_v3[iter]
  kmeans_ada2_pred<-combined_pred_all$kmeans_ada2[iter]
  minilm_pred <-combined_pred_all$euclidean_dist_MiniLM_L12_v2[iter]
  
  table_pred<-as.data.frame(table(c(kmeans_v3_pred,kmeans_ada2_pred,minilm_pred)))
  most_frequent <- as.character(table_pred$Var1[table_pred$Freq== max(table_pred$Freq)])
  
  if(length(most_frequent)==1){
    combined_pred_all$combined_predictions[iter]=most_frequent
  }else{
    combined_pred_all$combined_predictions[iter]=combined_pred_all$euclidean_dist_MiniLM_L12_v2[iter]
  }
  
  if(combined_pred_all$combined_predictions[iter] %in% ground_truths){
    combined_pred_all$valid_combined_predictions[iter]=1
  }else{
    combined_pred_all$valid_combined_predictions[iter]=0
  }
  
}

print(sum(combined_pred_all$valid_combined_predictions)/nrow(combined_pred_all))


# Convert MI matrix to long format for visualization
mi_melted_all <- melt(MI_mat_all)

# Plot heatmap
ggplot(mi_melted_all, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c()+
  theme_minimal() +
  labs(title = "Heatmap of Pairwise Mutual Information for WHO all editions", 
       x = "Methods", y = "Methods", fill = "MI")+
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),  # Reduce x-axis tick size and rotate
    axis.text.y = element_text(size = 8)  # Reduce y-axis tick size
  )