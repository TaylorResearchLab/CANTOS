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


#tumor_all_gt<-read.csv(paste(result_dir,"/tumor_sample_df_gt_annotated_all.csv",sep = ""))
#tumor_5thed_gt<-read.csv(paste(result_dir_5th,"/tumor_sample_df_gt_annotated_5th.csv",sep = ""))

tumor_all_gt<-read.csv(paste(result_dir,"/tumor_manually_validated_all.csv",sep = ""))
tumor_5thed_gt<-read.csv(paste(result_dir_5th,"/tumor_manually_validated_5th.csv",sep = ""))

tumor_all_gt<-tumor_all_gt[,c(-1)]
tumor_5thed_gt<-tumor_5thed_gt[,c(-1)]

tumor_all_gt<-tumor_all_gt%>%filter(ground_truth %in% c("G","MG","G-Manual"))
tumor_5thed_gt<-tumor_5thed_gt%>%filter(ground_truth %in% c("G","MG","G-Manual"))

data_5th<- tumor_5thed_gt[,seq(6,76,2)]
data_all<- tumor_all_gt[,seq(6,76,2)]

col_select <- c("valid_euclidean_dist_v3","valid_euclidean_dist_ada2","valid_af_ad2","valid_kmeans_v3",
                "valid_kmeans_ada2","valid_euclidean_dist_MiniLM_L6_v2","valid_euclidean_mpnet_base",
                "valid_euclidean_dist_e5_large","valid_euclidean_dist_MiniLM_L12_v2","valid_euclidean_dist_nomic")


data_5th<-data_5th%>%dplyr::select(col_select)
data_all<-data_all%>%dplyr::select(col_select)

cor_mat_5th = cor(data_5th)
cor_mat_all = cor(data_all)



which(cor_mat_5th==min(cor_mat_5th),arr.ind = TRUE)
which(cor_mat_all==min(cor_mat_all),arr.ind = TRUE)

combined_pred_MiniLM_L12_kmeans_v3_5th<-tumor_5thed_gt%>%dplyr::select(nct_id,Tumor_Names,ground_truth,ground_truth_val,
                                                                       kmeans_v3,valid_kmeans_v3,euclidean_dist_MiniLM_L12_v2,
                                                                       valid_euclidean_dist_MiniLM_L12_v2) 

combined_pred_MiniLM_L12_kmeans_v3_5th$combined_predictions<-NA
combined_pred_MiniLM_L12_kmeans_v3_5th$valid_combined_predictions<-NA

for (iter in 1:nrow(combined_pred_MiniLM_L12_kmeans_v3_5th)){
  if(combined_pred_MiniLM_L12_kmeans_v3_5th$valid_kmeans_v3[iter]==combined_pred_MiniLM_L12_kmeans_v3_5th$valid_euclidean_dist_MiniLM_L12_v2[iter]){
    combined_pred_MiniLM_L12_kmeans_v3_5th$combined_predictions[iter]<-combined_pred_MiniLM_L12_kmeans_v3_5th$euclidean_dist_MiniLM_L12_v2[iter]
  }else if(combined_pred_MiniLM_L12_kmeans_v3_5th$valid_kmeans_v3[iter]!=combined_pred_MiniLM_L12_kmeans_v3_5th$valid_euclidean_dist_MiniLM_L12_v2[iter]){
    if(combined_pred_MiniLM_L12_kmeans_v3_5th$valid_kmeans_v3[iter]==1){
      combined_pred_MiniLM_L12_kmeans_v3_5th$combined_predictions[iter]<-combined_pred_MiniLM_L12_kmeans_v3_5th$kmeans_v3[iter]
    }else if(combined_pred_MiniLM_L12_kmeans_v3_5th$valid_euclidean_dist_MiniLM_L12_v2[iter]==1){
      combined_pred_MiniLM_L12_kmeans_v3_5th$combined_predictions[iter]<-combined_pred_MiniLM_L12_kmeans_v3_5th$euclidean_dist_MiniLM_L12_v2[iter]
    }
  }
  ground_truths<- combined_pred_MiniLM_L12_kmeans_v3_5th$ground_truth_val[iter]
  ground_truths<-unique(unlist(strsplit(ground_truths,";")))
  
  if(combined_pred_MiniLM_L12_kmeans_v3_5th$combined_predictions[iter] %in% ground_truths){
    combined_pred_MiniLM_L12_kmeans_v3_5th$valid_combined_predictions[iter]=1
  }else{
    combined_pred_MiniLM_L12_kmeans_v3_5th$valid_combined_predictions[iter]=0
  }
  
  
}

print(sum(combined_pred_MiniLM_L12_kmeans_v3_5th$valid_combined_predictions)/1044)







######

combined_pred_MiniLM_L12_kmeans_v3_all<-tumor_all_gt%>%dplyr::select(nct_id,Tumor_Names,ground_truth,ground_truth_val,
                                                                       kmeans_v3,valid_kmeans_v3,euclidean_dist_MiniLM_L12_v2,
                                                                       valid_euclidean_dist_MiniLM_L12_v2)

combined_pred_MiniLM_L12_kmeans_v3_all$combined_predictions<-NA
combined_pred_MiniLM_L12_kmeans_v3_all$valid_combined_predictions<-NA

for (iter in 1:nrow(combined_pred_MiniLM_L12_kmeans_v3_all)){
  if(combined_pred_MiniLM_L12_kmeans_v3_all$valid_kmeans_v3[iter]==combined_pred_MiniLM_L12_kmeans_v3_all$valid_euclidean_dist_MiniLM_L12_v2[iter]){
    combined_pred_MiniLM_L12_kmeans_v3_all$combined_predictions[iter]<-combined_pred_MiniLM_L12_kmeans_v3_all$euclidean_dist_MiniLM_L12_v2[iter]
  }else if(combined_pred_MiniLM_L12_kmeans_v3_all$valid_kmeans_v3[iter]!=combined_pred_MiniLM_L12_kmeans_v3_all$valid_euclidean_dist_MiniLM_L12_v2[iter]){
    if(combined_pred_MiniLM_L12_kmeans_v3_all$valid_kmeans_v3[iter]==1){
      combined_pred_MiniLM_L12_kmeans_v3_all$combined_predictions[iter]<-combined_pred_MiniLM_L12_kmeans_v3_all$kmeans_v3[iter]
    }else if(combined_pred_MiniLM_L12_kmeans_v3_all$valid_euclidean_dist_MiniLM_L12_v2[iter]==1){
      combined_pred_MiniLM_L12_kmeans_v3_all$combined_predictions[iter]<-combined_pred_MiniLM_L12_kmeans_v3_all$euclidean_dist_MiniLM_L12_v2[iter]
    }
  }
  ground_truths<- combined_pred_MiniLM_L12_kmeans_v3_all$ground_truth_val[iter]
  ground_truths<-unique(unlist(strsplit(ground_truths,";")))
  
  if(combined_pred_MiniLM_L12_kmeans_v3_all$combined_predictions[iter] %in% ground_truths){
    combined_pred_MiniLM_L12_kmeans_v3_all$valid_combined_predictions[iter]=1
  }else{
    combined_pred_MiniLM_L12_kmeans_v3_all$valid_combined_predictions[iter]=0
  }
  
  
}

print(sum(combined_pred_MiniLM_L12_kmeans_v3_all$valid_combined_predictions)/1118)





### Repeat analysis with MI

MI_mat_5th <- matrix(0, nrow = 10, ncol = 10)

# Compute pairwise Mutual Information
for (i in 1:10) {
  for (j in i:10) {
    mi_score <- mutinformation(data_5th[[i]], data_5th[[j]])
    MI_mat_5th[i, j] <- mi_score
    MI_mat_5th[j, i] <- mi_score  # Fill both halves since it's symmetric
  }
}
colnames(MI_mat_5th) <- colnames(data_5th)
rownames(MI_mat_5th) <- colnames(data_5th)


MI_mat_all <- matrix(0, nrow = 10, ncol = 10)

# Compute pairwise Mutual Information
for (i in 1:10) {
  for (j in i:10) {
    mi_score <- mutinformation(data_all[[i]], data_all[[j]])
    MI_mat_all[i, j] <- mi_score
    MI_mat_all[j, i] <- mi_score  # Fill both halves since it's symmetric
  }
}
colnames(MI_mat_all) <- colnames(data_all)
rownames(MI_mat_all) <- colnames(data_all)


which(MI_mat_5th==min(MI_mat_5th),arr.ind = TRUE)
which(MI_mat_all==min(MI_mat_all),arr.ind = TRUE)  

# Same as MiniLM_L12_v2 and Kmean_v3