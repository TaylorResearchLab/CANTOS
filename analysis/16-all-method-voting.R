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


#### Measure for all possible vote
col_select <- c("valid_euclidean_dist_v3","valid_euclidean_dist_ada2","valid_af_ad2","valid_af_v3","valid_kmeans_v3",
                "valid_kmeans_ada2","valid_euclidean_dist_MiniLM_L6_v2","valid_euclidean_mpnet_base",
                "valid_euclidean_dist_e5_large","valid_euclidean_dist_MiniLM_L12_v2","valid_euclidean_dist_nomic")


tumor_all_gt<-read.csv(paste(result_dir,"/tumor_manually_validated_all.csv",sep = ""))
tumor_5thed_gt<-read.csv(paste(result_dir_5th,"/tumor_manually_validated_5th.csv",sep = ""))
tumor_all_gt<-tumor_all_gt[,c(-1)]
tumor_5thed_gt<-tumor_5thed_gt[,c(-1)]

col_select2<-col_select
col_select2 <- gsub("valid_","",col_select2)
col_select2[3]<-"af_ada2"  
tumor_5thed_gt<-tumor_5thed_gt%>%dplyr::select(nct_id,Tumor_Names,ground_truth,ground_truth_val,col_select,col_select2)

tumor_5thed_gt$combined_predictions<-NA
tumor_5thed_gt$valid_combined_predictions<-NA


tumor_all_gt<-tumor_all_gt%>%dplyr::select(nct_id,Tumor_Names,ground_truth,ground_truth_val,col_select,col_select2)

tumor_all_gt$combined_predictions<-NA
tumor_all_gt$valid_combined_predictions<-NA

for(iter in 1:nrow(tumor_5thed_gt)){
  ground_truths<- tumor_5thed_gt$ground_truth_val[iter]
  ground_truths<-unique(unlist(strsplit(ground_truths,";")))
  
  list_std<-as.character(tumor_5thed_gt[iter,16:26])
  
  table_pred<-as.data.frame(table(list_std))
  most_frequent <- as.character(table_pred$Var1[table_pred$Freq== max(table_pred$Freq)])
  
  if(length(most_frequent)==1){
    tumor_5thed_gt$combined_predictions[iter]=most_frequent
  }else{
    tumor_5thed_gt$combined_predictions[iter]=tumor_5thed_gt$euclidean_dist_MiniLM_L12_v2[iter]
  }
  
  if(tumor_5thed_gt$combined_predictions[iter] %in% ground_truths){
    tumor_5thed_gt$valid_combined_predictions[iter]=1
  }else{
    tumor_5thed_gt$valid_combined_predictions[iter]=0
  }
  
  
  
}

print(sum(tumor_5thed_gt$valid_combined_predictions)/nrow(tumor_5thed_gt))




####

for(iter in 1:nrow(tumor_all_gt)){
  ground_truths<- tumor_all_gt$ground_truth_val[iter]
  ground_truths<-unique(unlist(strsplit(ground_truths,";")))
  
  list_std<-as.character(tumor_all_gt[iter,16:26])
  
  table_pred<-as.data.frame(table(list_std))
  most_frequent <- as.character(table_pred$Var1[table_pred$Freq== max(table_pred$Freq)])
  
  if(length(most_frequent)==1){
    tumor_all_gt$combined_predictions[iter]=most_frequent
  }else{
    tumor_all_gt$combined_predictions[iter]=tumor_all_gt$euclidean_dist_v3[iter]
  }
  
  if(tumor_all_gt$combined_predictions[iter] %in% ground_truths){
    tumor_all_gt$valid_combined_predictions[iter]=1
  }else{
    tumor_all_gt$valid_combined_predictions[iter]=0
  }
  
  
  
}

print(sum(tumor_all_gt$valid_combined_predictions)/nrow(tumor_all_gt))







