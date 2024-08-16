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
  
})

# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
intermediate_dir <- file.path(analysis_dir,"intermediate")
result_dir <- file.path(analysis_dir,"results")
intermediate_dir_5th <- file.path(analysis_dir,"intermediate_5th")
result_dir_5th <- file.path(analysis_dir,"results_5th")

tumor_sample_df_all<-read.csv(paste(result_dir,"/tumor_sample_df_script10_validation.csv",sep = ""))
tumor_sample_df_5thed<-read.csv(paste(result_dir_5th,"/tumor_sample_df_script10_validation_5thed.csv",sep = ""))






# tumor_sample_df<-tumor_sample_df %>% filter(!is.na(valid_af_v3))
# 
# tumor_sample_df<-tumor_sample_df %>% filter(valid_af_v3==1|valid_af_ad2==1|valid_kmeans_v3==1|
#                                               valid_kmeans_ad2==1| valid_af_cosine==1| valid_af_jw==1|
#                                               valid_af_lv==1|valid_euclidean_dist_v3==1| valid_euclidean_dist_ada2==1|
#                                               valid_cosine_match==1|valid_jw_match==1|valid_jw_match==1)

tumor_sample_df_all$ground_truth <- NA
tumor_sample_df_5thed$ground_truth_val <- NA

for (iter in 1:dim(tumor_sample_df_all)[1]){
  ind<- which(tumor_sample_df_all[iter,seq(4,27,2)]==1)
  actual_ind<- seq(4,27,2)[ind]-1
  tumor_names <- unique(unlist(tumor_sample_df_all[iter,actual_ind]))
  if(length(ind)>0){
  if(length(tumor_names)>1){
    tumor_sample_df_all$ground_truth[iter]<-"MG"
    tumor_sample_df_all$ground_truth_val[iter]<-paste(tumor_names,collapse = ";")
    
  }else{
    tumor_sample_df_all$ground_truth[iter]<-"G"
    tumor_sample_df_all$ground_truth_val[iter]<-paste(tumor_names,collapse = ";")
  }
  }else{
    tumor_sample_df_all$ground_truth[iter]<-"U"
    tumor_sample_df_all$ground_truth_val[iter]<-"Unknown"
    }
}

tumor_sample_df_all<- tumor_sample_df_all[order(tumor_sample_df_all$ground_truth),]
rownames(tumor_sample_df_all)<-NULL


for (iter in 1:dim(tumor_sample_df_5thed)[1]){
  ind<- which(tumor_sample_df_5thed[iter,seq(4,27,2)]==1)
  actual_ind<- seq(4,27,2)[ind]-1
  tumor_names <- unique(unlist(tumor_sample_df_5thed[iter,actual_ind]))
  if(length(ind)>0){
    if(length(tumor_names)>1){
      tumor_sample_df_5thed$ground_truth[iter]<-"MG"
      tumor_sample_df_5thed$ground_truth_val[iter]<-paste(tumor_names,collapse = ";")
      
    }else{
      tumor_sample_df_5thed$ground_truth[iter]<-"G"
      tumor_sample_df_5thed$ground_truth_val[iter]<-paste(tumor_names,collapse = ";")
    }
  }else{
    tumor_sample_df_5thed$ground_truth[iter]<-"U"
    tumor_sample_df_5thed$ground_truth_val[iter]<-"Unknown"
  }
}
tumor_sample_df_5thed<- tumor_sample_df_5thed[order(tumor_sample_df_5thed$ground_truth),]
rownames(tumor_sample_df_5thed)<-NULL


write.csv(tumor_sample_df_all,paste(result_dir,"/tumor_sample_df_ground_truth_all.csv",sep = ""))
write.csv(tumor_sample_df_5thed,paste(result_dir_5th,"/tumor_sample_df_ground_truth_5thed.csv",sep = ""))
