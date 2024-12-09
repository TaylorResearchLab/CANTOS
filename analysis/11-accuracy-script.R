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

#tumor_sample_df_all<-read.csv(paste(result_dir,"/tumor_sample_df_script10_validation.csv",sep = ""))
#tumor_sample_df_5thed<-read.csv(paste(result_dir_5th,"/tumor_sample_df_script10_5thed_corrected_validated.csv",sep = ""))

tumor_sample_df_all<-read.csv(paste(result_dir,"/tumor_sample_df_script10_11sep.csv",sep = ""))
tumor_sample_df_5thed<-read.csv(paste(result_dir_5th,"/tumor_sample_df_script10_5thed_sep11.csv",sep = ""))




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

tumor_sample_df_5thed<-tumor_sample_df_5thed[,c(1,2,27,28,3:26)]
tumor_sample_df_all<-tumor_sample_df_all[,c(1,2,27,28,3:26)]


write.csv(tumor_sample_df_all,paste(result_dir,"/tumor_sample_df_gt_annotated_all.csv",sep = ""))
write.csv(tumor_sample_df_5thed,paste(result_dir_5th,"/tumor_sample_df_gt_annotated_5th.csv",sep = ""))


