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


tumor_sample_df<-read.csv(paste(result_dir,"/tumor_sample_df.csv",sep = ""))

tumor_sample_df<-tumor_sample_df %>% filter(!is.na(valid_af_v3))

tumor_sample_df<-tumor_sample_df %>% filter(valid_af_v3==1|valid_af_ad2==1|valid_kmeans_v3==1|
                                              valid_kmeans_ad2==1| valid_af_cosine==1| valid_af_jw==1|
                                              valid_af_lv==1|valid_euclidean_dist_v3==1| valid_euclidean_dist_ada2==1|
                                              valid_cosine_match==1|valid_jw_match==1|valid_jw_match==1)

tumor_sample_df$ground_truth <- NA
tumor_sample_df$ground_truth_val <- NA
tumor_sample_df$ID<-as.character(tumor_sample_df$ID)
for (iter in 1:dim(tumor_sample_df)[1]){
  ind<- which(tumor_sample_df[iter,seq(5,27,2)]==1)
  actual_ind<- seq(5,27,2)[ind]-1
  tumor_names <- unique(unlist(tumor_sample_df[iter,actual_ind]))
  if(length(ind)>0){
  if(length(tumor_names)>1){
    tumor_sample_df$ground_truth[iter]<-"MG"
    tumor_sample_df$ground_truth_val[iter]<-paste(tumor_names,collapse = ";")
    
  }else{
    tumor_sample_df$ground_truth[iter]<-"G"
    tumor_sample_df$ground_truth_val[iter]<-paste(tumor_names,collapse = ";")
  }
  }else{
    tumor_sample_df$ground_truth[iter]<-"U"
    tumor_sample_df$ground_truth_val[iter]<-"Unknown"
    }
}

tumor_sample_df <- tumor_sample_df %>% filter(ground_truth=="G")


#accuracy_df<- tumor_sample_df[,c(seq(5,27,2))]
#print(colSums(accuracy_df)/dim(accuracy_df)[1])
print(round(colSums(accuracy_df)/dim(accuracy_df)[1],2))

lv_df<-tumor_sample_df%>% dplyr::select(ID,nct_id,Tumor_Names,lv_match,valid_lv_match)