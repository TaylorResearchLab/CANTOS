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


tumor_all_gt<-read.csv(paste(result_dir,"/tumor_sample_df_gt_annotated_all_sep11.csv",sep = ""))
tumor_5thed_gt<-read.csv(paste(result_dir_5th,"/tumor_sample_df_gt_annotated_5th_sep11.csv",sep = ""))



tumor_all_gt<-tumor_all_gt[,c(-1)]
tumor_5thed_gt<-tumor_5thed_gt[,c(-1)]

tumor_all_gt<-tumor_all_gt%>%filter(ground_truth %in% c("G","MG","G-Manual"))
tumor_5thed_gt<-tumor_5thed_gt%>%filter(ground_truth %in% c("G","MG","G-Manual"))

all<-(colSums(tumor_all_gt[,c(seq(6,28,2))]))/1118
fifth<-(colSums(tumor_5thed_gt[,c(seq(6,28,2))]))/1033
all<-as.data.frame(all)
fifth<-as.data.frame(fifth)

all$method<-rownames(all)
fifth$method<-rownames(fifth)
rownames(all)<-NULL
rownames(fifth)<-NULL
all <- all %>%dplyr::select(method,all)
fifth<- fifth %>%dplyr::select(method,fifth)

all<-all[order(all$all,decreasing = TRUE),]
fifth<-fifth[order(fifth$fifth,decreasing = TRUE),]

for(iter in 1:12){
  if(all$method[iter]=="valid_euclidean_dist_v3" ){
    all$method[iter]="LTE-3 + Euclidean Dist"
  }else if(all$method[iter]=="valid_af_v3"){
    all$method[iter]="LTE-3 + AP"
    
  }else if(all$method[iter]=="valid_euclidean_dist_ada2"){
    all$method[iter]="ADA002 + Euclidean Dist"
    
  }else if(all$method[iter]=="valid_af_ad2"){
    all$method[iter]="ADA002 + AP"
    
  }else if(all$method[iter]=="valid_kmeans_v3"){
    all$method[iter]="LTE-3 + KMeans"
    
  }else if(all$method[iter]=="valid_kmeans_ada2"){
    all$method[iter]="ADA002 + KMeans"
    
  }else if(all$method[iter]=="valid_lv_match"){
    all$method[iter]="Levenshtein"
    
  }else if(all$method[iter]=="valid_af_lv"){
    all$method[iter]="Levenshtein + AP"
    
  }else if(all$method[iter]=="valid_jw_match"){
    all$method[iter]="Jarro Winkler"
  }else if(all$method[iter]=="valid_af_jw"){
    all$method[iter]="Jarro Winkler + AP"
    
  }else if(all$method[iter]=="valid_cosine_match"){
    all$method[iter]="Cosine"
    
  }else if(all$method[iter]=="valid_af_cosine"){
    all$method[iter]="Cosine + AP"
    
  }
  
}

for(iter in 1:12){
  if(fifth$method[iter]=="valid_euclidean_dist_v3" ){
    fifth$method[iter]="LTE-3 + Euclidean Dist"
  }else if(fifth$method[iter]=="valid_af_v3"){
    fifth$method[iter]="LTE-3 + AP"
    
  }else if(fifth$method[iter]=="valid_euclidean_dist_ada2"){
    fifth$method[iter]="ADA002 + Euclidean Dist"
    
  }else if(fifth$method[iter]=="valid_af_ad2"){
    fifth$method[iter]="ADA002 + AP"
    
  }else if(fifth$method[iter]=="valid_kmeans_v3"){
    fifth$method[iter]="LTE-3 + KMeans"
    
  }else if(fifth$method[iter]=="valid_kmeans_ada2"){
    fifth$method[iter]="ADA002 + KMeans"
    
  }else if(fifth$method[iter]=="valid_lv_match"){
    fifth$method[iter]="Levenshtein"
    
  }else if(fifth$method[iter]=="valid_af_lv"){
    fifth$method[iter]="Levenshtein + AP"
    
  }else if(fifth$method[iter]=="valid_jw_match"){
    fifth$method[iter]="Jarro Winkler"
  }else if(fifth$method[iter]=="valid_af_jw"){
    fifth$method[iter]="Jarro Winkler + AP"
    
  }else if(fifth$method[iter]=="valid_cosine_match"){
    fifth$method[iter]="Cosine"
    
  }else if(fifth$method[iter]=="valid_af_cosine"){
    fifth$method[iter]="Cosine + AP"
    
  }
  
}
rownames(all)<- NULL
rownames(fifth)<-NULL
print(all)
print(fifth)
