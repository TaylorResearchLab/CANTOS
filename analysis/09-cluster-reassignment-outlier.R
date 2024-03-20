suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(tidyverse)
  library(stringi)
  library(qdapRegex)
  library(ghql)
  library(readxl)
})

# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
intermediate_dir <- file.path(analysis_dir,"intermediate")

# Load data 
affinity_cluster_ADA2_df <- read.csv(paste(intermediate_dir,"/affinity_cluster_ADA2_df.csv",sep=""))
affinity_cluster_v3_df <- read.csv(paste(intermediate_dir,"/affinity_cluster_v3_df.csv",sep=""))

affinity_cluster_ADA2_df<-affinity_cluster_ADA2_df[,c(-1)]
affinity_cluster_v3_df<- affinity_cluster_v3_df[,c(-1)]


affinity_cluster_ADA2_reassigned_df<-affinity_cluster_ADA2_df %>% dplyr::select(Tumor_Names,Cluster_ID,Isolation_Outlier,LOF_Outlier,assigned_class,suggested_cluster_label)
affinity_cluster_v3_reassigned_df<-affinity_cluster_v3_df %>% dplyr::select(Tumor_Names,Cluster_ID,Isolation_Outlier,LOF_Outlier,assigned_class,suggested_cluster_label)

# Find all the tumors with outlier
outlier_tumor_ada2 <- affinity_cluster_ADA2_df%>%filter(Isolation_Outlier=="Yes" & LOF_Outlier=="Yes")%>%dplyr::select(Tumor_Names)
outlier_tumor_v3 <- affinity_cluster_v3_df%>%filter(Isolation_Outlier=="Yes" & LOF_Outlier=="Yes")%>%dplyr::select(Tumor_Names)

for (iter in 1: dim(outlier_tumor_ada2)[1]){
  current_tumor_name <- outlier_tumor_ada2$Tumor_Names[iter]
  current_tumor_details <-  affinity_cluster_ADA2_df %>% filter(Tumor_Names==current_tumor_name)
  
  current_tumor_assigned_class <- current_tumor_details$assigned_class
  
  cluster_matched_assigned_class<- affinity_cluster_ADA2_df$Cluster_ID[which(str_detect(string = affinity_cluster_ADA2_df$suggested_cluster_label,pattern=current_tumor_assigned_class))]
  cluster_matched_assigned_class<- unique(cluster_matched_assigned_class)
  
  if(length(cluster_matched_assigned_class)==1){ # No other cluster found
    affinity_cluster_ADA2_reassigned_df$suggested_cluster_label[which(affinity_cluster_ADA2_df$Tumor_Names==current_tumor_name)]<-current_tumor_assigned_class
    
  }
}

