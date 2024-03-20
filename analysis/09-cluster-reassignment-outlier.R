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


# Find all the indices with outlier
affinity_cluster_ADA2_df<-affinity_cluster_ADA2_df %>% dplyr::select(Tumor_Names,Cluster_ID,Isolation_Outlier,LOF_Outlier,assigned_class,suggested_cluster_label)
affinity_cluster_v3_df<-affinity_cluster_v3_df %>% dplyr::select(Tumor_Names,Cluster_ID,Isolation_Outlier,LOF_Outlier,assigned_class,suggested_cluster_label)