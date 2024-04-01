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

# Read file v3
affinity_cluster_v3_reassigned_df<-read.csv(paste(intermediate_dir,"/affinity_cluster_v3_reassigned_df.csv",sep=""))
affinity_cluster_v3_reassigned_df<-affinity_cluster_v3_reassigned_df[,c(-1,-8)]
colnames(affinity_cluster_v3_reassigned_df)[2]<-"Updated_Cluster_ID"

affinity_cluster_v3_reassigned_df<-affinity_cluster_v3_reassigned_df[order(affinity_cluster_v3_reassigned_df$Updated_Cluster_ID),]

affinity_cluster_v3_reassigned_df$ID <- seq.int(nrow(affinity_cluster_v3_reassigned_df))
affinity_cluster_v3_reassigned_df<-affinity_cluster_v3_reassigned_df[,c(7,1:6)]
