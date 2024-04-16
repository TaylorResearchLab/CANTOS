### DELETE THIS SCRIPT NOT REQUIRED

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
affinity_cluster_v3_reassigned_df<-affinity_cluster_v3_reassigned_df[,c(-1,-10)]
colnames(affinity_cluster_v3_reassigned_df)[2]<-"Updated_Cluster_ID"
affinity_cluster_v3_reassigned_df<-affinity_cluster_v3_reassigned_df[order(affinity_cluster_v3_reassigned_df$Updated_Cluster_ID),]
# 
affinity_cluster_v3_reassigned_df$ID <- seq.int(nrow(affinity_cluster_v3_reassigned_df))

# match location
affinity_cluster_v3_reassigned_df<-affinity_cluster_v3_reassigned_df %>% dplyr::mutate(exact_match=case_when(NCIT_distance==0~"Yes",
                                                                                                             WHO_distance==0~"Yes",
                                                                                                             TRUE~"No"))
affinity_cluster_v3_reassigned_df<-affinity_cluster_v3_reassigned_df[,c(9,1:8)]

# annotate file
affinity_cluster_v3_manual_annotation <- affinity_cluster_v3_reassigned_df %>% filter(exact_match=="No")

affinity_cluster_v3_manual_annotation_short <- affinity_cluster_v3_manual_annotation %>% dplyr::select(ID, Tumor_Names,Updated_Cluster_ID,who_cluster_label,ncit_cluster_label)


# write to gdrive
write.csv(affinity_cluster_v3_manual_annotation_short,paste(intermediate_dir,"/affinity_cluster_v3_manual_annotation_short.csv",sep=""))