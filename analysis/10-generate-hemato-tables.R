# Script to generate hematology tumors table
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



# 
#load(paste(intermediate_dir,"/affinity_cluster_df.RData",sep=""))
affinity_cluster_ADA2_df <- read.csv(paste(intermediate_dir,"/affinity_cluster_ADA2_reassigned_df.csv",sep=""))
affinity_cluster_V3_df <- read.csv(paste(intermediate_dir,"/affinity_cluster_v3_reassigned_df.csv",sep=""))

affinity_cluster_ADA2_df<-affinity_cluster_ADA2_df%>% dplyr::select(Tumor_Names,updated_ID,who_cluster_label,ncit_cluster_label)
affinity_cluster_V3_df<-affinity_cluster_V3_df%>% dplyr::select(Tumor_Names,updated_ID,who_cluster_label,ncit_cluster_label)

colnames(affinity_cluster_ADA2_df)[2]<-"Updated_Cluster_ID"
colnames(affinity_cluster_V3_df)[2]<-"Updated_Cluster_ID"

###### Lymphoma Luek analysis

lymphoma_leukemia_strings <- c("leukemia", "lymphoma", "leukemias", "lymphomas", "leukaemia", "leukaemias",
                               "leuk","hematologic tumors","hemato")

affinity_cluster_hema_ADA2_df <- affinity_cluster_ADA2_df %>% dplyr::filter(str_detect(Tumor_Names,paste(lymphoma_leukemia_strings, collapse = "|")))
hema_label_ADA2 <- unique(affinity_cluster_hema_ADA2_df$Updated_Cluster_ID)
affinity_cluster_hema_ADA2_df <- affinity_cluster_ADA2_df %>% dplyr::filter(Updated_Cluster_ID %in% hema_label_ADA2 )


affinity_cluster_hema_V3_df <- affinity_cluster_V3_df %>% dplyr::filter(str_detect(Tumor_Names,paste(lymphoma_leukemia_strings, collapse = "|")))
hema_label_V3 <- unique(affinity_cluster_hema_V3_df$Updated_Cluster_ID)
affinity_cluster_hema_V3_df <- affinity_cluster_V3_df %>% dplyr::filter(Updated_Cluster_ID %in% hema_label_V3 )



write.csv(affinity_cluster_hema_ADA2_df,paste(intermediate_dir,"/hemato_tumor_ADA2.csv",sep=""))
write.csv(affinity_cluster_hema_V3_df,paste(intermediate_dir,"/hemato_tumor_V3.csv",sep=""))

#save.image("script08.RData")