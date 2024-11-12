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
intermediate_dir_5th <- file.path(analysis_dir,"intermediate_5th")
result_dir <- file.path(analysis_dir,"results")
result_dir_5th <- file.path(analysis_dir,"results_5th")


# Load the annotations for 5th edition 
tumor_5th_edition<-read.csv(paste(result_dir_5th,"/tumor_sample_df_gt_annotated_5th.csv",sep = ""))
tumor_all_edition<-read.csv(paste(result_dir,"/tumor_sample_df_gt_annotated_all.csv",sep = ""))

# Pick only the Euclidean distance standardization with V3 embeddings
tumor_5th_edition<-tumor_5th_edition %>% dplyr::select(nct_id,Tumor_Names,ground_truth_val,ground_truth,
                                                       euclidean_dist_v3,valid_euclidean_dist_v3)
tumor_all_edition<-tumor_all_edition %>% dplyr::select(nct_id,Tumor_Names,ground_truth_val,ground_truth,
                                                       euclidean_dist_v3,valid_euclidean_dist_v3)


# Load 5th edition distance
affinity_cluster_v3_reassigned_5thed_df<-read.csv(paste(intermediate_dir_5th,"/affinity_cluster_v3_reassigned_df_5thed.csv",sep=""))
WHO_5th_edition<-affinity_cluster_v3_reassigned_5thed_df %>% dplyr::select(Tumor_Names,WHO_Matches,WHO_distance)


# Load all edition distance
affinity_cluster_v3_reassigned_df<-read.csv(paste(intermediate_dir,"/affinity_cluster_v3_reassigned_df.csv",sep=""))
WHO_all_edition<-affinity_cluster_v3_reassigned_df %>% dplyr::select(Tumor_Names,WHO_Matches,WHO_distance)


# Join 5th edition data
tumor_5th_edition<- tumor_5th_edition %>% dplyr::left_join(WHO_5th_edition,by="Tumor_Names")

#Join all edition data
tumor_all_edition<- tumor_all_edition %>% dplyr::left_join(WHO_all_edition,by="Tumor_Names")



# Remove cases where there were no ground truths (NF = Not Found)
tumor_5th_edition<- tumor_5th_edition %>%filter(ground_truth !="NF")
tumor_all_edition<- tumor_all_edition %>%filter(ground_truth !="NF")


# data table for 5th edition Euclidean V3 distance
