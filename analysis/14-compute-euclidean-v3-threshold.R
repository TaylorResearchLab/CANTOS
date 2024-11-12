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