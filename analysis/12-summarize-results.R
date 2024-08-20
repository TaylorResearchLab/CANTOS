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

tumor_all_gt<-read.csv(paste(result_dir,"/tumor_sample_df_ground_truth_all_edition_aug19.csv",sep = ""))
tumor_5thed_gt<-read.csv(paste(result_dir_5th,"/tumor_sample_df_script10_5thed_corrected_ground_truth_aug19.csv",sep = ""))

tumor_all_gt<-tumor_all_gt[,c(-1)]
tumor_5thed_gt<-tumor_5thed_gt[,c(-1)]

tumor_all_gt<-tumor_all_gt%>%filter(ground_truth %in% c("G","MG","G-Manual"))
tumor_5thed_gt<-tumor_5thed_gt%>%filter(ground_truth %in% c("G","MG","G-Manual"))

all<-(colSums(tumor_all_gt[,c(seq(6,28,2))]))/1118
fifth<-(colSums(tumor_all_gt[,c(seq(6,28,2))]))/1033
all<-as.data.frame(all)
fifth<-as.data.frame(fifth)