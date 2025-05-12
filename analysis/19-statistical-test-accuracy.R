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



tumor_all_gt<-read.csv(paste(result_dir,"/tumor_manually_validated_all.csv",sep = ""))
tumor_5thed_gt<-read.csv(paste(result_dir_5th,"/tumor_manually_validated_5th.csv",sep = ""))

tumor_all_gt<-tumor_all_gt[,c(-1)]
tumor_5thed_gt<-tumor_5thed_gt[,c(-1)]

tumor_all_gt<-tumor_all_gt%>%filter(ground_truth %in% c("G","MG","G-Manual"))
tumor_5thed_gt<-tumor_5thed_gt%>%filter(ground_truth %in% c("G","MG","G-Manual"))



correct_5thed_matrix<-tumor_5thed_gt[,seq(6,76,2)]
correct_all_matrix<-tumor_all_gt[,seq(6,76,2)]


cleaned_colnames_5thed <- gsub("^valid_", "", colnames(correct_5thed_matrix))
cleaned_colnames_all <- gsub("^valid_", "", colnames(correct_all_matrix))

colnames(correct_5thed_matrix)<-cleaned_colnames_5thed
colnames(correct_all_matrix)<-cleaned_colnames_all



