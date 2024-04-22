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


tumor_sample_df<-read.csv(paste(result_dir,"/tumor_sample_df.csv",sep = ""))

tumor_sample_df<-tumor_sample_df %>% filter(!is.na(valid_af_v3))

accuracy_df<- tumor_sample_df[,c(seq(5,27,2))]
print(colSums(accuracy_df)/dim(accuracy_df)[1])