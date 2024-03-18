suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(tidyverse)
  library(stringi)
  library(qdapRegex)
  library(jsonlite)
  library(httr)
  library(biomaRt)
  library(ghql)
  library(readxl)
  library(doParallel)
  library(foreach)
})

# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
intermediate_dir <- file.path(analysis_dir,"intermediate")


# Load affinity clusters for ADA2 and V3
affinity_cluster_ADA2_df <- read.csv(paste(intermediate_dir,"/affinity_cluster_ADA2_df.csv",sep=""))
affinity_cluster_V3_df <- read.csv(paste(intermediate_dir,"/affinity_cluster_v3_df.csv",sep=""))



# Load embeddings 
disease_transform_ADA2<- read.csv(paste(intermediate_dir,"/disease_transform_pca.csv",sep="") )
colnames(disease_transform_ADA2)[1]<-"Tumor_Names"
rownames(disease_transform_ADA2)<-disease_transform_ADA2$Tumor_Names 


disease_transform_V3<- read.csv(paste(intermediate_dir,"/disease_transform_pca_v3.csv",sep="") )
colnames(disease_transform_V3)[1]<-"Tumor_Names"
rownames(disease_transform_V3)<-disease_transform_ADA2$Tumor_Names 



