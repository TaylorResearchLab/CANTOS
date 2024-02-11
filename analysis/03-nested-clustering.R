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
})

# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
intermediate_dir <- file.path(analysis_dir,"intermediate")

# Load affinity Cluster
load(paste(intermediate_dir,"/affinity_cluster_annotation.RData",sep=""))
source(paste(util_dir,"/nested_affinity_cluster.R",sep=""))
source(paste(util_dir,"/cluster_label_assignment.R",sep=""))

# load disease_transfor
disease_transform <- read.csv(paste(intermediate_dir,"/disease_transform_pca.csv",sep=""))
colnames(disease_transform)[1]<-"Tumor_Name"

# Load data
ncit_match_df <- read.csv(paste(intermediate_dir,"/ncit_match_df.csv",sep=""))
who_match_df <- read.csv(paste(intermediate_dir,"/who_ct_distance_mat.csv",sep=""))


# Need for second run of clustering

affinity_cluster_nested <- affinity_cluster_annotation %>% dplyr::select(Tumor_Names,Pediatric_SubsetCluster_ID, SubsetCluster_IDs)
affinity_cluster_nested <- nested_affinity_cluster(n=3,affinity_cluster_nested,disease_transform)


affinity_cluster_nested<- affinity_cluster_nested %>% dplyr::left_join(ncit_match_df,by="Tumor_Names")
affinity_cluster_nested <- affinity_cluster_nested %>%dplyr::left_join(who_match_df,by="Tumor_Names")



affinity_cluster_nested <- affinity_cluster_nested %>% dplyr::mutate(assigned_class = case_when(ncit_distance < WHO_distance ~ NCIT_Matches,
                                                                                                ncit_distance > WHO_distance ~ WHO_Matches,
                                                                                                TRUE ~ "Both"))

# Cluster voting
affinity_cluster_nested<- cluster_label_assignment(affinity_cluster_nested)

disease_affinity_cluster_table<- affinity_cluster_nested %>% dplyr::select(Tumor_Names,cluster_label)

# Write
write.csv(affinity_cluster_nested,paste(intermediate_dir,"/affinity_cluster_nested.csv",sep=""))
