# Annotate Affinity cluster results of V3

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
intermediate_dir <- file.path(analysis_dir,"intermediate_5th")
source(paste(util_dir,"/cluster_label_assignment_refined.R",sep=""))


load(paste(intermediate_dir,"/affinity_cluster_v3_df_5thed.RData",sep=""))
load(paste(intermediate_dir,"/combined_embedding_v3_df.RData",sep=""))


WHO_embedding_df<-read.csv(paste(data_dir,"/WHO_Terms_All_V3.csv",sep=""))
NCIT_embedding_df<-read.csv(paste(data_dir,"/NCIT_Embeddings_V3.csv",sep=""))

WHO_embedding_df<-WHO_embedding_df[,c(-1)]
NCIT_embedding_df<-NCIT_embedding_df[,c(-1)]

cl <- makeCluster(50, outfile="")
registerDoParallel(cl)

CalculateEuclideanDistance <- function(vect1, vect2) sqrt(sum((vect1 - vect2)^2)) 

outer_who_final<-foreach(i = 1:dim(combined_embeddings_df)[1], .combine = rbind) %dopar% { #Two days
  print(i)
  embedding_pairwise<- as.matrix(rbind(combined_embeddings_df[i,],WHO_embedding_df[,2:3073]))
  euclidean_dist <- as.matrix(dist(embedding_pairwise,method = "euclidean"))
  d<-as.double(euclidean_dist[1,c(-1)])
}
colnames(outer_who_final)<-(WHO_embedding_df$Tumor_Names)
rownames(outer_who_final)<-rownames(combined_embeddings_df)

combined_embeddings_df
# 
start

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
