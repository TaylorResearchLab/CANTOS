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


# 
load(paste(intermediate_dir,"/affinity_cluster_df.RData",sep=""))
load(paste(intermediate_dir,"/combined_embedding_df.RData",sep=""))
NCIT_embedding_df <-read.csv(paste(data_dir,"/dt_input_file_6_dec/NCIT_Neoplasm_Core_terms_text-embedding-ada-002_embeddings.csv",sep=""))
WHO_embedding_df <-read.csv(paste(data_dir,"/dt_input_file_6_dec/WHO_Only_terms_text-embedding-ada-002_embeddings.csv",sep=""))
NCIT_embedding_df<-NCIT_embedding_df[c(-1),] # Remove the header (column name) embedding
WHO_embedding_df<-WHO_embedding_df[c(-1),] # Remove the header (column name) embedding

rownames(NCIT_embedding_df)<-NULL
rownames(WHO_embedding_df)<-NULL

#combined_embedding_df$Tumor_Names <- rownames(combined_embedding_df)
#CT_embedding<- combined_embedding_df %>% filter(!(Tumor_Names %in% unique(c(tolower(NCIT_embedding_df$Disease),tolower(WHO_embedding_df$Disease)) ) ))
#
tumor_distances_df <- as.data.frame(matrix(nrow=dim(combined_embedding_df)[1],ncol = 5))
colnames(tumor_distances_df)<- c("Tumor_Names","NCIT_Match","NCIT_Distance","WHO_Match","WHO_Distance")
tumor_distances_df$Tumor_Names<- rownames(combined_embedding_df)


#
cl <- makeCluster(5, outfile="")
registerDoParallel(cl)

CalculateEuclideanDistance <- function(vect1, vect2) sqrt(sum((vect1 - vect2)^2)) 


outer_who_final<-foreach(i = 1:100, .combine = rbind) %dopar% {
  print(i)
  #s <- apply(WHO_embedding_df[,2:1537],1,CalculateEuclideanDistance,vect2=combined_embedding_df[i,])
  embedding_pairwise<- as.matrix(rbind(combined_embedding_df[i,],WHO_embedding_df[,2:1537]))
  euclidean_dist <- as.matrix(dist(embedding_pairwise,method = "euclidean"))
  d<-as.double(euclidean_dist[1,c(-1)])
}
colnames(outer_who_final)<-(WHO_embedding_df$Disease)
rownames(outer_who_final)<-rownames(combined_embedding_df)