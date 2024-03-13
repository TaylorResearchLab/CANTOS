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
source(paste(util_dir,"/cluster_label_assignment.R",sep=""))


load(paste(intermediate_dir,"/affinity_cluster_v3_df.RData",sep=""))
embedding_v3_large <- read.csv(paste(data_dir,"/embedding_tumor_names_text-embedding-3-large_embeddings.csv",sep=""))
rownames(embedding_v3_large)<-embedding_v3_large$Disease

NCIT_embedding_df <-read.csv(paste(data_dir,"/dt_input_file_6_dec/NCIT_Neoplasm_Core_terms_text-embedding-ada-002_embeddings.csv",sep=""))
WHO_embedding_df <-read.csv(paste(data_dir,"/dt_input_file_6_dec/WHO_Only_terms_text-embedding-ada-002_embeddings.csv",sep=""))

NCIT_embedding_df<-NCIT_embedding_df[c(-1),]
WHO_embedding_df<-WHO_embedding_df[c(-1),]

NCIT_Tumors <- tolower(NCIT_embedding_df$Disease)
WHO_Tumors <- tolower(WHO_embedding_df$Disease)

rm(NCIT_embedding_df,WHO_embedding_df)


NCIT_embedding_df<- embedding_v3_large %>% dplyr::filter(Disease %in% NCIT_Tumors)
WHO_embedding_df<- embedding_v3_large %>% dplyr::filter(Disease %in% WHO_Tumors)
embedding_v3_large<-embedding_v3_large[,c(-1)]


tumor_distances_df <- as.data.frame(matrix(nrow=dim(embedding_v3_large)[1],ncol = 5))
colnames(tumor_distances_df)<- c("Tumor_Names","NCIT_Match","NCIT_Distance","WHO_Match","WHO_Distance")
tumor_distances_df$Tumor_Names<- rownames(embedding_v3_large)


cl <- makeCluster(6, outfile="")
registerDoParallel(cl)

CalculateEuclideanDistance <- function(vect1, vect2) sqrt(sum((vect1 - vect2)^2)) 

outer_who_final<-foreach(i = 1:dim(embedding_v3_large)[1], .combine = rbind) %dopar% { #3:34 pm -
  print(i)
  embedding_pairwise<- as.matrix(rbind(embedding_v3_large[i,],WHO_embedding_df[,2:3073]))
  euclidean_dist <- as.matrix(dist(embedding_pairwise,method = "euclidean"))
  d<-as.double(euclidean_dist[1,c(-1)])
}
colnames(outer_who_final)<-(WHO_embedding_df$Disease)
rownames(outer_who_final)<-rownames(embedding_v3_large)

stopCluster(cl)

save.image("script-7B.RData")


