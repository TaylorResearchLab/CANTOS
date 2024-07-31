# Annotate Affinity cluster results of ADA2
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

load(paste(intermediate_dir,"/affinity_cluster_df_ada2_5thed.RData",sep=""))
load(paste(intermediate_dir,"/combined_embedding_ada2_df_5thed.RData",sep=""))

NCIT_embedding_df <-read.csv(paste(data_dir,"/dt_input_file_6_dec/NCIT_Neoplasm_Core_terms_text-embedding-ada-002_embeddings.csv",sep=""))
NCIT_embedding_df<-NCIT_embedding_df[c(-1),] # Remove the header (column name) embedding


WHO_Terms_All <-readxl::read_xlsx(paste(data_dir,"/WHO_Tumors/result/WHO_Tumor_all_edition.xlsx",sep=""))
WHO_Terms_5th<-WHO_Terms_All%>%filter(edition_5th=="Yes")

WHO_embedding_df <- read.csv(paste(data_dir,"/WHO_Aggregate_ADA2.csv",sep="")) #
WHO_embedding_df<-WHO_embedding_df %>% dplyr::filter(Disease %in% WHO_Terms_5th$Tumor_Names)
WHO_embedding_df<-WHO_embedding_df[,c(-1)]
rownames(NCIT_embedding_df)<-NULL
rownames(WHO_embedding_df)<-NULL
colnames(WHO_embedding_df)<-colnames(NCIT_embedding_df)



cl <- makeCluster(6, outfile="")
registerDoParallel(cl)

CalculateEuclideanDistance <- function(vect1, vect2) sqrt(sum((vect1 - vect2)^2)) 

outer_who_final<-foreach(i = 1:dim(combined_embedding_df)[1], .combine = rbind) %dopar% { #12:10pm -
  print(i)
  #s <- apply(WHO_embedding_df[,2:1537],1,CalculateEuclideanDistance,vect2=combined_embedding_df[i,])
  embedding_pairwise<- as.matrix(rbind(combined_embedding_df[i,],WHO_embedding_df[,2:ncol(WHO_embedding_df)]))
  euclidean_dist <- as.matrix(dist(embedding_pairwise,method = "euclidean"))
  d<-as.double(euclidean_dist[1,c(-1)])
}
colnames(outer_who_final)<-(WHO_embedding_df$Disease)
rownames(outer_who_final)<-rownames(combined_embedding_df)

