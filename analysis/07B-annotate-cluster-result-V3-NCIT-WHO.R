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
intermediate_dir <- file.path(analysis_dir,"intermediate")
source(paste(util_dir,"/cluster_label_assignment_refined.R",sep=""))


load(paste(intermediate_dir,"/affinity_cluster_v3_df.RData",sep=""))
load(paste(intermediate_dir,"/combined_embedding_v3_df.RData",sep=""))


WHO_embedding_df<-read.csv(paste(data_dir,"/WHO_Terms_ALL_V3.csv",sep=""))
NCIT_embedding_df<-read.csv(paste(data_dir,"/NCIT_Embeddings_V3.csv",sep=""))

WHO_embedding_df<-WHO_embedding_df[,c(-1)]
NCIT_embedding_df<-NCIT_embedding_df[,c(-1)]


cl <- makeCluster(6, outfile="")
registerDoParallel(cl)

CalculateEuclideanDistance <- function(vect1, vect2) sqrt(sum((vect1 - vect2)^2)) 

combined_embedding<-combined_embeddings_df

outer_who_final<-foreach(i = 1:dim(combined_embeddings_df)[1], .combine = rbind) %dopar% { #Two days
  print(i)
  embedding_pairwise<- as.matrix(rbind(combined_embeddings_df[i,],WHO_embedding_df[,2:3073]))
  euclidean_dist <- as.matrix(dist(embedding_pairwise,method = "euclidean"))
  d<-as.double(euclidean_dist[1,c(-1)])
}
colnames(outer_who_final)<-(WHO_embedding_df$Tumor_Names)
rownames(outer_who_final)<-rownames(combined_embeddings_df)

save.image("script-7B.RData")


outer_NCIT_final<-foreach(i = 1:dim(combined_embeddings_df)[1], .combine = rbind) %dopar% { #2:03 pm Friday
  print(i)
  embedding_pairwise<- as.matrix(rbind(combined_embeddings_df[i,],NCIT_embedding_df[,2:3073]))
  euclidean_dist <- as.matrix(dist(embedding_pairwise,method = "euclidean"))
  d<-as.double(euclidean_dist[1,c(-1)])
}
colnames(outer_NCIT_final)<-(NCIT_embedding_df$Tumor_Names)
rownames(outer_NCIT_final)<-rownames(combined_embeddings_df)

stopCluster(cl)
save.image("script-7B.RData")




######################
index_min_who <- as.matrix(apply(outer_who_final, 1, which.min))

who_match_df <- cbind(rownames(outer_who_final))

colnames(who_match_df)<-"Tumor_Names"
who_match_df <-as.data.frame(who_match_df)

who_match_df$WHO_Matches<- NA
who_match_df$WHO_distance<-NA

for (iter in 1: dim(who_match_df)[1]){
  
  who_match_df$WHO_Matches[iter] <- colnames(outer_who_final)[index_min_who[iter]]
  who_match_df$WHO_distance[iter]<-outer_who_final[iter,index_min_who[iter]]
  
}



index_min_NCIT <- as.matrix(apply(outer_NCIT_final, 1, which.min))

NCIT_match_df <- cbind(rownames(outer_NCIT_final))

colnames(NCIT_match_df)<-"Tumor_Names"
NCIT_match_df <-as.data.frame(NCIT_match_df)

NCIT_match_df$NCIT_Matches<- NA
NCIT_match_df$NCIT_distance<-NA

for (iter in 1: dim(NCIT_match_df)[1]){
  
  NCIT_match_df$NCIT_Matches[iter] <- colnames(outer_NCIT_final)[index_min_NCIT[iter]]
  NCIT_match_df$NCIT_distance[iter]<-outer_NCIT_final[iter,index_min_NCIT[iter]]
  
}


########
affinity_cluster_v3_df<- affinity_cluster_v3_df %>% dplyr::left_join(who_match_df,by="Tumor_Names")
affinity_cluster_v3_df<- affinity_cluster_v3_df %>% dplyr::left_join(NCIT_match_df,by="Tumor_Names")

# affinity_cluster_v3_df <- affinity_cluster_v3_df %>% dplyr::mutate(assigned_class = case_when(NCIT_distance < WHO_distance ~ NCIT_Matches,
#                                                                                         NCIT_distance > WHO_distance ~ WHO_Matches,
#                                                                                         NCIT_distance==WHO_distance ~ WHO_Matches,
#                                                                                         TRUE~NA))

affinity_cluster_v3_df<- cluster_label_assignment_refined(affinity_cluster_v3_df)



save.image("script-7B.RData")
write.csv(affinity_cluster_v3_df,paste(intermediate_dir,"/affinity_cluster_v3_df.csv",sep=""))

