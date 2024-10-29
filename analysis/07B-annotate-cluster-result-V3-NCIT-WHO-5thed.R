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
load(paste(intermediate_dir,"/combined_embedding_v3_df_5thed.RData",sep=""))


WHO_embedding_df<-read.csv(paste(data_dir,"/WHO_Terms_All_V3.csv",sep=""))
NCIT_embedding_df<-read.csv(paste(data_dir,"/NCIT_Embeddings_V3.csv",sep=""))

WHO_embedding_df<-WHO_embedding_df[,c(-1)]
NCIT_embedding_df<-NCIT_embedding_df[,c(-1)]

WHO_Terms_All <-readxl::read_xlsx(paste(data_dir,"/WHO_Tumors/result/WHO_Tumor_all_edition.xlsx",sep=""))
WHO_Terms_5th<-WHO_Terms_All%>%filter(edition_5th=="Yes")

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

outer_who_final<-as.data.frame(outer_who_final)
outer_who_final<-outer_who_final %>% dplyr::select(any_of(WHO_Terms_5th$Tumor_Names))
outer_who_final<-as.matrix(outer_who_final)




combined_embedding<-combined_embeddings_df
# 


outer_NCIT_final<-foreach(i = 1:dim(combined_embeddings_df)[1], .combine = rbind) %dopar% { #2:03 pm Friday
  print(i)
  embedding_pairwise<- as.matrix(rbind(combined_embeddings_df[i,],NCIT_embedding_df[,2:3073]))
  euclidean_dist <- as.matrix(dist(embedding_pairwise,method = "euclidean"))
  d<-as.double(euclidean_dist[1,c(-1)])
}
colnames(outer_NCIT_final)<-(NCIT_embedding_df$Tumor_Names)
rownames(outer_NCIT_final)<-rownames(combined_embeddings_df)





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



affinity_cluster_v3_df<- cluster_label_assignment_refined(affinity_cluster_v3_df)
affinity_cluster_v3_df<-affinity_cluster_v3_df[,c(-1)]
tumor_id<- read.csv(paste(data_dir,"/Tumor_NCT_ID.csv",sep=""))
tumor_id<- tumor_id[,c(-1)]
affinity_cluster_v3_df<-affinity_cluster_v3_df%>%left_join(tumor_id,by="Tumor_Names")
affinity_cluster_v3_df<-affinity_cluster_v3_df[,c(9,1:8)]


##
write.csv(affinity_cluster_v3_df,paste(intermediate_dir,"/affinity_cluster_v3_df_5thed.csv",sep=""))
#save.image(file = "script7b-5thed_aug4.RData")
save.image(file = "script7b-5thed_aug18.RData")
