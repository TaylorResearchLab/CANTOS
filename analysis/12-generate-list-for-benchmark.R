suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(tidyverse)
  library(stringi)
  library(qdapRegex)
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
source(paste(util_dir,"/cluster_label_assignment.R",sep=""))

# Read file v3
affinity_cluster_v3_reassigned_df<-read.csv(paste(intermediate_dir,"/affinity_cluster_v3_reassigned_df.csv",sep=""))
affinity_cluster_v3_reassigned_df<-affinity_cluster_v3_reassigned_df[,c(-1)]
# Read ADA2 
affinity_cluster_ADA2_reassigned_df<-read.csv(paste(intermediate_dir,"/affinity_cluster_ADA2_reassigned_df.csv",sep=""))
affinity_cluster_ADA2_reassigned_df<-affinity_cluster_ADA2_reassigned_df[,c(-1)]
#Read Kmeans 
kmeans_clust_result_embedding_ADA2 <- read_csv("analysis/results/kmeans_clust_result_embedding.csv")
kmeans_clust_result_embedding_ADA2<-kmeans_clust_result_embedding_ADA2[,c(-1)]
kmeans_clust_result_embedding_V3 <- read_csv("analysis/results/kmeans_clust_result_embedding_v3.csv")
kmeans_clust_result_embedding_V3<-kmeans_clust_result_embedding_V3[,c(-1)]

# Read edit distance cluster
nested_affinity_cluster_cosine <- read_csv("analysis/results/nested_affinity_cluster_cosine.csv")
nested_affinity_cluster_jw <- read_csv("analysis/results/nested_affinity_cluster_jw.csv")
nested_affinity_cluster_lv <- read_csv("analysis/results/nested_affinity_cluster_lv.csv")

nested_affinity_cluster_cosine<-nested_affinity_cluster_cosine[,c(-1)]
nested_affinity_cluster_jw<-nested_affinity_cluster_jw[,c(-1)]
nested_affinity_cluster_lv<-nested_affinity_cluster_lv[,c(-1)]

# Find the WHO and NCIT closest matches 
who_ncit_match_ADA2 <- affinity_cluster_ADA2_reassigned_df %>% dplyr::select(Tumor_Names,WHO_Matches,WHO_distance,NCIT_Matches,NCIT_distance) 
who_ncit_match_v3 <- affinity_cluster_v3_reassigned_df %>% dplyr::select(Tumor_Names,WHO_Matches,WHO_distance,NCIT_Matches,NCIT_distance) 

# Join the matches to Kmeans 
kmeans_clust_result_embedding_ADA2<-kmeans_clust_result_embedding_ADA2 %>% dplyr::select(Tumors,cluster)
kmeans_clust_result_embedding_V3<-kmeans_clust_result_embedding_V3 %>% dplyr::select(Tumors,cluster)

colnames(kmeans_clust_result_embedding_ADA2)<-c("Tumor_Names","Cluster_ID")
colnames(kmeans_clust_result_embedding_V3)<-c("Tumor_Names","Cluster_ID")

kmeans_clust_result_embedding_ADA2<-kmeans_clust_result_embedding_ADA2 %>% dplyr::left_join(who_ncit_match_ADA2,by="Tumor_Names")
kmeans_clust_result_embedding_V3<-kmeans_clust_result_embedding_V3 %>% dplyr::left_join(who_ncit_match_v3,by="Tumor_Names")

kmeans_clust_result_embedding_ADA2<- cluster_label_assignment(kmeans_clust_result_embedding_ADA2)
kmeans_clust_result_embedding_V3<- cluster_label_assignment(kmeans_clust_result_embedding_V3)


# Join the matches to affinity
nested_affinity_cluster_cosine<-nested_affinity_cluster_cosine %>% dplyr::select(Tumor_Names,Cluster_ID)
nested_affinity_cluster_jw<-nested_affinity_cluster_jw %>% dplyr::select(Tumor_Names,Cluster_ID)
nested_affinity_cluster_lv<-nested_affinity_cluster_lv %>% dplyr::select(Tumor_Names,Cluster_ID)

dissimilarity_matrix_cosine<-read.csv(paste(intermediate_dir,"/dissimilarity_matrix_cosine.csv",sep=""))
dissimilarity_matrix_jw<-read.csv(paste(intermediate_dir,"/dissimilarity_matrix_jw.csv",sep=""))
dissimilarity_matrix_lv<-read.csv(paste(intermediate_dir,"/dissimilarity_matrix_lv.csv",sep=""))

dissimilarity_matrix_cosine<-as.data.frame(dissimilarity_matrix_cosine)
dissimilarity_matrix_jw<-as.data.frame(dissimilarity_matrix_jw)
dissimilarity_matrix_lv<-as.data.frame(dissimilarity_matrix_lv)

colnames(dissimilarity_matrix_cosine)[1]<-"Tumor_Names"
colnames(dissimilarity_matrix_jw)[1]<-"Tumor_Names"
colnames(dissimilarity_matrix_lv)[1]<-"Tumor_Names"

colnames(dissimilarity_matrix_cosine)[2:16197]<-dissimilarity_matrix_cosine$Tumor_Names
colnames(dissimilarity_matrix_jw)[2:16197]<-dissimilarity_matrix_jw$Tumor_Names
colnames(dissimilarity_matrix_lv)[2:16197]<-dissimilarity_matrix_lv$Tumor_Names

NCIT_embedding_df <-read.csv(paste(data_dir,"/dt_input_file_6_dec/NCIT_Neoplasm_Core_terms_text-embedding-ada-002_embeddings.csv",sep=""))
WHO_embedding_df <-read.csv(paste(data_dir,"/dt_input_file_6_dec/WHO_Only_terms_text-embedding-ada-002_embeddings.csv",sep=""))
NCIT_embedding_df<-NCIT_embedding_df[c(-1),] # Remove the header (column name) embedding
WHO_embedding_df<-WHO_embedding_df[c(-1),] # Remove the header (column name) embedding

rownames(NCIT_embedding_df)<-NULL
rownames(WHO_embedding_df)<-NULL

NCIT_Tumors<-tolower(NCIT_embedding_df$Disease)
WHO_Tumors<-tolower(WHO_embedding_df$Disease)
rm(NCIT_embedding_df,WHO_embedding_df)

dissimilarity_matrix_cosine_who <- dissimilarity_matrix_cosine %>% dplyr::select(.dots = c("Tumor_Names",WHO_Tumors))
dissimilarity_matrix_jw_who <- dissimilarity_matrix_jw %>% dplyr::select(.dots = c("Tumor_Names",WHO_Tumors))
dissimilarity_matrix_lv_who <- dissimilarity_matrix_lv %>% dplyr::select(.dots = c("Tumor_Names",WHO_Tumors))


dissimilarity_matrix_cosine_ncit <- dissimilarity_matrix_cosine %>% dplyr::select(.dots = c("Tumor_Names",NCIT_Tumors))
dissimilarity_matrix_jw_ncit <- dissimilarity_matrix_jw %>% dplyr::select(.dots = c("Tumor_Names",NCIT_Tumors))
dissimilarity_matrix_lv_ncit <- dissimilarity_matrix_lv %>% dplyr::select(.dots = c("Tumor_Names",NCIT_Tumors))


