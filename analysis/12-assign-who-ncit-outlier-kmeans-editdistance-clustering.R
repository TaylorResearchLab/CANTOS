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
source(paste(util_dir,"/cluster_label_assignment.R",sep=""))
source(paste(util_dir,"/outlier_detection_edit_dist.R",sep=""))
source(paste(util_dir,"/edit_distance_cluster_reassignment.R",sep=""))


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

rownames(dissimilarity_matrix_cosine)<-dissimilarity_matrix_cosine$Tumor_Names
rownames(dissimilarity_matrix_jw)<-dissimilarity_matrix_jw$Tumor_Names
rownames(dissimilarity_matrix_lv)<-dissimilarity_matrix_lv$Tumor_Names

dissimilarity_matrix_cosine<-dissimilarity_matrix_cosine[,c(-1)]
dissimilarity_matrix_jw<-dissimilarity_matrix_jw[,c(-1)]
dissimilarity_matrix_lv<-dissimilarity_matrix_lv[,c(-1)]


dissimilarity_matrix_cosine_who <- dissimilarity_matrix_cosine %>% dplyr::select(one_of(WHO_Tumors))
dissimilarity_matrix_jw_who <- dissimilarity_matrix_jw %>% dplyr::select(one_of(WHO_Tumors))
dissimilarity_matrix_lv_who <- dissimilarity_matrix_lv %>% dplyr::select(one_of(WHO_Tumors))




dissimilarity_matrix_cosine_ncit <- dissimilarity_matrix_cosine %>% dplyr:::select(one_of(NCIT_Tumors))
dissimilarity_matrix_jw_ncit <- dissimilarity_matrix_jw %>% dplyr:::select(one_of(NCIT_Tumors))
dissimilarity_matrix_lv_ncit <- dissimilarity_matrix_lv %>% dplyr:::select(one_of(NCIT_Tumors))

######## Now find the. WHO matches
index_min_who_cosine <- as.matrix(apply(dissimilarity_matrix_cosine_who, 1, which.min))
index_min_who_jw <- as.matrix(apply(dissimilarity_matrix_jw_who, 1, which.min))
index_min_who_lv <- as.matrix(apply(dissimilarity_matrix_lv_who, 1, which.min))



who_match_cosine_df <- cbind(rownames(dissimilarity_matrix_cosine_who))
who_match_jw_df <- cbind(rownames(dissimilarity_matrix_jw_who))
who_match_lv_df <- cbind(rownames(dissimilarity_matrix_lv_who))

colnames(who_match_cosine_df)<-"Tumor_Names"
colnames(who_match_jw_df)<-"Tumor_Names"
colnames(who_match_lv_df)<-"Tumor_Names"

who_match_cosine_df <-as.data.frame(who_match_cosine_df)
who_match_jw_df <-as.data.frame(who_match_jw_df)
who_match_lv_df <-as.data.frame(who_match_lv_df)

who_match_cosine_df$WHO_Matches<- NA
who_match_cosine_df$WHO_distance<-NA

who_match_jw_df$WHO_Matches<- NA
who_match_jw_df$WHO_distance<-NA

who_match_lv_df$WHO_Matches<- NA
who_match_lv_df$WHO_distance<-NA

for (iter in 1: dim(who_match_cosine_df)[1]){
  
  who_match_cosine_df$WHO_Matches[iter] <- colnames(dissimilarity_matrix_cosine_who)[index_min_who_cosine[iter]]
  who_match_cosine_df$WHO_distance[iter]<-dissimilarity_matrix_cosine_who[iter,index_min_who_cosine[iter]]
  
  who_match_jw_df$WHO_Matches[iter] <- colnames(dissimilarity_matrix_jw_who)[index_min_who_jw[iter]]
  who_match_jw_df$WHO_distance[iter]<-dissimilarity_matrix_jw_who[iter,index_min_who_jw[iter]]
  
  who_match_lv_df$WHO_Matches[iter] <- colnames(dissimilarity_matrix_lv_who)[index_min_who_lv[iter]]
  who_match_lv_df$WHO_distance[iter]<-dissimilarity_matrix_lv_who[iter,index_min_who_lv[iter]]
  
}


######## Now find the. NCIT matches
index_min_ncit_cosine <- as.matrix(apply(dissimilarity_matrix_cosine_ncit, 1, which.min))
index_min_ncit_jw <- as.matrix(apply(dissimilarity_matrix_jw_ncit, 1, which.min))
index_min_ncit_lv <- as.matrix(apply(dissimilarity_matrix_lv_ncit, 1, which.min))



ncit_match_cosine_df <- cbind(rownames(dissimilarity_matrix_cosine_ncit))
ncit_match_jw_df <- cbind(rownames(dissimilarity_matrix_jw_ncit))
ncit_match_lv_df <- cbind(rownames(dissimilarity_matrix_lv_ncit))

colnames(ncit_match_cosine_df)<-"Tumor_Names"
colnames(ncit_match_jw_df)<-"Tumor_Names"
colnames(ncit_match_lv_df)<-"Tumor_Names"

ncit_match_cosine_df <-as.data.frame(ncit_match_cosine_df)
ncit_match_jw_df <-as.data.frame(ncit_match_jw_df)
ncit_match_lv_df <-as.data.frame(ncit_match_lv_df)

ncit_match_cosine_df$NCIT_Matches<- NA
ncit_match_cosine_df$NCIT_distance<-NA

ncit_match_jw_df$NCIT_Matches<- NA
ncit_match_jw_df$NCIT_distance<-NA

ncit_match_lv_df$NCIT_Matches<- NA
ncit_match_lv_df$NCIT_distance<-NA

for (iter in 1: dim(ncit_match_cosine_df)[1]){
  
  ncit_match_cosine_df$NCIT_Matches[iter] <- colnames(dissimilarity_matrix_cosine_ncit)[index_min_ncit_cosine[iter]]
  ncit_match_cosine_df$NCIT_distance[iter]<-dissimilarity_matrix_cosine_ncit[iter,index_min_ncit_cosine[iter]]
  
  ncit_match_jw_df$NCIT_Matches[iter] <- colnames(dissimilarity_matrix_jw_ncit)[index_min_ncit_jw[iter]]
  ncit_match_jw_df$NCIT_distance[iter]<-dissimilarity_matrix_jw_ncit[iter,index_min_ncit_jw[iter]]
  
  ncit_match_lv_df$NCIT_Matches[iter] <- colnames(dissimilarity_matrix_lv_ncit)[index_min_ncit_lv[iter]]
  ncit_match_lv_df$NCIT_distance[iter]<-dissimilarity_matrix_lv_ncit[iter,index_min_ncit_lv[iter]]
  
}

nested_affinity_cluster_cosine<- nested_affinity_cluster_cosine %>% dplyr::left_join(who_match_cosine_df,by="Tumor_Names")
nested_affinity_cluster_cosine<- nested_affinity_cluster_cosine %>% dplyr::left_join(ncit_match_cosine_df,by="Tumor_Names")

nested_affinity_cluster_jw<- nested_affinity_cluster_jw %>% dplyr::left_join(who_match_jw_df,by="Tumor_Names")
nested_affinity_cluster_jw<- nested_affinity_cluster_jw %>% dplyr::left_join(ncit_match_jw_df,by="Tumor_Names")


nested_affinity_cluster_lv<- nested_affinity_cluster_lv %>% dplyr::left_join(who_match_lv_df,by="Tumor_Names")
nested_affinity_cluster_lv<- nested_affinity_cluster_lv %>% dplyr::left_join(ncit_match_lv_df,by="Tumor_Names")


## Cluster Label assignment
nested_affinity_cluster_cosine<- cluster_label_assignment(nested_affinity_cluster_cosine)
nested_affinity_cluster_jw<- cluster_label_assignment(nested_affinity_cluster_jw)
nested_affinity_cluster_lv<- cluster_label_assignment(nested_affinity_cluster_lv)


# Compute isolation forest for embedding based Kmeans

# Load embeddings 
disease_transform_ADA2<- read.csv(paste(intermediate_dir,"/disease_transform_pca.csv",sep="") )
colnames(disease_transform_ADA2)[1]<-"Tumor_Names"
rownames(disease_transform_ADA2)<-disease_transform_ADA2$Tumor_Names 


disease_transform_V3<- read.csv(paste(intermediate_dir,"/disease_transform_pca_v3.csv",sep="") )
colnames(disease_transform_V3)[1]<-"Tumor_Names"
rownames(disease_transform_V3)<-disease_transform_V3$Tumor_Names 


set.seed(13)

embedding_ADA2<-kmeans_clust_result_embedding_ADA2 %>% dplyr::left_join(disease_transform_ADA2,by="Tumor_Names")
kmeans_clust_result_embedding_ADA2$isolation_outlier_score<-NA


cluster_labels_ADA2 <- unique(embedding_ADA2$Cluster_ID)
for(iter in 1:length(cluster_labels_ADA2)){
  cluster_label_current <- cluster_labels_ADA2[iter]
  embedding_subset <- embedding_ADA2 %>% dplyr::filter(Cluster_ID==cluster_label_current)
  if(dim(embedding_subset)[1]>2){ # Need at least 2 data points to run isolation forest
    model <- isolation.forest(embedding_subset[1:nrow(embedding_subset),9:ncol(embedding_subset)], ndim=3, ntrees=100, nthreads=1) # ntrees 50 initially
    scores <- predict(model, embedding_subset[1:nrow(embedding_subset),9:ncol(embedding_subset)], type="score")
    ind_clust <- which(kmeans_clust_result_embedding_ADA2$Cluster_ID==cluster_label_current)
    kmeans_clust_result_embedding_ADA2$isolation_outlier_score[ind_clust]<-scores
  }else{
    ind_clust <- which(kmeans_clust_result_embedding_ADA2$Cluster_ID==cluster_label_current)
    kmeans_clust_result_embedding_ADA2$isolation_outlier_score[ind_clust]<-0
  }
}
kmeans_clust_result_embedding_ADA2<- kmeans_clust_result_embedding_ADA2 %>% dplyr::mutate(Isolation_Outlier = case_when(isolation_outlier_score>0.5 ~ "Yes", TRUE ~ "No"))
# LOF 
kmeans_clust_result_embedding_ADA2$LOF_Scores<-NA
lof_scores_minpts_list<-list()

for(iter in 1:length(cluster_labels_ADA2)){
  cluster_label_current <- cluster_labels_ADA2[iter]
  ind_clust <- which(kmeans_clust_result_embedding_ADA2$Cluster_ID==cluster_label_current)
  lof_scores_minpts_list<-list()
  
  embedding_subset <- embedding_ADA2 %>% dplyr::filter(Cluster_ID==cluster_label_current)
  if(dim(embedding_subset)[1]>2){ # Need at least 2 data points to run isolation forest
    min_pts<- 2:(dim(embedding_subset)[1]-1)
    for(iter_pts in min_pts){
      lof_scores_minpts <- lof(embedding_subset[,9:ncol(embedding_subset)],iter_pts)
      lof_scores_minpts_list[[as.character(iter_pts)]]<-lof_scores_minpts
    }
    lof_scores_minpts_list<- t(as.data.frame(lof_scores_minpts_list))
    lof_scores_minpts_list_median<-apply(lof_scores_minpts_list,2,median)
    kmeans_clust_result_embedding_ADA2$LOF_Scores[ind_clust]<-lof_scores_minpts_list_median
    
  }else{
    kmeans_clust_result_embedding_ADA2$LOF_Scores[ind_clust]<-0
  }
  
  
}

kmeans_clust_result_embedding_ADA2<- kmeans_clust_result_embedding_ADA2 %>% dplyr::mutate(LOF_Outlier = case_when(LOF_Scores>1 ~ "Yes", TRUE ~ "No"))


### Compute isolation forest for V3 embedding Kmeans

# Compute isolation forest

embedding_V3<-kmeans_clust_result_embedding_V3 %>% dplyr::left_join(disease_transform_V3,by="Tumor_Names")
kmeans_clust_result_embedding_V3$isolation_outlier_score<-NA

cluster_labels_V3 <- unique(embedding_V3$Cluster_ID)
for(iter in 1:length(cluster_labels_V3)){
  cluster_label_current <- cluster_labels_V3[iter]
  embedding_subset <- embedding_V3 %>% dplyr::filter(Cluster_ID==cluster_label_current)
  if(dim(embedding_subset)[1]>2){ # Need at least 2 data points to run isolation forest
    model <- isolation.forest(embedding_subset[1:nrow(embedding_subset),9:ncol(embedding_subset)], ndim=3, ntrees=100, nthreads=1) # ntrees 50 initially
    scores <- predict(model, embedding_subset[1:nrow(embedding_subset),9:ncol(embedding_subset)], type="score")
    ind_clust <- which(kmeans_clust_result_embedding_V3$Cluster_ID==cluster_label_current)
    kmeans_clust_result_embedding_V3$isolation_outlier_score[ind_clust]<-scores
  }else{
    ind_clust <- which(kmeans_clust_result_embedding_V3$Cluster_ID==cluster_label_current)
    kmeans_clust_result_embedding_V3$isolation_outlier_score[ind_clust]<-0
  }
}
kmeans_clust_result_embedding_V3<- kmeans_clust_result_embedding_V3 %>% dplyr::mutate(Isolation_Outlier = case_when(isolation_outlier_score>0.5 ~ "Yes", TRUE ~ "No"))


# Compute LOF 
kmeans_clust_result_embedding_V3$LOF_Scores<-NA
lof_scores_minpts_list<-list()

for(iter in 1:length(cluster_labels_V3)){
  cluster_label_current <- cluster_labels_V3[iter]
  ind_clust <- which(kmeans_clust_result_embedding_V3$Cluster_ID==cluster_label_current)
  lof_scores_minpts_list<-list()
  
  embedding_subset <- embedding_V3 %>% dplyr::filter(Cluster_ID==cluster_label_current)
  if(dim(embedding_subset)[1]>2){ # Need at least 2 data points to run isolation forest
    min_pts<- 2:(dim(embedding_subset)[1]-1)
    for(iter_pts in min_pts){
      lof_scores_minpts <- lof(embedding_subset[,9:ncol(embedding_subset)],iter_pts)
      lof_scores_minpts_list[[as.character(iter_pts)]]<-lof_scores_minpts
    }
    lof_scores_minpts_list<- t(as.data.frame(lof_scores_minpts_list))
    lof_scores_minpts_list_median<-apply(lof_scores_minpts_list,2,median)
    kmeans_clust_result_embedding_V3$LOF_Scores[ind_clust]<-lof_scores_minpts_list_median
    
  }else{
    kmeans_clust_result_embedding_V3$LOF_Scores[ind_clust]<-0
  }
  
  
}

kmeans_clust_result_embedding_V3<- kmeans_clust_result_embedding_V3 %>% dplyr::mutate(LOF_Outlier = case_when(LOF_Scores>1 ~ "Yes", TRUE ~ "No"))


#  Isolation forest and LOF analysis
nested_affinity_cluster_cosine<-outlier_detection_edit_dist(nested_affinity_cluster_cosine,dissimilarity_matrix_cosine)
nested_affinity_cluster_jw<-outlier_detection_edit_dist(nested_affinity_cluster_jw,dissimilarity_matrix_jw)
nested_affinity_cluster_lv<-outlier_detection_edit_dist(nested_affinity_cluster_lv,dissimilarity_matrix_lv)

# Cluster Reassignment 
nested_affinity_cluster_cosine_reassigned<-edit_distance_cluster_reassignment(nested_affinity_cluster_cosine)
nested_affinity_cluster_jw_reassigned<-edit_distance_cluster_reassignment(nested_affinity_cluster_jw)
nested_affinity_cluster_lv_reassigned<-edit_distance_cluster_reassignment(nested_affinity_cluster_lv)


# Sample 
affinity_cluster_ADA2_reassigned_df_short <- affinity_cluster_ADA2_reassigned_df %>% dplyr::select(Tumor_Names,who_cluster_label)
affinity_cluster_v3_reassigned_df_short <- affinity_cluster_v3_reassigned_df %>% dplyr::select(Tumor_Names,who_cluster_label)
kmeans_clust_result_embedding_ADA2_short <- kmeans_clust_result_embedding_ADA2 %>% dplyr::select(Tumor_Names,who_cluster_label)
kmeans_clust_result_embedding_V3_short <- kmeans_clust_result_embedding_V3 %>% dplyr::select(Tumor_Names,who_cluster_label)
nested_affinity_cluster_cosine_reassigned_short <- nested_affinity_cluster_cosine_reassigned %>% dplyr::select(Tumor_Names,who_cluster_label)
nested_affinity_cluster_jw_reassigned_short <- nested_affinity_cluster_jw_reassigned %>% dplyr::select(Tumor_Names,who_cluster_label)
nested_affinity_cluster_lv_reassigned_short <- nested_affinity_cluster_lv_reassigned %>% dplyr::select(Tumor_Names,who_cluster_label)

affinity_cluster_ADA2_dist_short<-affinity_cluster_ADA2_reassigned_df %>% dplyr::select(Tumor_Names,WHO_Matches)
affinity_cluster_v3_dist_short<-affinity_cluster_v3_reassigned_df %>% dplyr::select(Tumor_Names,WHO_Matches)


colnames(affinity_cluster_v3_reassigned_df_short)[2]<-"af_v3"
colnames(affinity_cluster_ADA2_reassigned_df_short)[2]<-"af_ada2"

colnames(kmeans_clust_result_embedding_V3_short)[2]<-"kmeans_v3"
colnames(kmeans_clust_result_embedding_ADA2_short)[2]<-"kmeans_ad2"

colnames(nested_affinity_cluster_cosine_reassigned_short)[2]<-"af_cosine"
colnames(nested_affinity_cluster_jw_reassigned_short)[2]<-"af_jw"
colnames(nested_affinity_cluster_lv_reassigned_short)[2]<-"af_lv"

colnames(affinity_cluster_v3_dist_short)[2]<-"euclidean_dist_v3"
colnames(affinity_cluster_ADA2_dist_short)[2]<-"euclidean_dist_ada2"

# Random samples
set.seed(13)
tumor_sample_df<-sample_n(affinity_cluster_v3_reassigned_df_short, 1000)
tumor_sample_df<-tumor_sample_df %>% dplyr::left_join(affinity_cluster_ADA2_reassigned_df_short,by="Tumor_Names")%>%
  dplyr::left_join(kmeans_clust_result_embedding_V3_short,by="Tumor_Names") %>%
  dplyr::left_join(kmeans_clust_result_embedding_ADA2_short,by="Tumor_Names") %>%
  dplyr::left_join(nested_affinity_cluster_cosine_reassigned_short,by="Tumor_Names") %>%
  dplyr::left_join(nested_affinity_cluster_jw_reassigned_short,by="Tumor_Names") %>%
  dplyr::left_join(nested_affinity_cluster_lv_reassigned_short,by="Tumor_Names") %>%
  dplyr::left_join(affinity_cluster_v3_dist_short,by="Tumor_Names") %>%
  dplyr::left_join(affinity_cluster_ADA2_dist_short,by="Tumor_Names")
  
  