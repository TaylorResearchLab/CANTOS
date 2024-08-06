# 
# Performs outlier detection and cluster label assignment for Kmeans V3 and ADA2
# Assigns the closest matching label to edit distance affinity clusters
# Samples 1600 tumors and generates file for validation
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(tidyverse)
  library(stringi)
  library(readxl)
  library(dbscan)
  library(isotree)
  library(qdapRegex)
  library(ghql)
})

setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
intermediate_dir <- file.path(analysis_dir,"intermediate_5th")
result_dir <- file.path(analysis_dir,"results_5th")

source(paste(util_dir,"/cluster_label_assignment_refined.R",sep=""))
source(paste(util_dir,"/outlier_detection_edit_dist.R",sep=""))
source(paste(util_dir,"/edit_distance_cluster_reassignment.R",sep=""))


# Read file v3
affinity_cluster_v3_reassigned_df<-read.csv(paste(intermediate_dir,"/affinity_cluster_v3_reassigned_df_5thed.csv",sep=""))
affinity_cluster_v3_reassigned_df<-affinity_cluster_v3_reassigned_df[,c(-1)]
# Read ADA2 
affinity_cluster_ADA2_reassigned_df<-read.csv(paste(intermediate_dir,"/affinity_cluster_ADA2_reassigned_df_5thed.csv",sep=""))
affinity_cluster_ADA2_reassigned_df<-affinity_cluster_ADA2_reassigned_df[,c(-1)]
#Read Kmeans 
kmeans_clust_result_embedding_ADA2 <- read_csv(paste(result_dir,"/kmeans_clust_result_embedding_ada2_5thed.csv",sep=""))
kmeans_clust_result_embedding_ADA2<-kmeans_clust_result_embedding_ADA2[,c(-1)]
kmeans_clust_result_embedding_V3 <- read_csv(paste(result_dir,"/kmeans_clust_result_embedding_v3_5thed.csv",sep=""))
kmeans_clust_result_embedding_V3<-kmeans_clust_result_embedding_V3[,c(-1)]

# Read edit distance cluster
nested_affinity_cluster_cosine <- read_csv(paste(result_dir,"/nested_affinity_cluster_cosine.csv",sep=""))
nested_affinity_cluster_jw <- read_csv(paste(result_dir,"/nested_affinity_cluster_jw.csv",sep=""))
nested_affinity_cluster_lv <- read_csv(paste(result_dir,"/nested_affinity_cluster_lv.csv",sep=""))



nested_affinity_cluster_cosine<-nested_affinity_cluster_cosine[,c(-1,-8)]
nested_affinity_cluster_jw<-nested_affinity_cluster_jw[,c(-1,-8)]
nested_affinity_cluster_lv<-nested_affinity_cluster_lv[,c(-1,-8)]


# Adjust NCT_IDs
affinity_cluster_v3_reassigned_df<-affinity_cluster_v3_reassigned_df[,c(-1)]
affinity_cluster_ADA2_reassigned_df<-affinity_cluster_ADA2_reassigned_df[,c(-1)]
nested_affinity_cluster_cosine<-nested_affinity_cluster_cosine[,c(-1)]
nested_affinity_cluster_jw<-nested_affinity_cluster_jw[,c(-1)]
nested_affinity_cluster_lv<-nested_affinity_cluster_lv[,c(-1)]

tumor_nct_map <- kmeans_clust_result_embedding_ADA2 %>% dplyr::select(nct_id,Tumor_Names)

affinity_cluster_v3_reassigned_df<-affinity_cluster_v3_reassigned_df %>% left_join(tumor_nct_map,by="Tumor_Names")
affinity_cluster_ADA2_reassigned_df<-affinity_cluster_ADA2_reassigned_df %>% left_join(tumor_nct_map,by="Tumor_Names")
nested_affinity_cluster_cosine<-nested_affinity_cluster_cosine %>% left_join(tumor_nct_map,by="Tumor_Names")
nested_affinity_cluster_jw<-nested_affinity_cluster_jw %>% left_join(tumor_nct_map,by="Tumor_Names")
nested_affinity_cluster_lv<-nested_affinity_cluster_lv %>% left_join(tumor_nct_map,by="Tumor_Names")

affinity_cluster_v3_reassigned_df<-affinity_cluster_v3_reassigned_df[,c(10,1:9)]
affinity_cluster_ADA2_reassigned_df<-affinity_cluster_ADA2_reassigned_df[,c(10,1:9)]

nested_affinity_cluster_cosine<-nested_affinity_cluster_cosine[,c(6,1:5)]
nested_affinity_cluster_jw<-nested_affinity_cluster_jw[,c(6,1:5)]
nested_affinity_cluster_lv<-nested_affinity_cluster_lv[,c(6,1:5)]

# Find the WHO and NCIT closest matches 
who_ncit_match_ADA2 <- affinity_cluster_ADA2_reassigned_df %>% dplyr::select(Tumor_Names,WHO_Matches,WHO_distance,NCIT_Matches,NCIT_distance) 
who_ncit_match_v3 <- affinity_cluster_v3_reassigned_df %>% dplyr::select(Tumor_Names,WHO_Matches,WHO_distance,NCIT_Matches,NCIT_distance) 

# Join the matches to Kmeans 
kmeans_clust_result_embedding_ADA2<-kmeans_clust_result_embedding_ADA2 %>% dplyr::select(nct_id,Tumor_Names,cluster)
kmeans_clust_result_embedding_V3<-kmeans_clust_result_embedding_V3 %>% dplyr::select(nct_id,Tumor_Names,cluster)

colnames(kmeans_clust_result_embedding_ADA2)[3]<-c("Cluster_ID")
colnames(kmeans_clust_result_embedding_V3)[3]<-c("Cluster_ID")

kmeans_clust_result_embedding_ADA2<-kmeans_clust_result_embedding_ADA2 %>% dplyr::left_join(who_ncit_match_ADA2,by="Tumor_Names")
kmeans_clust_result_embedding_V3<-kmeans_clust_result_embedding_V3 %>% dplyr::left_join(who_ncit_match_v3,by="Tumor_Names")

kmeans_clust_result_embedding_ADA2<- cluster_label_assignment_refined(kmeans_clust_result_embedding_ADA2)
kmeans_clust_result_embedding_V3<- cluster_label_assignment_refined(kmeans_clust_result_embedding_V3)


# Join the matches to affinity
nested_affinity_cluster_cosine<-nested_affinity_cluster_cosine %>% dplyr::select(nct_id,Tumor_Names,Cluster_ID)
nested_affinity_cluster_jw<-nested_affinity_cluster_jw %>% dplyr::select(nct_id,Tumor_Names,Cluster_ID)
nested_affinity_cluster_lv<-nested_affinity_cluster_lv %>% dplyr::select(nct_id,Tumor_Names,Cluster_ID)

dissimilarity_matrix_cosine<-load(paste(intermediate_dir,"/dissimilarity_matrix_cosine.RData",sep=""))
dissimilarity_matrix_jw<-load(paste(intermediate_dir,"/dissimilarity_matrix_jw.RData",sep=""))
dissimilarity_matrix_lv<-load(paste(intermediate_dir,"/dissimilarity_matrix_lv.RData",sep=""))


dissimilarity_matrix_cosine<-as.data.frame(dissimilarity_matrix_cosine)
dissimilarity_matrix_jw<-as.data.frame(dissimilarity_matrix_jw)
dissimilarity_matrix_lv<-as.data.frame(dissimilarity_matrix_lv)




NCIT_embedding_df <-read.csv(paste(data_dir,"/dt_input_file_6_dec/NCIT_Neoplasm_Core_terms_text-embedding-ada-002_embeddings.csv",sep=""))
NCIT_embedding_df<-NCIT_embedding_df[c(-1),] # Remove the header (column name) embedding

rownames(NCIT_embedding_df)<-NULL

NCIT_Tumors<-tolower(NCIT_embedding_df$Disease)
rm(NCIT_embedding_df)

WHO_Terms_All <-readxl::read_xlsx(paste(data_dir,"/WHO_Tumors/result/WHO_Tumor_all_edition.xlsx",sep=""))
WHO_Terms_5th<-WHO_Terms_All%>%filter(edition_5th=="Yes")


dissimilarity_matrix_cosine_who <- dissimilarity_matrix_cosine %>% dplyr::select(one_of(WHO_Terms_5th$Tumor_Names))
dissimilarity_matrix_jw_who <- dissimilarity_matrix_jw %>% dplyr::select(one_of(WHO_Terms_5th$Tumor_Names))
dissimilarity_matrix_lv_who <- dissimilarity_matrix_lv %>% dplyr::select(one_of(WHO_Terms_5th$Tumor_Names))

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
nested_affinity_cluster_cosine<- cluster_label_assignment_refined(nested_affinity_cluster_cosine)
nested_affinity_cluster_jw<- cluster_label_assignment_refined(nested_affinity_cluster_jw)
nested_affinity_cluster_lv<- cluster_label_assignment_refined(nested_affinity_cluster_lv)

# Compute isolation forest for embedding based Kmeans

# Load embeddings 
disease_transform_ADA2<- read.csv(paste(intermediate_dir,"/disease_transform_pca_ada2_5thed.csv",sep="") )
colnames(disease_transform_ADA2)[1]<-"Tumor_Names"
rownames(disease_transform_ADA2)<-disease_transform_ADA2$Tumor_Names 


disease_transform_V3<- read.csv(paste(intermediate_dir,"/disease_transform_pca_v3_5thed.csv",sep="") )
colnames(disease_transform_V3)[1]<-"Tumor_Names"
rownames(disease_transform_V3)<-disease_transform_V3$Tumor_Names 


set.seed(13)

embedding_ADA2<-kmeans_clust_result_embedding_ADA2 %>% dplyr::left_join(disease_transform_ADA2,by="Tumor_Names")
kmeans_clust_result_embedding_ADA2$isolation_outlier_score<-NA
idx_ada2<-which(colnames(embedding_ADA2)=="PC1")


cluster_labels_ADA2 <- unique(embedding_ADA2$Cluster_ID)
for(iter in 1:length(cluster_labels_ADA2)){
  cluster_label_current <- cluster_labels_ADA2[iter]
  embedding_subset <- embedding_ADA2 %>% dplyr::filter(Cluster_ID==cluster_label_current)
  if(dim(embedding_subset)[1]>2){ # Need at least 2 data points to run isolation forest
    model <- isolation.forest(embedding_subset[1:nrow(embedding_subset),idx_ada2:ncol(embedding_subset)], ndim=3, ntrees=100, nthreads=1) # ntrees 50 initially
    scores <- predict(model, embedding_subset[1:nrow(embedding_subset),idx_ada2:ncol(embedding_subset)], type="score")
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
      lof_scores_minpts <- lof(embedding_subset[,idx_ada2:ncol(embedding_subset)],iter_pts)
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
idx_v3<-which(colnames(embedding_V3)=="PC1")

cluster_labels_V3 <- unique(embedding_V3$Cluster_ID)
for(iter in 1:length(cluster_labels_V3)){
  cluster_label_current <- cluster_labels_V3[iter]
  embedding_subset <- embedding_V3 %>% dplyr::filter(Cluster_ID==cluster_label_current)
  if(dim(embedding_subset)[1]>2){ # Need at least 2 data points to run isolation forest
    model <- isolation.forest(embedding_subset[1:nrow(embedding_subset),idx_v3:ncol(embedding_subset)], ndim=3, ntrees=100, nthreads=1) # ntrees 50 initially
    scores <- predict(model, embedding_subset[1:nrow(embedding_subset),idx_v3:ncol(embedding_subset)], type="score")
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
      lof_scores_minpts <- lof(embedding_subset[,idx_v3:ncol(embedding_subset)],iter_pts)
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

nested_affinity_cluster_cosine_reassigned<-nested_affinity_cluster_cosine_reassigned %>% left_join(tumor_nct_map,by="Tumor_Names")
nested_affinity_cluster_jw_reassigned<-nested_affinity_cluster_jw_reassigned %>% left_join(tumor_nct_map,by="Tumor_Names")
nested_affinity_cluster_lv_reassigned<-nested_affinity_cluster_lv_reassigned %>% left_join(tumor_nct_map,by="Tumor_Names")

# Sample 
affinity_cluster_ADA2_reassigned_df_short <- affinity_cluster_ADA2_reassigned_df %>% dplyr::select(nct_id,Tumor_Names,who_cluster_label)
affinity_cluster_v3_reassigned_df_short <- affinity_cluster_v3_reassigned_df %>% dplyr::select(nct_id,Tumor_Names,who_cluster_label)
kmeans_clust_result_embedding_ADA2_short <- kmeans_clust_result_embedding_ADA2 %>% dplyr::select(nct_id,Tumor_Names,who_cluster_label)
kmeans_clust_result_embedding_V3_short <- kmeans_clust_result_embedding_V3 %>% dplyr::select(nct_id,Tumor_Names,who_cluster_label)
nested_affinity_cluster_cosine_reassigned_short <- nested_affinity_cluster_cosine_reassigned %>% dplyr::select(nct_id,Tumor_Names,who_cluster_label)
nested_affinity_cluster_jw_reassigned_short <- nested_affinity_cluster_jw_reassigned %>% dplyr::select(nct_id,Tumor_Names,who_cluster_label)
nested_affinity_cluster_lv_reassigned_short <- nested_affinity_cluster_lv_reassigned %>% dplyr::select(nct_id,Tumor_Names,who_cluster_label)

affinity_cluster_ADA2_dist_short<-affinity_cluster_ADA2_reassigned_df %>% dplyr::select(nct_id,Tumor_Names,WHO_Matches)
affinity_cluster_v3_dist_short<-affinity_cluster_v3_reassigned_df %>% dplyr::select(nct_id,Tumor_Names,WHO_Matches)

colnames(affinity_cluster_v3_reassigned_df_short)[3]<-"af_v3"
colnames(affinity_cluster_ADA2_reassigned_df_short)[3]<-"af_ada2"

colnames(kmeans_clust_result_embedding_V3_short)[3]<-"kmeans_v3"
colnames(kmeans_clust_result_embedding_ADA2_short)[3]<-"kmeans_ada2"

colnames(nested_affinity_cluster_cosine_reassigned_short)[3]<-"af_cosine"
colnames(nested_affinity_cluster_jw_reassigned_short)[3]<-"af_jw"
colnames(nested_affinity_cluster_lv_reassigned_short)[3]<-"af_lv"

colnames(affinity_cluster_v3_dist_short)[3]<-"euclidean_dist_v3"
colnames(affinity_cluster_ADA2_dist_short)[3]<-"euclidean_dist_ada2"

# Closest Cosine , LV, JW
min_dist_matches<- as.data.frame(affinity_cluster_ADA2_reassigned_df$Tumor_Names)
min_dist_matches$cosine_match<-NA
min_dist_matches$jw_match<-NA
min_dist_matches$lv_match<-NA

colnames(min_dist_matches)[1]<-"Tumor_Names"
for (iter in 1: dim(min_dist_matches)[1]){
  print(iter)
  ind_row_cosine <- which(rownames(dissimilarity_matrix_cosine_who)==min_dist_matches$Tumor_Names[iter])
  ind_col_cosine <- which(dissimilarity_matrix_cosine_who[ind_row_cosine,]==min(dissimilarity_matrix_cosine_who[ind_row_cosine,]))
  
  
  
  ind_row_jw <- which(rownames(dissimilarity_matrix_jw_who)==min_dist_matches$Tumor_Names[iter])
  ind_col_jw <- which(dissimilarity_matrix_jw_who[ind_row_jw,]==min(dissimilarity_matrix_jw_who[ind_row_jw,]))
  
  ind_row_lv <- which(rownames(dissimilarity_matrix_lv_who)==min_dist_matches$Tumor_Names[iter])
  ind_col_lv<- which(dissimilarity_matrix_lv_who[ind_row_lv,]==min(dissimilarity_matrix_lv_who[ind_row_lv,]))
  
  if(length(ind_col_cosine)>1){
    min_dist_matches$cosine_match[iter]<- paste(unique(colnames(dissimilarity_matrix_cosine_who)[ind_col_cosine]), collapse =";") 
    
  }else{
    min_dist_matches$cosine_match[iter]<- colnames(dissimilarity_matrix_cosine_who)[ind_col_cosine]
    
  }
  
  if(length(ind_col_jw)>1){
    min_dist_matches$jw_match[iter]<- paste(unique(colnames(dissimilarity_matrix_jw_who)[ind_col_jw]), collapse =";") 
    
  }else{
    min_dist_matches$jw_match[iter]<- colnames(dissimilarity_matrix_jw_who)[ind_col_jw]
    
  }
  
  if(length(ind_col_lv)>1){
    min_dist_matches$lv_match[iter]<- paste(unique(colnames(dissimilarity_matrix_lv_who)[ind_col_lv]), collapse =";") 
    
  }else{
    min_dist_matches$lv_match[iter]<- colnames(dissimilarity_matrix_lv_who)[ind_col_lv]
    
  }
  
  
  
}

affinity_cluster_v3_reassigned_df_short<-affinity_cluster_v3_reassigned_df_short%>%dplyr::filter(nct_id!="NA")

tumor_sample_df_all_editions<- read.csv(paste(analysis_dir,"/results/tumor_sample_df_script10.csv",sep=""))

tumor_sample_df<-affinity_cluster_v3_reassigned_df_short %>% filter(Tumor_Names %in% tumor_sample_df_all_editions$Tumor_Names)


tumor_sample_df<-tumor_sample_df %>% dplyr::left_join(affinity_cluster_ADA2_reassigned_df_short,by="Tumor_Names")%>%
  dplyr::left_join(kmeans_clust_result_embedding_V3_short,by="Tumor_Names") %>%
  dplyr::left_join(kmeans_clust_result_embedding_ADA2_short,by="Tumor_Names") %>%
  dplyr::left_join(nested_affinity_cluster_cosine_reassigned_short,by="Tumor_Names") %>%
  dplyr::left_join(nested_affinity_cluster_jw_reassigned_short,by="Tumor_Names") %>%
  dplyr::left_join(nested_affinity_cluster_lv_reassigned_short,by="Tumor_Names") %>%
  dplyr::left_join(affinity_cluster_v3_dist_short,by="Tumor_Names") %>%
  dplyr::left_join(affinity_cluster_ADA2_dist_short,by="Tumor_Names")%>%
  dplyr::left_join(min_dist_matches,by="Tumor_Names") %>% 
  dplyr::select(nct_id,Tumor_Names,af_v3,af_ada2,kmeans_v3,kmeans_ada2,af_cosine,af_jw,af_lv,
                euclidean_dist_v3,euclidean_dist_ada2,cosine_match,jw_match,lv_match)


tumor_sample_df<- add_column(tumor_sample_df,valid_af_v3="", .after = "af_v3")
tumor_sample_df<- add_column(tumor_sample_df,valid_af_ad2="", .after = "af_ada2")
tumor_sample_df<- add_column(tumor_sample_df,valid_kmeans_v3="", .after = "kmeans_v3")
tumor_sample_df<- add_column(tumor_sample_df,valid_kmeans_ada2="", .after = "kmeans_ada2")
tumor_sample_df<- add_column(tumor_sample_df,valid_af_cosine="", .after = "af_cosine")
tumor_sample_df<- add_column(tumor_sample_df,valid_af_jw="", .after = "af_jw")
tumor_sample_df<- add_column(tumor_sample_df,valid_af_lv="", .after = "af_lv")
tumor_sample_df<- add_column(tumor_sample_df,valid_euclidean_dist_v3="", .after = "euclidean_dist_v3")
tumor_sample_df<- add_column(tumor_sample_df,valid_euclidean_dist_ada2="", .after = "euclidean_dist_ada2")
tumor_sample_df<- add_column(tumor_sample_df,valid_cosine_match="", .after = "cosine_match")
tumor_sample_df<- add_column(tumor_sample_df,valid_jw_match="", .after = "jw_match")
tumor_sample_df<- add_column(tumor_sample_df,valid_lv_match="", .after = "lv_match")

tumor_sample_df<-tumor_sample_df[order(tumor_sample_df$Tumor_Names),]
rownames(tumor_sample_df)<-NULL
save.image("script10_aug6_5thed.RData")
