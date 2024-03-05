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
cl <- makeCluster(6, outfile="")
registerDoParallel(cl)

CalculateEuclideanDistance <- function(vect1, vect2) sqrt(sum((vect1 - vect2)^2)) 
#stopCluster(cl)



###### Lymphoma Luek analysis

lymphoma_leukemia_strings <- c("leukemia", "lymphoma", "leukemias", "lymphomas", "leukaemia", "leukaemias",
                               "leuk","hematologic tumors","hemato")

affinity_cluster_hema_df <- affinity_cluster_df %>% dplyr::filter(str_detect(Tumor_Names,paste(lymphoma_leukemia_strings, collapse = "|")))

hema_label <- unique(affinity_cluster_hema_df$Cluster_ID)

affinity_cluster_hema_df <- affinity_cluster_df %>% dplyr::filter(Cluster_ID %in% hema_label )


#sanity check 
#t1 <-t1 %>% filter(Var1 %in% hema_label)
#t2<-as.data.frame(table(affinity_cluster_hema_df$Cluster_ID))
#t3<- t1 %>% dplyr::left_join(t2, by="Var1")
#print(identical(t3[['Freq.x']],t3[['Freq.y']]))

# compute distances to hemato 
combined_embedding_hema_df <- combined_embedding_df %>% filter(rownames(combined_embedding_df) %in% affinity_cluster_hema_df$Tumor_Names)


outer_who_final_hema<-foreach(i = 1:dim(combined_embedding_hema_df)[1], .combine = rbind) %dopar% { #03:53pm - 6:58 pm
  print(i)
  #s <- apply(WHO_embedding_df[,2:1537],1,CalculateEuclideanDistance,vect2=combined_embedding_df[i,])
  embedding_pairwise<- as.matrix(rbind(combined_embedding_hema_df[i,],WHO_embedding_df[,2:1537]))
  euclidean_dist <- as.matrix(dist(embedding_pairwise,method = "euclidean"))
  d<-as.double(euclidean_dist[1,c(-1)])
}
colnames(outer_who_final_hema)<-(WHO_embedding_df$Disease)
rownames(outer_who_final_hema)<-rownames(combined_embedding_hema_df)


index_min_who_hema <- as.matrix(apply(outer_who_final_hema, 1, which.min))

who_match_hema_df <- cbind(rownames(outer_who_final_hema))

colnames(who_match_hema_df)<-"Tumor_Names"
who_match_hema_df <-as.data.frame(who_match_hema_df)

who_match_hema_df$WHO_Matches<- NA
who_match_hema_df$WHO_distance<-NA


for (iter in 1: dim(who_match_hema_df)[1]){
  
  who_match_hema_df$WHO_Matches[iter] <- colnames(outer_who_final_hema)[index_min_who_hema[iter]]
  who_match_hema_df$WHO_distance[iter]<-outer_who_final_hema[iter,index_min_who_hema[iter]]
  
}
############ NCIT 

outer_NCIT_final_hema<-foreach(i = 1:dim(combined_embedding_hema_df)[1], .combine = rbind) %dopar% { #7:20pm - 
  print(i)
  embedding_pairwise<- as.matrix(rbind(combined_embedding_hema_df[i,],NCIT_embedding_df[,2:1537]))
  euclidean_dist <- as.matrix(dist(embedding_pairwise,method = "euclidean"))
  d<-as.double(euclidean_dist[1,c(-1)])
}
colnames(outer_NCIT_final_hema)<-(NCIT_embedding_df$Disease)
rownames(outer_NCIT_final_hema)<-rownames(combined_embedding_hema_df)


index_min_NCIT_hema <- as.matrix(apply(outer_NCIT_final_hema, 1, which.min))

NCIT_match_hema_df <- cbind(rownames(outer_NCIT_final_hema))

colnames(NCIT_match_hema_df)<-"Tumor_Names"
NCIT_match_hema_df <-as.data.frame(NCIT_match_hema_df)

NCIT_match_hema_df$NCIT_Matches<- NA
NCIT_match_hema_df$NCIT_distance<-NA

for (iter in 1: dim(NCIT_match_hema_df)[1]){
  
  NCIT_match_hema_df$NCIT_Matches[iter] <- colnames(outer_NCIT_final_hema)[index_min_NCIT_hema[iter]]
  NCIT_match_hema_df$NCIT_distance[iter]<-outer_NCIT_final_hema[iter,index_min_NCIT_hema[iter]]
  
}
NCIT_match_hema_df$NCIT_Matches<-tolower(NCIT_match_hema_df$NCIT_Matches)


affinity_cluster_hema_df<- affinity_cluster_hema_df %>% dplyr::left_join(who_match_hema_df,by="Tumor_Names")
affinity_cluster_hema_df<- affinity_cluster_hema_df %>% dplyr::left_join(NCIT_match_hema_df,by="Tumor_Names")

affinity_cluster_hema_df <- affinity_cluster_hema_df %>% dplyr::mutate(assigned_class = case_when(NCIT_distance < WHO_distance ~ NCIT_Matches,
                                                                                                  NCIT_distance > WHO_distance ~ WHO_Matches,
                                                                                                  TRUE ~ "Both"))

affinity_cluster_hema_df<- cluster_label_assignment(affinity_cluster_hema_df)

affinity_cluster_hema_df <- affinity_cluster_hema_df %>% dplyr::select(Tumor_Names,Cluster_ID,WHO_Matches,NCIT_Matches,suggested_cluster_label)


#### ISOLATION
affinity_cluster_outlier<-affinity_cluster_hema_df%>%dplyr::select(Tumor_Names,Cluster_ID)

disease_transform<- read.csv(paste(intermediate_dir,"/disease_transform_pca.csv",sep="") )
colnames(disease_transform)[1]<-"Tumor_Names"
rownames(disease_transform)<-disease_transform$Tumor_Names 


hemato_embedding<-affinity_cluster_outlier %>% dplyr::left_join(disease_transform,by="Tumor_Names")

affinity_cluster_outlier$isolation_outlier_score<-NA

hemato_cluster_labels <- unique(hemato_embedding$Cluster_ID)
for(iter in 1:length(hemato_cluster_labels)){
  cluster_label_current <- hemato_cluster_labels[iter]
  set.seed(13)
  hemato_embedding_subset <- hemato_embedding %>% dplyr::filter(Cluster_ID==cluster_label_current)
  if(dim(hemato_embedding_subset)[1]>2){ # Need at least 2 data points to run isolation forest
  model <- isolation.forest(hemato_embedding_subset[1:nrow(hemato_embedding_subset),3:ncol(hemato_embedding_subset)], ndim=3, ntrees=50, nthreads=1)
  scores <- predict(model, hemato_embedding_subset[1:nrow(hemato_embedding_subset),3:ncol(hemato_embedding_subset)], type="score")
  ind_clust <- which(affinity_cluster_outlier$Cluster_ID==cluster_label_current)
  affinity_cluster_outlier$isolation_outlier_score[ind_clust]<-scores
  }else{
    ind_clust <- which(affinity_cluster_outlier$Cluster_ID==cluster_label_current)
    affinity_cluster_outlier$isolation_outlier_score[ind_clust]<-0
  }
}
affinity_cluster_outlier<- affinity_cluster_outlier %>% dplyr::mutate(Isolation_Outlier = case_when(isolation_outlier_score>0.5 ~ "Yes", TRUE ~ "No"))


affinity_cluster_outlier$LOF_Scores<-NA
lof_scores_minpts_list<-list()

for(iter in 1:length(hemato_cluster_labels)){
  cluster_label_current <- hemato_cluster_labels[iter]
  ind_clust <- which(affinity_cluster_outlier$Cluster_ID==cluster_label_current)
  lof_scores_minpts_list<-list()
  
  set.seed(13)
  hemato_embedding_subset <- hemato_embedding %>% dplyr::filter(Cluster_ID==cluster_label_current)
  if(dim(hemato_embedding_subset)[1]>2){ # Need at least 2 data points to run isolation forest
    min_pts<- 2:(dim(hemato_embedding_subset)[1]-1)
    for(iter_pts in min_pts){
      lof_scores_minpts <- lof(hemato_embedding_subset[,3:137],iter_pts)
      lof_scores_minpts_list[[as.character(iter_pts)]]<-lof_scores_minpts
    }
    lof_scores_minpts_list<- t(as.data.frame(lof_scores_minpts_list))
    lof_scores_minpts_list_median<-apply(lof_scores_minpts_list,2,median)
    affinity_cluster_outlier$LOF_Scores[ind_clust]<-lof_scores_minpts_list_median
    
  }else{
    affinity_cluster_outlier$LOF_Scores[ind_clust]<-0
  }
  

}

affinity_cluster_outlier<- affinity_cluster_outlier %>% dplyr::mutate(LOF_Outlier = case_when(LOF_Scores>1 ~ "Yes", TRUE ~ "No"))

affinity_cluster_hema_df2<- affinity_cluster_hema_df %>% dplyr::select(Tumor_Names,assigned_class,Cluster_ID)
affinity_cluster_hema_df2 <- affinity_cluster_hema_df2 %>% filter(!Tumor_Names %in% unique(c(NCIT_embedding_df$Disease,WHO_embedding_df$Disease)))
affinity_cluster_hema_df2 <- affinity_cluster_hema_df2 %>% dplyr::left_join(affinity_cluster_outlier[,c(1)],by="Tumor_Names")





##### ISOLATION FOR HEMA DF2
affinity_cluster_outlier2<-affinity_cluster_hema_df2%>%dplyr::select(Tumor_Names,assigned_class)




hemato_embedding<-affinity_cluster_outlier2 %>% dplyr::left_join(disease_transform,by="Tumor_Names")

affinity_cluster_hema_df2$isolation_outlier_score<-NA

hemato_cluster_labels <- unique(hemato_embedding$Cluster_ID)
for(iter in 1:length(hemato_cluster_labels)){
  cluster_label_current <- hemato_cluster_labels[iter]
  set.seed(13)
  hemato_embedding_subset <- hemato_embedding %>% dplyr::filter(Cluster_ID==cluster_label_current)
  if(dim(hemato_embedding_subset)[1]>2){ # Need at least 2 data points to run isolation forest
    model <- isolation.forest(hemato_embedding_subset[1:nrow(hemato_embedding_subset),3:ncol(hemato_embedding_subset)], ndim=3, ntrees=50, nthreads=1)
    scores <- predict(model, hemato_embedding_subset[1:nrow(hemato_embedding_subset),3:ncol(hemato_embedding_subset)], type="score")
    ind_clust <- which(affinity_cluster_hema_df2$Cluster_ID==cluster_label_current)
    affinity_cluster_hema_df2$isolation_outlier_score[ind_clust]<-scores
  }else{
    ind_clust <- which(affinity_cluster_hema_df2$Cluster_ID==cluster_label_current)
    affinity_cluster_outlier2$isolation_outlier_score[ind_clust]<-0
  }
}
affinity_cluster_outlier2<- affinity_cluster_outlier2 %>% dplyr::mutate(Isolation_Outlier = case_when(isolation_outlier_score>0.5 ~ "Yes", TRUE ~ "No"))


affinity_cluster_outlier$LOF_Scores<-NA
lof_scores_minpts_list<-list()

for(iter in 1:length(hemato_cluster_labels)){
  cluster_label_current <- hemato_cluster_labels[iter]
  ind_clust <- which(affinity_cluster_outlier$Cluster_ID==cluster_label_current)
  lof_scores_minpts_list<-list()
  
  set.seed(13)
  hemato_embedding_subset <- hemato_embedding %>% dplyr::filter(Cluster_ID==cluster_label_current)
  if(dim(hemato_embedding_subset)[1]>2){ # Need at least 2 data points to run isolation forest
    min_pts<- 2:(dim(hemato_embedding_subset)[1]-1)
    for(iter_pts in min_pts){
      lof_scores_minpts <- lof(hemato_embedding_subset[,3:137],iter_pts)
      lof_scores_minpts_list[[as.character(iter_pts)]]<-lof_scores_minpts
    }
    lof_scores_minpts_list<- t(as.data.frame(lof_scores_minpts_list))
    lof_scores_minpts_list_median<-apply(lof_scores_minpts_list,2,median)
    affinity_cluster_outlier$LOF_Scores[ind_clust]<-lof_scores_minpts_list_median
    
  }else{
    affinity_cluster_outlier$LOF_Scores[ind_clust]<-0
  }
  
  
}






write.csv(affinity_cluster_hema_df,"hemato_tumor.csv")

stopCluster(cl)
save.image("script07.RData")