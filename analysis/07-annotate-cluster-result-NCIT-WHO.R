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

outer_who_final<-foreach(i = 1:dim(combined_embedding_df)[1], .combine = rbind) %dopar% { #12:10pm -
  print(i)
  #s <- apply(WHO_embedding_df[,2:1537],1,CalculateEuclideanDistance,vect2=combined_embedding_df[i,])
  embedding_pairwise<- as.matrix(rbind(combined_embedding_df[i,],WHO_embedding_df[,2:1537]))
  euclidean_dist <- as.matrix(dist(embedding_pairwise,method = "euclidean"))
  d<-as.double(euclidean_dist[1,c(-1)])
}
colnames(outer_who_final)<-(WHO_embedding_df$Disease)
rownames(outer_who_final)<-rownames(combined_embedding_df)


outer_NCIT_final<-foreach(i = 1:dim(combined_embedding_df)[1], .combine = rbind) %dopar% { #7:20pm - 
  print(i)
  embedding_pairwise<- as.matrix(rbind(combined_embedding_df[i,],NCIT_embedding_df[,2:1537]))
  euclidean_dist <- as.matrix(dist(embedding_pairwise,method = "euclidean"))
  d<-as.double(euclidean_dist[1,c(-1)])
}

colnames(outer_NCIT_final)<-(tolower(NCIT_embedding_df$Disease))
rownames(outer_NCIT_final)<-rownames(combined_embedding_df)




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
##########################
affinity_cluster_df<- affinity_cluster_df %>% dplyr::left_join(who_match_df,by="Tumor_Names")
affinity_cluster_df<- affinity_cluster_df %>% dplyr::left_join(NCIT_match_df,by="Tumor_Names")

affinity_cluster_df <- affinity_cluster_df %>% dplyr::mutate(assigned_class = case_when(NCIT_distance < WHO_distance ~ NCIT_Matches,
                                                                                                  NCIT_distance > WHO_distance ~ WHO_Matches,
                                                                                                  TRUE ~ "Both"))

affinity_cluster_df<- cluster_label_assignment(affinity_cluster_df)


write.csv(affinity_cluster_df,"affinity_cluster_assignment.csv")



###### Lymphoma Luek analysis

lymphoma_leukemia_strings <- c("leukemia", "lymphoma", "leukemias", "lymphomas", "leukaemia", "leukaemias",
                               "leuk","hematologic tumors","hemato")

affinity_cluster_hema_df <- affinity_cluster_df %>% dplyr::filter(str_detect(Tumor_Names,paste(lymphoma_leukemia_strings, collapse = "|")))

hema_label <- unique(affinity_cluster_hema_df$Cluster_ID)

affinity_cluster_hema_df <- affinity_cluster_df %>% dplyr::filter(Cluster_ID %in% hema_label )


#sanity check 
t1 <-t1 %>% filter(Var1 %in% hema_label)
t2<-as.data.frame(table(affinity_cluster_hema_df$Cluster_ID))
t3<- t1 %>% dplyr::left_join(t2, by="Var1")
print(identical(t3[['Freq.x']],t3[['Freq.y']]))

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

write.csv(affinity_cluster_hema_df,"hemato_tumor.csv")

stopCluster(cl)
save.image("script-NCIT-WHO.RData")