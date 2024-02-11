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
})

# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
intermediate_dir <- file.path(analysis_dir,"intermediate")

#load combined embedding
combined_embedding_df <- read.csv(paste(intermediate_dir,"/combined_embedding_df.csv",sep="") )

source(paste(util_dir,"/nearest_ncit.R",sep=""))


# Compute distances 

tumor_distances_df <- as.data.frame(matrix(nrow=dim(combined_embedding_df)[1],ncol = 5))
colnames(tumor_distances_df)<- c("Tumor_Names","NCIT_Match","NCIT_Distance","WHO_Match","WHO_Distance")



tumor_distances_df$Tumor_Names<- rownames(combined_embedding_df)

for(iter in 1:dim(tumor_distances_df)[1]){ #dim(tumor_distances_df)[1] 2445 
  print(iter)
  ncit_info <- nearest_ncit(combined_embedding_df[iter,1:1536],NCIT_embedding_df) 
  tumor_distances_df$NCIT_Match[iter]<- ncit_info[[1]]
  tumor_distances_df$NCIT_Distance[iter]<- ncit_info[[2]]
  
}


cl <- makeCluster(5, outfile="")
registerDoParallel(cl)
iterations <- 18104
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


stopCluster(cl)


distance_ncit <- as.data.frame(matrix(nrow=dim(combined_embedding_df)[1], ncol=dim(NCIT_embedding_df)[1]))
rownames(distance_ncit)<-rownames(combined_embedding_df)
colnames(distance_ncit)<-(NCIT_embedding_df$Disease)

CalculateEuclideanDistance <- function(vect1, vect2) sqrt(sum((vect1 - vect2)^2)) 

for(iter in 1:dim(distance_ncit)[1]){ # Stopped at 1895
  distance_ncit[iter,]<-apply(NCIT_embedding_df[,2:1537],1,CalculateEuclideanDistance,vect2=combined_embedding_df[iter,])
  print(iter)
}



outer_ncit_all<-invisible(foreach(i = 1:18014, .combine = rbind, .options.snow = opts, .verbose = T) %dopar% {
  s <- apply(NCIT_embedding_df[,2:1537],1,CalculateEuclideanDistance,vect2=combined_embedding_df[i,])
  setTxtProgressBar(pb, i)
  return(s)
})


outer_ncit_all_df<-outer_ncit_all[1:18014,]
rownames(outer_ncit_all_df) <- rownames(combined_embedding_df)
colnames(outer_ncit_all_df) <- NCIT_embedding_df$Disease

index_min <- as.matrix(apply(outer_ncit_all_df, 1, which.min))
ncit_match_df <- cbind(rownames(outer_ncit_all_df))

colnames(ncit_match_df)<-"Tumor_Names"
ncit_match_df <-as.data.frame(ncit_match_df)

ncit_match_df$NCIT_Matches<- NA
ncit_match_df$ncit_distance<-NA

for (iter in 1: dim(ncit_match_df)[1]){
  ncit_match_df$NCIT_Matches[iter] <- colnames(outer_ncit_all_df)[index_min[iter]]
  ncit_match_df$ncit_distance[iter]<-outer_ncit_all_df[iter,index_min[iter]]
}

affinity_cluster_annotation2<-affinity_cluster_annotation 
affinity_cluster_annotation2<- affinity_cluster_annotation2 %>% dplyr::left_join(ncit_match_df,by="Tumor_Names")

###################### WHO MATCHES
outer_who_final<-invisible(foreach(i = 1:18014, .combine = rbind, .options.snow = opts, .verbose = T) %dopar% {
  s <- apply(WHO_embedding_df[,2:1537],1,CalculateEuclideanDistance,vect2=combined_embedding_df[i,])
  setTxtProgressBar(pb, i)
  return(s)
})

outer_who_final_df<-outer_who_final[1:18014,]
rownames(outer_who_final_df) <- rownames(combined_embedding_df)
colnames(outer_who_final_df) <- WHO_embedding_df$Disease 

index_min_who <- as.matrix(apply(outer_who_final_df, 1, which.min))
who_match_df <- cbind(rownames(outer_who_final_df))

colnames(who_match_df)<-"Tumor_Names"
who_match_df <-as.data.frame(who_match_df)

who_match_df$WHO_Matches<- NA
who_match_df$WHO_distance<-NA


for (iter in 1: dim(who_match_df)[1]){
  who_match_df$WHO_Matches[iter] <- colnames(outer_who_final_df)[index_min_who[iter]]
  who_match_df$WHO_distance[iter]<-outer_who_final_df[iter,index_min_who[iter]]
  
}

affinity_cluster_annotation2<- affinity_cluster_annotation2 %>% dplyr::left_join(who_match_df,by="Tumor_Names")

affinity_cluster_annotation2 <- affinity_cluster_annotation2 %>% dplyr::mutate(assigned_class = case_when(ncit_distance < WHO_distance ~ NCIT_Matches,
                                                                                                          ncit_distance > WHO_distance ~ WHO_Matches,
                                                                                                          TRUE ~ "Both"))
affinity_cluster_annotation3 <- affinity_cluster_annotation2 %>% dplyr::select(Tumor_Names,Pediatric_SubsetCluster_ID,SubsetCluster_IDs,NCIT_Tumor,WHO_Tumor,assigned_class)




#
write.csv(who_match_df,paste(intermediate_dir,"/who_ct_distance_mat.csv",sep=""))
write.csv(ncit_match_df,paste(intermediate_dir,"/ncit_match_df.csv",sep=""))
write.csv(affinity_cluster_annotation3,paste(intermediate_dir,"/affinity_cluster_annotation3.csv",sep=""))