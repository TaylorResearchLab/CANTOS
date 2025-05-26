# Annotate Affinity cluster results of V3
start_time_1<-Sys.time()

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
run_time_analysis<-file.path(analysis_dir,"run-time-analysis")
intermediate_dir <- file.path(run_time_analysis,"intermediate_5th")
results_dir <- file.path(run_time_analysis,"results_5th")

source(paste(util_dir,"/cluster_label_assignment_refined.R",sep=""))

print("CPT1")
load(paste(intermediate_dir,"/affinity_cluster_v3_df_5thed.RData",sep=""))
load(paste(intermediate_dir,"/combined_embedding_v3_df_5thed.RData",sep=""))


WHO_embedding_df<-read.csv(paste(data_dir,"/WHO_Terms_All_V3.csv",sep=""))
NCIT_embedding_df<-read.csv(paste(data_dir,"/NCIT_Embeddings_V3.csv",sep=""))

WHO_embedding_df<-WHO_embedding_df[,c(-1)]
NCIT_embedding_df<-NCIT_embedding_df[,c(-1)]

WHO_Terms_All <-readxl::read_xlsx(paste(data_dir,"/WHO_Tumors/result/WHO_Tumor_all_edition.xlsx",sep=""))
WHO_Terms_5th<-WHO_Terms_All%>%filter(edition_5th=="Yes")


#######
cl <- makeCluster(5, outfile="")
registerDoParallel(cl)
who_mat <- as.matrix(WHO_embedding_df[, 2:ncol(WHO_embedding_df)])
query_mat <- as.matrix(combined_embeddings_df)
print("CPT2")

outer_who_final <- foreach(i = 1:nrow(query_mat), .combine = rbind) %dopar% {
  print(i)
  v <- query_mat[i, ]
  
  diff_sq <- sweep(who_mat, 2, v)^2
  dists <- sqrt(rowSums(diff_sq))
  
  matrix(dists, nrow = 1)
}

colnames(outer_who_final)<-(WHO_embedding_df$Tumor_Names)
rownames(outer_who_final)<-rownames(combined_embeddings_df)


outer_who_final<-as.data.frame(outer_who_final)
outer_who_final<-outer_who_final %>% dplyr::select(any_of(WHO_Terms_5th$Tumor_Names))
outer_who_final<-as.matrix(outer_who_final)

stop_time_1<-Sys.time()
save.image(paste(run_time_analysis,"/LTE3_AP/07B-LTE_5th.RData",sep=""))

start_time_2<-Sys.time()


# 

print("CPT3")

ncit_mat <- as.matrix(NCIT_embedding_df[, 2:ncol(NCIT_embedding_df)])
query_mat <- as.matrix(combined_embeddings_df)

outer_NCIT_final <- foreach(i = 1:nrow(query_mat), .combine = rbind) %dopar% {
  print(i)
  v <- query_mat[i, ]
  
  diff_sq <- sweep(ncit_mat, 2, v)^2
  dists <- sqrt(rowSums(diff_sq))
  
  matrix(dists, nrow = 1)
}
colnames(outer_NCIT_final)<-tolower(NCIT_embedding_df$Tumor_Names)
rownames(outer_NCIT_final)<-rownames(combined_embeddings_df)

stop_time_2<-Sys.time()
save.image(paste(run_time_analysis,"/LTE3_AP/07B-LTE_5th.RData",sep=""))


start_time_3<-Sys.time()
print("CPT4")

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
#affinity_cluster_v3_df<-affinity_cluster_v3_df[,c(-1)]
tumor_id<- read.csv(paste(data_dir,"/Tumor_NCT_ID.csv",sep=""))
tumor_id<- tumor_id[,c(-1)]
affinity_cluster_v3_df<-affinity_cluster_v3_df%>%left_join(tumor_id,by="Tumor_Names")
affinity_cluster_v3_df<-affinity_cluster_v3_df[,c(9,1:8)]

stop_time_3<-Sys.time()
save.image(paste(run_time_analysis,"/LTE3_AP/07B-LTE_5th.RData",sep=""))

elapsed_time_1<- as.numeric(difftime(stop_time_1, start_time_1, units = "secs"))
elapsed_time_2<- as.numeric(difftime(stop_time_2, start_time_2, units = "secs"))
elapsed_time_3<- as.numeric(difftime(stop_time_3, start_time_3, units = "secs"))

elapsed_time<-elapsed_time_3+elapsed_time_2+elapsed_time_1

objs <- ls(envir = .GlobalEnv)

# Get the size of each object
sizes <- sapply(objs, function(x) object.size(get(x, envir = .GlobalEnv)))

# Convert sizes to readable format and summarize
sizes_df <- data.frame(
  Object = objs,
  Size_MB = round(sizes / (1024^2), 3)
)

# Sort by size descending
sizes_df <- sizes_df[order(-sizes_df$Size_MB), ]

# Print the total memory used
cat("Total memory used in Global Environment:", round(sum(sizes) / (1024^2), 3), "MB\n")

Runtime_7B_5th = round(elapsed_time, 2)
Memory_gb_7b_5th=round(sum(sizes) / (1024^2), 3)/1000

print(Runtime_7B_5th)
print(Memory_gb_7b_5th)






save.image(paste(run_time_analysis,"/LTE3_AP/07B-LTE_5th.RData",sep=""))

##
write.csv(affinity_cluster_v3_df,paste(intermediate_dir,"/affinity_cluster_v3_df_5thed.csv",sep=""))

