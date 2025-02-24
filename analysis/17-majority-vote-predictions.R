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
  library(infotheo)
  library(combinat)
  library(ggplot2)
  library(reshape2)
})

# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
intermediate_dir <- file.path(analysis_dir,"intermediate")
intermediate_dir_5th <- file.path(analysis_dir,"intermediate_5th")
result_dir <- file.path(analysis_dir,"results")
result_dir_5th <- file.path(analysis_dir,"results_5th")


# Load the Result files
WHO_Results_all_ose <- read.csv(paste(result_dir,"/WHO_Results_all_ose.csv",sep=""))
WHO_Results_all_ose<-WHO_Results_all_ose[,c(-1)]
WHO_Results_5thed_ose <- read.csv(paste(result_dir_5th,"/WHO_Results_5thed_ose.csv",sep=""))
WHO_Results_5thed_ose<-WHO_Results_5thed_ose[,c(-1)]

NCIT_Results_all_ose <- read.csv(paste(result_dir,"/NCIT_Results_all_ose.csv",sep=""))
NCIT_Results_all_ose<-NCIT_Results_all_ose[,c(-1)]

NCIT_Results_5thed_ose <- read.csv(paste(result_dir_5th,"/NCIT_Results_5thed_ose.csv",sep=""))
NCIT_Results_5thed_ose<-NCIT_Results_5thed_ose[,c(-1)]

WHO_Results_5thed_ose$majority_vote_predictions<-NA

## WHO 5th edition

for(iter in 1:nrow(WHO_Results_5thed_ose)){

  kmeans_v3_pred<- WHO_Results_5thed_ose$kmeans_v3[iter]
  af_ada2_pred<-WHO_Results_5thed_ose$af_ada2[iter]
  minilm_pred <-WHO_Results_5thed_ose$euclidean_dist_MiniLM_L12_v2[iter]
  
  table_pred<-as.data.frame(table(c(kmeans_v3_pred,af_ada2_pred,minilm_pred)))
  most_frequent <- as.character(table_pred$Var1[table_pred$Freq== max(table_pred$Freq)])
  
  if(length(most_frequent)==1){
    WHO_Results_5thed_ose$majority_vote_predictions[iter]=most_frequent
  }else{
    WHO_Results_5thed_ose$majority_vote_predictions[iter]=WHO_Results_5thed_ose$euclidean_dist_MiniLM_L12_v2[iter]
  }
  
  
}


## WHO all edition
WHO_Results_all_ose$majority_vote_predictions<-NA

for(iter in 1:nrow(WHO_Results_all_ose)){
  
  kmeans_v3_pred<- WHO_Results_all_ose$kmeans_v3[iter]
  kmeans_ada2_pred<-WHO_Results_all_ose$kmeans_ada2[iter]
  minilm_pred <-WHO_Results_all_ose$euclidean_dist_MiniLM_L12_v2[iter]
  
  table_pred<-as.data.frame(table(c(kmeans_v3_pred,kmeans_ada2_pred,minilm_pred)))
  most_frequent <- as.character(table_pred$Var1[table_pred$Freq== max(table_pred$Freq)])
  
  if(length(most_frequent)==1){
    WHO_Results_all_ose$majority_vote_predictions[iter]=most_frequent
  }else{
    WHO_Results_all_ose$majority_vote_predictions[iter]=WHO_Results_all_ose$euclidean_dist_MiniLM_L12_v2[iter]
  }
  
  
}



# Repeat the analysis for NCIT
NCIT_Results_5thed_ose$majority_vote_predictions<-NA

for(iter in 1:nrow(NCIT_Results_5thed_ose)){
  
  kmeans_v3_pred<- NCIT_Results_5thed_ose$kmeans_v3[iter]
  af_ada2_pred<-NCIT_Results_5thed_ose$af_ada2[iter]
  minilm_pred <-NCIT_Results_5thed_ose$euclidean_dist_MiniLM_L12_v2[iter]
  
  table_pred<-as.data.frame(table(c(kmeans_v3_pred,af_ada2_pred,minilm_pred)))
  most_frequent <- as.character(table_pred$Var1[table_pred$Freq== max(table_pred$Freq)])
  
  
  if(length(most_frequent)==1){
    NCIT_Results_5thed_ose$majority_vote_predictions[iter]=most_frequent
  }else{
    NCIT_Results_5thed_ose$majority_vote_predictions[iter]=NCIT_Results_5thed_ose$euclidean_dist_MiniLM_L12_v2[iter]
  }
  
  
}


# Repeat the analysis for NCIT
NCIT_Results_all_ose$majority_vote_predictions<-NA

for(iter in 1:nrow(NCIT_Results_all_ose)){
  
  kmeans_v3_pred<- NCIT_Results_all_ose$kmeans_v3[iter]
  kmeans_ada2_pred<-NCIT_Results_all_ose$kmeans_ada2[iter]
  minilm_pred <-NCIT_Results_all_ose$euclidean_dist_MiniLM_L12_v2[iter]
  
  table_pred<-as.data.frame(table(c(kmeans_v3_pred,kmeans_ada2_pred,minilm_pred)))
  most_frequent <- as.character(table_pred$Var1[table_pred$Freq== max(table_pred$Freq)])
  
  if(length(most_frequent)==1){
    NCIT_Results_all_ose$majority_vote_predictions[iter]=most_frequent
  }else{
    NCIT_Results_all_ose$majority_vote_predictions[iter]=NCIT_Results_all_ose$euclidean_dist_MiniLM_L12_v2[iter]
  }
  
}


# Write the results
write.csv(WHO_Results_5thed_ose,paste(result_dir_5th,"/WHO_Results_5thed_ose.csv",sep = ""))
write.csv(WHO_Results_all_ose,paste(result_dir,"/WHO_Results_all_ose.csv",sep = ""))

write.csv(NCIT_Results_5thed_ose,paste(result_dir_5th,"/NCIT_Results_5thed_ose.csv",sep = ""))
write.csv(NCIT_Results_all_ose,paste(result_dir,"/NCIT_Results_all_ose.csv",sep = ""))



