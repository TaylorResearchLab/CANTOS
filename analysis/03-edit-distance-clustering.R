# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(doParallel)
  library(foreach)
  library(apcluster)
  library(stringdist)
})

# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
intermediate_dir <- file.path(analysis_dir,"intermediate")
results_dir <- file.path(analysis_dir,"results")

source(paste(util_dir,"/string_normalizing.R",sep=""))
# Load data

load(paste(intermediate_dir,"/dissimilarity_matrix_lv.RData",sep=""))
load(paste(intermediate_dir,"/dissimilarity_matrix_jw.RData",sep=""))
load(paste(intermediate_dir,"/dissimilarity_matrix_cosine.RData",sep=""))


cluster_results_lv<- read.csv(paste(results_dir,"/cluster_lv.csv",sep=""))
cluster_results_jw<- read.csv(paste(results_dir,"/cluster_jw.csv",sep=""))
cluster_results_cosine<- read.csv(paste(results_dir,"/cluster_cosine.csv",sep=""))


# Kmeans 

simmilarity_matrix_cosine = 1 - dissimilarity_matrix_cosine
simmilarity_matrix_jw = 1-dissimilarity_matrix_jw



cl <- makeCluster(5, outfile="")
registerDoParallel(cl)
df_tumor_names<- colnames(dissimilarity_matrix_lv)
normalizing_matrix_lv <- as.data.frame(matrix(nrow=length(df_tumor_names),ncol=length(df_tumor_names)))
rownames(normalizing_matrix_lv)<-df_tumor_names
colnames(normalizing_matrix_lv)<-df_tumor_names


normalizing_matrix_lv<-foreach(iter=1:length(df_tumor_names),.combine=rbind) %dopar% {
  print(iter)
  disease_name <- colnames(dissimilarity_matrix_lv)[iter]
  norm_factors<-unlist(lapply(df_tumor_names,string_normalzing,S2=disease_name))
}
rownames(normalizing_matrix_lv) <- df_tumor_names
colnames(normalizing_matrix_lv) <- df_tumor_names

simmilarity_matrix_lv <- 1- (dissimilarity_matrix_lv/normalizing_matrix_lv)






apres1b <- apcluster(1-dissimilarity_matrix_jw[1:100,1:100])
affinity_cluster_df<-as.data.frame(matrix(nrow=1,ncol=2))
colnames(affinity_cluster_df)<-c("Tumor_Names","Cluster_ID")
for (iter in 1: length(apres1b@clusters)){
  affinity_cluster_df[iter,1] <- paste(names(unlist(apres1b@clusters[iter])),collapse = "@")
  affinity_cluster_df[iter,2] <- iter
}
affinity_cluster_df<- affinity_cluster_df %>% separate_rows(Tumor_Names, sep = '@')
