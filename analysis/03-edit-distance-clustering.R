# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(doParallel)
  library(foreach)
  library(apcluster)
  library(stringdist)
  library(tidyverse)
  library(magrittr)
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
source(paste(util_dir,"/nested_clust_edit_dist.R",sep=""))
source(paste(util_dir,"/compute_silhouette.R",sep=""))



# Load data

load(paste(intermediate_dir,"/dissimilarity_matrix_lv.RData",sep=""))
load(paste(intermediate_dir,"/dissimilarity_matrix_jw.RData",sep=""))
load(paste(intermediate_dir,"/dissimilarity_matrix_cosine.RData",sep=""))


cluster_results_lv<- read.csv(paste(results_dir,"/cluster_lv.csv",sep=""))
cluster_results_jw<- read.csv(paste(results_dir,"/cluster_jw.csv",sep=""))
cluster_results_cosine<- read.csv(paste(results_dir,"/cluster_cosine.csv",sep=""))


# Compute Simmilarity matrix for each edit distance
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


stopCluster(cl)

######### Cluster with LV ########
apclust_lv <- apcluster(simmilarity_matrix_lv) #10:26pm -12:45 pm
affinity_cluster_lv_df<-as.data.frame(matrix(nrow=1,ncol=2))
colnames(affinity_cluster_lv_df)<-c("Tumor_Names","Cluster_ID")
for (iter in 1: length(apclust_lv@clusters)){
  affinity_cluster_lv_df[iter,1] <- paste(names(unlist(apclust_lv@clusters[iter])),collapse = "@")
  affinity_cluster_lv_df[iter,2] <- iter
}
affinity_cluster_lv_df<- affinity_cluster_lv_df %>% separate_rows(Tumor_Names, sep = '@')

######### Cluster with jw ########
apclust_jw <- apcluster(simmilarity_matrix_jw) #8:32 pm -10:06 pm
affinity_cluster_jw_df<-as.data.frame(matrix(nrow=1,ncol=2))
colnames(affinity_cluster_jw_df)<-c("Tumor_Names","Cluster_ID")
for (iter in 1: length(apclust_jw@clusters)){
  affinity_cluster_jw_df[iter,1] <- paste(names(unlist(apclust_jw@clusters[iter])),collapse = "@")
  affinity_cluster_jw_df[iter,2] <- iter
}
affinity_cluster_jw_df<- affinity_cluster_jw_df %>% separate_rows(Tumor_Names, sep = '@')

######### Cluster with cosine ########
apclust_cosine <- apcluster(simmilarity_matrix_cosine)#11:07 am - 1:56 pm
affinity_cluster_cosine_df<-as.data.frame(matrix(nrow=1,ncol=2))
colnames(affinity_cluster_cosine_df)<-c("Tumor_Names","Cluster_ID")
for (iter in 1: length(apclust_cosine@clusters)){
  affinity_cluster_cosine_df[iter,1] <- paste(names(unlist(apclust_cosine@clusters[iter])),collapse = "@")
  affinity_cluster_cosine_df[iter,2] <- iter
}
affinity_cluster_cosine_df<- affinity_cluster_cosine_df %>% separate_rows(Tumor_Names, sep = '@')

save.image(file = "editdistancecluster.RData")

### Nested LV Clustering

nested_affinity_cluster_lv <- nested_clust_edit_dist(n=3,affinity_cluster_df = affinity_cluster_lv_df,
                                                     dist_mat = simmilarity_matrix_lv)

### Nested JW Clustering

nested_affinity_cluster_jw <- nested_clust_edit_dist(n=3,affinity_cluster_df = affinity_cluster_jw_df,
                                                     dist_mat = simmilarity_matrix_jw)
### Nested Cosine Clustering

nested_affinity_cluster_cosine <- nested_clust_edit_dist(n=3,affinity_cluster_df = affinity_cluster_cosine_df,
                                                     dist_mat = simmilarity_matrix_cosine)

# Compute Silo Dist
nested_affinity_cluster_lv<-compute_silhouette(cluster_df = nested_affinity_cluster_lv,dist_mat = dissimilarity_matrix_lv)
nested_affinity_cluster_jw<-compute_silhouette(cluster_df = nested_affinity_cluster_jw,dist_mat = dissimilarity_matrix_jw)
nested_affinity_cluster_cosine<-compute_silhouette(cluster_df = nested_affinity_cluster_cosine,dist_mat = dissimilarity_matrix_cosine)

# Select benchmarks 
benchmark_tumors <- c("b cell lymphoma", "neuroblastoma", "triple negative breast cancer",
                      "unresectable lung carcinoma", "liposarcoma","cancer of the liver",
                      "smoldering myeloma")

subcluster_lv <- nested_affinity_cluster_lv$SubsetCluster_IDs[nested_affinity_cluster_lv$Tumor_Names %in% benchmark_tumors]
benchmark_aff_clust_lv <- nested_affinity_cluster_lv %>% dplyr::filter(SubsetCluster_IDs %in% subcluster_lv)

subcluster_jw <- nested_affinity_cluster_jw$SubsetCluster_IDs[nested_affinity_cluster_jw$Tumor_Names %in% benchmark_tumors]
benchmark_aff_clust_jw <- nested_affinity_cluster_jw %>% dplyr::filter(SubsetCluster_IDs %in% subcluster_jw)

subcluster_cosine <- nested_affinity_cluster_cosine$SubsetCluster_IDs[nested_affinity_cluster_cosine$Tumor_Names %in% benchmark_tumors]
benchmark_aff_clust_cosine <- nested_affinity_cluster_cosine %>% dplyr::filter(SubsetCluster_IDs %in% subcluster_cosine)


display_table_benchmark <- display_table_benchmark %>% dplyr::select(Tumors,cluster_lv,cluster_jw,cluster_cosine)



## Write workspace
save.image(file = "editdistancecluster.RData")

write.csv(benchmark_aff_clust_lv,paste(results_dir,"/benchmark_aff_clust_lv.csv",sep=""))
write.csv(benchmark_aff_clust_jw,paste(results_dir,"/benchmark_aff_clust_jw.csv",sep=""))
write.csv(benchmark_aff_clust_cosine,paste(results_dir,"/benchmark_aff_clust_cosine.csv",sep=""))


write.csv(nested_affinity_cluster_lv,paste(results_dir,"/nested_affinity_cluster_lv.csv",sep=""))
write.csv(nested_affinity_cluster_jw,paste(results_dir,"/nested_affinity_cluster_jw.csv",sep=""))
write.csv(nested_affinity_cluster_cosine,paste(results_dir,"/nested_affinity_cluster_cosine.csv",sep=""))