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

# 
disease_transform <- read_csv(paste(intermediate_dir,"/disease_transform_pca.csv",sep=""))
colnames(disease_transform)[1]<-"Tumor_Name"
# PAM
pam_results<-fviz_nbclust_verbose(
  disease_transform,
  FUNcluster = pam,
  method = c("silhouette"),
  diss = NULL,
  k.max = 10,
  nboot = 100,
  verbose = interactive(),
  barfill = "steelblue",
  barcolor = "steelblue",
  linecolor = "steelblue",
  print.summary = TRUE,
)

pam_index_opt_clust<- which(pam_results$data$y==max(pam_results$data$y))
pam_opt_clust_size<- as.integer(pam_results$data$clusters[pam_index_opt_clust])

pam_results_cluster_vec<-fviz_nbclust_verbose(
  disease_transform,
  FUNcluster = pam,
  method = c("silhouette"),
  diss = NULL,
  k.max = 10,
  cluster_vec = c(1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000,13000),
  nboot = 100,
  verbose = interactive(),
  barfill = "steelblue",
  barcolor = "steelblue",
  linecolor = "steelblue",
  print.summary = TRUE,
)


# BIC
set.seed(13)

mclust_result <- Mclust(disease_transform, G=c(75,1599,2999,4899,12999))
mclust_best <- dim(mclust_result$z)[2]

cl <- makeCluster(4)
registerDoParallel(cl)
start_time <- Sys.time()
d_clust2_parallel <- foreach(i = 1:10,.combine = "cbind") %dopar% Mclust(disease_transform,G=i)





#hierarchical 
disease_transform_pca_scaled <- scale(disease_transform)
distance_pca_transformed <- dist(disease_transform_pca_scaled, method = "euclidean")
diseases_cluster_hierarchical <- hclust(distance_pca_transformed, method = "complete" )
diseases_clusterCut <- cutree(diseases_cluster_hierarchical, 4800)
diseases_clusterCut<-as.data.frame(diseases_clusterCut)


diseases_cluster_nb <- NbClust(data = disease_transform, diss = distance_pca_transformed,
                               distance=NULL, method="single", min.nc = 2, max.nc = 15)



# Gap Stats
gap.stat <- clusGap_verbose(disease_transform, FUNcluster = kmeans, K.max = 15)

gap_stat_hcut <- clusGap(disease_transform_pca_scaled, FUN = hcut, K.max = 2, B = 10)
k <- maxSE(hcluster$Tab[, "gap"], hcluster$Tab[, "SE.sim"], method="Tibs2001SEmax")

# Dunn Index
fviz_dunn <- function(data) {
  cluster_vec <- c(10,100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000,13000)
  dunnin <- c()
  for (i in 1:length(cluster_vec)) {
    print(cluster_vec[i])
    dunnin[i] <- dunn(distance = dist(data), clusters = kmeans(data, cluster_vec[i])$cluster)
  }
  plot(cluster_vec, dunnin, xlab =  "Cluster number k",
       ylab = "Dunn Index",
       main = "Dunn Plot", cex.main=1,
       col = "dodgerblue1", cex = 0.9 ,
       lty=1 , type="o" , lwd=1, pch=4,
       bty = "l",
       las = 1, cex.axis = 0.8, tcl  = -0.2)
  abline(v=which(dunnin==max(dunnin)) + 1, lwd=1, col="red", lty="dashed")
}

dunn_result <- fviz_dunn(disease_transform)


