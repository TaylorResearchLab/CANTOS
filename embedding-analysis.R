# Load libraries
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
  library(factoextra)
  library(NbClust)
  library(dbscan)
  library(apcluster)
  library(doParallel)
  library(foreach)
  library(doSNOW)
  library(tcltk)
})

# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
source(paste(util_dir,"/nearest_ncit.R",sep=""))

# Read CT embedding file 
embedding_df <- read.csv(paste(analyses_dir,"/embedding-analysis-dt/disease_embeddings.csv",sep=""))
embedding_df<-embedding_df[order(embedding_df$Disease),]

# Read missing CT embeddings 
CT_missing_embedding_df <-read.csv(paste(analyses_dir,"/embedding-analysis-dt/dt_input_file_6_dec/missed_CT_result_embeddings.csv",sep=""))

# Combine all embeddings 
#colnames(embedding_df_agg)<-colnames(CT_missing_embedding_df)
#rownames(embedding_df_agg)<-NULL
CT_embedding_df <- rbind(embedding_df,CT_missing_embedding_df)
CT_embedding_df <- unique(CT_embedding_df) # STILL CONTAINS DISEASES WITH EXACT STRING NAME BUT SLIGHTLY VARYING EMBEDDINGS

CT_embedding_df <- CT_embedding_df %>% group_by(Disease) %>% summarise(across(everything(), list(mean)))

# Read CT diseases to remove disease names corrupted added due to csv 
CT_only_df <- read.csv(paste(data_dir,"/who_classification/CT_Only.csv",sep=""))
CT_embedding_agg_df <- CT_only_df %>% left_join(CT_embedding_df, by=c("DISEASE_NAMES"="Disease"))
CT_embedding_agg_df<-CT_embedding_agg_df[order(CT_embedding_agg_df$DISEASE_NAMES),] # Order by names so the rows with missing values are next to their matched rows 

ind_missing_commas <- which(is.na(CT_embedding_agg_df[,2:1537]), arr.ind=TRUE)
ind_missing_commas <- unique(ind_missing_commas[,c(-2)])
ind_missing_commas_iteration_2 <-list()

for (iter in ind_missing_commas){
  print(iter)
  disease_name <- CT_embedding_agg_df$DISEASE_NAMES[iter]
  disease_name_comma_removed <- gsub(",","",disease_name)
  ind_replacement <- which(CT_missing_embedding_df$Disease==disease_name_comma_removed)
  if(length(ind_replacement)>0){
    CT_embedding_agg_df[iter,2:1537] <- CT_missing_embedding_df[ind_replacement,2:1537]
  }
  else{
    disease_name_comma_removed <- gsub(","," ",disease_name)
    ind_replacement <- which(CT_missing_embedding_df$Disease==disease_name_comma_removed)
    if(length(ind_replacement)>0){
      CT_embedding_agg_df[iter,2:1537] <- CT_missing_embedding_df[ind_replacement,2:1537]
    }else{
      ind_missing_commas_iteration_2 <- c(ind_missing_commas_iteration_2,iter)
    }
    
  }
}
ind_missing_commas_iteration_2 <- unlist(ind_missing_commas_iteration_2)
tumor_list <- CT_missing_embedding_df$Disease
tumor_list <- gsub(",","",tumor_list)
tumor_list <- gsub(" ","",tumor_list)

for(iter in ind_missing_commas_iteration_2){
  disease_name <- CT_embedding_agg_df$DISEASE_NAMES[iter]
  disease_name_comma_removed <- gsub(",","",disease_name)
  disease_name_comma_removed<- gsub(" ","",disease_name_comma_removed)
  ind_replacement <- which(tumor_list==disease_name_comma_removed)
  if(length(ind_replacement)>0){
    CT_embedding_agg_df[iter,2:1537]  <-CT_missing_embedding_df[ind_replacement,2:1537]
  }
}


ind_missing_commas_iteration_3 <-  which(is.na(CT_embedding_agg_df[,2:1537]), arr.ind=TRUE)
ind_missing_commas_iteration_3 <- unique(ind_missing_commas_iteration_3[,1])
CT_embedding_agg_df[ind_missing_commas_iteration_3[1],2:1537]<-CT_missing_embedding_df[1085,2:1537]
CT_embedding_agg_df[ind_missing_commas_iteration_3[2],2:1537]<-CT_missing_embedding_df[1228,2:1537]

print(length(which(is.na(CT_embedding_agg_df[,2:1537]), arr.ind=TRUE)))

CT_embedding_agg_df<-unique(CT_embedding_agg_df) # HERE START Again

# Read NCIT Terms and WHO Terms with embedding
NCIT_embedding_df <-read.csv(paste(analyses_dir,"/embedding-analysis-dt/dt_input_file_6_dec/NCIT_Neoplasm_Core_terms_text-embedding-ada-002_embeddings.csv",sep=""))
WHO_embedding_df <-read.csv(paste(analyses_dir,"/embedding-analysis-dt/dt_input_file_6_dec/WHO_Only_terms_text-embedding-ada-002_embeddings.csv",sep=""))

NCIT_embedding_df<-NCIT_embedding_df[c(-1),] # Remove the header (column name) embedding
WHO_embedding_df<-WHO_embedding_df[c(-1),] # Remove the header (column name) embedding

rownames(NCIT_embedding_df)<-NULL
rownames(WHO_embedding_df)<-NULL

# Combined Embedding File for PC and Cluster Analysis
colnames(CT_embedding_agg_df)<-colnames(NCIT_embedding_df)
combined_embedding_df <- rbind(CT_embedding_agg_df,NCIT_embedding_df,WHO_embedding_df)
combined_embedding_df <- as.data.frame(combined_embedding_df %>% group_by(Disease) %>% summarise(across(everything(), list(mean))))
colnames(combined_embedding_df)<-colnames(NCIT_embedding_df)
rownames(combined_embedding_df)<-combined_embedding_df$Disease
combined_embedding_df<-combined_embedding_df[,c(-1)]
# 
# # Read embedding file
# embedding_df <- read.csv(paste(analyses_dir,"/embedding-analysis-dt/disease_embeddings.csv",sep=""))
# embedding_df<-embedding_df[order(embedding_df$Disease),]
# 
# # Remove repeats
# embedding_df_agg<-embedding_df %>% group_by(Disease) %>% summarise(across(everything(), list(mean)))
# rownames(embedding_df_agg)<-embedding_df_agg$Disease
# embedding_df_agg<-as.data.frame(embedding_df_agg)
# rownames(embedding_df_agg)<-embedding_df_agg$Disease
# embedding_df_agg<-embedding_df_agg[,c(-1)]
# colnames(embedding_df_agg)<-colnames(embedding_df)[2:1537]

# Compute PCs
results <- prcomp(combined_embedding_df, scale = TRUE)
eigs <- results$sdev^2
ff<-rbind(SD = sqrt(eigs),Proportion = eigs/sum(eigs), Cumulative = cumsum(eigs)/sum(eigs))

# Check the number of components needed to capture 80% variance at least
#print(sum(ff[2,1:125]))
print(sum(ff[2,1:138]))

# Select the top 138 PCs
disease_transform = as.data.frame(-results$x[,1:138])

# Find optimal number of Clusters using KMeans Silhouette 
cluster_results<-fviz_nbclust(disease_transform, kmeans, method = 'silhouette',  k.max = 5000,iter.max=50)
cluster_results_verbose<-fviz_nbclust_verbose(disease_transform, kmeans, method = 'silhouette',  k.max = 13433,iter.max=50)
index_opt_clust<- which(cluster_results$data$y==max(cluster_results$data$y))
opt_clust_size<- as.integer(cluster_results$data$clusters[index_opt_clust]) # 4800
kmeans_disease = kmeans(disease_transform, centers = opt_clust_size, nstart = 100)
diseases_cluster <- as.data.frame(kmeans_disease$cluster)


## CHI Index
CH_Results<-CHCriterion(disease_transform_pca_scaled, kmax=13434,clustermethod="hclust", method = "average")
CH_scroes <- as.data.frame(CH_Results$data$CHIndex) # ratio of the between-cluster variance and the within-cluster variance
WSS_scores<- as.data.frame(CH_Results$data$wss)
CH_Results<-CHCriterion(disease_transform_pca_scaled, kmax=13434,clustermethod="hclust", method = "average")

# DB Scan
disease_transform_subset <- disease_transform[1:100,]
hdb <- hdbscan(disease_transform_subset, minPts = 2)
#opt <- optics(disease_transform_subset, eps = 1, minPts = 4)
#opt <- extractDBSCAN(opt, eps_cl = 0.4)
#plot(opt)
db <- dbscan(disease_transform_subset, eps = 0.4, minPts = 4)


#affinity cluster
set.seed(13)
d.apclus2 <- apcluster(negDistMat(r=2), disease_transform)
cat("affinity propogation optimal number of clusters:", length(d.apclus2@clusters), "\n") #1255

affinity_cluster_df<-as.data.frame(matrix(nrow=1,ncol=2))
colnames(affinity_cluster_df)<-c("Tumor_Names","Cluster_ID")
for (iter in 1: length(d.apclus2@clusters)){
  affinity_cluster_df[iter,1] <- paste(names(unlist(d.apclus2@clusters[iter])),collapse = "@")
  affinity_cluster_df[iter,2] <- iter
}
affinity_cluster_df<- affinity_cluster_df %>% separate_rows(Tumor_Names, sep = '@')


# Find cluster membership
#affinity_cluster_df$Cluster_Total_Members <- NA
count_table <- as.data.frame(table(affinity_cluster_df$Cluster_ID))
colnames(count_table)<- c("Cluster_ID","Primary_Cluster_Frequency")
count_table$Cluster_ID<- as.numeric(count_table$Cluster_ID)
affinity_cluster_df <- affinity_cluster_df %>% dplyr::left_join(count_table,by="Cluster_ID")


# Find children and pediatric cluster
ind_ped <- c(which(str_detect(affinity_cluster_df$Tumor_Names, "childhood")),
             which(str_detect(affinity_cluster_df$Tumor_Names, "children")),
             which(str_detect(affinity_cluster_df$Tumor_Names, "child")),
             which(str_detect(affinity_cluster_df$Tumor_Names, "pediatric")),
             which(str_detect(affinity_cluster_df$Tumor_Names, "paediatric")))

affinity_cluster_df$contains_pediatric_string <- "No"

affinity_cluster_df$contains_pediatric_string[ind_ped]<-"Yes"

pediatric_cluster_ID<-0

affinity_cluster_df$Pediatric_SubsetCluster_ID <- NA  
affinity_cluster_df$Pediatric_SubsetCluster_ID[ind_ped]<-pediatric_cluster_ID


#Nested affinity clusters
barplot(height=affinity_cluster_df$Primary_Cluster_Frequency, names=affinity_cluster_df$Cluster_ID, 
        xlab="cluster_id", 
        ylab="frequency", 
        main="Cluster_Frequency", 
        ylim=c(0,50)
)

large_cluster_index <- which(affinity_cluster_df$Primary_Cluster_Frequency > median(affinity_cluster_df$Primary_Cluster_Frequency))

affinity_cluster_df$High_Primary_Cluster_Membership<-"No"
affinity_cluster_df$High_Primary_Cluster_Membership[large_cluster_index]<-"Yes"

affinity_cluster_df <- affinity_cluster_df %>% dplyr::select(Tumor_Names,Cluster_ID,Primary_Cluster_Frequency,
                                                             High_Primary_Cluster_Membership,contains_pediatric_string,
                                                             Pediatric_SubsetCluster_ID)
Clusters_Names <- unique(affinity_cluster_df$Cluster_ID)

disease_transform$Tumor_Name<-rownames(disease_transform)

affinity_cluster_df$SubsetCluster_IDs <- NA



for (iter in 1:length(Clusters_Names)) {

  print(iter)
  is_large_cluster <- unique(affinity_cluster_df$High_Primary_Cluster_Membership[which(affinity_cluster_df$Cluster_ID==
                                                                                  Clusters_Names[iter])])
  
  if(is_large_cluster=="Yes"){
    
    subset_embedding_df <- as.data.frame(affinity_cluster_df$Tumor_Names[affinity_cluster_df$Cluster_ID==Clusters_Names[iter]])
    colnames(subset_embedding_df)<-"Tumor_Name"
    rownames(subset_embedding_df)<-subset_embedding_df$Tumor_Name
    subset_embedding_df<- subset_embedding_df %>% dplyr::left_join(disease_transform,by="Tumor_Name")
    rownames(subset_embedding_df)<-subset_embedding_df$Tumor_Name
    subset_embedding_df<-subset_embedding_df[,c(-1)]
    
    affinity_subset <- apcluster(negDistMat(r=2), subset_embedding_df)
    cat("affinity propogation optimal number of clusters:", length(affinity_subset@clusters), "\n")
    
    subset_affinity_df<-as.data.frame(affinity_subset@idx)
    subset_affinity_df<-as.data.frame(matrix(nrow=1,ncol=2))
    colnames(subset_affinity_df)<-c("Tumor_Names","SubCluster_ID")
    for (iter_subset in 1: length(affinity_subset@clusters)){
      subset_affinity_df[iter_subset,1] <- paste(names(unlist(affinity_subset@clusters[iter_subset])),collapse = "@")
      subset_affinity_df[iter_subset,2] <- iter_subset
    }
    subset_affinity_df<- subset_affinity_df %>% separate_rows(Tumor_Names, sep = '@')
    largest_cluster_id <- max(subset_affinity_df$SubCluster_ID)
    power_of_ten <- floor(log10(largest_cluster_id)) + 1
    subset_affinity_df$SubCluster_ID<-subset_affinity_df$SubCluster_ID/(10^power_of_ten)+Clusters_Names[iter]
    
    for (iter_affinity_cluser in 1: dim(subset_affinity_df)[1]){
      ind_location <- which (affinity_cluster_df$Tumor_Names==subset_affinity_df$Tumor_Names[iter_affinity_cluser])
      affinity_cluster_df$SubsetCluster_IDs[ind_location]<-subset_affinity_df$SubCluster_ID[iter_affinity_cluser]
    }
    
    
  }
}

ind_subcluster_na <- which(is.na(affinity_cluster_df$SubsetCluster_IDs),arr.ind = TRUE)

affinity_cluster_df$SubsetCluster_IDs[ind_subcluster_na] <- affinity_cluster_df$Cluster_ID[ind_subcluster_na]

count_subset_table <- as.data.frame(table(affinity_cluster_df$SubsetCluster_IDs))
colnames(count_subset_table) <- c("SubsetCluster_IDs", "Subset_Cluster_Freq")
count_subset_table$SubsetCluster_IDs <- as.numeric(as.character(sub("," , ".", count_subset_table$SubsetCluster_IDs)))

#count_subset_table$SubCluster_IDs<- as.numeric(as.character(count_subset_table$SubsetCluster_IDs))



affinity_cluster_df <- affinity_cluster_df %>% dplyr::left_join(count_subset_table,by ="SubsetCluster_IDs")


affinity_cluster_annotation <- affinity_cluster_df %>% dplyr::select(Tumor_Names,Pediatric_SubsetCluster_ID,SubsetCluster_IDs)

####### CHECK IF NCIT OR WHO 

affinity_cluster_annotation$NCIT_Tumor<-"No"
affinity_cluster_annotation$WHO_Tumor<-"No"

for (iter in 1:dim(affinity_cluster_annotation)[1]){
  if(affinity_cluster_annotation$Tumor_Names[iter] %in% NCIT_embedding_df$Disease){
    affinity_cluster_annotation$NCIT_Tumor[iter] <- "Yes"
  }else if(affinity_cluster_annotation$Tumor_Names[iter] %in% WHO_embedding_df$Disease){
    affinity_cluster_annotation$WHO_Tumor[iter] <- "Yes"
  }
}


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

ncit_parallel<-foreach(i=1:dim(tumor_distances_df)[1],.combine=rbind) %dopar% {
  print(iter)
  ncit_info <- nearest_ncit(combined_embedding_df[iter,1:1536],NCIT_embedding_df)
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


distance_ncit2 <- as.data.frame(matrix(nrow=dim(combined_embedding_df)[1], ncol=dim(NCIT_embedding_df)[1]))
rownames(distance_ncit2)<-rownames(combined_embedding_df)
colnames(distance_ncit2)<-(NCIT_embedding_df$Disease)

output_foreach<-foreach(i=1:dim(distance_ncit2)[1],.combine=rbind,.options.snow=opts) %dopar% {
apply(NCIT_embedding_df[,2:1537],1,CalculateEuclideanDistance,vect2=combined_embedding_df[i,])
  setTxtProgressBar(pb, i) 
}

 4

 outer_ncit_all<-invisible(foreach(i = 1:18014, .combine = rbind, .options.snow = opts, .verbose = T) %dopar% {
   s <- apply(NCIT_embedding_df[,2:1537],1,CalculateEuclideanDistance,vect2=combined_embedding_df[i,])
   setTxtProgressBar(pb, i)
   return(s)
 })


 outer_ncit_all_df<-outer_ncit_all_df[1:18014,]
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
 
 affinity_cluster_annotation2 <- affinity_cluster_annotation2 %>% dplyr::mutate(assigned_class = case_when(distance_ncit < WHO_distance ~ "NCIT",
                                                                                                    distance_ncit > WHO_distance ~ "WHO",
                                                                                                    TRUE ~ "Both"))
 affinity_cluster_annotation3 <- affinity_cluster_annotation2 %>% select(Tumor_Names,Pediatric_SubsetCluster_ID,SubsetCluster_IDs,NCIT_Tumor,WHO_Tumor,assigned_class)
 # Cluster voting
 
 affinity_cluster_annotation2$cluster_label <- NA
 subcluster_id_list<- unique(affinity_cluster_annotation2$SubsetCluster_IDs)
 
 for (iter in 1:length(subcluster_id_list)){
   
   cluster_id_affinity <- subcluster_id_list[iter]
   index_affinity <- which (affinity_cluster_annotation2$SubsetCluster_IDs==cluster_id_affinity)
   totaling_table <- as.data.frame(table(affinity_cluster_annotation2$assigned_class[index_affinity]))
   index_max<- which(totaling_table$Freq==max(totaling_table$Freq))
   assigned_cluster_labels <- unique(affinity_cluster_annotation2$assigned_class[index_max])
   affinity_cluster_annotation2$cluster_label[index_affinity]<- paste(assigned_cluster_labels,collapse = ";")
 }
 
 
 
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










# Read NCIT Codes
NCIT_df <- read.csv(paste(analyses_dir,"/embedding-analysis-dt/Neoplasm_Core.csv",sep=""))
NCIT_df$Preferred.Term <- tolower(NCIT_df$Preferred.Term )
NCIT_df$Synonyms <- tolower(NCIT_df$Synonyms )

df_tumor_leven$NCIT_Match <- NA
df_tumor_leven <- df_tumor_leven %>% dplyr::mutate(NCIT_Match=case_when(diseases %in% NCIT_df$Preferred.Term ~ "Yes", TRUE~"No" ))
df_tumor_leven <- df_tumor_leven %>% dplyr::mutate(NCIT_Code=case_when(diseases %in% NCIT_df$Preferred.Term ~ "Yes", TRUE~"No" ))
df_tumor_leven$NCIT_Code <- NA

ind_ncit <- which (df_tumor_leven$NCIT_Match=="Yes")

for (iter in ind_ncit){
  print(iter)
  disease_name <- df_tumor_leven$diseases[iter]
  ncit_index_preffered <- which(NCIT_df$Preferred.Term==disease_name)
  ncit_index_synonym <- which(NCIT_df$Synonyms==disease_name)
  
  
  if(length(ncit_index_preffered)>0){
    df_tumor_leven$NCIT_Code[iter]<-NCIT_df$Code[ncit_index_preffered]
  } else if(length(ncit_index_synonym)>0){
    df_tumor_leven$NCIT_Code[iter]<-NCIT_df$Code[ncit_index_synonym]
  }
  
}




# Find diseases that were not included in embeddings:
df <- read_excel((paste(analyses_dir,"/embedding-analysis-dt/diseases.xlsx",sep="")))
df<- df[,-1]

df2_embedding <- as.data.frame(rownames(embedding_df_agg))  
colnames(df2_embedding)<-"diseases"

ind_missed<- which(!(df$diseases %in% df2_embedding$diseases))
missed_CT_embeddings <- as.data.frame(df$diseases[ind_missed])
colnames(missed_CT_embeddings)<- "DISEASES"
write.csv(missed_CT_embeddings,"/Users/lahiria/Desktop/MTP_Paper/PMTL_paper/analyses/embedding-analysis-dt/missed_CT_embeddings.csv")


#results$rotation <- -1*results$rotation
#results$x <- -1*results$x
#disease_pca <- results$x
save(diseases_cluster,file=paste(analyses_dir,"/embedding-analysis-dt/diseases_cluster.Rdata",sep=""))
save(embedding_df_agg,file=paste(analyses_dir,"/embedding-analysis-dt/embedding_df_agg.Rdata",sep=""))

save.image(file='/Users/lahiria/Desktop/MTP_Paper/PMTL_paper/analyses/embedding-analysis-dt/workspace.RData')
load('/Users/lahiria/Desktop/MTP_Paper/PMTL_paper/analyses/embedding-analysis-dt/server_chop_download/workspace.RData')

write.csv(affinity_cluster_annotation,"/Users/lahiria/Desktop/MTP_Paper/PMTL_paper/analyses/embedding-analysis-dt/affinity_annotation.csv")

r <- mclapply(1:5, function(i) {rnorm(3)}, mc.cores = 5)

