#Load the manually annotated disease file with pediatric and adult cancer annotation.
# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(doParallel)
  library(foreach)
  library(stringdist)
  #library(DescTools)
  library(stringr)
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
results_dir <- file.path(analysis_dir,"results")

#Load Functions
source(paste(util_dir,"/string_dissimilarity.R",sep = ""))
source(paste(util_dir,"/distance_clusters.R",sep=""))
source(paste(util_dir,"/string_normalizing.R",sep=""))

# Read the annotated file
ct_disease_annot_adult_ped_df<-read.csv(paste(input_dir,"/tumor_annotated_adult_ped.csv",sep=""))
ct_tumor_df<-ct_disease_annot_adult_ped_df%>%filter(validated_cancer_tumor=="Yes")
ct_tumor_df<-ct_tumor_df[,c(-1)]
# Add NCIT and WHO Tumors 
# Read NCIT Terms and WHO Terms with embedding and join them to the rest of the embedding list 

# NEW addition : 
WHO_Terms_All <-readxl::read_xlsx(paste(data_dir,"/WHO_Tumors/result/WHO_Tumor_all_edition.xlsx",sep=""))


NCIT_Terms <-read.csv(paste(data_dir,"/dt_input_file_6_dec/NCIT_Neoplasm_Core_terms_text-embedding-ada-002_embeddings.csv",sep=""))[,1]

NCIT_Terms<-tolower(NCIT_Terms[c(-1)])

# Levenstein distance between tumors 
df_tumor_combined<-as.data.frame(unique(c(ct_tumor_df$diseases,NCIT_Terms,WHO_Terms_All$Tumor_Names)))
colnames(df_tumor_combined)[1]<-"Tumor"


df_tumor_names<-unique(df_tumor_combined$Tumor)

dissimilarity_matrix_lv <- as.data.frame(matrix(nrow=length(df_tumor_names),ncol=length(df_tumor_names)))
rownames(dissimilarity_matrix_lv)<-df_tumor_names
colnames(dissimilarity_matrix_lv)<-df_tumor_names


cl <- makeCluster(25, outfile="")
registerDoParallel(cl)


dissimilarity_matrix_lv<-foreach(iter=1:length(df_tumor_names),.combine=rbind) %dopar% {
  print(iter)
  disease_name <- colnames(dissimilarity_matrix_lv)[iter]
  distances<-unlist(lapply(df_tumor_names,string_dissimilarity,S2=disease_name,meth="lv"))
}
rownames(dissimilarity_matrix_lv) <- df_tumor_names
colnames(dissimilarity_matrix_lv) <- df_tumor_names

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
dissimilarity_matrix_lv<-dissimilarity_matrix_lv/normalizing_matrix_lv

save(dissimilarity_matrix_lv,file=paste(intermediate_dir,"/dissimilarity_matrix_lv.RData",sep=""))

stopCluster(cl)



# Jarro Winkler Distance
cl <- makeCluster(25, outfile="")
registerDoParallel(cl)
dissimilarity_matrix_jw <- as.data.frame(matrix(nrow=length(df_tumor_names),ncol=length(df_tumor_names)))
rownames(dissimilarity_matrix_jw) <- df_tumor_names
colnames(dissimilarity_matrix_jw) <- df_tumor_names

dissimilarity_matrix_jw<-foreach(iter=1:length(df_tumor_names),.combine=rbind) %dopar% {
  print(iter)
  disease_name <- colnames(dissimilarity_matrix_jw)[iter]
  distances<-unlist(lapply(df_tumor_names,string_dissimilarity,S2=disease_name,meth="jw"))
  
}
rownames(dissimilarity_matrix_jw) <- df_tumor_names
colnames(dissimilarity_matrix_jw) <- df_tumor_names
stopCluster(cl)
save(dissimilarity_matrix_jw,file=paste(intermediate_dir,"/dissimilarity_matrix_jw.RData",sep=""))

# Cosine Distance
cl <- makeCluster(25, outfile="")
registerDoParallel(cl)
dissimilarity_matrix_cosine <- as.data.frame(matrix(nrow=length(df_tumor_names),ncol=length(df_tumor_names)))
rownames(dissimilarity_matrix_cosine) <- df_tumor_names
colnames(dissimilarity_matrix_cosine) <- df_tumor_names

dissimilarity_matrix_cosine<-foreach(iter=1:length(df_tumor_names),.combine=rbind) %dopar% {
  print(iter)
  disease_name <- colnames(dissimilarity_matrix_cosine)[iter]
  distances<-unlist(lapply(df_tumor_names,string_dissimilarity,S2=disease_name,meth="cosine"))
  
}
rownames(dissimilarity_matrix_cosine) <- df_tumor_names
colnames(dissimilarity_matrix_cosine) <- df_tumor_names
stopCluster(cl)

save(dissimilarity_matrix_cosine,file=paste(intermediate_dir,"/dissimilarity_matrix_cosine.RData",sep=""))


# Find clusters of 
cl <- makeCluster(25, outfile="")
registerDoParallel(cl)
cluster_results_lv<-foreach(iter=1:length(df_tumor_names),.combine=rbind) %do% {
  print(iter)
  distance_clusters(dissimilarity_matrix_lv[iter,],cutoff = 0.0005,df_tumor_names)
}

cluster_results_lv<-as.data.frame(cluster_results_lv)
rownames(cluster_results_lv)<-NULL
cluster_results_lv <-cbind(cluster_results_lv,df_tumor_names)


colnames(cluster_results_lv) <- c("cluster_lv","Tumors")
cluster_results_lv <- cluster_results_lv %>% dplyr::select(Tumors,cluster_lv)

### JW

cluster_results_jw<-foreach(iter=1:length(df_tumor_names),.combine=rbind) %do% {
  print(iter)
  distance_clusters(dissimilarity_matrix_jw[iter,],cutoff = 0.0005,df_tumor_names)
}

cluster_results_jw<-as.data.frame(cluster_results_jw)
rownames(cluster_results_jw)<-NULL
cluster_results_jw <-cbind(cluster_results_jw,df_tumor_names)


colnames(cluster_results_jw) <- c("cluster_jw","Tumors")
cluster_results_jw <- cluster_results_jw %>%  dplyr::select(Tumors,cluster_jw)

# Cosine

cluster_results_cosine<-foreach(iter=1:length(df_tumor_names),.combine=rbind) %do% {
  print(iter)
  distance_clusters(dissimilarity_matrix_cosine[iter,],cutoff = 0.0005,df_tumor_names)
}

cluster_results_cosine<-as.data.frame(cluster_results_cosine)
rownames(cluster_results_cosine)<-NULL
cluster_results_cosine <-cbind(cluster_results_cosine,df_tumor_names)


colnames(cluster_results_cosine) <- c("cluster_cosine","Tumors")
cluster_results_cosine <- cluster_results_cosine %>%  dplyr::select(Tumors,cluster_cosine)


# Number of clusters
cluster_results_lv<-cluster_results_lv %>% mutate(cluster_members = str_count(cluster_lv, ";")+1)
cluster_results_jw<-cluster_results_jw %>% mutate(cluster_members = str_count(cluster_jw, ";")+1)
cluster_results_cosine<-cluster_results_cosine %>% mutate(cluster_members = str_count(cluster_cosine, ";")+1)  
  
# filter by ct tumors
cluster_results_lv_ct_filtered<- cluster_results_lv %>% filter(Tumors %in% ct_tumor_df$diseases)
cluster_results_jw_ct_filtered<- cluster_results_jw %>% filter(Tumors %in% ct_tumor_df$diseases)
cluster_results_cosine_ct_filtered<- cluster_results_cosine %>% filter(Tumors %in% ct_tumor_df$diseases)

# add nct_id
ct_tumor_df_disease_short <- ct_tumor_df %>% dplyr::select(nct_id,diseases)
cluster_results_lv_ct_filtered <- cluster_results_lv_ct_filtered %>% left_join(ct_tumor_df_disease_short,by=c("Tumors"="diseases"))
cluster_results_lv_ct_filtered <- cluster_results_lv_ct_filtered[,c(4,1:3)]


cluster_results_jw_ct_filtered <- cluster_results_jw_ct_filtered %>% left_join(ct_tumor_df_disease_short,by=c("Tumors"="diseases"))
cluster_results_jw_ct_filtered <- cluster_results_jw_ct_filtered[,c(4,1:3)]


cluster_results_cosine_ct_filtered <- cluster_results_cosine_ct_filtered %>% left_join(ct_tumor_df_disease_short,by=c("Tumors"="diseases"))
cluster_results_cosine_ct_filtered <- cluster_results_cosine_ct_filtered[,c(4,1:3)]

# organize by alphabetical order
cluster_results_lv_ct_filtered <- cluster_results_lv_ct_filtered[order(cluster_results_lv_ct_filtered$Tumors),]
cluster_results_jw_ct_filtered <- cluster_results_jw_ct_filtered[order(cluster_results_jw_ct_filtered$Tumors),]
cluster_results_cosine_ct_filtered <- cluster_results_cosine_ct_filtered[order(cluster_results_cosine_ct_filtered$Tumors),]



# Null the row names
rownames(cluster_results_lv_ct_filtered)<-NULL
rownames(cluster_results_jw_ct_filtered)<-NULL
rownames(cluster_results_cosine_ct_filtered)<-NULL


# Show results with following tumors 
benchmark_tumors <- c("b cell lymphoma", "neuroblastoma", "triple negative breast cancer",
                      "unresectable lung carcinoma", "liposarcoma","cancer of the liver",
                      "smoldering myeloma")

display_table_benchmark <- cbind(cluster_results_lv_ct_filtered %>% filter(Tumors %in% benchmark_tumors) %>% dplyr::select(nct_id,Tumors,cluster_lv),
                                 cluster_results_jw_ct_filtered %>% filter(Tumors %in%  benchmark_tumors) %>% dplyr::select(cluster_jw),
                                 cluster_results_cosine_ct_filtered %>% filter(Tumors %in% benchmark_tumors) %>% dplyr::select(cluster_cosine)
                                 )
                                 

display_table_benchmark <- display_table_benchmark %>% dplyr::select(nct_id,Tumors,cluster_lv,cluster_jw,cluster_cosine)

# To organize from good to worse
benchmark_tumors<- as.data.frame(benchmark_tumors)
colnames(benchmark_tumors)<-"Tumors"
display_table_benchmark <- benchmark_tumors %>% dplyr::left_join(display_table_benchmark,by="Tumors")

display_table_benchmark <- display_table_benchmark %>% dplyr::select(nct_id,Tumors,cluster_lv,cluster_jw,cluster_cosine)


# Stop the cluster
stopCluster(cl)


# Write Results of Clusters
write.csv(cluster_results_lv_ct_filtered,paste(results_dir,"/cluster_lv.csv",sep=""))
write.csv(cluster_results_jw_ct_filtered,paste(results_dir,"/cluster_jw.csv",sep=""))
write.csv(cluster_results_cosine_ct_filtered,paste(results_dir,"/cluster_cosine.csv",sep=""))
write.csv(display_table_benchmark,paste(results_dir,"/edit_distance_bench_mark.csv",sep=""))

save.image(file = "script2.RData")
