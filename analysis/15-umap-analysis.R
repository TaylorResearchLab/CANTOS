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
  library(rollama)
  library(doParallel)
  library(foreach)
  library(stringdist)
  library(rdist)
  library(readr)
  library(umap)
  
})
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"data","input")
analyses_dir <- file.path(root_dir,"analyses")

# Load embeddings ADA2
CT_ada2<- read.csv(paste(data_dir,"/CT_Embeddings_ADA2.csv",sep=""))
CT_ada2<-CT_ada2[,c(-1)]

WHO_ada2<- read.csv(paste(data_dir,"/WHO_Aggregate_ADA2.csv",sep=""))
WHO_ada2<-WHO_ada2[,c(-1)]

NCIT_ada2 <-read.csv(paste(data_dir,"/dt_input_file_6_dec/NCIT_Neoplasm_Core_terms_text-embedding-ada-002_embeddings.csv",sep=""))
NCIT_ada2<-NCIT_ada2[c(-1),] # Remove the header (column name) embedding
colnames(NCIT_ada2)[1]<-"Tumor_Names"

colnames(WHO_ada2)<-colnames(NCIT_ada2)
colnames(CT_ada2)<-colnames(NCIT_ada2)

ADA_002_Embeddings<-rbind(CT_ada2,WHO_ada2,NCIT_ada2)
ADA_002_Embeddings$Tumor_Names<-tolower(ADA_002_Embeddings$Tumor_Names)
ADA_002_Embeddings<-ADA_002_Embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")

# Load embeddings V3
CT_V3<- read.csv(paste(data_dir,"/CT_Embeddings_V3.csv",sep=""))
WHO_V3<- read.csv(paste(data_dir,"/WHO_Terms_All_V3.csv",sep=""))
NCIT_V3<-read.csv(paste(data_dir,"/NCIT_Embeddings_V3.csv",sep=""))
LTE_V3_Embeddings<-rbind(CT_V3,WHO_V3,NCIT_V3)
LTE_V3_Embeddings<-LTE_V3_Embeddings[,c(-1)]
LTE_V3_Embeddings<-LTE_V3_Embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")


# Load medllama 
medllama_Embeddings <- read.csv(paste(data_dir,"/Embeddings/medllama-13b.csv",sep=""))
# Load llama 
llama_Embeddings <- read.csv(paste(data_dir,"/Embeddings/embeddings_llama.csv",sep=""))
llama_Embeddings<-llama_Embeddings[,c(-1)]
# Load biobert
biobert_Embeddings <- read.csv(paste(data_dir,"/Embeddings/biobert_embedding.csv",sep=""))
# Load pubmedbert
pubmedbert_Embeddings <- read.csv(paste(data_dir,"/Embeddings/pubmedbert-base-embeddings.csv",sep=""))

medllama_Embeddings<-medllama_Embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
llama_Embeddings<-llama_Embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
biobert_Embeddings<-biobert_Embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
pubmedbert_Embeddings<-pubmedbert_Embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")

save.image("15-umap-analysis.RData")

WHO_Terms_All <-readxl::read_xlsx(paste(data_dir,"/WHO_Tumors/result/WHO_Tumor_all_edition.xlsx",sep=""))
WHO_Terms_5th<-WHO_Terms_All%>%filter(edition_5th=="Yes")


# Scale the embeddings
scale_data <- function(embeddings_data){
  Tumors<- embeddings_data$Tumor_Names
  embeddings_data<-embeddings_data[,c(-1)]
  embeddings_data<-scale(embeddings_data)
  embeddings_data<-cbind(Tumor_Names =Tumors , as.data.frame(embeddings_data))
  #apply(embeddings_data, 2, sd)
  #colMeans(embeddings_data)
  return(embeddings_data)
}



ADA_002_Embeddings_Scaled<-scale_data(ADA_002_Embeddings)
LTE_V3_Embeddings_Scaled<-scale_data(LTE_V3_Embeddings)
medllama_Embeddings_Scaled<-scale_data(medllama_Embeddings)
llama_Embeddings_Scaled<-scale_data(llama_Embeddings)
biobert_Embeddings_Scaled<-scale_data(biobert_Embeddings)
pubmedbert_Embeddings_Scaled<-scale_data(pubmedbert_Embeddings)

set.seed(44)

UMAP_data <- function(Scaled_Data){
  Tumors<-Scaled_Data$Tumor_Names
  Scaled_Data<- Scaled_Data[,c(-1)]
  # Perform UMAP to reduce to 3 dimensions
  umap_config <- umap.defaults
  umap_config$n_components <- 3  # Focus on 3 dimensions
  umap_result <- umap(Scaled_Data, config = umap_config)
  # Extract the UMAP reduced data
  reduced_data <- as.data.frame(umap_result$layout)
  colnames(reduced_data) <- c("UMAP1", "UMAP2", "UMAP3")  # Rename columns
  reduced_data$Tumor_Names <- Tumors  # Add tumor names back to the reduced data
  reduced_data<-reduced_data[,c(4,1,2,3)]
  return(reduced_data)
}

ADA_002_Umap<-UMAP_data(ADA_002_Embeddings_Scaled)
LTE_V3_Umap<-UMAP_data(LTE_V3_Embeddings_Scaled)
medllama_Umap<-UMAP_data(medllama_Embeddings_Scaled)
llama_Umap<-UMAP_data(llama_Embeddings_Scaled)
biobert_Umap<-UMAP_data(biobert_Embeddings_Scaled)
pubmedbert_Umap<-UMAP_data(pubmedbert_Embeddings_Scaled)


perform_kmeans <- function(embedding_data){
  set.seed(44)
  rownames(embedding_data)<-embedding_data$Tumor_Names
  embedding_data<-embedding_data[,c(-1)]
  km.res <- eclust(embedding_data, "kmeans", k = 6200,nstart = 25, graph = FALSE)
  kmeans_clust_result <- as.data.frame(km.res$cluster)
  kmeans_clust_result$Tumor_Names<-rownames(kmeans_clust_result)
  colnames(kmeans_clust_result)[1]<-"cluster"
  kmeans_clust_result <- kmeans_clust_result %>% dplyr::select(Tumor_Names,cluster)
  
}

ADA002_KMeans<-perform_kmeans(ADA_002_Umap)
LTE_V3_KMeans<-perform_kmeans(LTE_V3_Umap)
medllama_KMeans<-perform_kmeans(medllama_Umap)
llama_KMeans<-perform_kmeans(llama_Umap)
biobert_KMeans<-perform_kmeans(biobert_Umap)
pubmedbert_KMeans<-perform_kmeans(pubmedbert_Umap)



# Entropy analysis
calculate_entropy <- function(data) {
  class_proportions <- table(data$class_label) / nrow(data)
  entropy <- -sum(class_proportions * log2(class_proportions))
  return(entropy)
}

ADA_002_Umap$class_label <- ADA002_KMeans$cluster
cluster_entropy_ada <- ADA_002_Umap %>% group_by(cluster) %>%summarise(entropy = calculate_entropy(cur_data()))

