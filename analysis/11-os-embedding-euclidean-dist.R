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
  library(rollama)
  library(doParallel)
  library(foreach)
  library(stringdist)
  library(rdist)
  library(readr)
})
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"data","input")
analysis_dir <- file.path(root_dir,"analysis")

source("~/Desktop/MTP_Paper/CT-Embedding-Paper/util/compute_embedding_distance.R")
# Read the embedding terms from the LTE-3 files
ct_terms <- read.csv(paste(data_dir,"/CT_Embeddings_V3.csv",sep = ""))
ncit_terms <- read.csv(paste(data_dir,"/NCIT_Embeddings_V3.csv",sep = ""))
who_terms <- read.csv(paste(data_dir,"/WHO_Terms_All_V3.csv",sep = ""))


ct_terms<- ct_terms %>% dplyr::select(Tumor_Names)
ncit_terms<- ncit_terms %>% dplyr::select(Tumor_Names)
who_terms<- who_terms %>% dplyr::select(Tumor_Names)

ct_terms2<-ct_terms %>% dplyr::pull(Tumor_Names)
who_terms2<-who_terms %>% dplyr::pull(Tumor_Names)
ncit_terms2<-ncit_terms %>% dplyr::pull(Tumor_Names)


tumor_terms<-c(ct_terms2,who_terms2,ncit_terms2)

#llama_embeddings<- embed_text(tumor_terms,model = "llama3",server = NULL,model_params = NULL,
#                             verbose = getOption("rollama_verbose", default = interactive()))

llama_embeddings<- read.csv(paste(data_dir,"/Embeddings/embeddings_llama.csv",sep="")) 
llama_embeddings<-llama_embeddings[,c(-1)]
biobert_embeddings <- read.csv(paste(data_dir,"/Embeddings/biobert_embedding.csv",sep=""))
medllama_embeddings <- read.csv(paste(data_dir,"/Embeddings/medllama-13b.csv",sep="")) #5121 dimension
medllama_embeddings_7b<-read.csv(paste(data_dir,"/Embeddings/medllama-7b.csv",sep=""))
pubmedbert_embeddings <-read.csv(paste(data_dir,"/Embeddings/pubmedbert-base-embeddings.csv",sep="")) #768 dimension
modernbert_embeddings<-read.csv(paste(data_dir,"/Embeddings/mordernbert_embeddings.csv",sep=""))
colnames(modernbert_embeddings)[1]<-"Tumor_Names"

llama_32_3b_embeddings<-read.csv(paste(data_dir,"/Embeddings/llama32_3B.csv",sep="")) 
llama_32_3b_embeddings<-llama_32_3b_embeddings[,c(-1)]
phi4_embeddings <- read.csv(paste(data_dir,"/Embeddings/phi4.csv",sep="")) 
phi4_embeddings<-phi4_embeddings[,c(-1)]



llama_embeddings<-llama_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
biobert_embeddings<-biobert_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
medllama_embeddings<-medllama_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
pubmedbert_embeddings<-pubmedbert_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
modernbert_embeddings<-modernbert_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")

medllama_embeddings_7b<-cbind(medllama_embeddings$Tumor_Names,medllama_embeddings_7b)
colnames(medllama_embeddings_7b)[1]<-"Tumor_Names"
medllama_embeddings_7b<-medllama_embeddings_7b %>% group_by(Tumor_Names) %>% summarise_all("mean")

llama_32_3b_embeddings<-llama_32_3b_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
phi4_embeddings<-phi4_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")




WHO_Terms_All <-readxl::read_xlsx(paste(data_dir,"/WHO_Tumors/result/WHO_Tumor_all_edition.xlsx",sep=""))
WHO_Terms_5th<-WHO_Terms_All%>%filter(edition_5th=="Yes")


distance_matrix_llama<- compute_embedding_distance(llama_all[1:nrow(llama_embeddings),2:ncol(llama_embeddings)],"euclidean")
rownames(distance_matrix_llama)<-llama_embeddings$Tumor_Names
colnames(distance_matrix_llama)<-llama_embeddings$Tumor_Names


distance_matrix_biobert<-compute_embedding_distance(biobert_embeddings[1:nrow(biobert_embeddings),2:ncol(biobert_embeddings)],"euclidean")
rownames(distance_matrix_biobert)<-biobert_embeddings$Tumor_Names
colnames(distance_matrix_biobert)<-biobert_embeddings$Tumor_Names

distance_matrix_medllama<-compute_embedding_distance(medllama_embeddings[1:nrow(medllama_embeddings),2:ncol(medllama_embeddings)],"euclidean")
rownames(distance_matrix_medllama)<-medllama_embeddings$Tumor_Names
colnames(distance_matrix_medllama)<-medllama_embeddings$Tumor_Names


distance_matrix_pubmedbert<-compute_embedding_distance(pubmedbert_embeddings[1:nrow(pubmedbert_embeddings),2:ncol(pubmedbert_embeddings)],"euclidean")
rownames(distance_matrix_pubmedbert)<-pubmedbert_embeddings$Tumor_Names
colnames(distance_matrix_pubmedbert)<-pubmedbert_embeddings$Tumor_Names


distance_matrix_modernbert<-compute_embedding_distance(modernbert_embeddings[1:nrow(modernbert_embeddings),2:ncol(modernbert_embeddings)],"euclidean")
rownames(distance_matrix_modernbert)<-modernbert_embeddings$Tumor_Names
colnames(distance_matrix_modernbert)<-modernbert_embeddings$Tumor_Names


distance_matrix_medllama_7b<- compute_embedding_distance(medllama_embeddings_7b[1:nrow(medllama_embeddings_7b),2:ncol(medllama_embeddings_7b)],"euclidean")
rownames(distance_matrix_medllama_7b)<-medllama_embeddings_7b$Tumor_Names
colnames(distance_matrix_medllama_7b)<-medllama_embeddings_7b$Tumor_Names

distance_matrix_llama32<- compute_embedding_distance(llama_32_3b_embeddings[1:nrow(llama_32_3b_embeddings),2:ncol(llama_32_3b_embeddings)],"euclidean")
rownames(distance_matrix_llama32)<-llama_32_3b_embeddings$Tumor_Names
colnames(distance_matrix_llama32)<-llama_32_3b_embeddings$Tumor_Names

distance_matrix_phi<- compute_embedding_distance(phi4_embeddings[1:nrow(phi4_embeddings),2:ncol(phi4_embeddings)],"euclidean")
rownames(distance_matrix_phi)<-phi4_embeddings$Tumor_Names
colnames(distance_matrix_phi)<-phi4_embeddings$Tumor_Names




distance_matrix_llama<-as.data.frame(distance_matrix_llama) # EXECUTED
distance_matrix_biobert<-as.data.frame(distance_matrix_biobert) # EXECUTED
distance_matrix_medllama<-as.data.frame(distance_matrix_medllama) # EXECUTED
distance_matrix_pubmedbert<-as.data.frame(distance_matrix_pubmedbert) # EXECUTED
distance_matrix_modernbert<-as.data.frame(distance_matrix_modernbert) # EXECUTED
distance_matrix_medllama_7b<-as.data.frame(distance_matrix_medllama_7b) # EXECUTED


# WHO Matrices All edition 
distance_matrix_llama_all <- distance_matrix_llama %>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_biobert_all <- distance_matrix_biobert %>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_medllama_all <- distance_matrix_medllama %>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_pubmedbert_all <- distance_matrix_pubmedbert %>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_modernbert_all <- distance_matrix_modernbert %>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_medllama_7b_all <- distance_matrix_medllama_7b %>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))

# WHO Matrices 5th edition 
distance_matrix_llama_5th <- distance_matrix_llama %>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))
distance_matrix_biobert_5th <- distance_matrix_biobert %>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))
distance_matrix_medllama_5th <- distance_matrix_medllama %>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))
distance_matrix_pubmedbert_5th <- distance_matrix_pubmedbert %>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))
distance_matrix_modernbert_5th <- distance_matrix_modernbert %>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))
distance_matrix_medllama_7b_5th <- distance_matrix_medllama_7b %>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))


###
llama_match_all<- nearest_match_embeddings(distance_matrix_llama_all,"llama") # Executed
biobert_match_all<- nearest_match_embeddings(distance_matrix_biobert_all,"biobert")
medllama_match_all<- nearest_match_embeddings(distance_matrix_medllama_all,"medllama")
pubmedbert_match_all<- nearest_match_embeddings(distance_matrix_pubmedbert_all,"pubmedbert")
modernbert_match_all<- nearest_match_embeddings(distance_matrix_modernbert_all,"modernbert")
medllama_7b_match_all<- nearest_match_embeddings(distance_matrix_medllama_7b_all,"medllama_7b")


embedding_nearest_all <- llama_match_all %>%
  full_join(biobert_match_all, by = "Tumor_Names") %>%
  full_join(medllama_match_all, by = "Tumor_Names") %>%
  full_join(pubmedbert_match_all, by = "Tumor_Names") %>%
  full_join(modernbert_match_all, by = "Tumor_Names") %>%
  full_join(medllama_7b_match_all, by = "Tumor_Names")

####
llama_match_5th<- nearest_match_embeddings(distance_matrix_llama_5th,"llama") # Executed
biobert_match_5th<- nearest_match_embeddings(distance_matrix_biobert_5th,"biobert")
medllama_match_5th<- nearest_match_embeddings(distance_matrix_medllama_5th,"medllama")
pubmedbert_match_5th<- nearest_match_embeddings(distance_matrix_pubmedbert_5th,"pubmedbert")
modernbert_match_5th<- nearest_match_embeddings(distance_matrix_modernbert_5th,"modernbert")
medllama_7b_match_5th<- nearest_match_embeddings(distance_matrix_medllama_7b_5th,"medllama_7b")


embedding_nearest_5th <- llama_match_5th %>%
  full_join(biobert_match_5th, by = "Tumor_Names") %>%
  full_join(medllama_match_5th, by = "Tumor_Names") %>%
  full_join(pubmedbert_match_5th, by = "Tumor_Names") %>%
  full_join(modernbert_match_5th, by = "Tumor_Names")%>%
  full_join(medllama_7b_match_5th, by = "Tumor_Names")


#### NCIT 
distance_matrix_llama_ncit <- distance_matrix_llama %>% dplyr::select(all_of(ncit_terms$Tumor_Names))
distance_matrix_biobert_ncit <- distance_matrix_biobert %>% dplyr::select(all_of(ncit_terms$Tumor_Names))
distance_matrix_medllama_ncit <- distance_matrix_medllama %>% dplyr::select(all_of(ncit_terms$Tumor_Names))
distance_matrix_pubmedbert_ncit <- distance_matrix_pubmedbert %>% dplyr::select(all_of(ncit_terms$Tumor_Names))
distance_matrix_modernbert_ncit <- distance_matrix_modernbert %>% dplyr::select(all_of(ncit_terms$Tumor_Names))
distance_matrix_medllama_7b_ncit <- distance_matrix_medllama_7b %>% dplyr::select(all_of(ncit_terms$Tumor_Names))

####
llama_match_ncit<- nearest_match_embeddings(distance_matrix_llama_ncit,"llama") # Executed
biobert_match_ncit<- nearest_match_embeddings(distance_matrix_biobert_ncit,"biobert")
medllama_match_ncit<- nearest_match_embeddings(distance_matrix_medllama_ncit,"medllama")
pubmedbert_match_ncit<- nearest_match_embeddings(distance_matrix_pubmedbert_ncit,"pubmedbert")
modernbert_match_ncit<- nearest_match_embeddings(distance_matrix_modernbert_ncit,"modernbert")
medllama_7b_match_ncit<- nearest_match_embeddings(distance_matrix_medllama_7b_ncit,"medllama_7b")


embedding_nearest_ncit <- llama_match_ncit %>%
  full_join(biobert_match_ncit, by = "Tumor_Names") %>%
  full_join(medllama_match_ncit, by = "Tumor_Names") %>%
  full_join(pubmedbert_match_ncit, by = "Tumor_Names") %>%
  full_join(modernbert_match_ncit, by = "Tumor_Names") %>%
  full_join(medllama_7b_match_ncit, by = "Tumor_Names")

# Reconcile Results

NCIT_Results_all <- read_csv("analysis/results/NCIT_Results_all.csv")
NCIT_Results_5th <- read_csv("analysis/results_5th/NCIT_Results_5thed.csv")
WHO_Results_all <- read_csv("analysis/results/WHO_Results_all.csv")
WHO_Results_5th <- read_csv("analysis/results_5th/WHO_Results_5thed.csv")

NCIT_Results_all<-NCIT_Results_all[,c(-1)]
NCIT_Results_5th<-NCIT_Results_5th[,c(-1)]
WHO_Results_all<-WHO_Results_all[,c(-1)]
WHO_Results_5th<-WHO_Results_5th[,c(-1)]

colnames(embedding_nearest_ncit) <- c("Tumor_Names","euclidean_dist_llama","euclidean_dist_biobert",
                                      "euclidean_dist_medllama","euclidean_dist_pubmedbert","euclidean_dist_modernbert","euclidean_dist_medllama_7b")
colnames(embedding_nearest_all) <- c("Tumor_Names","euclidean_dist_llama","euclidean_dist_biobert",
                                      "euclidean_dist_medllama","euclidean_dist_pubmedbert","euclidean_dist_modernbert","euclidean_dist_medllama_7b")
colnames(embedding_nearest_5th) <- c("Tumor_Names","euclidean_dist_llama","euclidean_dist_biobert",
                                      "euclidean_dist_medllama","euclidean_dist_pubmedbert","euclidean_dist_modernbert","euclidean_dist_medllama_7b")

NCIT_Results_all <- NCIT_Results_all %>% left_join(embedding_nearest_ncit , by="Tumor_Names")
NCIT_Results_5th <- NCIT_Results_5th %>% left_join(embedding_nearest_ncit , by="Tumor_Names")

WHO_Results_all <- WHO_Results_all %>% left_join(embedding_nearest_all , by="Tumor_Names")
WHO_Results_5th <- WHO_Results_5th %>% left_join(embedding_nearest_5th , by="Tumor_Names")

write.csv(NCIT_Results_all,paste(analysis_dir,"/results/NCIT_Results_all_ose.csv",sep=""))
write.csv(WHO_Results_all,paste(analysis_dir,"/results/WHO_Results_all_ose.csv",sep=""))

write.csv(NCIT_Results_5th,paste(analysis_dir,"/results_5th/NCIT_Results_5thed_ose.csv",sep=""))
write.csv(WHO_Results_5th,paste(analysis_dir,"/results_5th/WHO_Results_5thed_ose.csv",sep=""))



#save.image("12-prior.RData")


