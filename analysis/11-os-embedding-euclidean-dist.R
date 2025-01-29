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

source(paste(util_dir,"/compute_embedding_distance.R",sep=""))
source(paste(util_dir,"/nearest_match_embeddings.R",sep=""))

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

llama_33_70b_embeddings<-read.csv(paste(data_dir,"/Embeddings/all_embedding_llama_33_70b.csv",sep="")) 
colnames(llama_33_70b_embeddings)[1]<-"Tumor_Names"

all_MiniLM_L6_v2_embeddings<-read.csv(paste(data_dir,"/Embeddings/all_MiniLM_L6_v2.csv",sep="")) 
all_mpnet_base_v2_embeddings<-read.csv(paste(data_dir,"/Embeddings/all_mpnet_base_v2.csv",sep=""))

e5_large_embeddings<-read.csv(paste(data_dir,"/Embeddings/e5-large.csv",sep="")) 
gtr_t5_large_embeddings<-read.csv(paste(data_dir,"/Embeddings/gtr-t5-large.csv",sep="")) 

roberta_embeddings<-read.csv(paste(data_dir,"/Embeddings/roberta.csv",sep="")) 
MiniLM_L12_v2_embeddings<-read.csv(paste(data_dir,"/Embeddings/all-MiniLM-L12-v2.csv",sep=""))


labse_embeddings<-read.csv(paste(data_dir,"/Embeddings/tumor_embeddings_labse.csv",sep="")) 
sciBERT_embeddings<-read.csv(paste(data_dir,"/Embeddings/tumor_embeddings_scibert.csv",sep="")) 
sapBERT_embeddings<-read.csv(paste(data_dir,"/Embeddings/tumor_embeddings_sapbert.csv",sep="")) 
cohere_english_v2_embeddings<-read.csv(paste(data_dir,"/Embeddings/cohere_embeddings_embed_english_v2.csv",sep="")) 

deepseek_8b_embeddings<-read.csv(paste(data_dir,"/Embeddings/deepseek_8b.csv",sep = ""))


BioGPT_embeddings<-read.csv(paste(data_dir,"/Embeddings/output_tumor_embeddings_biogpt.csv",sep="")) 
clincalBERT_embeddings<-read.csv(paste(data_dir,"/Embeddings/output_tumor_embeddings_clinicalBERT.csv",sep="")) 

#############
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
llama_33_70b_embeddings<-llama_33_70b_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
all_MiniLM_L6_v2_embeddings<-all_MiniLM_L6_v2_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
all_mpnet_base_v2_embeddings<-all_mpnet_base_v2_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")

e5_large_embeddings<-e5_large_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
gtr_t5_large_embeddings<-gtr_t5_large_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")

roberta_embeddings<-roberta_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
MiniLM_L12_v2_embeddings<-MiniLM_L12_v2_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")



labse_embeddings<-labse_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
sciBERT_embeddings<-sciBERT_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
sapBERT_embeddings<-sapBERT_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
cohere_english_v2_embeddings<-cohere_english_v2_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")

deepseek_8b_embeddings<-deepseek_8b_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")

BioGPT_embeddings<-BioGPT_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
clincalBERT_embeddings<-clincalBERT_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")



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

distance_matrix_llama33_70b<- compute_embedding_distance(llama_33_70b_embeddings[1:nrow(llama_33_70b_embeddings),2:ncol(llama_33_70b_embeddings)],"euclidean")
rownames(distance_matrix_llama33_70b)<-llama_33_70b_embeddings$Tumor_Names
colnames(distance_matrix_llama33_70b)<-llama_33_70b_embeddings$Tumor_Names

distance_matrix_all_MiniLM_L6_v2<- compute_embedding_distance(all_MiniLM_L6_v2_embeddings[1:nrow(all_MiniLM_L6_v2_embeddings),2:ncol(all_MiniLM_L6_v2_embeddings)],"euclidean")
rownames(distance_matrix_all_MiniLM_L6_v2)<-all_MiniLM_L6_v2_embeddings$Tumor_Names
colnames(distance_matrix_all_MiniLM_L6_v2)<-all_MiniLM_L6_v2_embeddings$Tumor_Names

distance_matrix_all_mpnet_base_v2<- compute_embedding_distance(all_mpnet_base_v2_embeddings[1:nrow(all_mpnet_base_v2_embeddings),2:ncol(all_mpnet_base_v2_embeddings)],"euclidean")
rownames(distance_matrix_all_mpnet_base_v2)<-all_mpnet_base_v2_embeddings$Tumor_Names
colnames(distance_matrix_all_mpnet_base_v2)<-all_mpnet_base_v2_embeddings$Tumor_Names

distance_matrix_e5_large<- compute_embedding_distance(e5_large_embeddings[1:nrow(e5_large_embeddings),2:ncol(e5_large_embeddings)],"euclidean")
rownames(distance_matrix_e5_large)<-e5_large_embeddings$Tumor_Names
colnames(distance_matrix_e5_large)<-e5_large_embeddings$Tumor_Names

distance_matrix_gtr_t5_large<- compute_embedding_distance(gtr_t5_large_embeddings[1:nrow(gtr_t5_large_embeddings),2:ncol(gtr_t5_large_embeddings)],"euclidean")
rownames(distance_matrix_gtr_t5_large)<-gtr_t5_large_embeddings$Tumor_Names
colnames(distance_matrix_gtr_t5_large)<-gtr_t5_large_embeddings$Tumor_Names

distance_matrix_roberta<- compute_embedding_distance(roberta_embeddings[1:nrow(roberta_embeddings),2:ncol(roberta_embeddings)],"euclidean")
rownames(distance_matrix_roberta)<-roberta_embeddings$Tumor_Names
colnames(distance_matrix_roberta)<-roberta_embeddings$Tumor_Names

distance_matrix_MiniLM_L12_v2<- compute_embedding_distance(MiniLM_L12_v2_embeddings[1:nrow(MiniLM_L12_v2_embeddings),2:ncol(MiniLM_L12_v2_embeddings)],"euclidean")
rownames(distance_matrix_MiniLM_L12_v2)<-MiniLM_L12_v2_embeddings$Tumor_Names
colnames(distance_matrix_MiniLM_L12_v2)<-MiniLM_L12_v2_embeddings$Tumor_Names



distance_matrix_labse<- compute_embedding_distance(labse_embeddings[1:nrow(labse_embeddings),2:ncol(labse_embeddings)],"euclidean")
rownames(distance_matrix_labse)<-labse_embeddings$Tumor_Names
colnames(distance_matrix_labse)<-labse_embeddings$Tumor_Names

distance_matrix_sciBERT<- compute_embedding_distance(sciBERT_embeddings[1:nrow(sciBERT_embeddings),2:ncol(sciBERT_embeddings)],"euclidean")
rownames(distance_matrix_sciBERT)<-sciBERT_embeddings$Tumor_Names
colnames(distance_matrix_sciBERT)<-sciBERT_embeddings$Tumor_Names

distance_matrix_sapBERT<- compute_embedding_distance(sapBERT_embeddings[1:nrow(sapBERT_embeddings),2:ncol(sapBERT_embeddings)],"euclidean")
rownames(distance_matrix_sapBERT)<-sapBERT_embeddings$Tumor_Names
colnames(distance_matrix_sapBERT)<-sapBERT_embeddings$Tumor_Names

distance_matrix_cohere<- compute_embedding_distance(cohere_english_v2_embeddings[1:nrow(cohere_english_v2_embeddings),2:ncol(cohere_english_v2_embeddings)],"euclidean")
rownames(distance_matrix_cohere)<-cohere_english_v2_embeddings$Tumor_Names
colnames(distance_matrix_cohere)<-cohere_english_v2_embeddings$Tumor_Names


distance_matrix_deepseek<- compute_embedding_distance(deepseek_8b_embeddings[1:nrow(deepseek_8b_embeddings),2:ncol(deepseek_8b_embeddings)],"euclidean")
rownames(distance_matrix_deepseek)<-deepseek_8b_embeddings$Tumor_Names
colnames(distance_matrix_deepseek)<-deepseek_8b_embeddings$Tumor_Names

distance_matrix_BioGPT<- compute_embedding_distance(BioGPT_embeddings[1:nrow(BioGPT_embeddings),2:ncol(BioGPT_embeddings)],"euclidean")
rownames(distance_matrix_deepseek)<-BioGPT_embeddings$Tumor_Names
colnames(distance_matrix_deepseek)<-BioGPT_embeddings$Tumor_Names

distance_matrix_clinicalBERT<- compute_embedding_distance(clincalBERT_embeddings[1:nrow(clincalBERT_embeddings),2:ncol(clincalBERT_embeddings)],"euclidean")
rownames(distance_matrix_deepseek)<-clincalBERT_embeddings$Tumor_Names
colnames(distance_matrix_deepseek)<-clincalBERT_embeddings$Tumor_Names




distance_matrix_llama<-as.data.frame(distance_matrix_llama) # EXECUTED
distance_matrix_biobert<-as.data.frame(distance_matrix_biobert) # EXECUTED
distance_matrix_medllama<-as.data.frame(distance_matrix_medllama) # EXECUTED
distance_matrix_pubmedbert<-as.data.frame(distance_matrix_pubmedbert) # EXECUTED
distance_matrix_modernbert<-as.data.frame(distance_matrix_modernbert) # EXECUTED
distance_matrix_medllama_7b<-as.data.frame(distance_matrix_medllama_7b) # EXECUTED
distance_matrix_llama32<-as.data.frame(distance_matrix_llama32) # EXECUTED
distance_matrix_phi<-as.data.frame(distance_matrix_phi) # EXECUTED
distance_matrix_llama33_70b<-as.data.frame(distance_matrix_llama33_70b) # EXECUTED
distance_matrix_all_MiniLM_L6_v2<-as.data.frame(distance_matrix_all_MiniLM_L6_v2)# EXECUTED
distance_matrix_all_mpnet_base_v2<-as.data.frame(distance_matrix_all_mpnet_base_v2)# EXECUTED
distance_matrix_e5_large<-as.data.frame(distance_matrix_e5_large)
distance_matrix_gtr_t5_large<-as.data.frame(distance_matrix_gtr_t5_large)
distance_matrix_roberta<-as.data.frame(distance_matrix_roberta)
distance_matrix_MiniLM_L12_v2<-as.data.frame(distance_matrix_MiniLM_L12_v2)
distance_matrix_labse<-as.data.frame(distance_matrix_labse)
distance_matrix_sciBERT<-as.data.frame(distance_matrix_sciBERT)
distance_matrix_sapBERT<-as.data.frame(distance_matrix_sapBERT)
distance_matrix_cohere<-as.data.frame(distance_matrix_cohere)
distance_matrix_deepseek<-as.data.frame(distance_matrix_deepseek)

# WHO Matrices All edition 
distance_matrix_llama_all <- distance_matrix_llama %>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_biobert_all <- distance_matrix_biobert %>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_medllama_all <- distance_matrix_medllama %>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_pubmedbert_all <- distance_matrix_pubmedbert %>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_modernbert_all <- distance_matrix_modernbert %>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_medllama_7b_all <- distance_matrix_medllama_7b %>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_llama32_all <- distance_matrix_llama32 %>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_phi_all <- distance_matrix_phi %>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_llama33_70b_all <- distance_matrix_llama33_70b %>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_MiniLM_L6_v2_all <- distance_matrix_all_MiniLM_L6_v2 %>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_mpnet_base_all <- distance_matrix_all_mpnet_base_v2 %>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_e5_large_all<- distance_matrix_e5_large %>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_gtr_t5_large_all<-distance_matrix_gtr_t5_large %>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_roberta_all<-distance_matrix_roberta%>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_MiniLM_L12_v2_all<-distance_matrix_MiniLM_L12_v2%>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))

distance_matrix_labse_all<-distance_matrix_labse%>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_sciBERT_all<-distance_matrix_sciBERT%>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_sapBERT_all<-distance_matrix_sapBERT%>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_cohere_all<-distance_matrix_cohere%>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))

distance_matrix_deepseek_all<-distance_matrix_deepseek%>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))

# WHO Matrices 5th edition 
distance_matrix_llama_5th <- distance_matrix_llama %>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))
distance_matrix_biobert_5th <- distance_matrix_biobert %>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))
distance_matrix_medllama_5th <- distance_matrix_medllama %>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))
distance_matrix_pubmedbert_5th <- distance_matrix_pubmedbert %>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))
distance_matrix_modernbert_5th <- distance_matrix_modernbert %>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))
distance_matrix_medllama_7b_5th <- distance_matrix_medllama_7b %>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))
distance_matrix_llama32_5th <- distance_matrix_llama32 %>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))
distance_matrix_phi_5th <- distance_matrix_phi %>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))
distance_matrix_llama33_70b_5th <- distance_matrix_llama33_70b %>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))

distance_matrix_MiniLM_L6_v2_5th <- distance_matrix_all_MiniLM_L6_v2 %>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))
distance_matrix_mpnet_base_5th <- distance_matrix_all_mpnet_base_v2 %>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))

distance_matrix_e5_large_5th<- distance_matrix_e5_large %>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))
distance_matrix_gtr_t5_large_5th<<-distance_matrix_gtr_t5_large %>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))

distance_matrix_roberta_5th<-distance_matrix_roberta%>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))
distance_matrix_MiniLM_L12_v2_5th<-distance_matrix_MiniLM_L12_v2%>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))

distance_matrix_labse_5th<-distance_matrix_labse%>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))
distance_matrix_sciBERT_5th<-distance_matrix_sciBERT%>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))
distance_matrix_sapBERT_5th<-distance_matrix_sapBERT%>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))
distance_matrix_cohere_5th<-distance_matrix_cohere%>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))

distance_matrix_deepseek_5th<-distance_matrix_deepseek%>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))

###
llama_match_all<- nearest_match_embeddings(distance_matrix_llama_all,"llama") # Executed
biobert_match_all<- nearest_match_embeddings(distance_matrix_biobert_all,"biobert")
medllama_match_all<- nearest_match_embeddings(distance_matrix_medllama_all,"medllama")
pubmedbert_match_all<- nearest_match_embeddings(distance_matrix_pubmedbert_all,"pubmedbert")
modernbert_match_all<- nearest_match_embeddings(distance_matrix_modernbert_all,"modernbert")
medllama_7b_match_all<- nearest_match_embeddings(distance_matrix_medllama_7b_all,"medllama_7b")
llama_32_3b_match_all<- nearest_match_embeddings(distance_matrix_llama32_all,"llama_32_3b")
phi4_match_all<- nearest_match_embeddings(distance_matrix_phi_all,"phi4")
llama33_70b_all<-nearest_match_embeddings(distance_matrix_llama33_70b_all,"llama_33_70b")

MiniLM_L6_v2_all<-nearest_match_embeddings(distance_matrix_MiniLM_L6_v2_all,"MiniLM_L6_v2")
mpnet_base_all<-nearest_match_embeddings(distance_matrix_mpnet_base_all,"mpnet_base")

e5_large_all<-nearest_match_embeddings(distance_matrix_e5_large_all,"e5_large")
gtr_t5_large_all<-nearest_match_embeddings(distance_matrix_gtr_t5_large_all,"gtr_t5_large")

roberta_all<-nearest_match_embeddings(distance_matrix_roberta_all,"roberta")
MiniLM_L12_v2_all<-nearest_match_embeddings(distance_matrix_MiniLM_L12_v2_all,"MiniLM_L12_v2")

labse_all<-nearest_match_embeddings(distance_matrix_labse_all,"labase")
sciBERT_all<-nearest_match_embeddings(distance_matrix_sciBERT_all,"sciBERT")
sapBERT_all<-nearest_match_embeddings(distance_matrix_sapBERT_all,"sapBERT")
cohere_all<-nearest_match_embeddings(distance_matrix_cohere_all,"cohere")

deepseek_all<-nearest_match_embeddings(distance_matrix_deepseek_all,"deepseek")


embedding_nearest_all <- llama_match_all %>%
  full_join(biobert_match_all, by = "Tumor_Names") %>%
  full_join(medllama_match_all, by = "Tumor_Names") %>%
  full_join(pubmedbert_match_all, by = "Tumor_Names") %>%
  full_join(modernbert_match_all, by = "Tumor_Names") %>%
  full_join(medllama_7b_match_all, by = "Tumor_Names") %>%
  full_join(llama_32_3b_match_all, by = "Tumor_Names") %>%
  full_join(phi4_match_all, by = "Tumor_Names")%>%
  full_join(llama_match_all,by = "Tumor_Names")%>%
  full_join(MiniLM_L6_v2_all,by="Tumor_Names") %>%
  full_join(mpnet_base_all,by="Tumor_Names")%>%
  full_join(roberta_all,by="Tumor_Names")%>%
  full_join(MiniLM_L12_v2_all,by="Tumor_Names")%>%
  full_join(e5_large_all,by="Tumor_Names")%>%
  full_join(gtr_t5_large_all,by="Tumor_Names")%>%
  full_join(labse_all,by="Tumor_Names")%>%
  full_join(sciBERT_all,by="Tumor_Names")%>%
  full_join(sapBERT_all,by="Tumor_Names")%>%
  full_join(cohere_all,by="Tumor_Names")%>%
  full_join(deepseek_all,by="Tumor_Names")
####
llama_match_5th<- nearest_match_embeddings(distance_matrix_llama_5th,"llama") # Executed
biobert_match_5th<- nearest_match_embeddings(distance_matrix_biobert_5th,"biobert")
medllama_match_5th<- nearest_match_embeddings(distance_matrix_medllama_5th,"medllama")
pubmedbert_match_5th<- nearest_match_embeddings(distance_matrix_pubmedbert_5th,"pubmedbert")
modernbert_match_5th<- nearest_match_embeddings(distance_matrix_modernbert_5th,"modernbert")
medllama_7b_match_5th<- nearest_match_embeddings(distance_matrix_medllama_7b_5th,"medllama_7b")
llama_32_3b_match_5th<- nearest_match_embeddings(distance_matrix_llama32_5th,"llama_32_3b")
phi4_match_5th<- nearest_match_embeddings(distance_matrix_phi_5th,"phi4")
llama33_70b_5th<-nearest_match_embeddings(distance_matrix_llama33_70b_5th,"llama_33_70b")

MiniLM_L6_v2_5th<-nearest_match_embeddings(distance_matrix_MiniLM_L6_v2_5th,"MiniLM_L6_v2")
mpnet_base_5th<-nearest_match_embeddings(distance_matrix_mpnet_base_5th,"mpnet_base")

e5_large_5th<-nearest_match_embeddings(distance_matrix_e5_large_5th,"e5_large")
gtr_t5_large_5th<-nearest_match_embeddings(distance_matrix_gtr_t5_large_5th,"gtr_t5_large")

roberta_5th<-nearest_match_embeddings(distance_matrix_roberta_5th,"roberta")
MiniLM_L12_v2_5th<-nearest_match_embeddings(distance_matrix_MiniLM_L12_v2_5th,"MiniLM_L12_v2")

labse_5th<-nearest_match_embeddings(distance_matrix_labse_5th,"labase")
sciBERT_5th<-nearest_match_embeddings(distance_matrix_sciBERT_5th,"sciBERT")
sapBERT_5th<-nearest_match_embeddings(distance_matrix_sapBERT_5th,"sapBERT")
cohere_5th<-nearest_match_embeddings(distance_matrix_cohere_5th,"cohere")

deepseek_5th<-nearest_match_embeddings(distance_matrix_deepseek_5th,"deepseek")



embedding_nearest_5th <- llama_match_5th %>%
  full_join(biobert_match_5th, by = "Tumor_Names") %>%
  full_join(medllama_match_5th, by = "Tumor_Names") %>%
  full_join(pubmedbert_match_5th, by = "Tumor_Names") %>%
  full_join(modernbert_match_5th, by = "Tumor_Names")%>%
  full_join(medllama_7b_match_5th, by = "Tumor_Names") %>%
  full_join(llama_32_3b_match_5th, by = "Tumor_Names") %>%
  full_join(phi4_match_5th, by = "Tumor_Names")%>%
  full_join(llama_match_5th,by = "Tumor_Names")%>%
  full_join(MiniLM_L6_v2_5th,by="Tumor_Names") %>%
  full_join(mpnet_base_5th,by="Tumor_Names")%>%
  full_join(roberta_5th,by="Tumor_Names")%>%
  full_join(MiniLM_L12_v2_5th,by="Tumor_Names")%>%
  full_join(e5_large_5th,by="Tumor_Names")%>%
  full_join(gtr_t5_large_5th,by="Tumor_Names")%>%
  full_join(labse_5th,by="Tumor_Names")%>%
  full_join(sciBERT_5th,by="Tumor_Names")%>%
  full_join(sapBERT_5th,by="Tumor_Names")%>%
  full_join(cohere_5th,by="Tumor_Names")%>%
  full_join(deepseek_5th,by="Tumor_Names")

#### NCIT 
distance_matrix_llama_ncit <- distance_matrix_llama %>% dplyr::select(all_of(ncit_terms$Tumor_Names))
distance_matrix_biobert_ncit <- distance_matrix_biobert %>% dplyr::select(all_of(ncit_terms$Tumor_Names))
distance_matrix_medllama_ncit <- distance_matrix_medllama %>% dplyr::select(all_of(ncit_terms$Tumor_Names))
distance_matrix_pubmedbert_ncit <- distance_matrix_pubmedbert %>% dplyr::select(all_of(ncit_terms$Tumor_Names))
distance_matrix_modernbert_ncit <- distance_matrix_modernbert %>% dplyr::select(all_of(ncit_terms$Tumor_Names))
distance_matrix_medllama_7b_ncit <- distance_matrix_medllama_7b %>% dplyr::select(all_of(ncit_terms$Tumor_Names))
distance_matrix_llama32_ncit <- distance_matrix_llama32 %>% dplyr::select(all_of(ncit_terms$Tumor_Names))
distance_matrix_phi_ncit <- distance_matrix_phi %>% dplyr::select(all_of(ncit_terms$Tumor_Names))
distance_matrix_llama33_70b_ncit <- distance_matrix_llama33_70b %>% dplyr::select(all_of(ncit_terms$Tumor_Names))
distance_matrix_MiniLM_L6_v2_ncit <- distance_matrix_all_MiniLM_L6_v2 %>% dplyr::select(all_of(ncit_terms$Tumor_Names))
distance_matrix_mpnet_base_ncit <- distance_matrix_all_mpnet_base_v2 %>% dplyr::select(all_of(ncit_terms$Tumor_Names))
distance_matrix_e5_large_ncit<- distance_matrix_e5_large %>% dplyr::select(all_of(ncit_terms$Tumor_Names))
distance_matrix_gtr_t5_large_ncit<-distance_matrix_gtr_t5_large %>% dplyr::select(all_of(ncit_terms$Tumor_Names))
distance_matrix_roberta_ncit<-distance_matrix_roberta%>% dplyr::select(all_of(ncit_terms$Tumor_Names))
distance_matrix_MiniLM_L12_v2_ncit<-distance_matrix_MiniLM_L12_v2%>% dplyr::select(all_of(ncit_terms$Tumor_Names))

distance_matrix_labse_ncit<-distance_matrix_labse%>% dplyr::select(all_of(ncit_terms$Tumor_Names))
distance_matrix_sciBERT_ncit<-distance_matrix_sciBERT%>% dplyr::select(all_of(ncit_terms$Tumor_Names))
distance_matrix_sapBERT_ncit<-distance_matrix_sapBERT%>% dplyr::select(all_of(ncit_terms$Tumor_Names))
distance_matrix_cohere_ncit<-distance_matrix_cohere%>% dplyr::select(all_of(ncit_terms$Tumor_Names))

distance_matrix_deepseek_ncit<-distance_matrix_deepseek%>% dplyr::select(all_of(ncit_terms$Tumor_Names))

####
llama_match_ncit<- nearest_match_embeddings(distance_matrix_llama_ncit,"llama") # Executed
biobert_match_ncit<- nearest_match_embeddings(distance_matrix_biobert_ncit,"biobert")
medllama_match_ncit<- nearest_match_embeddings(distance_matrix_medllama_ncit,"medllama")
pubmedbert_match_ncit<- nearest_match_embeddings(distance_matrix_pubmedbert_ncit,"pubmedbert")
modernbert_match_ncit<- nearest_match_embeddings(distance_matrix_modernbert_ncit,"modernbert")
medllama_7b_match_ncit<- nearest_match_embeddings(distance_matrix_medllama_7b_ncit,"medllama_7b")
llama_32_3b_match_ncit<- nearest_match_embeddings(distance_matrix_llama32_ncit,"llama_32_3b")
phi4_match_ncit<- nearest_match_embeddings(distance_matrix_phi_ncit,"phi4")
llama33_70b_ncit<-nearest_match_embeddings(distance_matrix_llama33_70b_ncit,"llama_33_70b")

MiniLM_L6_v2_ncit<-nearest_match_embeddings(distance_matrix_MiniLM_L6_v2_ncit,"MiniLM_L6_v2")
mpnet_base_ncit<-nearest_match_embeddings(distance_matrix_mpnet_base_ncit,"mpnet_base")

e5_large_ncit<-nearest_match_embeddings(distance_matrix_e5_large_ncit,"e5_large")
gtr_t5_large_ncit<-nearest_match_embeddings(distance_matrix_gtr_t5_large_ncit,"gtr_t5_large")
roberta_ncit<-nearest_match_embeddings(distance_matrix_roberta_ncit,"roberta")
MiniLM_L12_v2_ncit<-nearest_match_embeddings(distance_matrix_MiniLM_L12_v2_ncit,"MiniLM_L12_v2")

labse_ncit<-nearest_match_embeddings(distance_matrix_labse_ncit,"labase")
sciBERT_ncit<-nearest_match_embeddings(distance_matrix_sciBERT_ncit,"sciBERT")
sapBERT_ncit<-nearest_match_embeddings(distance_matrix_sapBERT_ncit,"sapBERT")
cohere_ncit<-nearest_match_embeddings(distance_matrix_cohere_ncit,"cohere")

deepseek_ncit<-nearest_match_embeddings(distance_matrix_deepseek_ncit,"deepseek")


embedding_nearest_ncit <- llama_match_ncit %>%
  full_join(biobert_match_ncit, by = "Tumor_Names") %>%
  full_join(medllama_match_ncit, by = "Tumor_Names") %>%
  full_join(pubmedbert_match_ncit, by = "Tumor_Names") %>%
  full_join(modernbert_match_ncit, by = "Tumor_Names") %>%
  full_join(medllama_7b_match_ncit, by = "Tumor_Names")%>%
  full_join(llama_32_3b_match_ncit, by = "Tumor_Names") %>%
  full_join(phi4_match_ncit, by = "Tumor_Names")%>%
  full_join(llama33_70b_ncit, by = "Tumor_Names")%>%
  full_join(MiniLM_L6_v2_ncit, by = "Tumor_Names")%>%
  full_join(mpnet_base_ncit, by = "Tumor_Names")%>%
  full_join(roberta_ncit,by="Tumor_Names")%>%
  full_join(MiniLM_L12_v2_ncit,by="Tumor_Names")%>%
  full_join(e5_large_ncit,by="Tumor_Names")%>%
  full_join(gtr_t5_large_ncit,by="Tumor_Names")%>%
  full_join(labse_ncit,by="Tumor_Names")%>%
  full_join(sciBERT_ncit,by="Tumor_Names")%>%
  full_join(sapBERT_ncit,by="Tumor_Names")%>%
  full_join(cohere_ncit,by="Tumor_Names")%>%
  full_join(deepseek_ncit,by="Tumor_Names")


# Reconcile Results

NCIT_Results_all <- read_csv("analysis/results/NCIT_Results_all.csv")
NCIT_Results_5th <- read_csv("analysis/results_5th/NCIT_Results_5thed.csv")
WHO_Results_all <- read_csv("analysis/results/WHO_Results_all.csv")
WHO_Results_5th <- read_csv("analysis/results_5th/WHO_Results_5thed.csv")

NCIT_Results_all<-NCIT_Results_all[,c(-1)]
NCIT_Results_5th<-NCIT_Results_5th[,c(-1)]
WHO_Results_all<-WHO_Results_all[,c(-1)]
WHO_Results_5th<-WHO_Results_5th[,c(-1)]

names_col <- c("Tumor_Names","euclidean_dist_llama","euclidean_dist_biobert","euclidean_dist_medllama", "euclidean_dist_pubmedbert",
               "euclidean_dist_modernbert","euclidean_dist_medllama_7b",   "euclidean_dist_llama_32_3b","euclidean_dist_phi4",
               "euclidean_dist_llama_33_70b","euclidean_dist_MiniLM_L6_v2", "euclidean_mpnet_base", "euclidean_dist_roberta",
               "euclidean_dist_MiniLM_L12_v2","euclidean_dist_e5_large","euclidean_dist_gtr_t5_large", "euclidean_dist_labse",
               "euclidean_dist_sciBERT","euclidean_dist_sapBERT","euclidean_dist_cohere","euclidean_dist_deepseek")


colnames(embedding_nearest_ncit) <- names_col
colnames(embedding_nearest_all) <- names_col
colnames(embedding_nearest_5th) <- names_col
  
NCIT_Results_all <- NCIT_Results_all %>% left_join(embedding_nearest_ncit , by="Tumor_Names")
NCIT_Results_5th <- NCIT_Results_5th %>% left_join(embedding_nearest_ncit , by="Tumor_Names")

WHO_Results_all <- WHO_Results_all %>% left_join(embedding_nearest_all , by="Tumor_Names")
WHO_Results_5th <- WHO_Results_5th %>% left_join(embedding_nearest_5th , by="Tumor_Names")

write.csv(NCIT_Results_all,paste(analysis_dir,"/results/NCIT_Results_all_ose.csv",sep=""))
write.csv(WHO_Results_all,paste(analysis_dir,"/results/WHO_Results_all_ose.csv",sep=""))

write.csv(NCIT_Results_5th,paste(analysis_dir,"/results_5th/NCIT_Results_5thed_ose.csv",sep=""))
write.csv(WHO_Results_5th,paste(analysis_dir,"/results_5th/WHO_Results_5thed_ose.csv",sep=""))



#save.image("12-prior.RData")

#save.image("12-prior-jan15-onlyllama32-phi4.RData")
NCIT_Results_all <- read_csv(paste(analysis_dir,"/results/NCIT_Results_all_ose.csv",sep=""))
NCIT_Results_5th <- read_csv(paste(analysis_dir,"/results_5th/NCIT_Results_5thed_ose.csv",sep=""))
WHO_Results_all <- read_csv(paste(analysis_dir,"/results/WHO_Results_all_ose.csv",sep=""))
WHO_Results_5th <- read_csv(paste(analysis_dir,"/results_5th/WHO_Results_5thed_ose.csv",sep=""))
NCIT_Results_all<-NCIT_Results_all[,c(-1)]
NCIT_Results_5th<-NCIT_Results_5th[,c(-1)]
WHO_Results_all<-WHO_Results_all[,c(-1)]
WHO_Results_5th<-WHO_Results_5th[,c(-1)]

# NCIT_Results_all<-NCIT_Results_all%>%left_join(llama_32_3b_match_ncit, by = "Tumor_Names")
# NCIT_Results_all<-NCIT_Results_all%>%left_join(phi4_match_ncit, by = "Tumor_Names")
# NCIT_Results_all<-NCIT_Results_all%>%left_join(llama33_70b_ncit, by = "Tumor_Names")
# NCIT_Results_all<-NCIT_Results_all%>%left_join(MiniLM_L6_v2_ncit, by = "Tumor_Names")
# NCIT_Results_all<-NCIT_Results_all%>%left_join(mpnet_base_ncit, by = "Tumor_Names")
# NCIT_Results_all<-NCIT_Results_all%>%left_join(e5_large_ncit, by = "Tumor_Names")
#NCIT_Results_all<-NCIT_Results_all%>%left_join(gtr_t5_large_ncit, by = "Tumor_Names")

# NCIT_Results_5th<-NCIT_Results_5th%>%left_join(llama_32_3b_match_ncit, by = "Tumor_Names")
# NCIT_Results_5th<-NCIT_Results_5th%>%left_join(phi4_match_ncit, by = "Tumor_Names")
# NCIT_Results_5th<-NCIT_Results_5th%>%left_join(llama33_70b_ncit, by = "Tumor_Names")
#NCIT_Results_5th<-NCIT_Results_5th%>%left_join(MiniLM_L6_v2_ncit, by = "Tumor_Names")
#NCIT_Results_5th<-NCIT_Results_5th%>%left_join(mpnet_base_ncit, by = "Tumor_Names")
#NCIT_Results_5th<-NCIT_Results_5th%>%left_join(e5_large_ncit, by = "Tumor_Names")
#NCIT_Results_5th<-NCIT_Results_5th%>%left_join(gtr_t5_large_ncit, by = "Tumor_Names")

# WHO_Results_5th<-WHO_Results_5th%>%left_join(llama_32_3b_match_5th, by = "Tumor_Names")
# WHO_Results_5th<-WHO_Results_5th%>%left_join(phi4_match_5th, by = "Tumor_Names")
# WHO_Results_5th<-WHO_Results_5th%>%left_join(llama33_70b_5th, by = "Tumor_Names")
# WHO_Results_5th<-WHO_Results_5th%>%left_join(MiniLM_L6_v2_5th, by = "Tumor_Names")
#WHO_Results_5th<-WHO_Results_5th%>%left_join(mpnet_base_5th, by = "Tumor_Names")
#WHO_Results_5th<-WHO_Results_5th%>%left_join(e5_large_5th, by = "Tumor_Names")
#WHO_Results_5th<-WHO_Results_5th%>%left_join(gtr_t5_large_5th, by = "Tumor_Names")

# WHO_Results_all<-WHO_Results_all%>%left_join(llama_32_3b_match_all, by = "Tumor_Names")
# WHO_Results_all<-WHO_Results_all%>%left_join(phi4_match_all, by = "Tumor_Names")
# WHO_Results_all<-WHO_Results_all%>%left_join(llama33_70b_all, by = "Tumor_Names")
# WHO_Results_all<-WHO_Results_all%>%left_join(MiniLM_L6_v2_all, by = "Tumor_Names")
#WHO_Results_all<-WHO_Results_all%>%left_join(mpnet_base_all, by = "Tumor_Names")
#WHO_Results_all<-WHO_Results_all%>%left_join(e5_large_all, by = "Tumor_Names")
#WHO_Results_all<-WHO_Results_all%>%left_join(gtr_t5_large_all, by = "Tumor_Names")

# NCIT_Results_all<-NCIT_Results_all%>%left_join(roberta_ncit, by = "Tumor_Names")
#  NCIT_Results_all<-NCIT_Results_all%>%left_join(MiniLM_L12_v2_ncit, by = "Tumor_Names")
#  NCIT_Results_5th<-NCIT_Results_5th%>%left_join(roberta_ncit, by = "Tumor_Names")
# > NCIT_Results_5th<-NCIT_Results_5th%>%left_join(MiniLM_L12_v2_ncit, by = "Tumor_Names")
# > WHO_Results_5th<-WHO_Results_5th%>%left_join(roberta_5th,by="Tumor_Names")
# > WHO_Results_5th<-WHO_Results_5th%>%left_join(MiniLM_L12_v2_5th,by="Tumor_Names")
# > WHO_Results_all<-WHO_Results_all%>%left_join(roberta_all,by="Tumor_Names")
# > WHO_Results_all<-WHO_Results_all%>%left_join(MiniLM_L12_v2_all,by="Tumor_Names")
# > WHO_Results_all<-WHO_Results_all%>%left_join(labse_all, by = "Tumor_Names")
# > WHO_Results_all<-WHO_Results_all%>%left_join(sciBERT_all, by = "Tumor_Names")
# > WHO_Results_all<-WHO_Results_all%>%left_join(sapBERT_all, by = "Tumor_Names")
# > WHO_Results_all<-WHO_Results_all%>%left_join(cohere_all, by = "Tumor_Names")
# > 
#   > WHO_Results_5th<-WHO_Results_5th%>%left_join(labse_5th, by = "Tumor_Names")
# > WHO_Results_5th<-WHO_Results_5th%>%left_join(sciBERT_5th, by = "Tumor_Names")
# > WHO_Results_5th<-WHO_Results_5th%>%left_join(sapBERT_5th, by = "Tumor_Names")
# > WHO_Results_5th<-WHO_Results_5th%>%left_join(cohere_5th, by = "Tumor_Names")

# colnames(NCIT_Results_all)[c(28,29)]<-c("euclidean_dist_roberta", "euclidean_dist_MiniLM_L12_v2")
# colnames(NCIT_Results_5th)[c(28,29)]<-c("euclidean_dist_roberta", "euclidean_dist_MiniLM_L12_v2")
# colnames(WHO_Results_all)[c(28,29)]<-c("euclidean_dist_roberta", "euclidean_dist_MiniLM_L12_v2")
# colnames(WHO_Results_5th)[c(28,29)]<-c("euclidean_dist_roberta", "euclidean_dist_MiniLM_L12_v2")