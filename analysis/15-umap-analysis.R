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
# Load biobert
biobert_Embeddings <- read.csv(paste(data_dir,"/Embeddings/biobert_embedding.csv",sep=""))
# Load pubmedbert
pubmedbert_Embeddings <- read.csv(paste(data_dir,"/Embeddings/pubmedbert-base-embeddings.csv",sep=""))

medllama_Embeddings<-medllama_Embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
llama_Embeddings<-llama_Embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
biobert_Embeddings<-biobert_Embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
pubmedbert_Embeddings<-pubmedbert_Embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")




