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
analyses_dir <- file.path(root_dir,"analyses")


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

#embeddings_lama<- embed_text(tumor_terms,model = "llama3",server = NULL,model_params = NULL,
#                             verbose = getOption("rollama_verbose", default = interactive()))

embeddings_lama<- read.csv(paste(data_dir,"/embeddings_llama.csv",sep=""))
biobert_embeddings <- read.csv(paste(data_dir,"/biobert_embedding.csv",sep=""))
medllama_embeddings <- read.csv(paste(data_dir,"/medllama-13b.csv",sep="")) #5121 dimension
pubmedbert_embeddings <-read.csv(paste(data_dir,"/pubmedbert-base-embeddings.csv",sep="")) #768 dimension

embeddings_lama<-embeddings_lama %>% group_by(Tumor_Names) %>% summarise_all("mean")
biobert_embeddings<-biobert_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
medllama_embeddings<-medllama_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
pubmedbert_embeddings<-pubmedbert_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
