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

# Load embeddings
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
ADA_002_Embeddings<-ADA_002_Embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")



