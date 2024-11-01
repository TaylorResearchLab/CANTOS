# Loads embeddings ADA2 for CT, WHO, NCIT  Tumors and then performs PCA. 
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
intermediate_dir <- file.path(analysis_dir,"intermediate_5th")

#
WHO_Terms_All <-readxl::read_xlsx(paste(data_dir,"/WHO_Tumors/result/WHO_Tumor_all_edition.xlsx",sep=""))
WHO_Terms_5th<-WHO_Terms_All%>%filter(edition_5th=="Yes")


CT_embedding_df <- read.csv(paste(data_dir,"/CT_Embeddings_ADA2.csv",sep=""))
CT_embedding_df<-CT_embedding_df[,c(-1)]
# Read NCIT Terms and WHO Terms with embedding and join them to the rest of the embedding list 
NCIT_embedding_df <-read.csv(paste(data_dir,"/dt_input_file_6_dec/NCIT_Neoplasm_Core_terms_text-embedding-ada-002_embeddings.csv",sep=""))
NCIT_embedding_df<-NCIT_embedding_df[c(-1),] # Remove the header (column name) embedding

#WHO_embedding_df <-read.csv(paste(data_dir,"/dt_input_file_6_dec/WHO_Only_terms_text-embedding-ada-002_embeddings.csv",sep=""))
WHO_embedding_df <- read.csv(paste(data_dir,"/WHO_Aggregate_ADA2.csv",sep="")) #
WHO_embedding_df<-WHO_embedding_df[,c(-1)]
WHO_embedding_df<-WHO_embedding_df %>% dplyr::filter(Disease %in% WHO_Terms_5th$Tumor_Names)

rownames(NCIT_embedding_df)<-NULL
rownames(WHO_embedding_df)<-NULL

# Combined Embedding File for PC and Cluster Analysis
colnames(CT_embedding_df)<-colnames(NCIT_embedding_df) # CT embedding columns need to be fixed
colnames(WHO_embedding_df)<-colnames(NCIT_embedding_df)

combined_embedding_df <- rbind(CT_embedding_df,NCIT_embedding_df,WHO_embedding_df)
combined_embedding_df$Disease<- tolower(combined_embedding_df$Disease)
combined_embedding_df <- as.data.frame(combined_embedding_df %>% group_by(Disease) %>% summarise(across(everything(), list(mean))))
colnames(combined_embedding_df)<-colnames(NCIT_embedding_df)
rownames(combined_embedding_df)<-combined_embedding_df$Disease
combined_embedding_df<-combined_embedding_df[,c(-1)]

# Now perform PCA 
results <- prcomp(combined_embedding_df, scale = TRUE)
eigs <- results$sdev^2
ff<-rbind(SD = sqrt(eigs),Proportion = eigs/sum(eigs), Cumulative = cumsum(eigs)/sum(eigs))
print(sum(ff[2,1:136]))
# Select the top 136 PCs
disease_transform = as.data.frame(-results$x[,1:136])

write.csv(disease_transform,file =paste(intermediate_dir,"/disease_transform_pca_ada2_5thed.csv",sep="") )

#save.image(file = "script4a_5thed.RData")
save(combined_embedding_df,file = paste(intermediate_dir,"/combined_embedding_ada2_df_5thed.RData",sep=""))

