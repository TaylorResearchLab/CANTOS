# Loads embeddings V3 for CT, WHO, NCIT  Tumors and then performs PCA. 
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


# Read CT embedding file 
ct_tumor_embeddings_v3_df<-read.csv(paste(data_dir,"/CT_Embeddings_V3.csv",sep=""))
ct_tumor_embeddings_v3_df<-ct_tumor_embeddings_v3_df[,c(-1)]

WHO_Terms_V3 <- read.csv(paste(data_dir,"/WHO_Terms_All_V3.csv",sep="")) #
WHO_Terms_V3<-WHO_Terms_V3[,c(-1)]

NCIT_embedding_v3<-read.csv(paste(data_dir,"/NCIT_Embeddings_V3.csv",sep="")) #1395
NCIT_embedding_v3<-NCIT_embedding_v3[,c(-1)]

combined_embeddings_df<-rbind(ct_tumor_embeddings_v3_df,WHO_Terms_V3,NCIT_embedding_v3)#19345
combined_embeddings_df <- combined_embeddings_df%>%group_by(Tumor_Names) %>% summarise(across(everything(), list(mean)))#18284
combined_embeddings_df<-as.data.frame(combined_embeddings_df)
rownames(combined_embeddings_df)<-combined_embeddings_df$Tumor_Names
colnames(combined_embeddings_df)<-colnames(NCIT_embedding_v3)

combined_embeddings_df<-combined_embeddings_df[,c(-1)]
# PCA for v3 

results_v3 <- prcomp(combined_embeddings_df, scale = TRUE)
eigs_v3 <- results_v3$sdev^2
ff_v3<-rbind(SD = sqrt(eigs_v3),Proportion = eigs_v3/sum(eigs_v3), Cumulative = cumsum(eigs_v3)/sum(eigs_v3))
# Check the number of components needed to capture 80% variance at least
print(sum(ff_v3[2,1:185]))
disease_transform_v3<-as.data.frame(-results_v3$x[,1:185])
# Save this file "

write.csv(disease_transform_v3,file =paste(intermediate_dir,"/disease_transform_pca_v3.csv",sep="") )
save.image(file = "script4b.RData")
save(combined_embeddings_df,file = paste(intermediate_dir,"/combined_embedding_v3_df.RData",sep=""))
