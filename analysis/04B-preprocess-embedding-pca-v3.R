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
embedding_v3_large <- read.csv(paste(data_dir,"/embedding_tumor_names_text-embedding-3-large_embeddings.csv",sep="")) #16139
colnames(embedding_v3_large)[1]<-"Tumor_Names"

missing_v3_tumors <- read.csv(paste(data_dir,"/missing_V3_tumors.csv",sep="")) # 57
colnames(missing_v3_tumors)[1]<-"Tumor_Names"

embedding_v3_large<-rbind(embedding_v3_large,missing_v3_tumors)#16196








# Read the annotated file
ct_disease_df <- read.csv(paste(input_dir,"/CT-Aug22-2023-Disease-File - clinical_trial_disease_aug_22_2023.csv",sep=""))
ct_tumor_df<- ct_disease_df %>% filter(validated_cancer_tumor=="Yes") %>%dplyr::select(diseases)
colnames(ct_tumor_df)<-"Tumor_Names"
ct_tumor_embeddings_df<- ct_tumor_df %>% left_join(embedding_v3_large,by="Tumor_Names")


WHO_3rd_4th_5th_text_embedding_3_large_embeddings <- read_csv(paste(data_dir,"/WHO_3rd_4th_5th_text-embedding-3-large_embeddings.csv",sep="")) #2330
colnames(WHO_3rd_4th_5th_text_embedding_3_large_embeddings)<- colnames(embedding_v3_large)

NCIT_embedding_df <-read.csv(paste(data_dir,"/dt_input_file_6_dec/NCIT_Neoplasm_Core_terms_text-embedding-ada-002_embeddings.csv",sep=""))
NCIT_embedding_df<-NCIT_embedding_df[c(-1),] # Remove the header (column name) embedding
NCIT_embedding_v3 <- as.data.frame(NCIT_embedding_df$Disease)
colnames(NCIT_embedding_v3)<-"Tumor_Names"
NCIT_embedding_v3$Tumor_Names<-tolower(NCIT_embedding_v3$Tumor_Names)
NCIT_embedding_v3<- NCIT_embedding_v3 %>% left_join(embedding_v3_large,by="Tumor_Names")

combined_embeddings_df<-rbind(ct_tumor_embeddings_df,WHO_3rd_4th_5th_text_embedding_3_large_embeddings,NCIT_embedding_v3)#17054
combined_embeddings_df <- combined_embeddings_df%>%group_by(Tumor_Names) %>% summarise(across(everything(), list(mean)))#18284


# PCA for v3 
rownames(combined_embeddings_df)<-combined_embeddings_df$Tumor_Names
combined_embeddings_df<-combined_embeddings_df[,c(-1)]
results_v3 <- prcomp(combined_embeddings_df, scale = TRUE)
eigs_v3 <- results_v3$sdev^2
ff_v3<-rbind(SD = sqrt(eigs_v3),Proportion = eigs_v3/sum(eigs_v3), Cumulative = cumsum(eigs_v3)/sum(eigs_v3))
# Check the number of components needed to capture 80% variance at least
print(sum(ff_v3[2,1:187]))
disease_transform_v3<-as.data.frame(-results_v3$x[,1:187])

# Save this file "
write.csv(disease_transform,file =paste(intermediate_dir,"/disease_transform_pca.csv",sep="") )
save(combined_embedding_df,file = paste(intermediate_dir,"/combined_embedding_df.RData",sep=""))
write.csv(disease_transform_v3,file =paste(intermediate_dir,"/disease_transform_pca_v3.csv",sep="") )
