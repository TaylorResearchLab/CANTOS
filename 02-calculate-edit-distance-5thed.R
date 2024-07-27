suppressPackageStartupMessages({
  library(dplyr)
  library(doParallel)
  library(foreach)
  library(stringdist)
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

# Load dissimalirity matrices
load("/Users/lahiria/Desktop/MTP_Paper/temp/CT-Large-File-June21/dissimilarity_matrix_lv.RData")
load("/Users/lahiria/Desktop/MTP_Paper/temp/CT-Large-File-June21/dissimilarity_matrix_cosine.RData")
load("/Users/lahiria/Desktop/MTP_Paper/temp/CT-Large-File-June21/dissimilarity_matrix_jw.RData")
source(paste(util_dir,"/distance_clusters.R",sep=""))

dissimilarity_matrix_cosine<-as.data.frame(dissimilarity_matrix_cosine)
dissimilarity_matrix_lv<-as.data.frame(dissimilarity_matrix_lv)
dissimilarity_matrix_jw<-as.data.frame(dissimilarity_matrix_jw)


WHO_Terms_All <-readxl::read_xlsx(paste(data_dir,"/WHO_Tumors/result/WHO_Tumor_all_edition.xlsx",sep=""))
WHO_Terms_5th<-WHO_Terms_All%>%filter(edition_5th=="Yes")

NCIT_Terms <-read.csv(paste(data_dir,"/dt_input_file_6_dec/NCIT_Neoplasm_Core_terms_text-embedding-ada-002_embeddings.csv",sep=""))[,1]
NCIT_Terms<-tolower(NCIT_Terms[c(-1)])

ct_disease_annot_adult_ped_df<-read.csv(paste(input_dir,"/tumor_annotated_adult_ped.csv",sep=""))
ct_tumor_df<-ct_disease_annot_adult_ped_df%>%filter(validated_cancer_tumor=="Yes")
ct_tumor_df<-ct_tumor_df[,c(-1)]

tumor_names<- unique(c(ct_tumor_df$diseases,NCIT_Terms,WHO_Terms_5th$Tumor_Names))

dissimilarity_matrix_cosine<- dissimilarity_matrix_cosine %>% dplyr::filter(rownames(dissimilarity_matrix_cosine) %in% 
                                                                              tumor_names )
dissimilarity_matrix_lv<- dissimilarity_matrix_lv %>% dplyr::filter(rownames(dissimilarity_matrix_lv) %in% 
                                                                              tumor_names )
dissimilarity_matrix_jw<- dissimilarity_matrix_jw %>% dplyr::filter(rownames(dissimilarity_matrix_jw) %in% 
                                                                              tumor_names )


dissimilarity_matrix_cosine<- dissimilarity_matrix_cosine %>%dplyr::select(any_of(tumor_names))
dissimilarity_matrix_lv<- dissimilarity_matrix_lv %>%dplyr::select(any_of(tumor_names))
dissimilarity_matrix_jw<- dissimilarity_matrix_jw %>% dplyr::select(any_of(tumor_names))


df_tumor_names<- colnames(dissimilarity_matrix_lv)



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
