method_name <- "Levenshtein"  
start_time <- Sys.time()

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


source(paste(util_dir,"/string_dissimilarity.R",sep = ""))
source(paste(util_dir,"/distance_clusters.R",sep=""))
source(paste(util_dir,"/string_normalizing.R",sep=""))

# Read the annotated file
ct_disease_annot_adult_ped_df<-read.csv(paste(input_dir,"/tumor_annotated_adult_ped.csv",sep=""))
ct_tumor_df<-ct_disease_annot_adult_ped_df%>%filter(validated_cancer_tumor=="Yes")
ct_tumor_df<-ct_tumor_df[,c(-1)]
WHO_Terms_All <-readxl::read_xlsx(paste(data_dir,"/WHO_Tumors/result/WHO_Tumor_all_edition.xlsx",sep=""))
NCIT_Terms <-read.csv(paste(data_dir,"/dt_input_file_6_dec/NCIT_Neoplasm_Core_terms_text-embedding-ada-002_embeddings.csv",sep=""))[,1]
NCIT_Terms<-tolower(NCIT_Terms[c(-1)])
WHO_Terms_5th<-WHO_Terms_All%>%filter(edition_5th=="Yes")






# Levenstein distance between tumors 
df_tumor_combined<-as.data.frame(unique(c(ct_tumor_df$diseases,NCIT_Terms,WHO_Terms_All$Tumor_Names)))
colnames(df_tumor_combined)[1]<-"Tumor"

df_tumor_names<-unique(df_tumor_combined$Tumor)
dissimilarity_matrix_lv <- as.data.frame(matrix(nrow=length(df_tumor_names),ncol=length(df_tumor_names)))
rownames(dissimilarity_matrix_lv)<-df_tumor_names
colnames(dissimilarity_matrix_lv)<-df_tumor_names


cl <- makeCluster(5, outfile="")
registerDoParallel(cl)


dissimilarity_matrix_lv<-foreach(iter=1:length(df_tumor_names),.combine=rbind) %dopar% {
#dissimilarity_matrix_lv<-foreach(iter=1:20,.combine=rbind) %dopar% {
    print(iter)
  disease_name <- colnames(dissimilarity_matrix_lv)[iter]
  distances<-unlist(lapply(df_tumor_names,string_dissimilarity,S2=disease_name,meth="lv"))
}

#rownames(dissimilarity_matrix_lv) <- df_tumor_names[1:20]
rownames(dissimilarity_matrix_lv) <- df_tumor_names
colnames(dissimilarity_matrix_lv) <- df_tumor_names
stopCluster(cl)


df_tumor_names<- colnames(dissimilarity_matrix_lv)
normalizing_matrix_lv <- as.data.frame(matrix(nrow=length(df_tumor_names),ncol=length(df_tumor_names)))
#normalizing_matrix_lv <- as.data.frame(matrix(nrow=20,ncol=length(df_tumor_names)))
#rownames(normalizing_matrix_lv)<-df_tumor_names[1:20]
rownames(normalizing_matrix_lv)<-df_tumor_names
colnames(normalizing_matrix_lv)<-df_tumor_names

cl <- makeCluster(5, outfile="")
registerDoParallel(cl)
normalizing_matrix_lv<-foreach(iter=1:length(df_tumor_names),.combine=rbind) %dopar% {
#normalizing_matrix_lv<-foreach(iter=1:20,.combine=rbind) %dopar% {
  print(iter)
  disease_name <- colnames(dissimilarity_matrix_lv)[iter]
  norm_factors<-unlist(lapply(df_tumor_names,string_normalzing,S2=disease_name))
}
stopCluster(cl)


#rownames(normalizing_matrix_lv) <- df_tumor_names[1:20]
rownames(normalizing_matrix_lv) <- df_tumor_names
colnames(normalizing_matrix_lv) <- df_tumor_names
dissimilarity_matrix_lv<-dissimilarity_matrix_lv/normalizing_matrix_lv


dissimilarity_matrix_lv<-as.data.frame(dissimilarity_matrix_lv)

dissimilarity_matrix_lv_who_all <- dissimilarity_matrix_lv %>% dplyr::select(one_of(WHO_Terms_All$Tumor_Names))
dissimilarity_matrix_lv_who_5th <- dissimilarity_matrix_lv %>% dplyr::select(one_of(WHO_Terms_5th$Tumor_Names))
dissimilarity_matrix_lv_ncit <- dissimilarity_matrix_lv %>% dplyr:::select(one_of(NCIT_Terms))

index_min_who_all_lv <- as.matrix(apply(dissimilarity_matrix_lv_who_all, 1, which.min))
who_all_match_lv_df <- as.data.frame(cbind(rownames(dissimilarity_matrix_lv_who_all)))

index_min_who_5th_lv <- as.matrix(apply(dissimilarity_matrix_lv_who_5th, 1, which.min))
who_5th_match_lv_df <- as.data.frame(cbind(rownames(dissimilarity_matrix_lv_who_5th)))

index_min_ncit_lv <- as.matrix(apply(dissimilarity_matrix_lv_ncit, 1, which.min))
ncit_match_lv_df <- as.data.frame(cbind(rownames(dissimilarity_matrix_lv_ncit)))


who_all_match_lv_df$WHO_Matches<- NA
who_all_match_lv_df$WHO_distance<-NA

who_5th_match_lv_df$WHO_Matches<- NA
who_5th_match_lv_df$WHO_distance<-NA

ncit_match_lv_df$WHO_Matches<- NA
ncit_match_lv_df$WHO_distance<-NA

for (iter in 1: dim(who_all_match_lv_df)[1]){
  
  who_all_match_lv_df$WHO_Matches[iter] <- colnames(dissimilarity_matrix_lv_who_all)[index_min_who_all_lv[iter]]
  who_all_match_lv_df$WHO_distance[iter]<-dissimilarity_matrix_lv_who_all[iter,index_min_who_all_lv[iter]]
  
  who_5th_match_lv_df$WHO_Matches[iter] <- colnames(dissimilarity_matrix_lv_who_5th)[index_min_who_5th_lv[iter]]
  who_5th_match_lv_df$WHO_distance[iter]<-dissimilarity_matrix_lv_who_5th[iter,index_min_who_5th_lv[iter]]
  
  ncit_match_lv_df$ncit_Matches[iter] <- colnames(dissimilarity_matrix_lv_ncit)[index_min_ncit_lv[iter]]
  ncit_match_lv_df$ncit_distance[iter]<-dissimilarity_matrix_lv_ncit[iter,index_min_ncit_lv[iter]]
  
}

end_time <- Sys.time()
elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

benchmark <- data.frame(
  Method = method_name,
  Runtime_sec = round(elapsed_time, 2),
  stringsAsFactors = FALSE
)

output_file <- paste(analysis_dir,"/run-time-analysis/runtime_benchmarks_all_methods.csv",sep="")
if (!file.exists(output_file)) {
  write.csv(benchmark, output_file, row.names = FALSE)
} else {
  write.table(benchmark, output_file, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
}
save.image(paste(analysis_dir,"/run-time-analysis/lv_match.RData",sep=""))