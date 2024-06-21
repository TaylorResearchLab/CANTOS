#Load the manually annotated disease file with pediatric and adult cancer annotation.


# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(doParallel)
  library(foreach)
  library(stringdist)
  #library(DescTools)
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

#Load Functions
source(paste(util_dir,"/string_dissimilarity.R",sep = ""))
source(paste(util_dir,"/distance_clusters.R",sep=""))
# Read the annotated file
#ct_disease_df <- read.csv(paste(input_dir,"/CT-Aug22-2023-Disease-File - clinical_trial_disease_aug_22_2023.csv",sep=""))
#ct_disease_df <- read.csv(paste(input_dir,"/cancer_annotated_file_ammended.csv",sep=""))
#ct_tumor_df<- ct_disease_df %>% filter(validated_cancer_tumor=="Yes")
#ct_tumor_df<-read_xlsx(paste(input_dir,"/cancer_annotated_file_ammended.xlsx",sep=""))
#ct_tumor_df<- ct_disease_df %>% filter(validated_cancer_tumor=="Yes")
ct_disease_annot_adult_ped_df<-read.csv(paste(input_dir,"/tumor_annotated_adult_ped.csv",sep=""))
ct_tumor_df<-ct_disease_annot_adult_ped_df%>%filter(validated_cancer_tumor=="Yes")
ct_tumor_df<-ct_tumor_df[,c(-1)]
# Add NCIT and WHO Tumors 
# Read NCIT Terms and WHO Terms with embedding and join them to the rest of the embedding list 

# NEW addition : 
WHO_Terms_All <-readxl::read_xlsx(paste(data_dir,"/WHO_Tumors/result/WHO_Tumor_all_edition.xlsx",sep=""))


NCIT_Terms <-read.csv(paste(data_dir,"/dt_input_file_6_dec/NCIT_Neoplasm_Core_terms_text-embedding-ada-002_embeddings.csv",sep=""))[,1]
# WHO_Terms <-read.csv(paste(data_dir,"/dt_input_file_6_dec/WHO_Only_terms_text-embedding-ada-002_embeddings.csv",sep=""))[,1]
# 
NCIT_Terms<-tolower(NCIT_Terms[c(-1)])
# WHO_Terms<-tolower(WHO_Terms[c(-1)])



# show which were removed due to excessive typos ct_tumor_df$diseases[which(! ct_tumor_df$diseases %in% CT_embedding_agg_df$DISEASE_NAMES)]
# [1] "diffuse large b cell lymphomaÔºådlbcl"                                                                                      
# [2] "kaposi¬¥s sarcoma"                                                                                                          
# [3] "hodgkin¬¥s lymphoma"                                                                                                        
# [4] "lymphocyte predominant hodgkin¬¥s lymphoma (lphd)"                                                                          
# [5] "follicular non-hodgkin¬¥s lymphoma"                                                                                         
# [6] "neoplasmsÔºånon-small cell lung cancer"                                                                                     
# [7] "waldenstr√∂m macroglobulinemia"                                                                                             
# [8] "angiogenesis inhibitorsÔºåovarian neoplasms"                                                                                
# [9] "non-small cell lung cancer stage ‚Ö≤a"                                                                                      
# [10] "follicular non-hodgking¬¥s lymphoma refractory or relapsed after treatment with r-chemotherapy in first line."              
# [11] "waldenstr√∂m's macroglobulinemia"                                                                                           
# [12] "chemotherapyÔºõadvanced gastric cancerÔºõcisplatinÔºõdisulfiram"                                                            
# [13] "non-small cell lung cancer stage ‚Ö±"                                                                                       
# [14] "follicular lymphoma grade iii (fl iii¬∞)"                                                                                   
# [15] "stage-‚Ö± colorectal cancer"                                                                                                
# [16] "mycosis fungoides and s√©zary syndrome"                                                                                     
# [17] "non-metastatic, hormone na√Øve prostate cancer"                                                                             
# [18] "carcinomaÔºånon-small-cell lung"                                                                                            
# [19] "childhood non-hodgkin lymphoma"                                                                                             
# [20] "rituximab, lenalidomide, zebutinib Ôºåmantle cell lymphoma"                                                                 
# [21] "treatment-na—óve mantle cell lymphoma"                                                                                      
# [22] "effectivenessÔºåsafetyÔºåthymic cancer"                                                                                     
# [23] "hif-2Œ± mutated cancers"                                                                                                    
# [24] "hormone receptor positive,human epidermal receptor 2 negative, node-positive, high risk, early stageÔºåfemale breast cancer"
# [25] "recurrent, or metastatic cervical cancer with pd-l1 positive (cps‚â•1)"                                                     
# [26] "acute lymphocytic leukemiaÔºå b-cell"                                                                                       
# [27] "advanced digestive system neuroendocrine neoplasm"                                                                          
# [28] "metastatic melanoma (stage iiic non-r√©s√©cable or no surgically curable or stage iv with classification ajcc)"             
# [29] "transformed follicular lymphoma with ‚â• 50% diffuse large cell component"                                                  
# [30] "hormone receptor positiveÔºåher2-negative breast cancer"                                                                    
# [31] "prostate cancer with ‚â§10 bone metastases"                                                                                 
# [32] "lymphomaÔºåmalignant"                                                                                                       
# [33] "cervical cancer ‚â• figo iib and or lymph node metastases"                                                                  
# [34] "waldenstr√∂m macroglobulinemia (wm)"                                                                                        
# [35] "relapsedÔºèrefractory b-cell lymphoma"   

# Levenstein distance between tumors 
#df_tumor_combined<-as.data.frame(unique(c(ct_tumor_df$diseases,NCIT_Terms,WHO_Terms)))
df_tumor_combined<-as.data.frame(unique(c(ct_tumor_df$diseases,NCIT_Terms,WHO_Terms_All$Tumor_Names)))
colnames(df_tumor_combined)[1]<-"Tumor"


df_tumor_names<-unique(df_tumor_combined$Tumor)

dissimilarity_matrix_lv <- as.data.frame(matrix(nrow=length(df_tumor_names),ncol=length(df_tumor_names)))
rownames(dissimilarity_matrix_lv)<-df_tumor_names
colnames(dissimilarity_matrix_lv)<-df_tumor_names




# for (iter in 1:dim(dissimilarity_matrix_lv)[1]){
#   print(iter)
#   disease_name <- colnames(dissimilarity_matrix_lv)[iter]
#   distances<-unlist(lapply(df_tumor_names,string_dissimilarity,S2=disease_name,meth="lv"))
#   dissimilarity_matrix_lv[iter,]<-distances
# }


cl <- makeCluster(25, outfile="")
registerDoParallel(cl)


dissimilarity_matrix_lv<-foreach(iter=1:length(df_tumor_names),.combine=rbind) %dopar% {
  print(iter)
  disease_name <- colnames(dissimilarity_matrix_lv)[iter]
  distances<-unlist(lapply(df_tumor_names,string_dissimilarity,S2=disease_name,meth="lv"))
}
rownames(dissimilarity_matrix_lv) <- df_tumor_names
colnames(dissimilarity_matrix_lv) <- df_tumor_names


stopCluster(cl)

save(dissimilarity_matrix_lv,file=paste(intermediate_dir,"/dissimilarity_matrix_lv.RData",sep=""))




# Jarro Winkler Distance
cl <- makeCluster(25, outfile="")
registerDoParallel(cl)
dissimilarity_matrix_jw <- as.data.frame(matrix(nrow=length(df_tumor_names),ncol=length(df_tumor_names)))
rownames(dissimilarity_matrix_jw) <- df_tumor_names
colnames(dissimilarity_matrix_jw) <- df_tumor_names

dissimilarity_matrix_jw<-foreach(iter=1:length(df_tumor_names),.combine=rbind) %dopar% {
  print(iter)
  disease_name <- colnames(dissimilarity_matrix_jw)[iter]
  distances<-unlist(lapply(df_tumor_names,string_dissimilarity,S2=disease_name,meth="jw"))
  
}
rownames(dissimilarity_matrix_jw) <- df_tumor_names
colnames(dissimilarity_matrix_jw) <- df_tumor_names
stopCluster(cl)
save(dissimilarity_matrix_jw,file=paste(intermediate_dir,"/dissimilarity_matrix_jw.RData",sep=""))

# Cosine Distance
cl <- makeCluster(25, outfile="")
registerDoParallel(cl)
dissimilarity_matrix_cosine <- as.data.frame(matrix(nrow=length(df_tumor_names),ncol=length(df_tumor_names)))
rownames(dissimilarity_matrix_cosine) <- df_tumor_names
colnames(dissimilarity_matrix_cosine) <- df_tumor_names

dissimilarity_matrix_cosine<-foreach(iter=1:length(df_tumor_names),.combine=rbind) %dopar% {
  print(iter)
  disease_name <- colnames(dissimilarity_matrix_cosine)[iter]
  distances<-unlist(lapply(df_tumor_names,string_dissimilarity,S2=disease_name,meth="cosine"))
  
}
rownames(dissimilarity_matrix_cosine) <- df_tumor_names
colnames(dissimilarity_matrix_cosine) <- df_tumor_names
stopCluster(cl)

save(dissimilarity_matrix_cosine,file=paste(intermediate_dir,"/dissimilarity_matrix_cosine.RData",sep=""))


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

### JW

cluster_results_jw<-foreach(iter=1:length(df_tumor_names),.combine=rbind) %do% {
  print(iter)
  distance_clusters(dissimilarity_matrix_jw[iter,],cutoff = 0.0005,df_tumor_names)
}

cluster_results_jw<-as.data.frame(cluster_results_jw)
rownames(cluster_results_jw)<-NULL
cluster_results_jw <-cbind(cluster_results_jw,df_tumor_names)


colnames(cluster_results_jw) <- c("cluster_jw","Tumors")
cluster_results_jw <- cluster_results_jw %>%  dplyr::select(Tumors,cluster_jw)

# Cosine

cluster_results_cosine<-foreach(iter=1:length(df_tumor_names),.combine=rbind) %do% {
  print(iter)
  distance_clusters(dissimilarity_matrix_cosine[iter,],cutoff = 0.0005,df_tumor_names)
}

cluster_results_cosine<-as.data.frame(cluster_results_cosine)
rownames(cluster_results_cosine)<-NULL
cluster_results_cosine <-cbind(cluster_results_cosine,df_tumor_names)


colnames(cluster_results_cosine) <- c("cluster_cosine","Tumors")
cluster_results_cosine <- cluster_results_cosine %>%  dplyr::select(Tumors,cluster_cosine)


# Number of clusters
cluster_results_lv<-cluster_results_lv %>% mutate(cluster_members = str_count(cluster_lv, ";")+1)
cluster_results_jw<-cluster_results_jw %>% mutate(cluster_members = str_count(cluster_jw, ";")+1)
cluster_results_cosine<-cluster_results_cosine %>% mutate(cluster_members = str_count(cluster_cosine, ";")+1)  
  
# filter by ct tumors
cluster_results_lv_ct_filtered<- cluster_results_lv %>% filter(Tumors %in% ct_tumor_df$diseases)
cluster_results_jw_ct_filtered<- cluster_results_jw %>% filter(Tumors %in% ct_tumor_df$diseases)
cluster_results_cosine_ct_filtered<- cluster_results_cosine %>% filter(Tumors %in% ct_tumor_df$diseases)

# add nct_id
ct_tumor_df_disease_short <- ct_tumor_df %>% dplyr::select(nct_id,diseases)
cluster_results_lv_ct_filtered <- cluster_results_lv_ct_filtered %>% left_join(ct_tumor_df_disease_short,by=c("Tumors"="diseases"))
cluster_results_lv_ct_filtered <- cluster_results_lv_ct_filtered[,c(4,1:3)]


cluster_results_jw_ct_filtered <- cluster_results_jw_ct_filtered %>% left_join(ct_tumor_df_disease_short,by=c("Tumors"="diseases"))
cluster_results_jw_ct_filtered <- cluster_results_jw_ct_filtered[,c(4,1:3)]


cluster_results_cosine_ct_filtered <- cluster_results_cosine_ct_filtered %>% left_join(ct_tumor_df_disease_short,by=c("Tumors"="diseases"))
cluster_results_cosine_ct_filtered <- cluster_results_cosine_ct_filtered[,c(4,1:3)]

# organize by alphabetical order
cluster_results_lv_ct_filtered <- cluster_results_lv_ct_filtered[order(cluster_results_lv_ct_filtered$Tumors),]
cluster_results_jw_ct_filtered <- cluster_results_jw_ct_filtered[order(cluster_results_jw_ct_filtered$Tumors),]
cluster_results_cosine_ct_filtered <- cluster_results_cosine_ct_filtered[order(cluster_results_cosine_ct_filtered$Tumors),]

stopCluster(cl)


# Null the row names
rownames(cluster_results_lv_ct_filtered)<-NULL
rownames(cluster_results_jw_ct_filtered)<-NULL
rownames(cluster_results_cosine_ct_filtered)<-NULL


# Show results with following tumors 
benchmark_tumors <- c("b cell lymphoma", "neuroblastoma", "triple negative breast cancer",
                      "unresectable lung carcinoma", "liposarcoma","cancer of the liver",
                      "smoldering myeloma")

display_table_benchmark <- cbind(cluster_results_lv_ct_filtered %>% filter(Tumors %in% benchmark_tumors) %>% dplyr::select(nct_id,Tumors,cluster_lv),
                                 cluster_results_jw_ct_filtered %>% filter(Tumors %in%  benchmark_tumors) %>% dplyr::select(cluster_jw),
                                 cluster_results_cosine_ct_filtered %>% filter(Tumors %in% benchmark_tumors) %>% dplyr::select(cluster_cosine)
                                 )
                                 

display_table_benchmark <- display_table_benchmark %>% dplyr::select(nct_id,Tumors,cluster_lv,cluster_jw,cluster_cosine)

# To organize from good to worse
benchmark_tumors<- as.data.frame(benchmark_tumors)
colnames(benchmark_tumors)<-"Tumors"
display_table_benchmark <- benchmark_tumors %>% dplyr::left_join(display_table_benchmark,by="Tumors")

display_table_benchmark <- display_table_benchmark %>% dplyr::select(nct_id,Tumors,cluster_lv,cluster_jw,cluster_cosine)


# Stop the cluster
stopCluster(cl)


# Write Results of Clusters
write.csv(cluster_results_lv,paste(results_dir,"/cluster_lv.csv",sep=""))
write.csv(cluster_results_jw,paste(results_dir,"/cluster_jw.csv",sep=""))
write.csv(cluster_results_cosine,paste(results_dir,"/cluster_cosine.csv",sep=""))
write.csv(display_table_benchmark,paste(results_dir,"/edit_distance_bench_mark.csv",sep=""))

 