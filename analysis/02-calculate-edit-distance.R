# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(doParallel)
  library(foreach)
  library(stringdist)
})


# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
intermediate_dir <- file.path(analysis_dir,"intermediate")

#Load Functions
source(paste(util_dir,"/string_dissimilarity.R",sep = ""))

# Read the annotated file
ct_disease_df <- read.csv(paste(input_dir,"/CT-Aug22-2023-Disease-File - clinical_trial_disease_aug_22_2023.csv",sep=""))
ct_tumor_df<- ct_disease_df %>% filter(validated_cancer_tumor=="Yes")

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
df_tumor_leven<-as.data.frame(unique(ct_tumor_df$diseases))
colnames(df_tumor_leven)[1]<-"Tumor"


df_tumor_names<-unique(df_tumor_leven$Tumor)

dissimilarity_matrix_lv <- as.data.frame(matrix(nrow=length(df_tumor_names),ncol=length(df_tumor_names)))
rownames(dissimilarity_matrix_lv)<-df_tumor_names
colnames(dissimilarity_matrix_lv)<-df_tumor_names


# for (iter in 1:dim(dissimilarity_matrix_lv)[1]){
#   print(iter)
#   disease_name <- colnames(dissimilarity_matrix_lv)[iter]
#   distances<-unlist(lapply(df_tumor_names,string_dissimilarity,S2=disease_name,meth="lv"))
#   dissimilarity_matrix_lv[iter,]<-distances
# }


cl <- makeCluster(5, outfile="")
registerDoParallel(cl)

dissimilarity_matrix_lv<-foreach(iter=1:length(df_tumor_names),.combine=rbind) %dopar% {
  print(iter)
  disease_name <- colnames(dissimilarity_matrix_lv)[iter]
  distances<-unlist(lapply(df_tumor_names,string_dissimilarity,S2=disease_name,meth="lv"))
  
}
rownames(dissimilarity_matrix_lv) <- df_tumor_names
colnames(dissimilarity_matrix_lv) <- df_tumor_names

# Jarro Winkler Distance
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

# Cosine Distance
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

