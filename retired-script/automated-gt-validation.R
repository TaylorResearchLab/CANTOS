tumor_all_review_gdrive<- tumor_all_review_gdrive[,1:28]
tumor_5thed_review_gdrive<-tumor_5thed_review_gdrive[,1:28]


load("12-prior.RData")

colnames(embedding_nearest_all)<- c("Tumor_Names","euclidean_dist_llama","euclidean_dist_biobert",
                                    "euclidean_dist_medllama","euclidean_dist_pubmedbert","euclidean_dist_modernbert")

colnames(embedding_nearest_5th)<- c("Tumor_Names","euclidean_dist_llama","euclidean_dist_biobert",
                                    "euclidean_dist_medllama","euclidean_dist_pubmedbert","euclidean_dist_modernbert")


embedding_nearest_all$valid_euclidean_dist_llama<-NA
embedding_nearest_all$valid_euclidean_dist_biobert<-NA
embedding_nearest_all$valid_euclidean_dist_medllama<-NA
embedding_nearest_all$valid_euclidean_dist_pubmedbert<-NA
embedding_nearest_all$valid_euclidean_dist_modernbert<-NA
embedding_nearest_all<-embedding_nearest_all[,c(1,2,7,3,8,4,9,5,10,6,11)]

embedding_nearest_5th$valid_euclidean_dist_llama<-NA
embedding_nearest_5th$valid_euclidean_dist_biobert<-NA
embedding_nearest_5th$valid_euclidean_dist_medllama<-NA
embedding_nearest_5th$valid_euclidean_dist_pubmedbert<-NA
embedding_nearest_5th$valid_euclidean_dist_modernbert<-NA
embedding_nearest_5th<-embedding_nearest_5th[,c(1,2,7,3,8,4,9,5,10,6,11)]

tumor_all_review_gdrive <- tumor_all_review_gdrive %>% left_join(embedding_nearest_all, by="Tumor_Names")
tumor_5thed_review_gdrive <- tumor_5thed_review_gdrive %>% left_join(embedding_nearest_5th, by="Tumor_Names")


for(iter in 1:1600){
  ground_truths<- tumor_all_review_gdrive$ground_truth_val[iter]
  ground_truths<-unique(unlist(strsplit(ground_truths,";")))
  
  
  if(tumor_all_review_gdrive$euclidean_dist_llama[iter] %in% ground_truths){
    tumor_all_review_gdrive$valid_euclidean_dist_llama[iter]=1
  }else{
    tumor_all_review_gdrive$valid_euclidean_dist_llama[iter]=0
  }
  
  if(tumor_all_review_gdrive$euclidean_dist_biobert[iter] %in% ground_truths ){
    tumor_all_review_gdrive$valid_euclidean_dist_biobert[iter]=1
  }else{
    tumor_all_review_gdrive$valid_euclidean_dist_biobert[iter]=0
  }
  
  if(tumor_all_review_gdrive$euclidean_dist_medllama[iter] %in% ground_truths){
    tumor_all_review_gdrive$valid_euclidean_dist_medllama[iter]=1
  }else{
    tumor_all_review_gdrive$valid_euclidean_dist_medllama[iter]=0
  }
  
  if(tumor_all_review_gdrive$euclidean_dist_pubmedbert[iter] %in% ground_truths){
    tumor_all_review_gdrive$valid_euclidean_dist_pubmedbert[iter]=1
  }else{
    tumor_all_review_gdrive$valid_euclidean_dist_pubmedbert[iter]=0
  }
  if(tumor_all_review_gdrive$euclidean_dist_modernbert[iter] %in% ground_truths){
    tumor_all_review_gdrive$valid_euclidean_dist_modernbert[iter]=1
  }else{
    tumor_all_review_gdrive$valid_euclidean_dist_modernbert[iter]=0
  }
}





for(iter in 1:1600){
  ground_truths<- tumor_5thed_review_gdrive$ground_truth_val[iter]
  ground_truths<-unique(unlist(strsplit(ground_truths,";")))
  
  
  if(tumor_5thed_review_gdrive$euclidean_dist_llama[iter] %in% ground_truths){
    tumor_5thed_review_gdrive$valid_euclidean_dist_llama[iter]=1
  }else{
    tumor_5thed_review_gdrive$valid_euclidean_dist_llama[iter]=0
  }
  
  if(tumor_5thed_review_gdrive$euclidean_dist_biobert[iter] %in% ground_truths ){
    tumor_5thed_review_gdrive$valid_euclidean_dist_biobert[iter]=1
  }else{
    tumor_5thed_review_gdrive$valid_euclidean_dist_biobert[iter]=0
  }
  
  if(tumor_5thed_review_gdrive$euclidean_dist_medllama[iter] %in% ground_truths){
    tumor_5thed_review_gdrive$valid_euclidean_dist_medllama[iter]=1
  }else{
    tumor_5thed_review_gdrive$valid_euclidean_dist_medllama[iter]=0
  }
  
  if(tumor_5thed_review_gdrive$euclidean_dist_pubmedbert[iter] %in% ground_truths){
    tumor_5thed_review_gdrive$valid_euclidean_dist_pubmedbert[iter]=1
  }else{
    tumor_5thed_review_gdrive$valid_euclidean_dist_pubmedbert[iter]=0
  }
  if(tumor_5thed_review_gdrive$euclidean_dist_modernbert[iter] %in% ground_truths){
    tumor_5thed_review_gdrive$valid_euclidean_dist_modernbert[iter]=1
  }else{
    tumor_5thed_review_gdrive$valid_euclidean_dist_modernbert[iter]=0
  }
}

write.csv(tumor_all_review_gdrive,"analysis/results/tumor_manually_validated_tmp.csv")
write.csv(tumor_5thed_review_gdrive,"analysis/results_5th/tumor_manually_validated_tmp_5th.csv")