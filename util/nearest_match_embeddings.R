nearest_match_embeddings <- function(distance_matrix_df,col_name,standardization_type){
  
  embedding_df<-as.data.frame(matrix(nrow = nrow(distance_matrix_df),ncol=2))
  colnames(embedding_df)<-c("Tumor_Names",col_name)
  embedding_df$Tumor_Names<-rownames(distance_matrix_df)

  for(iter in 1:nrow(embedding_df)){
    ind_min <- which(distance_matrix_df[iter,]==min(distance_matrix_df[iter,]),arr.ind = TRUE)[2]
    names_extracted <- colnames(distance_matrix_df)[ind_min]
    if(length(names_extracted)>1){
      names_extracted<-paste(names_extracted,collapse = "*;*")
    }
    embedding_df[[col_name]][iter]<-names_extracted
    #print(iter)
  }
  return(embedding_df)  

}