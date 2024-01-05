nearest_ncit<-function (disease_embedding, embedding_df_NCIT){
  
  df<-as.data.frame(matrix(nrow=dim(embedding_df_NCIT)[1], ncol=1))
  rownames(df)<- embedding_df_NCIT$Disease
  for(iter in 1:dim(embedding_df_NCIT)[1]){
    df$V1[iter]<- dist(rbind(disease_embedding, embedding_df_NCIT[iter,2:1537]))[1]
  }
  ind_min <- which(df==min(df),arr.ind = TRUE)[1]
  ncit_info <- list(rownames(df)[ind_min], df$V1[ind_min])
  return(ncit_info)
}