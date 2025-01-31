find_euclidean_match <- function(result_df,embedding_df,tumor_data,col_names){
  
  colnames(tumor_data)[2]<-col_names[1]
  
  result_df<-result_df%>%left_join(tumor_data,by="Tumor_Names")
  result_df[col_names[2]]<-NA
  
  for(iter in 1:nrow(result_df)){
    CTR_name<-result_df$Tumor_Names[iter]
    match_name<-result_df[col_names[1]][iter,1]
  
      
    dist_vector<- embedding_df %>%filter(Tumor_Names %in% c(CTR_name,match_name))
    
    if(nrow(dist_vector)==1){
      distance<-0
    }else{
    
    distance<-dist(dist_vector[,2:ncol(dist_vector)])[1]
    }
    result_df[col_names[2]][iter,1]<-distance
  }
  return(result_df)
}