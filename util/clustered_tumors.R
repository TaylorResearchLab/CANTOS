clustered_tumors<- function (target_word, candidate_words,df_tumor_leven){
  candidate_words <- unlist(strsplit(candidate_words,split=";"))
  cluster_list <- list()
  for(iter in 1: length(candidate_words)){
    disease_term<-candidate_words[iter]
    status<- contains_target_word(disease_term,df_tumor_leven,target_word)
    if(status){
     cluster_list<-append(cluster_list,disease_term)  
    }
  }
  cluster_list<-unlist(paste(cluster_list,collapse = ";"))
  return(cluster_list)
}
  