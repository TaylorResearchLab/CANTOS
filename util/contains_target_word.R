contains_target_word<- function (disease_term,df_tumor_leven,target_word){
  
  ind_candidate <- which(df_tumor_leven$diseases==disease_term)
  
  disease_names <- unlist(strsplit(df_tumor_leven$all_word_cluster[ind_candidate],split=";"))
  status <- target_word %in% disease_names
  
  # ind_target<-which(df_tumor$diseases==target_word)
  # 
  # if(status){
  #    if(df_tumor$PedCanTumor[ind_target]=="Yes" & df_tumor$PedCanTumor[ind_candidate]=="Yes"){
  #      if(df_tumor$Ped_Evidence[ind])
  #    }
  #   
  # }
  
  
  return(status)
}
  